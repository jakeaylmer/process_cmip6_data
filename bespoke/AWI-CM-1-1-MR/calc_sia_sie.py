import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_raw_data as lrd,
    metadata as md,
    netcdf as nf,
    script_tools
)

from ice_edge_latitude.utilities import regions as iel_regions


# Basic diagnostic names (details such as time averaging, zonal
# mean, etc., and the corresponding netCDF variable name, are
# added automatically as specified in the nf module):
diag_area_name   = "sia"
diag_extent_name = "sie"

nc_extent_long_name = (
    "{}ern hemisphere sea ice extent (area of "
    f"ocean where siconc >= {nf.nc_siconc_threshold_name})"
    + "{}")

nc_area_long_name = "{}ern hemisphere sea ice area{}"

nc_var_extent_attrs = {
    "standard_name": "sea_ice_extent",
    "cell_methods" : f"{nf.nc_time_name}: mean",
    "coordinates"  : nf.nc_siconc_threshold_name,  # scalar coordinate variable
    "units"        : nf.field_units["seaiceextent"]
}

nc_var_area_attrs = {
    "standard_name": "sea_ice_area",
    "cell_methods" : f"{nf.nc_time_name}: mean",
    "units"        : nf.field_units["seaicearea"]
}


def main():
    
    prsr = script_tools.default_cmd_argument_parser()
    prsr.add_argument("-t", "--threshold", type=float,
                      default=md.default_sie_threshold)
    cmd = prsr.parse_args()
    cmd.model = "AWI-CM-1-1-MR"
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    nt = yr_e - yr_s + 1
    n_ens = len(ens_members)
    
    # Prepare siconc here rather than using the
    # general routine lrd.prepare_siconc()
    
    # Use coordinates from areacello rather than load for each
    # ensemble member:
    lon, lat, areacello = lrd.areacello(cmd.model)
    
    
    # Can also determine the lake masks here to speed up;
    # lake_mask = 1 if allowed, nan if masked:
    print("Determining lake mask")
    lon = lon % 360.0
    lake_mask = np.ones(np.shape(areacello), dtype=np.float64)
    
    for reg in ["lakes", "baltic_sea", "black_sea"]:
        
        reg_def = getattr(iel_regions, reg)
        
        for sub_reg in reg_def.keys():
            lake_mask *= np.where(
                (lon >= reg_def[sub_reg][0]) &
                (lon <= reg_def[sub_reg][1]) &
                (lat >= reg_def[sub_reg][2]) &
                (lat <= reg_def[sub_reg][3]),
                np.nan, 1.0)
    
    load_kw = {
        "model_id": cmd.model,
        "experiment_id": cmd.experiment
    }
    
    sia_n = np.zeros((nt*12, n_ens), dtype=np.float64)
    sia_s = np.zeros((nt*12, n_ens), dtype=np.float64)
    sie_n = np.zeros((nt*12, n_ens), dtype=np.float64)
    sie_s = np.zeros((nt*12, n_ens), dtype=np.float64)
    
    lat_mask_n = np.where(lat > 0.0, 1.0, 0.0)
    lat_mask_s = np.where(lat < 0.0, 1.0, 0.0)
    
    for m in range(n_ens):
        
        print(f"Processing member: {ens_members[m]} "
            + f"({m+1} of {n_ens})")
        
        # This model has siconc as percentage:
        siconc_m = lrd.field_2D("siconc",
            member_id=ens_members[m], **load_kw)[-1]/100.0
        siconc_m = np.maximum(0.0, siconc_m)
        siconc_m = np.minimum(1.0, siconc_m)
        siconc_m *= lake_mask
        
        print(f"Calculating sea ice extent and area")
        
        sie_mask = np.where(siconc_m >= cmd.threshold,
            1.0, 0.0)
        sia_mask = siconc_m
        
        sie_n[:,m] = np.nansum(
            sie_mask*areacello[np.newaxis,:]*lat_mask_n[np.newaxis,:],
            axis=1)
        
        sie_s[:,m] = np.nansum(
            sie_mask*areacello[np.newaxis,:]*lat_mask_s[np.newaxis,:],
            axis=1)
        
        sia_n[:,m] = np.nansum(
            sia_mask*areacello[np.newaxis,:]*lat_mask_n[np.newaxis,:],
            axis=1)
        
        sia_s[:,m] = np.nansum(
            sia_mask*areacello[np.newaxis,:]*lat_mask_s[np.newaxis,:],
            axis=1)
    
    sie_n /= 1.0E12  # put in 1e6 km2
    sie_s /= 1.0E12
    sia_n /= 1.0E12
    sia_s /= 1.0E12
    
    sie_n_yr = diags.year_mean_time_series_multi_member(sie_n)
    sie_s_yr = diags.year_mean_time_series_multi_member(sie_s)
    sia_n_yr = diags.year_mean_time_series_multi_member(sia_n)
    sia_s_yr = diags.year_mean_time_series_multi_member(sia_s)
    
    # ------------------------------------------------------- #
    
    print("Saving to NetCDF...")
    
    save_nc_kw = {
        "model_id": cmd.model,
        "member_ids": ens_members,
        "experiment_id": cmd.experiment,
        "year_range": (yr_s, yr_e)
    }
    
    sie_coord_var_attrs = {
        "name": nf.nc_siconc_threshold_name,
        "dtype": np.float64,
        "value": cmd.threshold,
        "attrs": {
            "standard_name": "sea_ice_area_fraction",
            "units": "1"
        }
    }
    
    diag_kw = {
        "name": diag_area_name,
        "time_methods": nf.diag_nq_monthly
    }
    
    nc_var_kw = {
        "name": diag_area_name,
        "time_methods": nf.nc_var_nq_monthly
    }
    
    nf.save_time_series_by_hemi(sia_n, sia_s,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(hemi="n", **nc_var_kw),
        nf.nc_var_name(hemi="s", **nc_var_kw),
        nc_field_type=sia_n.dtype,
        nc_field_attrs_n={
            "long_name": nc_area_long_name.format("North", ""),
            **nc_var_area_attrs},
        nc_field_attrs_s={
            "long_name": nc_area_long_name.format("South", ""),
            **nc_var_area_attrs},
        monthly=True,
        **save_nc_kw
    )
    
    diag_kw["time_methods"] = nf.diag_nq_yearly
    nc_var_kw["time_methods"] = nf.nc_var_nq_yearly
    
    nf.save_time_series_by_hemi(sia_n_yr, sia_s_yr,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(hemi="n", **nc_var_kw),
        nf.nc_var_name(hemi="s", **nc_var_kw),
        nc_field_type=sia_n.dtype,
        nc_field_attrs_n={
            "long_name": nc_area_long_name.format("North",
                ", annually averaged"),
            **nc_var_area_attrs},
        nc_field_attrs_s={
            "long_name": nc_area_long_name.format("South",
                ", annually averaged"),
            **nc_var_area_attrs},
        **save_nc_kw
    )
    
    diag_kw = {
        "name": diag_extent_name,
        "time_methods": nf.diag_nq_monthly
    }
    
    nc_var_kw = {
        "name": diag_extent_name,
        "time_methods": nf.nc_var_nq_monthly
    }
    
    nf.save_time_series_by_hemi(sie_n, sie_s,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(hemi="n", **nc_var_kw),
        nf.nc_var_name(hemi="s", **nc_var_kw),
        nc_field_type=sie_n.dtype,
        nc_field_attrs_n={
            "long_name": nc_extent_long_name.format("North",""),
            **nc_var_extent_attrs},
        nc_field_attrs_s={
            "long_name": nc_extent_long_name.format("South",""),
            **nc_var_extent_attrs},
        scalar_coord_var_n_kw=sie_coord_var_attrs,
        monthly=True,
        **save_nc_kw
    )
    
    diag_kw["time_methods"] = nf.diag_nq_yearly
    nc_var_kw["time_methods"] = nf.nc_var_nq_yearly
    
    nf.save_time_series_by_hemi(sie_n_yr, sie_s_yr,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(hemi="n", **nc_var_kw),
        nf.nc_var_name(hemi="s", **nc_var_kw),
        nc_field_type=sie_n.dtype,
        nc_field_attrs_n={
            "long_name": nc_extent_long_name.format("North",
                ", annually averaged"),
            **nc_var_extent_attrs},
        nc_field_attrs_s={
            "long_name": nc_extent_long_name.format("South",
                ", annually averaged"),
            **nc_var_extent_attrs},
        scalar_coord_var_n_kw=sie_coord_var_attrs,
        **save_nc_kw
    )
    



if __name__ == "__main__":
    main()
