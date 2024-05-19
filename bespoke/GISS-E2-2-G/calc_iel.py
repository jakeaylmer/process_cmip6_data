from pathlib import Path
import numpy as np

from process_cmip6_data.src.diagnostics import (
    year_mean_time_series_multi_member as ymtsmm)

from process_cmip6_data.src import (
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools,
    utils
)

from ice_edge_latitude.diagnostics import ice_edge_latitude as iel
from ice_edge_latitude.utilities import regions as iel_regions

# Basic diagnostic names (details such as time averaging, zonal
# mean, etc., and the corresponding nf variable name, are
# added automatically as specified in the nf module):
diag_name = "iel"

nc_long_name = (
    "{}ern hemisphere{} sea ice-edge latitude{}")

nc_var_attrs = {
    "cell_methods": f"{nf.nc_time_name}: mean {nf.nc_lon_1d_name}: mean",
    "coordinates": f"{nf.nc_siconc_threshold_name} {nf.nc_lon_1d_name}",
    "standard_name": "latitude",
    "units": nf.field_units["seaiceedge"]
}

nc_var_zm_attrs = {
    "cell_methods": f"{nf.nc_time_name}: mean",
    "coordinates": nf.nc_siconc_threshold_name,  # scalar coordinate variable
    "standard_name": "latitude",
    "units"      : nf.field_units["seaiceedge"]
}


def main():
    
    prsr = script_tools.default_cmd_argument_parser()
    prsr.add_argument("-t", "--threshold", type=float,
                      default=md.default_iel_threshold)
    # -------------------------------------------------------- #
    # The following arguments must be passed from the bash
    # script (as these settings have already been applied or
    # defined during remapping phase):
    # -------------------------------------------------------- #
    prsr.add_argument("--datadir", type=str,
        default="/storage/basic/cpom/gb919150/CMIP6/_swap")
    cmd = prsr.parse_args()
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    nt = yr_e - yr_s + 1  # number of years
    n_ens = len(ens_members)
    
    # No remapping done, so can load original (atmosphere grid)
    # coordinates in advance:
    lon, lat, _ = lpd.areacella(cmd.model)
    lon_bnds, _ = lpd.bndscella(cmd.model)
    n_lon = len(lon)
    n_lat = len(lat)
    
    # Need 2D coordinates for lake masking and 0-360 lon range
    # for lake masking:
    lon_m, lat_m = np.meshgrid(lon, lat)
    lon_m = lon_m % 360.0
    
    # Determine lake mask:
    lake_mask = np.ones((n_lat, n_lon), dtype=np.float64)
    for reg in ["lakes", "baltic_sea", "black_sea"]:
        
        reg_def = getattr(iel_regions, reg)
        
        for sub_reg in reg_def.keys():
            sub_reg_mask = np.where(
                (lon_m >= reg_def[sub_reg][0]) &
                (lon_m <= reg_def[sub_reg][1]) &
                (lat_m >= reg_def[sub_reg][2]) &
                (lat_m <= reg_def[sub_reg][3]),
                np.nan, 1.0)
            
            lake_mask *= sub_reg_mask
    
    # Pre-allocate ice edge as a function of longitude -- here
    # the number of longitudes is the same as the original grid:
    iel_n = np.zeros((nt*12, n_ens, n_lon))
    iel_s = np.zeros((nt*12, n_ens, n_lon))
    
    for m in range(n_ens):
        
        # Prepare siconc here rather than using the
        # general routine lpd.prepare_siconc()
        print(f"Preparing masked siconc ({cmd.model}, "
            + f"{cmd.experiment}, {ens_members[m]})")
        
        data_files_in = sorted(
            [str(x) for x in Path(cmd.datadir).glob(
                    f"masked_siconca_*{cmd.model}*"
                    + f"{cmd.experiment}*{ens_members[m]}*.nc"
                )
            ])
        
        siconc_m, = utils.nc_get_arrays(data_files_in, [],
                                        ["siconca"])
        
        siconc_m = np.where(
            siconc_m >= md.default_original_missing_value,
            np.nan, siconc_m)
        
        siconc_m /= 100.0  # [siconca] = % for GISS-E2-2-G
        
        # Restrict range to [0.0, 1.0]:
        siconc_m = np.maximum(0.0, siconc_m)
        siconc_m = np.minimum(1.0, siconc_m)
        
        siconc_m *= lake_mask[np.newaxis,:,:]
        
        print("Calculating sea ice edge latitudes")
        iel_n[:,m,:],iel_s[:,m,:] = iel.ice_edge(lat, siconc_m)
        
    iel_zm_n = np.nanmean(iel_n, axis=2)
    iel_zm_s = np.nanmean(iel_s, axis=2)
    
    iel_zm_n_yr = ymtsmm(iel_zm_n, keep_nan=False)
    iel_zm_s_yr = ymtsmm(iel_zm_s, keep_nan=False)
    
    # ------------------------------------------------------- #
    
    print("Saving to NetCDF...")
    
    save_nc_kw = {
        "model_id": cmd.model,
        "member_ids": ens_members,
        "experiment_id": cmd.experiment,
        "year_range": (yr_s, yr_e)
    }
    
    scalar_coord_var_attrs = {
        "name": nf.nc_siconc_threshold_name,
        "dtype": np.float64,
        "value": cmd.threshold,
        "attrs": {
            "standard_name": "sea_ice_area_fraction",
            "units": "1"
        }
    }
    
    # Applies to full and zonal-mean data:
    nc_var_attrs["comment"] = (
        "Computed on the native grid from siconca")
    
    nc_var_zm_attrs["comment"] = nc_var_attrs["comment"]
    
    diag_kw = {
        "name": diag_name,
        "time_methods": nf.diag_nq_monthly,
        "other_methods": nf.diag_nq_native
    }
    
    nc_var_kw = {
        "name": diag_name,
        "time_methods": nf.nc_var_nq_monthly,
        "other_methods": nf.nc_var_nq_native
    }
    
    
    # Ice edge as a function of longitude:
    nf.save_ice_edge_latitude_per_longitude(
        iel_n, iel_s, lon,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(hemi="n", **nc_var_kw),
        nf.nc_var_name(hemi="s", **nc_var_kw),
        nc_field_type=iel_n.dtype,
        lon_bnds=lon_bnds,
        monthly=True,
        nc_field_attrs_n={
            "long_name": nc_long_name.format("North", "", ""),
            **nc_var_attrs},
        nc_field_attrs_s={
            "long_name": nc_long_name.format("South", "", ""),
            **nc_var_attrs},
        scalar_coord_var_n_kw=scalar_coord_var_attrs,
        nc_title_str="sea ice-edge latitude",
        **save_nc_kw
    )
    
    diag_kw["space_methods"] = nf.diag_nq_zonal_mean
    nc_var_kw["space_methods"] = nf.nc_var_nq_zonal_mean
    
    nf.save_time_series_by_hemi(iel_zm_n, iel_zm_s,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(hemi="n", **nc_var_kw),
        nf.nc_var_name(hemi="s", **nc_var_kw),
        nc_field_type=iel_zm_n.dtype,
        nc_field_attrs_n={
            "long_name": nc_long_name.format("North",
                " zonal-mean", ""),
            **nc_var_zm_attrs},
        nc_field_attrs_s={
            "long_name": nc_long_name.format("South",
                " zonal-mean", ""),
            **nc_var_zm_attrs},
        scalar_coord_var_n_kw=scalar_coord_var_attrs,
        monthly=True,
        nc_title_str="sea ice-edge latitude zonal mean",
        **save_nc_kw
    )
    
    diag_kw["time_methods"] = nf.diag_nq_yearly
    nc_var_kw["time_methods"] = nf.nc_var_nq_yearly
    
    nf.save_time_series_by_hemi(iel_zm_n_yr, iel_zm_s_yr,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(hemi="n", **nc_var_kw),
        nf.nc_var_name(hemi="s", **nc_var_kw),
        nc_field_type=iel_zm_n.dtype,
        nc_field_attrs_n={
            "long_name": nc_long_name.format("North",
                " zonal-mean", ", annually averaged"),
            **nc_var_zm_attrs},
        nc_field_attrs_s={
            "long_name": nc_long_name.format("South",
                " zonal-mean", ", annually averaged"),
            **nc_var_zm_attrs},
        scalar_coord_var_n_kw=scalar_coord_var_attrs,
        monthly=False,
        nc_title_str="sea ice-edge latitude zonal mean",
        **save_nc_kw
    )



if __name__ == "__main__":
    main()
