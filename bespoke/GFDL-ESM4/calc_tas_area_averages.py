import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools
)

# Basic diagnostic name (details such as time averaging, zonal
# mean, etc., and the corresponding netCDF variable name, are
# added automatically as specified in the netcdf module):
diag_name = "tas"

nc_long_name = ("Near-surface air temperature, averaged "
    + "annually and everywhere {}ward of {}")

nc_var_attrs_n = {
    "units"        : nf.field_units["temperature"],
    "standard_name": "air_temperature",
    "cell_methods" : f"{nf.nc_time_name}: mean "
                     + f"{nf.nc_ref_lat_n_name}: mean "
                     + "(area-weighted)"
}

nc_var_attrs_s = nc_var_attrs_n.copy()
nc_var_attrs_s["cell_methods"] = (f"{nf.nc_time_name}: "
    + f"mean {nf.nc_ref_lat_s_name}: mean (area-weighted)")

# Short description added to netCDF "title" attribute (need
# not be completely accurate/detailed here):
nc_title_str = "near-surface air temperature polar-cap averages"


def main():
    
    prsr = script_tools.default_cmd_argument_parser()
    prsr.add_argument("-a", "--approx", action="store_true",
                      help="Do approximate version ("
                           + f"{nf.diag_nq_cell_center_approx})"
                           + " anyway")
    cmd = prsr.parse_args()
    cmd.model = "GFDL-ESM4"
    
    if cmd.experiment not in ["piControl", "historical"]:
        raise ValueError("This script only works for experiment"
            + f" piControl or historical, not {cmd.experiment}")
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    nt = yr_e - yr_s + 1
    n_ens = len(ens_members)
    ny, nx = md.grid_dims_atm[cmd.model]
    
    ref_lats_n = md.default_ref_lats_n
    ref_lats_s = md.default_ref_lats_s
    
    ref_lats_n_bnds = np.array([[ref_lats_n[k], 90.0]
        for k in range(len(ref_lats_n))])
    ref_lats_s_bnds = np.array([[-90.0, ref_lats_s[k]]
        for k in range(len(ref_lats_s))])
    
    lon, lat, areacella = lpd.areacella(cmd.model)
    lon_bnds, lat_bnds  = lpd.bndscella(cmd.model)
    
    # The following kwargs happen to be the same for the data
    # loaded and saved in this script -- but the latter needs
    # extra information (later):
    diag_kw = {"name": diag_name,
        "time_methods": nf.diag_nq_yearly}
    
    nc_var_kw = {"name": diag_name,
        "time_methods": nf.nc_var_nq_yearly}
    
    load_kw = {
        "diagnostic_id": nf.diag_name(**diag_kw),
        "nc_field_name": nf.nc_var_name(**nc_var_kw),
        "model_id"     : cmd.model,
        "experiment_id": cmd.experiment
    }
    
    tas = lpd.field_2D(**load_kw)[-1]
    
    save_nc_kw = {
        "model_id"     : cmd.model,
        "member_ids"   : ens_members,
        "experiment_id": cmd.experiment,
        "year_range"   : (yr_s, yr_e),
        "nc_global_attrs": {
            nf.nc_file_attrs_experiment_name: \
                f"esm-{cmd.experiment}"},
        "nc_title_str" : nc_title_str
    }
    
    # Above: overwrite global attribute experiment_id with esm
    # variant (but save with same directory/filename structure)
    
    diag_kw["space_methods"] = nf.diag_nq_area_mean
    nc_var_kw["space_methods"] = nf.nc_var_nq_area_mean
    
    if np.ndim(lat) == 1:
        
        true_lats_n, true_lats_n_bnds, true_tas_avg_n = \
            diags.integrate_horizontally_exact(tas, lat_bnds,
                areacella, hemi="n", normalise=True)
        
        true_lats_s, true_lats_s_bnds, true_tas_avg_s = \
            diags.integrate_horizontally_exact(tas, lat_bnds,
                areacella, hemi="s", normalise=True)
        
        interp_tas_avg_n = diags.interpolate_to_ref_latitudes(
            true_lats_n, true_tas_avg_n, ref_lats_n)
        interp_tas_avg_s = diags.interpolate_to_ref_latitudes(
            true_lats_s, true_tas_avg_s, ref_lats_s)
        
        # --------------------------------------------------- #
        
        print("Saving to NetCDF...")
        
        diag_kw["other_methods"] = nf.diag_nq_native
        nc_var_kw["other_methods"] = nf.nc_var_nq_native
        
        nf.save_yearly_ref_lat(
            true_tas_avg_n, true_tas_avg_s,
            true_lats_n, true_lats_s,
            nf.diag_name(**diag_kw),
            nf.nc_var_name(hemi="n", **nc_var_kw),
            nf.nc_var_name(hemi="s", **nc_var_kw),
            ref_lat_n_bnds=true_lats_n_bnds,
            ref_lat_s_bnds=true_lats_s_bnds,
            unlimited_ref_lat_dim=False,
            nc_field_attrs_n={
                "long_name": nc_long_name.format("north",
                    "native-grid latitude bounds"),
                **nc_var_attrs_n},
            nc_field_attrs_s={
                "long_name": nc_long_name.format("south",
                    "native-grid latitude bounds"),
                **nc_var_attrs_s},
            **save_nc_kw
        )
        
        diag_kw["other_methods"] = nf.diag_nq_native_interp
        nc_var_kw["other_methods"] = nf.nc_var_nq_native_interp
        
        nf.save_yearly_ref_lat(
            interp_tas_avg_n, interp_tas_avg_s,
            ref_lats_n, ref_lats_s,
            nf.diag_name(**diag_kw),
            nf.nc_var_name(hemi="n", **nc_var_kw),
            nf.nc_var_name(hemi="s", **nc_var_kw),
            ref_lat_n_bnds=ref_lats_n_bnds,
            ref_lat_s_bnds=ref_lats_s_bnds,
            nc_field_attrs_n={
                "long_name": nc_long_name.format("north",
                    "native-grid latitude bounds, interpolated "
                    + "to reference latitudes"),
                **nc_var_attrs_n},
            nc_field_attrs_s={
                "long_name": nc_long_name.format("south",
                    "native-grid latitude bounds, interpolated "
                    + "to reference latitudes"),
                **nc_var_attrs_s},
            **save_nc_kw
        )
        
    
    if np.ndim(lat) != 1 or cmd.approx:
        
        # Have to do it the "old fashioned" way
        # (integrate cells poleward of reference latitude,
        # based on cell center coordinates)
        # 
        # Either that, or want the approximate version anyway
        # in which case need to create mesh grid latitude if
        # it is 1D:
        if np.ndim(lat) == 1:
            lon, lat = np.meshgrid(lon, lat)
        
        # This function returns both the integral and mean,
        # so just take second return value:
        tas_avg_n = diags.integrate_horizontally_multi_members(
            tas, areacella, lat, ref_lats=ref_lats_n,
            hemi="n")[1]
        
        tas_avg_s = diags.integrate_horizontally_multi_members(
            tas, areacella, lat, ref_lats=ref_lats_s,
            hemi="s")[1]
        
        # --------------------------------------------------- #
        
        print("Saving to NetCDF...")
        
        diag_kw["other_methods"] = \
            nf.diag_nq_cell_center_approx
        nc_var_kw["other_methods"] = \
            nf.nc_var_nq_cell_center_approx
        
        nf.save_yearly_ref_lat(tas_avg_n, tas_avg_s,
            ref_lats_n, ref_lats_s,
            nf.diag_name(**diag_kw),
            nf.nc_var_name(hemi="n", **nc_var_kw),
            nf.nc_var_name(hemi="s", **nc_var_kw),
            ref_lat_n_bnds=ref_lats_n_bnds,
            ref_lat_s_bnds=ref_lats_s_bnds,
            nc_field_attrs_n={
                "long_name": nc_long_name.format("north",
                    "reference latitudes based on cell-center "
                    + "latitudes"),
                **nc_var_attrs_n},
            nc_field_attrs_s={
                "long_name": nc_long_name.format("south",
                    "reference latitudes based on cell-center "
                    + "latitudes"),
                **nc_var_attrs_s},
            **save_nc_kw
        )



if __name__ == '__main__':
    main()
