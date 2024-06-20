"""Calculate monthly and yearly-mean, zonal and zonally-averaged
sea ice-edge latitudes from sea ice concentration. There is a
companion bash script, calc_iel.sh, that first interpolates sea
ice concentration to a regular longitude-latitude grid of
appropriate resolution, which this python script and the
calculation of the ice edge requires. That bash script should be
run, not this script directly.

The ice-edge latitude is calculated using an external python
package, "ice_edge_latitude" (doi:10.5281/zenodo.5494523) that
must be available on the python PATH for this to work. The
default settings of that procedure are used.
"""

import numpy as np

from process_cmip6_data.src.diagnostics import (
    year_mean_time_series_multi_member as ymtsmm)

from process_cmip6_data.src import (
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools)

from ice_edge_latitude.diagnostics import ice_edge_latitude as iel

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
    "units": nf.field_units["seaiceedge"]}

nc_var_zm_attrs = {
    "cell_methods": f"{nf.nc_time_name}: mean",
    "coordinates": nf.nc_siconc_threshold_name,  # scalar coordinate variable
    "standard_name": "latitude",
    "units"      : nf.field_units["seaiceedge"]}

# Interpolation descriptions
# (these are followed by " interpolation"):
remap_methods = {
    "bil" : "bilinear",
    "con" : "first-order conservative",
    "con2": "second-order conservative",
    "dis" : "distance-weighted",
    "nn"  : "nearest-neighbour"}


def main():
    
    prsr = script_tools.default_cmd_argument_parser()
    prsr.add_argument("-t", "--threshold", type=float,
                      default=md.default_iel_threshold)
    # -------------------------------------------------------- #
    # The following arguments must be passed from the bash
    # script (as these settings have already been applied or
    # defined during remapping phase):
    # -------------------------------------------------------- #
    prsr.add_argument("-g", "--gridresolution", type=str,
                      default="1")
    prsr.add_argument("-r", "--remapmethod", type=str,
                      default="bil")
    prsr.add_argument("--datadir", type=str,
        default="/storage/basic/cpom/gb919150/CMIP6/_swap")
    cmd = prsr.parse_args()
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    nt = yr_e - yr_s + 1  # number of years
    n_ens = len(ens_members)
    
    load_kw = {
        "model_id": cmd.model,
        "experiment_id": cmd.experiment,
        "remapped_data_dir": cmd.datadir,
        "remap_method": cmd.remapmethod,
        "remap_res": cmd.gridresolution}
    
    # Ice edge as a function of longitude (determine number of
    # longitudes from the interpolated grid resolution):
    lon_res = float(cmd.gridresolution)
    n_lon = int(360 // lon_res)
    iel_n = np.zeros((nt*12, n_ens, n_lon))
    iel_s = np.zeros((nt*12, n_ens, n_lon))
    
    for m in range(n_ens):
        
        print(f"Processing member: {ens_members[m]} "
            + f"({m+1} of {n_ens})")
        
        lon_m, lat_m, siconc_m = lpd.prepare_siconc_remapped(
            member_id=ens_members[m], **load_kw)
        
        print("Calculating sea ice edge latitudes")
        iel_n[:,m,:],iel_s[:,m,:] = iel.ice_edge(lat_m,siconc_m)
        
    iel_zm_n = np.nanmean(iel_n, axis=2)
    iel_zm_s = np.nanmean(iel_s, axis=2)
    
    iel_zm_n_yr = ymtsmm(iel_zm_n, keep_nan=False)
    iel_zm_s_yr = ymtsmm(iel_zm_s, keep_nan=False)
    
    lon_bnds = np.array([
        [lon_m[k] - lon_res/2.0, lon_m[k] + lon_res/2.0]
        for k in range(n_lon)
    ])
    
    # ------------------------------------------------------- #
    
    print("Saving to NetCDF...")
    
    save_nc_kw = {
        "model_id": cmd.model,
        "member_ids": ens_members,
        "experiment_id": cmd.experiment,
        "year_range": (yr_s, yr_e)}
    
    scalar_coord_var_attrs = {
        "name": nf.nc_siconc_threshold_name,
        "dtype": np.float64,
        "value": cmd.threshold,
        "attrs": {
            "standard_name": "sea_ice_area_fraction",
            "units": "1"}
    }
    
    # Applies to full and zonal-mean data:
    nc_var_attrs["comment"] = (
        "Computed from siconc re-mapped to a "
        + f"{cmd.gridresolution} degree global rectilinear "
        + f"grid using {remap_methods[cmd.remapmethod]} "
        + "interpolation")
    
    nc_var_zm_attrs["comment"] = nc_var_attrs["comment"]
    
    def res_str(x):
        if x == "0.5":
            return "05deg"
        elif x == "0.25":
            return "025deg"
        else:
            return f"{x}deg"
    
    other_methods = getattr(nf, "diag_nq_"
        + f"{res_str(cmd.gridresolution)}_{cmd.remapmethod}")
    
    diag_kw = {
        "name": diag_name,
        "time_methods": nf.diag_nq_monthly,
        "other_methods": other_methods}
    
    # Currently set empty "" anyway:
    other_methods = getattr(nf, "nc_var_nq_"
        + f"{res_str(cmd.gridresolution)}_{cmd.remapmethod}")
    
    nc_var_kw = {
        "name": diag_name,
        "time_methods": nf.nc_var_nq_monthly,
        "other_methods": other_methods}
    
    
    # Ice edge as a function of longitude:
    nf.save_ice_edge_latitude_per_longitude(
        iel_n, iel_s, lon_m,
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
        **save_nc_kw)
    
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
        **save_nc_kw)
    
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
        **save_nc_kw)


if __name__ == "__main__":
    main()
