"""Calculate spatial integrals of ocean heat content tendency
between reference latitudes and the north or south pole, from
that already computed for the net downward heat flux into the
ocean (hfds) and native northward ocean heat transport data
(hfbasin). Data for the yearly-averaged and area-integrated
hfds, and hfbasin, should already be processed and saved.
"""

import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools)

# Diagnostic names:
diag_name = ("o{}temptend" + nf.diag_nq_vertical_integral
    + "_from_hfds_hfbasin")
nc_var_name = "o{}temptend" + nf.nc_var_nq_vertical_integral

nc_standard_name = ("tendency_of_sea_water_{}_temperature_"
                    + "expressed_as_heat_content")

nc_long_name = (
    "Ocean {} temperature expressed as heat content "
    + "(o{}temptend) integrated vertically and everywhere "
    + "{}ward of reference latitudes based on cell-center "
    + "coordinates, where the depth-integrated o{}temptend is "
    + "calculated as the residual of net flux into the ocean "
    + "column (hfds) and native-grid northward ocean heat "
    + "transport (hfbasin)")

nc_var_attrs_n = dict()
nc_var_attrs_n["cell_methods"] = (
    f"{nf.nc_time_name}: mean " +
    f"{nf.nc_ref_lat_n_name}: sum (area-weighted)")

nc_var_attrs_s = nc_var_attrs_n.copy()
nc_var_attrs_s["cell_methods"] = (
    f"{nf.nc_time_name}: mean " + 
    f"{nf.nc_ref_lat_s_name}: sum (area-weighted)")

# Short description added to netCDF "title" attribute (need
# not be completely accurate/detailed here):
nc_title_str = "ocean column heat content tendency polar-cap integrals"


def main():
    
    cmd = script_tools.default_cmd_args()
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    try:
        ftype = md.ocn_prognostic_temperature[cmd.model]
    except KeyError:
        ftype = md.ocn_prognostic_temperature_when_does_not_exist
    
    ftype_long = "potential" if ftype=="pot" else "conservative"
    
    diag_name_kw = {"name": "hfds",
        "time_methods": nf.diag_nq_yearly,
        "space_methods": nf.diag_nq_area_integral,
        "other_methods": nf.diag_nq_cell_center_approx}
    
    nc_var_name_kw = {"name": "hfds",
        "time_methods": nf.nc_var_nq_yearly,
        "space_methods": nf.nc_var_nq_area_integral,
        "other_methods": nf.nc_var_nq_cell_center_approx}
    
    ref_lats_n, ref_lats_s, _, _, _, _, hfds_n, hfds_s \
        = lpd.time_series_ref_lats(nf.diag_name(**diag_name_kw),
            nf.nc_var_name(hemi="n", **nc_var_name_kw),
            nf.nc_var_name(hemi="s", **nc_var_name_kw),
            cmd.model, cmd.experiment)
    
    ref_lats_n_bnds = np.array([[ref_lats_n[k], 90.0]
        for k in range(len(ref_lats_n))])
    ref_lats_s_bnds = np.array([[-90.0, ref_lats_s[k]]
        for k in range(len(ref_lats_s))])
    
    diag_name_kw = {"name": "oht_from_hfbasin",
        "time_methods": nf.diag_nq_yearly,
        "other_methods": nf.diag_nq_native_interp}
    nc_var_name_kw = {"name": "oht",
        "time_methods": nf.nc_var_nq_yearly,
        "other_methods": nf.nc_var_nq_native_interp}
    
    lat_hfbasin, _, _, _, _, hfbasin = \
        lpd.time_series_single_ref_lat(
            nf.diag_name(**diag_name_kw),
            nf.nc_var_name(hemi="", **nc_var_name_kw),
            cmd.model, cmd.experiment)
    
    nt, n_ens, n_ref_lats_n = np.shape(hfds_n)
    _, _, n_ref_lats_s = np.shape(hfds_s)
    
    if np.shape(hfbasin)[0] != nt:
        raise Exception(f"Time axis length of hfbasin "
            + f"({np.shape(hfbasin)[0]}) does not match that "
            + f"that of hfds ({nt})")
    
    if np.shape(hfbasin)[1] != n_ens:
        raise Exception(f"Member axis length of hfbasin "
            + f"({np.shape(hfbasin)[1]}) does not match that "
            + f"that of hfds ({n_ens})")
    
    # These are horizontal integrals of implied vertical
    # integrals (variable names are getting too long):
    otemptend_n = np.zeros((nt, n_ens, n_ref_lats_n))
    otemptend_s = np.zeros((nt, n_ens, n_ref_lats_s))
    
    for k in range(n_ref_lats_n):
        j = np.argmin(abs(lat_hfbasin - ref_lats_n[k]))
        otemptend_n[:,:,k] = hfds_n[:,:,k] + hfbasin[:,:,j]
    
    for k in range(n_ref_lats_s):
        j = np.argmin(abs(lat_hfbasin - ref_lats_s[k]))
        otemptend_s[:,:,k] = hfds_s[:,:,k] - hfbasin[:,:,j]
    
    # ------------------------------------------------------- #
    
    save_nc_kw = {
        "model_id"     : cmd.model,
        "member_ids"   : ens_members,
        "experiment_id": cmd.experiment,
        "year_range"   : (yr_s, yr_e),
        "nc_title_str" : nc_title_str}
    
    print("Saving to NetCDF...")
    
    diag_name_kw = {
        "name": diag_name.format(ftype),
        "time_methods": nf.diag_nq_yearly,
        "space_methods": nf.diag_nq_area_integral,
        "other_methods": nf.diag_nq_cell_center_approx}
    
    nc_var_name_kw = {
        "name": nc_var_name.format(ftype),
        "time_methods": nf.nc_var_nq_yearly,
        "space_methods": nf.nc_var_nq_area_integral,
        "other_methods": nf.nc_var_nq_cell_center_approx}
    
    nf.save_yearly_ref_lat(otemptend_n, otemptend_s,
        ref_lats_n, ref_lats_s,
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(hemi="n", **nc_var_name_kw),
        nf.nc_var_name(hemi="s", **nc_var_name_kw),
        ref_lat_n_bnds=ref_lats_n_bnds,
        ref_lat_s_bnds=ref_lats_s_bnds,
        nc_field_attrs_n={
            "units": nf.field_units["heattransport"],
            "long_name": nc_long_name.format(ftype_long, ftype,
                "north", ftype),
            "standard_name": nc_standard_name.format(ftype_long),
            **nc_var_attrs_n},
        nc_field_attrs_s={
            "units": nf.field_units["heattransport"],
            "long_name": nc_long_name.format(ftype_long, ftype,
                "south", ftype),
            "standard_name": nc_standard_name.format(ftype_long),
            **nc_var_attrs_s},
        **save_nc_kw)


if __name__ == "__main__":
    main()
