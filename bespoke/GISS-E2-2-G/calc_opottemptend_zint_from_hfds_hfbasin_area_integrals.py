import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools
)

# Diagnostic names:
diag_name = ("opottemptend" + nf.diag_nq_vertical_integral
    + "_from_hfds_hfbasin")
nc_var_name = "opottemptend" + nf.nc_var_nq_vertical_integral

nc_standard_name = ("tendency_of_sea_water_potential_temperatur"
                    + "e_expressed_as_heat_content")

nc_long_name = (
    "Ocean potential temperature expressed as heat content "
    + "(opottemptend) integrated vertically and everywhere "
    + "{}ward of {}, where the depth-integrated opottemptend "
    + "is calculated as the residual of net flux into the "
    + "ocean column (hfds) and native-grid northward ocean "
    + "heat transport (hfbasin)")

nc_var_attrs_n = dict()
nc_var_attrs_n["cell_methods"] = (
    f"{nf.nc_time_name}: mean " +
    f"{nf.nc_ref_lat_n_name}: sum (area-weighted)")

nc_var_attrs_s = nc_var_attrs_n.copy()
nc_var_attrs_s["cell_methods"] = (
    f"{nf.nc_time_name}: mean " + 
    f"{nf.nc_ref_lat_s_name}: sum (area-weighted)")


def main():
    
    cmd = script_tools.default_cmd_args()
    cmd.model = "GISS-E2-2-G"
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    # There are two different versions to do here:
    # 
    # (1) Calculated on the native grid and interpolated to
    #     standard reference latitudes
    # (2) Calculated approximately based on cell-center
    #     latitudes (as most other models)
    # 
    # Load hfbasin first (this is the same in both cases):
    # 
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
    
    # Load hfds native interp (for method 1):
    diag_name_kw = {"name": "hfds",
        "time_methods": nf.diag_nq_yearly,
        "space_methods": nf.diag_nq_area_integral,
        "other_methods": nf.diag_nq_native_interp}
    
    nc_var_name_kw = {"name": "hfds",
        "time_methods": nf.nc_var_nq_yearly,
        "space_methods": nf.nc_var_nq_area_integral,
        "other_methods": nf.nc_var_nq_native_interp}
    
    ref_lats_n, ref_lats_s, _, _, _, _, hfds_n, hfds_s \
        = lpd.time_series_ref_lats(nf.diag_name(**diag_name_kw),
            nf.nc_var_name(hemi="n", **nc_var_name_kw),
            nf.nc_var_name(hemi="s", **nc_var_name_kw),
            cmd.model, cmd.experiment)
    
    ref_lats_n_bnds = np.array([[ref_lats_n[k], 90.0]
        for k in range(len(ref_lats_n))])
    ref_lats_s_bnds = np.array([[-90.0, ref_lats_s[k]]
        for k in range(len(ref_lats_s))])
    
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
    
    # Compute opottemptend_zint (native, interp):
    # 
    # These are horizontal integrals of implied vertical
    # integrals (variable names are getting too long):
    dhdt_gn_n = np.zeros((nt, n_ens, n_ref_lats_n))
    dhdt_gn_s = np.zeros((nt, n_ens, n_ref_lats_s))
    
    for k in range(n_ref_lats_n):
        j = np.argmin(abs(lat_hfbasin - ref_lats_n[k]))
        dhdt_gn_n[:,:,k] = hfds_n[:,:,k] + hfbasin[:,:,j]
    
    for k in range(n_ref_lats_s):
        j = np.argmin(abs(lat_hfbasin - ref_lats_s[k]))
        dhdt_gn_s[:,:,k] = hfds_s[:,:,k] - hfbasin[:,:,j]
    
    
    # Repeat with the cell-center approx. version of hfds:
    diag_name_kw["other_methods"] = \
        nf.diag_nq_cell_center_approx
    nc_var_name_kw["other_methods"] = \
        nf.nc_var_nq_cell_center_approx
    
    _, _, _, _, _, _, hfds_n, hfds_s \
        = lpd.time_series_ref_lats(nf.diag_name(**diag_name_kw),
            nf.nc_var_name(hemi="n", **nc_var_name_kw),
            nf.nc_var_name(hemi="s", **nc_var_name_kw),
            cmd.model, cmd.experiment)
    
    dhdt_cc_n = np.zeros((nt, n_ens, n_ref_lats_n))
    dhdt_cc_s = np.zeros((nt, n_ens, n_ref_lats_s))
    
    for k in range(n_ref_lats_n):
        j = np.argmin(abs(lat_hfbasin - ref_lats_n[k]))
        dhdt_cc_n[:,:,k] = hfds_n[:,:,k] + hfbasin[:,:,j]
    
    for k in range(n_ref_lats_s):
        j = np.argmin(abs(lat_hfbasin - ref_lats_s[k]))
        dhdt_cc_s[:,:,k] = hfds_s[:,:,k] - hfbasin[:,:,j]
    
    # ------------------------------------------------------- #
    
    save_nc_kw = {
        "model_id"     : cmd.model,
        "member_ids"   : ens_members,
        "experiment_id": cmd.experiment,
        "year_range"   : (yr_s, yr_e)
    }
    
    print("Saving to NetCDF...")
    
    diag_name_kw = {
        "name": diag_name,
        "time_methods": nf.diag_nq_yearly,
        "space_methods": nf.diag_nq_area_integral,
        "other_methods": nf.diag_nq_native_interp}
    
    nc_var_name_kw = {
        "name": nc_var_name,
        "time_methods": nf.nc_var_nq_yearly,
        "space_methods": nf.nc_var_nq_area_integral,
        "other_methods": nf.nc_var_nq_native_interp}
    
    nf.save_yearly_ref_lat(dhdt_gn_n, dhdt_gn_s,
        ref_lats_n, ref_lats_s,
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(hemi="n", **nc_var_name_kw),
        nf.nc_var_name(hemi="s", **nc_var_name_kw),
        ref_lat_n_bnds=ref_lats_n_bnds,
        ref_lat_s_bnds=ref_lats_s_bnds,
        nc_field_attrs_n={
            "units": nf.field_units["heattransport"],
            "long_name": nc_long_name.format("north",
                "native-grid latitude bounds, interpolated "
                    + "to reference latitudes"),
            "standard_name": nc_standard_name,
                **nc_var_attrs_n},
        nc_field_attrs_s={
            "units": nf.field_units["heattransport"],
            "long_name": nc_long_name.format("south",
                "native-grid latitude bounds, interpolated "
                    + "to reference latitudes"),
            "standard_name": nc_standard_name,
                **nc_var_attrs_s},
        **save_nc_kw
    )
    
    diag_name_kw["other_methods"] = \
        nf.diag_nq_cell_center_approx
    
    nc_var_name_kw["other_methods"] = \
        nf.nc_var_nq_cell_center_approx
    
    nf.save_yearly_ref_lat(dhdt_cc_n, dhdt_cc_s,
        ref_lats_n, ref_lats_s,
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(hemi="n", **nc_var_name_kw),
        nf.nc_var_name(hemi="s", **nc_var_name_kw),
        ref_lat_n_bnds=ref_lats_n_bnds,
        ref_lat_s_bnds=ref_lats_s_bnds,
        nc_field_attrs_n={
            "units": nf.field_units["heattransport"],
            "long_name": nc_long_name.format("north",
                "reference latitudes based on cell-center "
                    + "coordinates"),
            "standard_name": nc_standard_name,
                **nc_var_attrs_n},
        nc_field_attrs_s={
            "units": nf.field_units["heattransport"],
            "long_name": nc_long_name.format("south",
                "reference latitudes based on cell-center "
                    + "coordinates"),
            "standard_name": nc_standard_name,
                **nc_var_attrs_s},
        **save_nc_kw
    )



if __name__ == '__main__':
    main()
