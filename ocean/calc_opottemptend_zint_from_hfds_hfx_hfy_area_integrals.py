import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools
)

# Diagnostic names:
diag_name = ("o{}temptend" + nf.diag_nq_vertical_integral
    + "_from_hfds_hfx_hfy")
nc_var_name = "o{}temptend" + nf.nc_var_nq_vertical_integral

nc_standard_name = ("tendency_of_sea_water_{}_temperature_"
                    + "expressed_as_heat_content")

nc_long_name = (
    "Ocean {} temperature expressed as heat content "
    + "(o{}temptend) integrated vertically and{} everywhere "
    + "{}ward of reference latitudes based on cell-center "
    + "coordinates")

nc_var_attrs_mean_n = dict()
nc_var_attrs_mean_n["cell_methods"] = (
    f"{nf.nc_time_name}: mean " +
    f"{nf.nc_ref_lat_n_name}: mean (area-weighted)")

nc_var_attrs_mean_s = nc_var_attrs_mean_n.copy()
nc_var_attrs_mean_s["cell_methods"] = (
    f"{nf.nc_time_name}: mean " + 
    f"{nf.nc_ref_lat_s_name}: mean (area-weighted)")

nc_var_attrs_int_n = nc_var_attrs_mean_n.copy()
nc_var_attrs_int_n["cell_methods"] = (
    f"{nf.nc_time_name}: mean {nf.nc_ref_lat_n_name}: sum "
    + "(area-weighted)")

nc_var_attrs_int_s = nc_var_attrs_int_n.copy()
nc_var_attrs_int_s["cell_methods"] = (
    f"{nf.nc_time_name}: mean {nf.nc_ref_lat_s_name}: sum "
    + "(area-weighted)")

def main():
    
    cmd = script_tools.default_cmd_args()
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    ref_lats_n = md.default_ref_lats_n
    ref_lats_s = md.default_ref_lats_s
    
    ref_lats_n_bnds = np.array([[ref_lats_n[k], 90.0]
        for k in range(len(ref_lats_n))])
    ref_lats_s_bnds = np.array([[-90.0, ref_lats_s[k]]
        for k in range(len(ref_lats_s))])
    
    try:
        ftype = md.ocn_prognostic_temperature[cmd.model]
    except KeyError:
        ftype = md.ocn_prognostic_temperature_when_does_not_exist
    
    ftype_long = "potential" if ftype=="pot" else "conservative"
    
    _, lat, areacello = lpd.areacello(cmd.model)
    
    diag_name_kw = {"name": "hfds",
        "time_methods": nf.diag_nq_yearly}
    
    nc_var_name_kw = {"name": "hfds",
        "time_methods": nf.nc_var_nq_yearly}
    
    hfds = lpd.field_2D(
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(**nc_var_name_kw),
        cmd.model, cmd.experiment)[-1]
    
    diag_name_kw["name"] = "ohtc_from_hfx_hfy"
    nc_var_name_kw["name"] = "ohtc"
    
    ohtc = lpd.field_2D(
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(**nc_var_name_kw),
        cmd.model, cmd.experiment)[-1]
    
    # -------------------------------------------------------- #
    # Set missing to zero in case missing values (NaN, land) do
    # not match up for whatever reason (necessary to close
    # budget properly with the other terms that are calculated
    # from the same fields but separately). In this case, for
    # the mean values in W m-2, need a separate land mask
    # otherwise normalisation will be over all grid cells:
    # 
    # While I've not as of yet checked this thoroughly, I assume
    # the reason for different location of NaN in the hfds and
    # ohtc-from-hfx-hfy data is because the latter will
    # potentially have rows of missing data at the north pole,
    # for example, depending on the underlying ocean grid and
    # calculation method.
    # 
    # Also, it therefore makes sense to use the hfds data to
    # determine the overall land mask for the normalised
    # diagnostics, rather than the ohtc data.
    # -------------------------------------------------------- #
    land_mask = np.where(np.all(np.isnan(hfds),
                                axis=(0, 1)),
                         0.0, 1.0)
    hfds = np.where(np.isnan(hfds), 0.0, hfds)
    ohtc = np.where(np.isnan(ohtc), 0.0, ohtc)
    
    opottemptend_zint = hfds + ohtc
    
    load_kw = {
        "model_id": cmd.model,
        "experiment_id":cmd.experiment
    }
    
    opottemptend_zint_hint_n, opottemptend_zint_mean_n = \
        diags.integrate_horizontally_multi_members(
            opottemptend_zint, areacello, lat,
            ref_lats=ref_lats_n, hemi="n", land_mask=land_mask)
    
    opottemptend_zint_hint_s, opottemptend_zint_mean_s = \
        diags.integrate_horizontally_multi_members(
            opottemptend_zint, areacello, lat,
            ref_lats=ref_lats_s, hemi="s", land_mask=land_mask)
    
    opottemptend_zint_hint_n /= 1.0E15  # in petawatts
    opottemptend_zint_hint_s /= 1.0E15
    
    # ------------------------------------------------------- #
    
    save_nc_kw = {
        "model_id"     : cmd.model,
        "member_ids"   : ens_members,
        "experiment_id": cmd.experiment,
        "year_range"   : (yr_s, yr_e)
    }
    
    print("Saving to NetCDF...")
    
    diag_name_kw["name"] = diag_name.format(ftype)
    diag_name_kw["space_methods"] = nf.diag_nq_area_mean
    diag_name_kw["other_methods"] = \
        nf.diag_nq_cell_center_approx
    
    nc_var_name_kw["name"] = nc_var_name.format(ftype)
    nc_var_name_kw["space_methods"] = nf.nc_var_nq_area_mean
    nc_var_name_kw["other_methods"] = \
        nf.nc_var_nq_cell_center_approx
    
    nf.save_yearly_ref_lat(
        opottemptend_zint_mean_n,
        opottemptend_zint_mean_s,
        ref_lats_n, ref_lats_s,
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(hemi="n", **nc_var_name_kw),
        nf.nc_var_name(hemi="s", **nc_var_name_kw),
        ref_lat_n_bnds=ref_lats_n_bnds,
        ref_lat_s_bnds=ref_lats_s_bnds,
        nc_field_attrs_n={
            "units": nf.field_units["heatflux"],
            "long_name": nc_long_name.format(ftype_long, ftype,
                " averaged", "north"),
            "standard_name": nc_standard_name.format(ftype_long),
            **nc_var_attrs_mean_n},
        nc_field_attrs_s={
            "units": nf.field_units["heatflux"],
            "long_name": nc_long_name.format(ftype_long, ftype,
                " averaged", "south"),
            "standard_name": nc_standard_name.format(ftype_long),
            **nc_var_attrs_mean_s},
        **save_nc_kw
    )
    
    diag_name_kw["space_methods"] = nf.diag_nq_area_integral
    nc_var_name_kw["space_methods"] = nf.nc_var_nq_area_integral
    
    nf.save_yearly_ref_lat(
        opottemptend_zint_hint_n,
        opottemptend_zint_hint_s,
        ref_lats_n, ref_lats_s,
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(hemi="n", **nc_var_name_kw),
        nf.nc_var_name(hemi="s", **nc_var_name_kw),
        ref_lat_n_bnds=ref_lats_n_bnds,
        ref_lat_s_bnds=ref_lats_s_bnds,
        nc_field_attrs_n={
            "units": nf.field_units["heattransport"],
            "long_name": nc_long_name.format(ftype_long, ftype,
                "", "north"),
            "standard_name": nc_standard_name.format(ftype_long),
            **nc_var_attrs_int_n},
        nc_field_attrs_s={
            "units": nf.field_units["heattransport"],
            "long_name": nc_long_name.format(ftype_long, ftype,
                "", "south"),
            "standard_name": nc_standard_name.format(ftype_long),
            **nc_var_attrs_int_s},
        **save_nc_kw
    )



if __name__ == '__main__':
    main()
