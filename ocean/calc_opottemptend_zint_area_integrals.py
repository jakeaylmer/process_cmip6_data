import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools
)

# Diagnostic names:
diag_name = "o{}temptend"

nc_var_name = "o{}temptend"

nc_standard_name = ("tendency_of_sea_water_{}_temperature_"
                    + "expressed_as_heat_content")

nc_long_name = (
    "Ocean {} temperature expressed as heat content "
    + "(o{}temptend) integrated (summed) vertically and{} "
    + "everywhere {}ward of reference latitudes based on cell-"
    + "center coordinates")

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
    
    if cmd.model == "MRI-ESM2-0" and cmd.experiment == "ssp370":
        print("WARNING: selecting only first ens. member in list")
        ens_members = ens_members[:1]
    
    ftype = md.ocn_prognostic_temperature[cmd.model]
    ftype_long = "potential" if ftype == "pot" else "conservative"
    
    _, lat, areacello = lpd.areacello(cmd.model)
    
    diag_name_kw = {"name": diag_name.format(ftype),
        "time_methods": nf.diag_nq_yearly,
        "space_methods": nf.diag_nq_vertical_integral}
    
    nc_var_name_kw = {"name": nc_var_name.format(ftype),
        "time_methods": nf.nc_var_nq_yearly,
        "space_methods": nf.nc_var_nq_vertical_integral}
    
    opottemptend_zint = lpd.field_2D(
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(**nc_var_name_kw),
        cmd.model, cmd.experiment)[-1]
    
    # Note: do not need to remove or otherwise worry about NaN
    # ----- (land masks) here as we only load one field. The
    #       integrate routine automatically accounts for NaN and
    #       normalises correctly.
    
    opottemptend_zint_hint_n, opottemptend_zint_mean_n = \
        diags.integrate_horizontally_multi_members(
            opottemptend_zint, areacello, lat,
            ref_lats=ref_lats_n, hemi="n")
    
    opottemptend_zint_hint_s, opottemptend_zint_mean_s = \
        diags.integrate_horizontally_multi_members(
            opottemptend_zint, areacello, lat,
            ref_lats=ref_lats_s, hemi="s")
    
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
    
    diag_name_kw["space_methods"] += nf.diag_nq_area_mean
    diag_name_kw["other_methods"] = \
        nf.diag_nq_cell_center_approx
    nc_var_name_kw["space_methods"] += nf.nc_var_nq_area_mean
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
                                             " averaged",
                                             "north"),
            "standard_name": nc_standard_name.format(ftype_long),
            **nc_var_attrs_mean_n},
        nc_field_attrs_s={
            "units": nf.field_units["heatflux"],
            "long_name": nc_long_name.format(ftype_long, ftype,
                                             " averaged",
                                             "south"),
            "standard_name": nc_standard_name.format(ftype_long),
            **nc_var_attrs_mean_s},
        **save_nc_kw
    )
    
    diag_name_kw["space_methods"] = \
        nf.diag_nq_vertical_integral + nf.diag_nq_area_integral
    nc_var_name_kw["space_methods"] = (
        nf.nc_var_nq_vertical_integral
        + nf.nc_var_nq_area_integral)
    
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


if __name__ == "__main__":
    main()
