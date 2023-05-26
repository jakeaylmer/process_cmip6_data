import numpy as np

from src import (
    diagnostics as diags,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools
)

# Diagnostic names:
diag_name = "hfds"
nc_var_name = "hfds"

nc_long_name = ("Surface downward heat flux in sea water "
    + "(hfds), averaged annually and{} everywhere {}ward of {}")

nc_var_attrs_mean_n = dict()
nc_var_attrs_mean_n["standard_name"] = \
    "surface_downward_heat_flux_in_sea_water"
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
    
    prsr = script_tools.default_cmd_argument_parser()
    prsr.add_argument("-a", "--approx", action="store_true",
                      help="Do approximate version ("
                           + f"{nf.diag_nq_cell_center_approx})"
                           + " anyway")
    cmd = prsr.parse_args()
    cmd.model = "GISS-E2-2-G"
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    nt = yr_e - yr_s + 1
    n_ens = len(ens_members)
    
    ref_lats_n = md.default_ref_lats_n
    ref_lats_s = md.default_ref_lats_s
    
    ref_lats_n_bnds = np.array([[ref_lats_n[k], 90.0]
        for k in range(len(ref_lats_n))])
    ref_lats_s_bnds = np.array([[-90.0, ref_lats_s[k]]
        for k in range(len(ref_lats_s))])
    
    # Yes, this should be areacella not areacello
    # (see README.txt):
    lon, lat, areacello = lpd.areacella(cmd.model)
    lon_bnds, lat_bnds  = lpd.bndscella(cmd.model)
    
    diag_kw = {"name": "hfds",
        "time_methods": nf.diag_nq_yearly}
    
    nc_var_kw = {"name": "hfds",
        "time_methods": nf.nc_var_nq_yearly}
    
    hfds = lpd.field_2D(
        nf.diag_name(**diag_kw),
        nf.nc_var_name(**nc_var_kw),
        cmd.model, cmd.experiment)[-1]
    
    # No land mask in areacella (need to add it for integrals):
    areacello *= np.where(np.all(np.isnan(hfds), axis=(0,1)),
                          np.nan, 1.0)
    
    save_nc_kw = {
        "model_id"     : cmd.model,
        "member_ids"   : ens_members,
        "experiment_id": cmd.experiment,
        "year_range"   : (yr_s, yr_e)
    }
    
    true_lats_n, true_lats_n_bnds, true_hfds_avg_n = \
        diags.integrate_horizontally_exact(hfds, lat_bnds,
            areacello, hemi="n", normalise=True)
    
    true_lats_s, true_lats_s_bnds, true_hfds_avg_s = \
        diags.integrate_horizontally_exact(hfds, lat_bnds,
            areacello, hemi="s", normalise=True)
    
    _, _, true_hfds_hint_n = \
        diags.integrate_horizontally_exact(hfds, lat_bnds,
            areacello, hemi="n", normalise=False)
    
    _, _, true_hfds_hint_s = \
        diags.integrate_horizontally_exact(hfds, lat_bnds,
            areacello, hemi="s", normalise=False)
    
    true_hfds_hint_n /= 1.0E15  # in PW
    true_hfds_hint_s /= 1.0E15
    
    interp_hfds_avg_n = diags.interpolate_to_ref_latitudes(
        true_lats_n, true_hfds_avg_n, ref_lats_n)
    
    interp_hfds_avg_s = diags.interpolate_to_ref_latitudes(
        true_lats_s, true_hfds_avg_s, ref_lats_s)
    
    interp_hfds_hint_n = diags.interpolate_to_ref_latitudes(
        true_lats_n, true_hfds_hint_n, ref_lats_n)
    
    interp_hfds_hint_s = diags.interpolate_to_ref_latitudes(
        true_lats_s, true_hfds_hint_s, ref_lats_s)
    
    diag_kw["space_methods"] = nf.diag_nq_area_mean
    diag_kw["other_methods"] = nf.diag_nq_native
    nc_var_kw["space_methods"] = nf.nc_var_nq_area_mean
    nc_var_kw["other_methods"] = nf.nc_var_nq_native
    
    nf.save_yearly_ref_lat(true_hfds_avg_n, true_hfds_avg_s,
        true_lats_n, true_lats_s,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(hemi="n", **nc_var_kw),
        nf.nc_var_name(hemi="s", **nc_var_kw),
        ref_lat_n_bnds=true_lats_n_bnds,
        ref_lat_s_bnds=true_lats_s_bnds,
        nc_field_attrs_n={
            "units": nf.field_units["heatflux"],
            "long_name": nc_long_name.format("", "north",
                "native-grid latitude bounds"),
            **nc_var_attrs_mean_n},
        nc_field_attrs_s={
            "units": nf.field_units["heatflux"],
            "long_name": nc_long_name.format("", "south",
                "native-grid latitude bounds"),
            **nc_var_attrs_mean_s},
        **save_nc_kw
    )
    
    diag_kw["space_methods"] = nf.diag_nq_area_integral
    nc_var_kw["space_methods"] = nf.nc_var_nq_area_integral
    
    nf.save_yearly_ref_lat(true_hfds_hint_n, true_hfds_hint_s,
        true_lats_n, true_lats_s,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(hemi="n", **nc_var_kw),
        nf.nc_var_name(hemi="s", **nc_var_kw),
        ref_lat_n_bnds=true_lats_n_bnds,
        ref_lat_s_bnds=true_lats_s_bnds,
        nc_field_attrs_n={
            "units": nf.field_units["heattransport"],
            "long_name": nc_long_name.format(" integrated",
                "north", "native-grid latitude bounds"),
            **nc_var_attrs_mean_n},
        nc_field_attrs_s={
            "units": nf.field_units["heattransport"],
            "long_name": nc_long_name.format(" integrated",
                "south", "native-grid latitude bounds"),
            **nc_var_attrs_mean_s},
        **save_nc_kw
    )
    
    diag_kw["space_methods"] = nf.diag_nq_area_mean
    diag_kw["other_methods"] = nf.diag_nq_native_interp
    nc_var_kw["space_methods"] = nf.nc_var_nq_area_mean
    nc_var_kw["other_methods"] = nf.nc_var_nq_native_interp
    
    nf.save_yearly_ref_lat(interp_hfds_avg_n, interp_hfds_avg_s,
        ref_lats_n, ref_lats_s,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(hemi="n", **nc_var_kw),
        nf.nc_var_name(hemi="s", **nc_var_kw),
        ref_lat_n_bnds=ref_lats_n_bnds,
        ref_lat_s_bnds=ref_lats_s_bnds,
        nc_field_attrs_n={
            "units": nf.field_units["heatflux"],
            "long_name": nc_long_name.format("", "north",
                    "native-grid latitude bounds, interpolated "
                    + "to reference latitudes"),
            **nc_var_attrs_mean_n},
        nc_field_attrs_s={
            "units": nf.field_units["heatflux"],
            "long_name": nc_long_name.format("", "south",
                "native-grid latitude bounds, interpolated "
                    + "to reference latitudes"),
            **nc_var_attrs_mean_s},
        **save_nc_kw
    )
    
    diag_kw["space_methods"] = nf.diag_nq_area_integral
    nc_var_kw["space_methods"] = nf.nc_var_nq_area_integral
    
    nf.save_yearly_ref_lat(interp_hfds_hint_n,
        interp_hfds_hint_s,
        ref_lats_n, ref_lats_s,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(hemi="n", **nc_var_kw),
        nf.nc_var_name(hemi="s", **nc_var_kw),
        ref_lat_n_bnds=ref_lats_n_bnds,
        ref_lat_s_bnds=ref_lats_s_bnds,
        nc_field_attrs_n={
            "units": nf.field_units["heattransport"],
            "long_name": nc_long_name.format(" integrated",
                "north", "native-grid latitude bounds, "
                    + "interpolated to reference latitudes"),
            **nc_var_attrs_mean_n},
        nc_field_attrs_s={
            "units": nf.field_units["heattransport"],
            "long_name": nc_long_name.format(" integrated",
                "south", "native-grid latitude bounds, "
                    + "interpolated to reference latitudes"),
            **nc_var_attrs_mean_s},
        **save_nc_kw
    )
    
    
    
    if cmd.approx:
        
        lon, lat = np.meshgrid(lon, lat)
        
        hfds_hint_n, hfds_mean_n = \
            diags.integrate_horizontally_multi_members(
                hfds, areacello, lat, ref_lats=ref_lats_n,
                hemi="n")
        
        hfds_hint_s, hfds_mean_s = \
            diags.integrate_horizontally_multi_members(
                hfds, areacello, lat, ref_lats=ref_lats_s,
                hemi="s")
        
        hfds_hint_n /= 1.0E15  # in petawatts
        hfds_hint_s /= 1.0E15
        
        # ------------------------------------------------------- #
        print("Saving to NetCDF...")
        
        diag_kw["space_methods"] = nf.diag_nq_area_mean
        diag_kw["other_methods"] = \
            nf.diag_nq_cell_center_approx
        
        nc_var_kw["space_methods"] = nf.nc_var_nq_area_mean
        nc_var_kw["other_methods"] = \
            nf.nc_var_nq_cell_center_approx
        
        nf.save_yearly_ref_lat(hfds_mean_n, hfds_mean_s,
            ref_lats_n, ref_lats_s,
            nf.diag_name(**diag_kw),
            nf.nc_var_name(hemi="n", **nc_var_kw),
            nf.nc_var_name(hemi="s", **nc_var_kw),
            ref_lat_n_bnds=ref_lats_n_bnds,
            ref_lat_s_bnds=ref_lats_s_bnds,
            nc_field_attrs_n={
                "units": nf.field_units["heatflux"],
                "long_name": nc_long_name.format("", "north",
                    "reference latitudes based on cell-center "
                        + "coordinates"),
                **nc_var_attrs_mean_n},
            nc_field_attrs_s={
                "units": nf.field_units["heatflux"],
                "long_name": nc_long_name.format("", "south",
                    "reference latitudes based on cell-center "
                        + "coordinates"),
                **nc_var_attrs_mean_s},
            **save_nc_kw
        )
        
        diag_kw["space_methods"] = nf.diag_nq_area_integral
        nc_var_kw["space_methods"] = \
            nf.nc_var_nq_area_integral
        
        nf.save_yearly_ref_lat(hfds_hint_n, hfds_hint_s,
            ref_lats_n, ref_lats_s,
            nf.diag_name(**diag_kw),
            nf.nc_var_name(hemi="n", **nc_var_kw),
            nf.nc_var_name(hemi="s", **nc_var_kw),
            ref_lat_n_bnds=ref_lats_n_bnds,
            ref_lat_s_bnds=ref_lats_s_bnds,
            nc_field_attrs_n={
                "units": nf.field_units["heattransport"],
                "long_name": nc_long_name.format(" integrated",
                    "north", "reference latitudes based on "
                        + "cell-center coordinates"),
                **nc_var_attrs_int_n},
            nc_field_attrs_s={
                "units": nf.field_units["heattransport"],
                "long_name": nc_long_name.format(" integrated",
                    "south", "reference latitudes based on "
                        + "cell-center coordinates"),
                **nc_var_attrs_int_s},
            **save_nc_kw
        )


if __name__ == '__main__':
    main()
