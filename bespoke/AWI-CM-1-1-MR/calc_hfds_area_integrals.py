import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_raw_data as lrd,
    metadata as md,
    netcdf as nf,
    script_tools
)

# Diagnostic names:
diag_name = "hfds"
nc_var_name = "hfds"

nc_long_name = ("Surface downward heat flux in sea water "
    + "(hfds), averaged annually and{} everywhere {}ward of "
    + "reference latitudes based on cell-center coordinates")

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

# Short description added to netCDF "title" attribute (need
# not be completely accurate/detailed here):
nc_title_str = "heat flux into ocean surface polar-cap "
nc_title_str_int = nc_title_str + "integrals"
nc_title_str_mean = nc_title_str + "averages"

def main():
    
    cmd = script_tools.default_cmd_args()
    cmd.model = "AWI-CM-1-1-MR"
    
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
    
    # For this model, these are 1D arrays:
    _, lat, areacello = lrd.areacello(cmd.model)
    n_cells = len(areacello)
    
    hfds = np.zeros((nt, n_ens, n_cells), dtype=np.float64)
    
    for m in range(n_ens):
        # Download monthly data, compute yearly averages
        # separately from this script before loading here [also
        # 1D (+ time -> 2D)]:
        hfds[:,m,:] = lrd.field_2D("hfds",
            member_id=ens_members[m], model_id=cmd.model,
            experiment_id=cmd.experiment)[-1]
    
    # Add dimension of size 1 so the diagnostic routines still
    # work:
    lat = np.expand_dims(lat, axis=1)
    areacello = np.expand_dims(areacello, axis=1)
    hfds = np.expand_dims(hfds, axis=3)
    
    hfds_hint_n, hfds_mean_n = \
        diags.integrate_horizontally_multi_members(
            hfds, areacello, lat, ref_lats=ref_lats_n, hemi="n"
        )
    
    hfds_hint_s, hfds_mean_s = \
        diags.integrate_horizontally_multi_members(
            hfds, areacello, lat, ref_lats=ref_lats_s, hemi="s"
        )
    
    hfds_hint_n /= 1.0E15  # in petawatts
    hfds_hint_s /= 1.0E15
    
    # ------------------------------------------------------- #
    
    save_nc_kw = {
        "model_id"     : cmd.model,
        "member_ids"   : ens_members,
        "experiment_id": cmd.experiment,
        "year_range"   : (yr_s, yr_e),
        "nc_title_str" : nc_title_str
    }
    
    print("Saving to NetCDF...")
    
    diag_name_kw = {
        "name": "hfds",
        "time_methods": nf.diag_nq_yearly,
        "space_methods": nf.diag_nq_area_mean,
        "other_methods": nf.diag_nq_cell_center_approx}
    
    nc_var_name_kw = {
        "name": "hfds",
        "time_methods": nf.nc_var_nq_yearly,
        "space_methods": nf.nc_var_nq_area_mean,
        "other_methods": nf.nc_var_nq_cell_center_approx}
    
    nf.save_yearly_ref_lat(hfds_mean_n, hfds_mean_s,
        ref_lats_n, ref_lats_s,
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(hemi="n", **nc_var_name_kw),
        nf.nc_var_name(hemi="s", **nc_var_name_kw),
        ref_lat_n_bnds=ref_lats_n_bnds,
        ref_lat_s_bnds=ref_lats_s_bnds,
        nc_field_attrs_n={
            "units": nf.field_units["heatflux"],
            "long_name": nc_long_name.format("", "north"),
            **nc_var_attrs_mean_n},
        nc_field_attrs_s={
            "units": nf.field_units["heatflux"],
            "long_name": nc_long_name.format("", "south"),
            **nc_var_attrs_mean_s},
        nc_title_str=nc_title_str_mean,
        **save_nc_kw
    )
    
    diag_name_kw["space_methods"] = nf.diag_nq_area_integral
    nc_var_name_kw["space_methods"] = nf.nc_var_nq_area_integral
    
    nf.save_yearly_ref_lat(hfds_hint_n, hfds_hint_s,
        ref_lats_n, ref_lats_s,
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(hemi="n", **nc_var_name_kw),
        nf.nc_var_name(hemi="s", **nc_var_name_kw),
        ref_lat_n_bnds=ref_lats_n_bnds,
        ref_lat_s_bnds=ref_lats_s_bnds,
        nc_field_attrs_n={
            "units": nf.field_units["heattransport"],
            "long_name": nc_long_name.format(" integrated",
                "north"),
            **nc_var_attrs_int_n},
        nc_field_attrs_s={
            "units": nf.field_units["heattransport"],
            "long_name": nc_long_name.format(" integrated",
                "south", "area integral"),
            **nc_var_attrs_int_s},
        nc_title_str=nc_title_str_int,
        **save_nc_kw
    )


if __name__ == '__main__':
    main()
