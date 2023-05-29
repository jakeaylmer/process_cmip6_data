from process_cmip6_data.src import (
    diagnostics as diags,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools
)

# Diagnostic names:
diag_name = "oht_from_hfx_hfy"
nc_var_name = "oht"

nc_long_name = ("Northward ocean heat transport "
    + "(OHT) calculated by integrating OHT convergence "
    + "(OHTC) {}ward of reference latitudes based on cell-"
    + "center latitudes, where OHTC is computed on "
    + "the native grid from hfx and hfy")

nc_var_attrs_n = dict()
nc_var_attrs_n["standard_name"] = \
    "northward_ocean_heat_transport"
nc_var_attrs_n["units"] = nf.field_units["heattransport"]
nc_var_attrs_n["cell_methods"] = (
    f"{nf.nc_time_name}: mean " +
    f"{nf.nc_ref_lat_n_name}: point")

nc_var_attrs_s = nc_var_attrs_n.copy()
nc_var_attrs_s["cell_methods"] = (
    f"{nf.nc_time_name}: mean " + 
    f"{nf.nc_ref_lat_s_name}: point")



def main():
    
    cmd = script_tools.default_cmd_args()
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    ref_lats_n = md.default_ref_lats_oht_n
    ref_lats_s = md.default_ref_lats_oht_s
    
    _, lat, areacello = lpd.areacello(cmd.model)
    
    diag_name_kw = {"name": "ohtc_from_hfx_hfy",
        "time_methods": nf.diag_nq_yearly}
    
    nc_var_name_kw = {"name": "ohtc",
        "time_methods": nf.nc_var_nq_yearly}
    
    ohtc = lpd.field_2D(
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(**nc_var_name_kw),
        cmd.model, cmd.experiment)[-1]
    
    oht_n = diags.integrate_horizontally_multi_members(
        ohtc, areacello, lat, ref_lats=ref_lats_n,
        hemi="n")[0] / 1.0E15
    
    oht_s = -diags.integrate_horizontally_multi_members(
        ohtc, areacello, lat, ref_lats=ref_lats_s,
        hemi="s")[0] / 1.0E15
    
    # ------------------------------------------------------- #
    
    print("Saving to NetCDF...")
    
    diag_name_kw["name"] = diag_name
    diag_name_kw["other_methods"] = \
        nf.diag_nq_cell_center_approx
    
    nc_var_name_kw["name"] = nc_var_name
    nc_var_name_kw["other_methods"] = \
        nf.nc_var_nq_cell_center_approx
    
    nf.save_yearly_ref_lat(oht_n, oht_s,
        ref_lats_n, ref_lats_s,
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(hemi="n", **nc_var_name_kw),
        nf.nc_var_name(hemi="s", **nc_var_name_kw),
        cmd.model, ens_members,
        experiment_id=cmd.experiment,
        year_range=(yr_s, yr_e),
        nc_field_attrs_n={
            "long_name": nc_long_name.format("north"),
            **nc_var_attrs_n},
        nc_field_attrs_s={
            "long_name": nc_long_name.format("south"),
            **nc_var_attrs_s},
    )



if __name__ == '__main__':
    main()
