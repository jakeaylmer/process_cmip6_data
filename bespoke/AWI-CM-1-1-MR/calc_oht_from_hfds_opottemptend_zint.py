import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_raw_data as lrd,
    metadata as md,
    netcdf as nf,
    script_tools
)

# Diagnostic names:
diag_name = "oht_from_hfds_opottemptend_" \
    + nf.diag_nq_vertical_integral
nc_var_name = "oht"

nc_long_name = ("Northward ocean heat transport "
    + "(OHT) calculated by integrating OHT convergence "
    + "(OHTC) {}ward of reference latitudes based on cell-"
    + "center latitudes, where OHTC is computed on "
    + "the native grid from the residual of net flux into "
    + "the ocean column (hfds) and depth-integrated ocean "
    + "heat content tendency (opottemptend)")

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

# Short description added to netCDF "title" attribute (need
# not be completely accurate/detailed here):
nc_title_str = "ocean heat transport"

def main():
    
    cmd = script_tools.default_cmd_args()
    cmd.model = "AWI-CM-1-1-MR"
    
    yr_s, yr_e  = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    nt = yr_e - yr_s + 1
    n_ens = len(ens_members)
    
    ref_lats_n = md.default_ref_lats_oht_n
    ref_lats_s = md.default_ref_lats_oht_s
    
    _, lat, areacello = lrd.areacello(cmd.model)
    n_cells = len(areacello)
    
    ohtc = np.zeros((nt, n_ens, n_cells), dtype=np.float64)
    
    for m in range(n_ens):
        
        print(f"Processing member {ens_members[m]} "
            + f"({m+1} of {n_ens})")
        
        # Download monthly data, compute yearly averages
        # separately from this script before loading here [also
        # 1D (+ time -> 2D)]:
        hfds_m = lrd.field_2D("hfds",
            member_id=ens_members[m], model_id=cmd.model,
            experiment_id=cmd.experiment)[-1]
        
        # Download yearly data, compute depth integrals (sum),
        # separately from this script, before loading here
        # [these fields should also be 1D (+ time -> 2D), after
        # vertical sum]:
        dhdt_m = lrd.field_2D("opottemptend",
            model_id=cmd.model, member_id=ens_members[m],
            experiment_id=cmd.experiment)[-1]
        
        # In case missing values (NaN, land) do not line up for
        # whatever reason (important for closing budget; here we
        # do not need to worry about land masks being removed as
        # we only integrate -- no normalisation versions of the
        # diagnostics):
        # 
        # (I don't think there are any missing values in the
        # data for AWI-CM-1-1-MR anyway, but this doesn't hurt)
        # 
        hfds_m = np.where(np.isnan(hfds_m), 0.0, hfds_m)
        dhdt_m = np.where(np.isnan(dhdt_m), 0.0, dhdt_m)
        
        ohtc[:,m,:] = dhdt_m - hfds_m
    
    # Add dimension of size 1 so the diagnostic routines still
    # work:
    lat = np.expand_dims(lat, axis=1)
    areacello = np.expand_dims(areacello, axis=1)
    ohtc = np.expand_dims(ohtc, axis=3)
    
    oht_n = diags.integrate_horizontally_multi_members(
        ohtc, areacello, lat, ref_lats=ref_lats_n,
        hemi='n')[0] / 1.0E15
    
    oht_s = -diags.integrate_horizontally_multi_members(
        ohtc, areacello, lat, ref_lats=ref_lats_s,
        hemi='s')[0] / 1.0E15
    
    # ------------------------------------------------------- #
    
    print("Saving to NetCDF...")
    
    save_nc_kw = {
        "model_id"     : cmd.model,
        "member_ids"   : ens_members,
        "experiment_id": cmd.experiment,
        "year_range"   : (yr_s, yr_e),
        "nc_title_str" : nc_title_str
    }
    
    diag_name_kw = {
        "name": diag_name,
        "time_methods": nf.diag_nq_yearly,
        "space_methods": "",
        "other_methods": nf.diag_nq_cell_center_approx}
    
    nc_var_name_kw = {
        "name": nc_var_name,
        "time_methods": nf.nc_var_nq_yearly,
        "space_methods": "",
        "other_methods": nf.nc_var_nq_cell_center_approx}
    
    nf.save_yearly_ref_lat(oht_n, oht_s,
        ref_lats_n, ref_lats_s,
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(hemi="n", **nc_var_name_kw),
        nf.nc_var_name(hemi="s", **nc_var_name_kw),
        nc_field_attrs_n={
            "long_name": nc_long_name.format("north"),
            **nc_var_attrs_n},
        nc_field_attrs_s={
            "long_name": nc_long_name.format("south"),
            **nc_var_attrs_s},
        **save_nc_kw
    )



if __name__ == '__main__':
    main()
