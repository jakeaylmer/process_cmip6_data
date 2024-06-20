"""Calculate atmospheric heat transport (AHT) from the net heat
flux into the atmospheric column. Data for the latter should 
already be saved; need to run save_areacella.py and
save_atm_flux_yearly_fields.py first.
"""

import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools)

# Basic diagnostic name (details such as time averaging, zonal
# mean, etc., and the corresponding netCDF variable name, are
# added automatically as specified in the nf module).
# 
# Here, this one needs a bit more information about about the
# calculation method next to the diagnostic that does not need
# to be in the variable name (so define a separate one):
diag_name = "aht_from_net_flux"
var_name  = "aht"

nc_long_name = ("Northward atmospheric heat (moist-static "
    + "energy) transport (AHT) calculated by integrating AHT "
    + "convergence (AHTC) {}ward of {}, where AHTC is computed "
    + "on the native grid from the net flux into the "
    + "atmospheric column")

nc_var_attrs_n = {
    "units"        : nf.field_units["heattransport"],
    "standard_name": "northward_atmosphere_heat_transport",
    "cell_methods" :  f"{nf.nc_time_name}: mean "
                      + f"{nf.nc_ref_lat_n_name}: point"
}

nc_var_attrs_s = nc_var_attrs_n.copy()
nc_var_attrs_s["cell_methods"] = (f"{nf.nc_time_name}: "
    + f"mean {nf.nc_ref_lat_s_name}: point")

# Short description added to netCDF "title" attribute (need
# not be completely accurate/detailed here):
nc_title_str = "atmospheric heat transport"


def main():
    
    prsr = script_tools.default_cmd_argument_parser()
    prsr.add_argument("-a", "--approx", action="store_true",
        help="Do approximate version ("
            + f"{nf.diag_nq_cell_center_approx}) anyway")
    cmd = prsr.parse_args()
    
    yr_s, yr_e  = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    nt = yr_e - yr_s + 1
    n_ens = len(ens_members)
    
    # For the interpolated AHT:
    ref_lats_n = md.default_ref_lats_aht_n
    ref_lats_s = md.default_ref_lats_aht_s
    
    lon, lat, areacella = lpd.areacella(cmd.model)
    lon_bnds, lat_bnds  = lpd.bndscella(cmd.model)
    
    ahtc = lpd.field_2D("ahtc_from_net_flux_yr", "ahtc",
        cmd.model, cmd.experiment)[-1]
    
    save_nc_kw = {
        "model_id"     : cmd.model,
        "member_ids"   : ens_members,
        "experiment_id": cmd.experiment,
        "year_range"   : (yr_s, yr_e),
        "nc_title_str" : nc_title_str}
    
    # Will need to add "hemi" and "other_methods" kw to these:
    diag_kw = {"name": diag_name,
        "time_methods": nf.diag_nq_yearly}
    nc_var_kw = {"name": var_name,
        "time_methods": nf.nc_var_nq_yearly}
    
    if np.ndim(lat) == 1:
        
        true_lats_n, true_lats_n_bnds, true_aht_n = \
            diags.integrate_horizontally_exact(ahtc, lat_bnds,
                areacella, hemi="n", normalise=False)
        
        true_lats_s, true_lats_s_bnds, true_aht_s = \
            diags.integrate_horizontally_exact(ahtc, lat_bnds,
                areacella, hemi="s", normalise=False)
        
        true_aht_n /= 1.0E15  # in petawatts
        true_aht_s /= -1.0E15  # make northward too
        
        interp_aht_n = diags.interpolate_to_ref_latitudes(
            true_lats_n, true_aht_n, ref_lats_n)
        interp_aht_s = diags.interpolate_to_ref_latitudes(
            true_lats_s, true_aht_s, ref_lats_s)
        
        # --------------------------------------------------- #
        
        print("Saving to NetCDF...")
        
        diag_kw["other_methods"] = nf.diag_nq_native
        nc_var_kw["other_methods"] = nf.nc_var_nq_native
        
        nf.save_yearly_ref_lat(true_aht_n, true_aht_s,
            true_lats_n, true_lats_s,
            nf.diag_name(**diag_kw),
            nf.nc_var_name(hemi="n", **nc_var_kw),
            nf.nc_var_name(hemi="s", **nc_var_kw),
            unlimited_ref_lat_dim=False,
            nc_field_attrs_n={
                "long_name": nc_long_name.format("north",
                    "native-grid latitude bounds"),
                **nc_var_attrs_n},
            nc_field_attrs_s={
                "long_name": nc_long_name.format("south",
                    "native-grid latitude bounds"),
                **nc_var_attrs_s},
            **save_nc_kw)
        
        diag_kw["other_methods"] = nf.diag_nq_native_interp
        nc_var_kw["other_methods"] = nf.nc_var_nq_native_interp
        
        nf.save_yearly_ref_lat(interp_aht_n, interp_aht_s,
            ref_lats_n, ref_lats_s,
            nf.diag_name(**diag_kw),
            nf.nc_var_name(hemi="n", **nc_var_kw),
            nf.nc_var_name(hemi="s", **nc_var_kw),
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
            **save_nc_kw)
        
    
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
        
        aht_n = diags.integrate_horizontally_multi_members(
            ahtc, areacella, lat, ref_lats=ref_lats_n,
            hemi='n')[0] / 1.0E15
        
        aht_s = -diags.integrate_horizontally_multi_members(
            ahtc, areacella, lat, ref_lats=ref_lats_s,
            hemi='s')[0] / 1.0E15
        
        # --------------------------------------------------- #
        
        diag_kw["other_methods"] = \
            nf.diag_nq_cell_center_approx
        nc_var_kw["other_methods"] = \
            nf.nc_var_nq_cell_center_approx
        
        nf.save_yearly_ref_lat(aht_n, aht_s,
            ref_lats_n, ref_lats_s,
            nf.diag_name(**diag_kw),
            nf.nc_var_name(hemi="n", **nc_var_kw),
            nf.nc_var_name(hemi="s", **nc_var_kw),
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
            **save_nc_kw)    


if __name__ == "__main__":
    main()
