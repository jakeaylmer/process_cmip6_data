import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools
)

# Basic diagnostic name (details such as time averaging, zonal
# mean, etc., and the corresponding netCDF variable name, are
# added automatically as specified in the nf module):
diag_names = ["f_down", "f_olr", "f_sw_surf", "f_sw_toa",
              "f_up"]

nc_long_names = [
    "Surface downwelling longwave radiation",
    "Top-of-atmosphere (TOA) outgoing longwave radiation",
    "Surface net downward shortwave radiation",
    "Top-of-atmosphere (TOA) net downward shortwave radiation",
    "Surface upwelling longwave radiation plus net turbulent heat flux upward"
]

nc_long_names = [
    x + ", averaged annually and everywhere {}ward of {}"
    for x in nc_long_names
]

nc_standard_names = [
    "downwelling_longwave_flux_in_air",
    "toa_outgoing_longwave_flux",
    "surface_net_downward_shortwave_flux",
    "toa_net_downward_shortwave_flux",
    None
]

nc_var_attrs_n = {
    "units"       : nf.field_units["heatflux"],
    "cell_methods": f"{nf.nc_time_name}: mean "
                    + f"{nf.nc_ref_lat_n_name}: mean "
                    + "(area-weighted)"
}

nc_var_attrs_s = nc_var_attrs_n.copy()
nc_var_attrs_s["cell_methods"] = (f"{nf.nc_time_name}: "
    + f"mean {nf.nc_ref_lat_s_name}: mean (area-weighted)")



def main():
    
    prsr = script_tools.default_cmd_argument_parser()
    prsr.add_argument("-a", "--approx", action="store_true",
        help="Do approximate version ("
            + f"{nf.diag_nq_cell_center_approx}) anyway"
    )
    cmd = prsr.parse_args()
    
    yr_s, yr_e  = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    nt = yr_e - yr_s + 1
    n_ens = len(ens_members)
    
    # For the interpolated values:
    ref_lats_n = md.default_ref_lats_n
    ref_lats_s = md.default_ref_lats_s
    
    ref_lats_n_bnds = np.array([[ref_lats_n[k], 90.0]
        for k in range(len(ref_lats_n))])
    ref_lats_s_bnds = np.array([[-90.0, ref_lats_s[k]]
        for k in range(len(ref_lats_s))])
    
    lon, lat, areacella = lpd.areacella(cmd.model)
    _, lat_bnds         = lpd.bndscella(cmd.model)
    
    load_nc_kw = {
        "model_id": cmd.model,
        "experiment_id": cmd.experiment
    }
    
    save_nc_kw = {
        "model_id"     : cmd.model,
        "member_ids"   : ens_members,
        "experiment_id": cmd.experiment,
        "year_range"   : (yr_s, yr_e)
    }
    
    # Will need to set "name", "hemi", "other_methods"
    # in loop:
    diag_kw = {
        "time_methods" : nf.diag_nq_yearly,
        "space_methods": nf.diag_nq_area_mean
    }
    nc_var_kw = {
        "time_methods" : nf.nc_var_nq_yearly,
        "space_methods": nf.nc_var_nq_area_mean
    }
    
    if np.ndim(lat) == 1:
        
        for dn, ln, sn in zip(
                diag_names, nc_long_names, nc_standard_names
            ):
            
            # Diagnostic (ldn) and netCDF variable (lvn)
            # names to load:
            ldn = nf.diag_name(dn,
                time_methods=nf.diag_nq_yearly)
            lvn = nf.nc_var_name(dn,
                time_methods=nf.nc_var_nq_yearly)
            
            data = lpd.field_2D(ldn, lvn, **load_nc_kw)[-1]
            
            true_lats_n, true_lats_n_bnds, true_data_mean_n = \
                diags.integrate_horizontally_exact(data,
                    lat_bnds, areacella, hemi="n",
                    normalise=True)
            
            true_lats_s, true_lats_s_bnds, true_data_mean_s = \
                diags.integrate_horizontally_exact(data,
                    lat_bnds, areacella, hemi="s",
                    normalise=True)
            
            interp_data_mean_n = \
                diags.interpolate_to_ref_latitudes(true_lats_n,
                    true_data_mean_n, ref_lats_n)
            
            interp_data_mean_s = \
                diags.interpolate_to_ref_latitudes(true_lats_s,
                    true_data_mean_s, ref_lats_s)
            
            # ----------------------------------------------- #
            
            print(f"Saving {dn} to NetCDF...")
            
            if sn is None:
                extra_nc_var_attrs = {}
            else:
                extra_nc_var_attrs = {"standard_name": sn}
            
            diag_kw["name"] = dn
            nc_var_kw["name"] = dn
            
            diag_kw["other_methods"] = nf.diag_nq_native
            nc_var_kw["other_methods"] = nf.nc_var_nq_native
            
            nf.save_yearly_ref_lat(true_data_mean_n,
                true_data_mean_s, true_lats_n, true_lats_s,
                nf.diag_name(**diag_kw),
                nf.nc_var_name(hemi="n", **nc_var_kw),
                nf.nc_var_name(hemi="s", **nc_var_kw),
                ref_lat_n_bnds=true_lats_n_bnds,
                ref_lat_s_bnds=true_lats_s_bnds,
                unlimited_ref_lat_dim=False,
                nc_field_attrs_n={
                    "long_name": ln.format("north",
                        "native-grid latitude bounds"),
                    **nc_var_attrs_n, **extra_nc_var_attrs},
                nc_field_attrs_s={
                    "long_name": ln.format("south",
                        "native-grid latitude bounds"),
                    **nc_var_attrs_s, **extra_nc_var_attrs},
                **save_nc_kw
            )
            
            diag_kw["other_methods"] = nf.diag_nq_native_interp
            nc_var_kw["other_methods"] = \
                nf.nc_var_nq_native_interp
            
            nf.save_yearly_ref_lat(interp_data_mean_n,
                interp_data_mean_s, ref_lats_n, ref_lats_s,
                nf.diag_name(**diag_kw),
                nf.nc_var_name(hemi="n", **nc_var_kw),
                nf.nc_var_name(hemi="s", **nc_var_kw),
                ref_lat_n_bnds=ref_lats_n_bnds,
                ref_lat_s_bnds=ref_lats_s_bnds,
                nc_field_attrs_n={
                    "long_name": ln.format("north",
                        "native-grid latitude bounds, interpol"
                        + "ated to reference latitudes"),
                    **nc_var_attrs_n, **extra_nc_var_attrs},
                nc_field_attrs_s={
                    "long_name": ln.format("south",
                        "native-grid latitude bounds, interpol"
                        + "ated to reference latitudes"),
                    **nc_var_attrs_s, **extra_nc_var_attrs},
                **save_nc_kw
            )
    
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
        
        for dn, ln, sn in zip(
                diag_names, nc_long_names, nc_standard_names
            ):
            
            # Diagnostic (ldn) and netCDF variable (lvn)
            # names to load:
            ldn = nf.diag_name(dn,
                time_methods=nf.diag_nq_yearly)
            lvn = nf.nc_var_name(dn,
                time_methods=nf.nc_var_nq_yearly)
            
            data = lpd.field_2D(ldn, lvn, **load_nc_kw)[-1]
            
            # This function returns both the integral and
            # mean, so just take second return value:
            data_mean_n = \
                diags.integrate_horizontally_multi_members(
                    data, areacella, lat, ref_lats=ref_lats_n,
                    hemi="n")[1]
            
            data_mean_s = \
                diags.integrate_horizontally_multi_members(
                    data, areacella, lat, ref_lats=ref_lats_s,
                    hemi="s")[1]
            
            # ----------------------------------------------- #
            
            print(f"Saving {dn} to NetCDF...")
            
            if sn is None:
                extra_nc_var_attrs = {}
            else:
                extra_nc_var_attrs = {"standard_name": sn}
            
            diag_kw["name"] = dn
            nc_var_kw["name"] = dn
            
            diag_kw["other_methods"] = \
                nf.diag_nq_cell_center_approx
            nc_var_kw["other_methods"] = \
                nf.nc_var_nq_cell_center_approx
            
            nf.save_yearly_ref_lat(data_mean_n,
                data_mean_s, ref_lats_n, ref_lats_s,
                nf.diag_name(**diag_kw),
                nf.nc_var_name(hemi="n", **nc_var_kw),
                nf.nc_var_name(hemi="s", **nc_var_kw),
                ref_lat_n_bnds=ref_lats_n_bnds,
                ref_lat_s_bnds=ref_lats_s_bnds,
                nc_field_attrs_n={
                    "long_name": ln.format("north",
                        "reference latitudes based on cell-"
                        "center latitudes"),
                    **nc_var_attrs_n, **extra_nc_var_attrs},
                nc_field_attrs_s={
                    "long_name": ln.format("south",
                        "reference latitudes based on cell-"
                        "center latitudes"),
                    **nc_var_attrs_s, **extra_nc_var_attrs},
                **save_nc_kw
            )



if __name__ == '__main__':
    main()
