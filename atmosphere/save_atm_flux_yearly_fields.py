"""Save various atmospheric vertical heat fluxes and the
atmospheric heat transport converence (AHTC) as yearly-averaged
fields, from input monthly fields. Need to run save_areacella.py
first (as the coordinates are obtained from there). 
"""

import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_processed_data as lpd,
    load_raw_data as lrd,
    metadata as md,
    netcdf as nf,
    script_tools)

# Basic diagnostic names (details such as time averaging,
# zonal mean, etc., and the corresponding netCDF variable
# name, are added automatically as specified in the nf
# module):
diag_name_ahtc     = "ahtc_from_net_flux"
diag_name_fup      = "f_up"
diag_name_fdn      = "f_down"
diag_name_folr     = "f_olr"
diag_name_fsw_toa  = "f_sw_toa"
diag_name_fsw_surf = "f_sw_surf"

# Differerent variable name for ahtc:
var_name_ahtc = "ahtc"

# Common attributes:
nc_var_attrs = {
    "units": nf.field_units["heatflux"],
    "cell_measures": "area: areacella",
    "cell_methods": f"area: mean {nf.nc_time_name}: mean"}


def main():
    
    cmd = script_tools.default_cmd_args()
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    nt = yr_e - yr_s + 1
    n_ens = len(ens_members)
    ny, nx = md.grid_dims_atm[cmd.model]
    
    lon, lat, _ = lpd.areacella(cmd.model)
    lon_bnds, lat_bnds = lpd.bndscella(cmd.model)
    
    ahtc = np.zeros((nt, n_ens, ny, nx))
    fup = np.zeros((nt, n_ens, ny, nx))
    fdn = np.zeros((nt, n_ens, ny, nx))
    folr = np.zeros((nt, n_ens, ny, nx))
    fsw_toa = np.zeros((nt, n_ens, ny, nx))
    fsw_surf = np.zeros((nt, n_ens, ny, nx))
    
    load_kw = {
        "model_id"      : cmd.model,
        "experiment_id" : cmd.experiment}
    
    # Raw data flux names:
    flux_names = ["hfls", "hfss", "rlds", "rlus", "rlut",
                  "rsds", "rsdt", "rsus", "rsut"]
    
    for m in range(n_ens):
        
        print(f"Processing member: {ens_members[m]} "
            + f"({m+1} of {n_ens})")
        
        fluxes_m = []
        
        for flux in flux_names:
        
            flux_dat = lrd.field_2D(flux,
                member_id=ens_members[m], **load_kw)[-1]
        
            fluxes_m.append(diags.year_mean_2D(flux_dat).copy())
        
        # Easier to think in terms of heat fluxes into the
        # atmospheric column (hence extra minus sign):
        ahtc[:,m,:,:] = -(
            fluxes_m[flux_names.index("hfls")]
            + fluxes_m[flux_names.index("hfss")]
            + fluxes_m[flux_names.index("rlus")]
            + fluxes_m[flux_names.index("rsus")]
            + fluxes_m[flux_names.index("rsdt")]
            - fluxes_m[flux_names.index("rlds")]
            - fluxes_m[flux_names.index("rlut")]
            - fluxes_m[flux_names.index("rsds")]
            - fluxes_m[flux_names.index("rsut")]
        )
        
        fup[:,m,:,:] = (
            fluxes_m[flux_names.index("hfls")]
            + fluxes_m[flux_names.index("hfss")]
            + fluxes_m[flux_names.index("rlus")])
        
        fdn[:,m,:,:] = fluxes_m[flux_names.index("rlds")]
        folr[:,m,:,:] = fluxes_m[flux_names.index("rlut")]
        
        fsw_toa[:,m,:,:] = (
            fluxes_m[flux_names.index("rsdt")]
            - fluxes_m[flux_names.index("rsut")])
        
        fsw_surf[:,m,:,:] = (
            fluxes_m[flux_names.index("rsds")]
            - fluxes_m[flux_names.index("rsus")])
    
    # ------------------------------------------------------- #
    
    print("Saving to NetCDF...")
    
    save_nc_kw = {
        "model_id": cmd.model,
        "member_ids": ens_members,
        "experiment_id": cmd.experiment,
        "year_range": (yr_s, yr_e),
        "longitude": lon,
        "latitude": lat,
        "longitude_bnds": lon_bnds,
        "latitude_bnds": lat_bnds,
        "nc_global_attrs": {"external_variables": "areacella"}}
    
    # Common properties for diagnostic and variable names
    # (just time here):
    diag_kw = {"time_methods": nf.diag_nq_yearly}
    nc_var_kw = {"time_methods": nf.nc_var_nq_yearly}
    
    long_name = ("Atmospheric heat (moist-static energy) "
        + "transport convergence calculated from net heat "
        + "flux into atmospheric column, annually averaged")
    
    nf.save_yearly(ahtc,
        nf.diag_name(name=diag_name_ahtc, **diag_kw),
        nf.nc_var_name(name=var_name_ahtc, **nc_var_kw),
        nc_field_type=ahtc.dtype,
        nc_field_attrs={"long_name":long_name, **nc_var_attrs},
        nc_title_str="atmospheric heat transport convergence",
        **save_nc_kw)
    
    long_name = ("Surface downwelling longwave radiation, "
        + "annually averaged")
    
    nf.save_yearly(fdn,
        nf.diag_name(name=diag_name_fdn, **diag_kw),
        nf.nc_var_name(name=diag_name_fdn, **nc_var_kw),
        nc_field_type=fdn.dtype,
        nc_field_attrs={"long_name":long_name,
            "standard_name": "downwelling_longwave_flux_in_air",
            **nc_var_attrs},
        nc_title_str="surface downwelling longwave radiation",
        **save_nc_kw)
    
    long_name = ("Top-of-atmosphere (TOA) outgoing longwave "
        + "radiation, annually averaged")
    
    nf.save_yearly(folr,
        nf.diag_name(name=diag_name_folr, **diag_kw),
        nf.nc_var_name(name=diag_name_folr, **nc_var_kw),
        nc_field_type=folr.dtype,
        nc_field_attrs={"long_name":long_name,
            "standard_name": "toa_outgoing_longwave_flux",
            **nc_var_attrs},
        nc_title_str="outgoing longwave radiation",
        **save_nc_kw)
    
    long_name = ("Surface net downward shortwave radiation, "
        + "annually averaged")
    
    nf.save_yearly(fsw_surf,
        nf.diag_name(name=diag_name_fsw_surf, **diag_kw),
        nf.nc_var_name(name=diag_name_fsw_surf, **nc_var_kw),
        nc_field_type=fsw_surf.dtype,
        nc_field_attrs={"long_name":long_name,
            "standard_name": "surface_net_downward_shortwave_flux",
            **nc_var_attrs},
        nc_title_str="surface net downward shortwave radiation",
        **save_nc_kw)
    
    long_name = ("Top-of-atmosphere (TOA) net downward shortwave "
        + "radiation, annually averaged")
    
    nf.save_yearly(fsw_toa,
        nf.diag_name(name=diag_name_fsw_toa, **diag_kw),
        nf.nc_var_name(name=diag_name_fsw_toa, **nc_var_kw),
        nc_field_type=fsw_toa.dtype,
        nc_field_attrs={"long_name":long_name,
            "standard_name": "toa_net_downward_shortwave_flux",
            **nc_var_attrs},
        nc_title_str="top of atmosphere net shortwave radiation",
        **save_nc_kw)
    
    long_name = ("Surface upwelling longwave radiation plus "
        + "net turbulent heat flux upward, annually averaged")
    
    nf.save_yearly(fup,
        nf.diag_name(name=diag_name_fup, **diag_kw),
        nf.nc_var_name(name=diag_name_fup, **nc_var_kw),
        nc_field_type=fup.dtype,
        nc_field_attrs={"long_name":long_name, **nc_var_attrs},
        nc_title_str="surface upward heat flux",
        **save_nc_kw)


if __name__ == "__main__":
    main()
