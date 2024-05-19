import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_raw_data as lrd,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools)

diag_name   = "hfds"
nc_var_name = "hfds"

nc_var_attrs = {
    "cell_measures": "area: areacella",
    "cell_methods" : f"area: mean {nf.nc_time_name}: mean",
    "long_name"    : "Net surface downward heat flux in sea "
                     + "water (hfds), annually averaged",
    "standard_name": "surface_downward_heat_flux_in_sea_water",
    "units"        : nf.field_units["heatflux"]
}

# Short description added to netCDF "title" attribute (need
# not be completely accurate/detailed here):
nc_title_str = "heat flux into ocean surface"

def main():
    
    cmd = script_tools.default_cmd_args()
    cmd.model = "GISS-E2-2-G"
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    nt = yr_e - yr_s + 1
    n_ens = len(ens_members)
    
    # Data for hfds is provided on the atmosphere grid (although
    # it is not labelled as such, nor obviously documented
    # anywhere, one has to assume this is the case given that
    # the sea ice concentration is explicitly provided on the
    # atmosphere grid and the dimensions of hfds match that of
    # the atmospheric fields):
    # 
    ny, nx = md.grid_dims_atm[cmd.model]
    lon, lat, _        = lpd.areacella(cmd.model)
    lon_bnds, lat_bnds = lpd.bndscella(cmd.model)
    
    hfds = np.zeros((nt, n_ens, ny, nx), dtype=np.float64)
    
    load_kw = {
        "model_id": cmd.model,
        "experiment_id": cmd.experiment
    }
    
    for m in range(n_ens):
        
        print(f"Processing member: {ens_members[m]} "
            + f"({m+1} of {n_ens})")
        
        hfds_data = lrd.field_2D("hfds",
            member_id=ens_members[m], **load_kw)[-1]
        
        hfds[:,m,:,:] = diags.year_mean_2D(hfds_data)
    
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
        "nc_global_attrs": {"external_variables": "areacella"},
        "nc_title_str": nc_title_str
    }
    
    diag_name_kw = {"name": diag_name,
       "time_methods": nf.diag_nq_yearly}
    
    nc_var_name_kw = {"name": nc_var_name,
        "time_methods": nf.nc_var_nq_yearly}
    
    nf.save_yearly(hfds,
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(hemi="", **nc_var_name_kw),
        nc_field_type=hfds.dtype, nc_field_attrs=nc_var_attrs,
        **save_nc_kw
    )



if __name__ == "__main__":
    main()
