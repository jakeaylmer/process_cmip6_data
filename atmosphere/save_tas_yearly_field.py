import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_raw_data as lrd,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools
)

# Basic diagnostic name (details such as time averaging,
# zonal mean, etc., and the corresponding nf variable
# name, are added automatically as specified in the nf
# module):
diag_name = "tas"

nc_var_attrs = {
    "cell_measures": "area: areacella",
    "cell_methods" : f"area: mean {nf.nc_time_name}: mean",
    "long_name"    : "Near-surface air temperature, annually "
                     + "averaged",
    "standard_name": "air_temperature",
    "units"        : nf.field_units["temperature"]
}


def main():
    
    cmd = script_tools.default_cmd_args()
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    nt = yr_e - yr_s + 1
    n_ens = len(ens_members)
    ny, nx = md.grid_dims_atm[cmd.model]
    
    lon, lat, _        = lpd.areacella(cmd.model)
    lon_bnds, lat_bnds = lpd.bndscella(cmd.model)
    
    tas = np.zeros((nt, n_ens, ny, nx)).astype(np.float64)
    
    load_kw = {
        "model_id"     : cmd.model,
        "experiment_id": cmd.experiment
    }
    
    for m in range(n_ens):
        
        print(f"Processing member: {ens_members[m]} "
            + f"({m+1} of {n_ens})")
        
        tas_data = lrd.field_2D("tas",
            member_id=ens_members[m], **load_kw)[-1]
    
        tas[:,m,:,:] = diags.year_mean_2D(tas_data)
    
    if nc_var_attrs["units"] == "K":
        if not all(tas[~np.isnan(tas)] > 100.0):
            # Safe bet (!) that we are in degrees_celsius
            tas = md.degree_C_to_K(tas)
    else:  # want degrees_C
        if all(tas[~np.isnan(tas)] > 100.0):
            # Safe bet (!) that we are in K
            tas = md.K_to_degree_C(tas)
            tas = md.K_to_degree_C(tas)
    
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
        "nc_global_attrs": {"external_variables": "areacella"}
    }
    
    diag_kw = {
        "name"        : diag_name,
        "time_methods": nf.diag_nq_yearly
    }
    
    nc_var_kw = {
        "name"        : diag_name,
        "time_methods": nf.nc_var_nq_yearly
    }
    
    nf.save_yearly(tas,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(**nc_var_kw),
        nc_field_type=tas.dtype, nc_field_attrs=nc_var_attrs,
        **save_nc_kw
    )



if __name__ == '__main__':
    main()
