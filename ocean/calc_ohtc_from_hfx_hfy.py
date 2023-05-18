import numpy as np

from src import (
    diagnostics as diags,
    load_raw_data as lrd,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools
)

diag_name   = "ohtc_from_hfx_hfy"
nc_var_name = "ohtc"

nc_var_attrs = {
    "cell_measures": "area: areacello",
    "cell_methods" : f"area: mean {nf.nc_time_name}: mean",
    "long_name"    : "Ocean heat transport convergence (OHTC) "
                     + "computed on the native grid from hfx "
                     + "and hfy",
    "units"        : nf.field_units["heatflux"]
}


def calc_ohtc_C_grid(hfx, hfy, areacello):
    """"""
    
    # Set missing to 0 for calculation:
    hfx = np.where(np.isnan(hfx), 0.0, hfx)
    hfy = np.where(np.isnan(hfy), 0.0, hfy)
    
    # Compute OHT divergence (to save inverting
    # all signs, which can be done at the end):
    ohtd = np.zeros(np.shape(hfx)).astype(np.float64)
    
    nbl = 1  # number of x-border cells, left hand side
    nbr = 1  # number of x-border cells, right hand side
    nbt = 2  # number of y-border cells, top
    
    # Set data for all rows (j) and columns (i) except first
    # and border cells):
    ohtd[:,1:-nbt,1+nbl:-nbr] = \
        hfx[:,1:-nbt,1+nbl:-nbr] - hfx[:,1:-nbt,nbl:-nbr-1] \
        + hfy[:,1:-nbt,1+nbl:-nbr] - hfy[:,:-1-nbt,1+nbl:-nbr]
    
    # Wrap around in x-direction, accounting for border cells:
    ohtd[:,1:-nbt,nbl] = \
        hfx[:,1:-nbt,nbl] - hfx[:,1:-nbt,-nbr-1] \
        + hfy[:,1:-nbt,nbl] - hfy[:,:-1-nbt,nbl]
    
    # Wrap around at north pole; grid cells are in two halfs
    # horizontally (NEMO grid only):
    
    
    # Determine grid-cell mean ohtd
    # (follows from Gauss' theorem):
    ohtd /= areacello[np.newaxis,:,:]
    
    return -ohtd


def process_member(member_id, lon, lat, areacello,
        model_id='CanESM5', experiment_id='historical'
    ):
    """"""
    
    kw = {
        "member_id": member_id,
        "model_id": model_id,
        "experiment_id": experiment_id
    }
    
    _, _, hfx, = lrd.field_2D("hfx", **kw)
    _, _, hfy, = lrd.field_2D("hfy", **kw)
    
    land = np.where(np.isnan(hfx[0,:,:]), np.nan, 1.0)
    
    ohtc = calc_ohtc_C_grid(hfx, hfy, areacello)
    
    ohtc_ym = diags.year_mean_2D(ohtc)*land[np.newaxis,:,:]
    
    return ohtc_ym


def main():
    
    cmd = script_tools.default_cmd_args()
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    nt = yr_e - yr_s + 1
    n_ens = len(ens_members)
    ny, nx = md.grid_dims_ocn[cmd.model]
    
    lon, lat, areacello = lpd.areacello(cmd.model)
    lon_bnds, lat_bnds  = lpd.bndscello(cmd.model)
    
    pm_kw = {
        'lon': lon,
        'lat': lat,
        'areacello': areacello,
        'model_id': cmd.model,
        'experiment_id': cmd.experiment
    }
    
    ohtc = np.zeros((nt, n_ens, ny, nx), dtype=np.float64)
    
    for m in range(n_ens):
        print(f"Processing member: {ens_members[m]} "
            + f"({m+1} of {n_ens})")
        ohtc[:,m,:,:] = process_member(ens_members[m], **pm_kw)
    
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
        "nc_global_attrs": {"external_variables": "areacello"}
    }
    
    diag_name_kw = {"name": diag_name,
        "time_methods": nf.diag_nq_yearly}
    
    nc_var_name_kw = {"name": nc_var_name,
        "time_methods": nf.nc_var_nq_yearly}
    
    nf.save_yearly(ohtc, nf.diag_name(**diag_name_kw),
        nf.nc_var_name(hemi="", **nc_var_name_kw),
        nc_field_type=ohtc.dtype, nc_field_attrs=nc_var_attrs,
        **save_nc_kw
    )



if __name__ == "__main__":
    main()
