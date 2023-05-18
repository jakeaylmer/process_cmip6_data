import numpy as np

from src import (
    diagnostics as diags,
    load_raw_data as lrd,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools
)

diag_name   = "o{}temptend"
nc_var_name = "o{}temptend"

nc_standard_name = ("tendency_of_sea_water_{}_temperature_"
                    + "expressed_as_heat_content")

nc_long_name = ("Ocean {} temperature expressed as heat content "
                + "(o{}temptend) integrated (summed) vertically")

nc_var_attrs = {
    "cell_measures": "area: areacello",
    "cell_methods" : f"area: mean {nf.nc_time_name}: mean",
    "units"        : nf.field_units["heatflux"]
}


def main():
    
    cmd = script_tools.default_cmd_args()
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    if cmd.model == "MRI-ESM2-0" and cmd.experiment == "ssp370":
        print("WARNING: selecting only first member in list")
        ens_members = ens_members[:1]
        n_ens = 1
    
    nt = yr_e - yr_s + 1
    n_ens = len(ens_members)
    ny, nx = md.grid_dims_ocn[cmd.model]
    
    lon, lat, _        = lpd.areacello(cmd.model)
    lon_bnds, lat_bnds = lpd.bndscello(cmd.model)
    
    ftype = md.ocn_prognostic_temperature[cmd.model]
    ftype_long = "potential" if ftype=="pot" else "conservative"
    
    load_kw = {
        "model_id": cmd.model,
        "experiment_id":cmd.experiment
    }
    
    opottemptend_zint = np.zeros((nt, n_ens, ny, nx),
                                 dtype=np.float64)
    
    for m in range(n_ens):
        
        print(f"Processing member: {ens_members[m]} "
              + f"({m+1} of {n_ens})")
        
        htend_data = lrd.field_3D(f"o{ftype}temptend",
                                  ens_members[m], **load_kw)[-1]
        
        if m == 0:
            # land mask (where missing for all time and depth):
            # (only need to calculate once)
            land_m = np.where(
                np.all(np.isnan(htend_data), axis=(0,1)),
                np.nan, 1.0
            )
        
        # Integrate vertically (sum, as dz already included
        # in the raw field o*temptend):
        opottemptend_zint[:,m,:,:] = \
            np.nansum(htend_data, axis=1)*land_m[np.newaxis,:,:]
    
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
    
    diag_name_kw = {"name": diag_name.format(ftype),
        "time_methods": nf.diag_nq_yearly,
        "space_methods": nf.diag_nq_vertical_integral}
    
    nc_var_name_kw = {"name": nc_var_name.format(ftype),
        "time_methods": nf.nc_var_nq_yearly,
        "space_methods": nf.nc_var_nq_vertical_integral}
    
    nf.save_yearly(opottemptend_zint,
        nf.diag_name(**diag_name_kw),
        nf.nc_var_name(hemi="", **nc_var_name_kw),
        nc_field_type=opottemptend_zint.dtype,
        nc_field_attrs={
            "standard_name": nc_standard_name.format(ftype_long),
            "long_name": nc_long_name.format(ftype_long, ftype),
            **nc_var_attrs},
        **save_nc_kw
    )


if __name__ == "__main__":
    main()
