import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_raw_data,
    metadata as md,
    netcdf as nf,
    script_tools
)

nc_var_attrs = {
    "standard_name": "cell_area",
    "long_name": "Atmosphere native-grid cell area",
    "units": nf.field_units["cellarea"]
}


def main():
    
    cmd = script_tools.default_cmd_args()
    
    cmd.experiment = md.areacella_file_kw[cmd.model][0]
    member_id      = md.areacella_file_kw[cmd.model][1]
    grid_label     = md.areacella_file_kw[cmd.model][2]
    
    lon, lat, areacella = load_raw_data.areacella(cmd.model)
    
    try:
        lon_bnds, lat_bnds  = load_raw_data.bndscella(cmd.model)
    
    except FileNotFoundError:
        
        print("WARNING: no atmosphere grid bounds detected for "
              + f"{cmd.model} -- assuming independent grid "
              + "coordinates and estimating bounds accordingly")
        
        lon_bnds, lat_bnds = diags.estimate_lon_lat_bnds_1D(
            lon, lat)
    
    # ------------------------------------------------------- #
    
    save_nc_kw = {
        "model_id": cmd.model,
        "member_id": member_id,
        "experiment_id": cmd.experiment,
        "grid_label": grid_label,
        "which": "a",
        "longitude": lon,
        "latitude": lat,
        "longitude_bnds": lon_bnds,
        "latitude_bnds": lat_bnds,
        "nc_field_attrs": nc_var_attrs,
        "nc_global_attrs": {"comment": "Climate model "
            + "diagnostics derived from the CMIP6 archive for "
            + "the analysis in the work, 'Modulation of future "
            + "sea ice loss by ocean heat transport'. This "
            + "dataset contains atmosphere grid data for one "
            + f"model ({nf.nc_file_attrs_model_name}), "
            + f"taken from one experiment "
            + f"({nf.nc_file_attrs_experiment_name}) and "
            + f"ensemble member ("
            + f"{nf.nc_file_attrs_member_name}).",
            "title": f"Atmosphere grid data for CMIP6 model: "
                + cmd.model
        }
    }
    
    nf.save_areacell(areacella, **save_nc_kw)



if __name__ == '__main__':
    main()
