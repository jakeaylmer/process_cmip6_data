import numpy as np

from src import (
    load_raw_data as lrd,
    metadata as md,
    netcdf as nf,
    script_tools
)

nc_var_attrs = {
    "standard_name": "cell_area",
    "long_name": "Ocean native-grid cell area",
    "units": nf.field_units["cellarea"]
}


def main():
    
    cmd = script_tools.default_cmd_args()
    
    cmd.experiment = md.areacello_file_kw[cmd.model][0]
    member_id      = md.areacello_file_kw[cmd.model][1]
    grid_label     = md.areacello_file_kw[cmd.model][2]
    
    lon, lat, areacello = lrd.areacello(cmd.model)
    lon_bnds, lat_bnds  = lrd.bndscello(cmd.model)
    
    areacello = np.where(
        areacello >= md.default_original_missing_value,
        md.default_new_missing_value, areacello)
    
    # ------------------------------------------------------- #
    
    save_nc_kw = {
        "model_id": cmd.model,
        "member_id": member_id,
        "experiment_id": cmd.experiment,
        "grid_label": grid_label,
        "which": "o",
        "longitude": lon,
        "latitude": lat,
        "longitude_bnds": lon_bnds,
        "latitude_bnds": lat_bnds,
        "nc_field_attrs": nc_var_attrs,
        "nc_global_attrs": {"comment": "Climate model "
            + "diagnostics derived from the CMIP6 archive for "
            + "the analysis in the work, 'Modulation of future "
            + "sea ice loss by ocean heat transport'. This "
            + "dataset contains ocean grid data for one "
            + f"model ({nf.nc_file_attrs_model_name}), "
            + f"taken from one experiment "
            + f"({nf.nc_file_attrs_experiment_name}) and "
            + f"ensemble member ("
            + f"{nf.nc_file_attrs_member_name}).",
            "title": f"Ocean grid data for CMIP6 model: "
                + cmd.model
        }
    }
    
    nf.save_areacell(areacello, **save_nc_kw)



if __name__ == '__main__':
    main()
