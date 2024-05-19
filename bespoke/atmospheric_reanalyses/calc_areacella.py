from argparse import ArgumentParser
from pathlib import Path
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
    
    prsr = ArgumentParser()
    
    prsr.add_argument("-r", "--reanalysis", type=str,
                      default=md.default_reanalysis,
                      choices=md.defined_reanalyses,
                      help="Reanalysis name (case sensitive)")
    
    cmd = prsr.parse_args()
    
    year_sel = 2021 if cmd.reanalysis == "CFSv2" else 1980
    
    lon, lat = load_raw_data._load_coordinate_arrays(
        Path(md.dir_raw_nc_data_reanalyses, cmd.reanalysis,
             md.reanalysis_nc_file_fmt["tas"].format(year_sel)),
        md.reanalysis_nc_coord_names[cmd.reanalysis])
    
    # Latitude should be increasing:
    if lat[1] < lat[0]:
        lat = lat[::-1]
    
    print(f"Estimating areacella for {cmd.reanalysis}")
    lon_bnds, lat_bnds = diags.estimate_lon_lat_bnds_1D(lon, lat)
    areacella = diags.grid_cell_area(lon_bnds, lat_bnds)
    
    # Sanity check:
    total_area = np.sum(areacella)*1.0E-6  # km2
    true_area = 4.0*np.pi*(6.371E3)**2  # km2
    perc_error = 100.0*(total_area - true_area) / true_area
    print(f"sum(areacella): {total_area:.2f} km2")
    print(f"Should be     : {true_area:.2f} km2")
    print(f"Perc. error   : {'+' if perc_error > 0 else '-'}"
          + f"{abs(perc_error):.6f} %")
    
    # ------------------------------------------------------- #
    
    save_nc_kw = {
        "model_id": cmd.reanalysis,
        "member_id": "N/A",
        "experiment_id": "reanalysis",
        "grid_label": "gr",
        "which": "a",
        "longitude": lon,
        "latitude": lat,
        "longitude_bnds": lon_bnds,
        "latitude_bnds": lat_bnds,
        "nc_field_attrs": nc_var_attrs,
        "nc_global_attrs": {
            "comment": "Atmospheric reanalysis diagnostics "
                       + "for the analysis presented in Aylmer "
                       + "et al. 2024 [1]. This dataset "
                       + "contains atmosphere grid data for "
                       + "one reanalysis "
                       + f"({nf.nc_file_attrs_model_name}) "
                       + "[2,3].",
            "title": "Atmosphere grid data for reanalysis: "
                     + cmd.reanalysis
        },
        "save_dir" : Path(md.dir_out_nc_data_reanalyses,
                          "areacella"),
        "file_name": f"areacella_{cmd.reanalysis}.nc"
    }
    
    nf.save_areacell(areacella, **save_nc_kw)



if __name__ == "__main__":
    main()
