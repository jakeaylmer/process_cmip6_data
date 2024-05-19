from argparse import ArgumentParser
from pathlib import Path
import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_raw_data as lrd,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools,
    utils
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

nc_title_str = "near-surface air temperature"


def main():
    
    prsr = ArgumentParser()
    
    prsr.add_argument("-r", "--reanalysis", type=str,
                      default=md.default_reanalysis,
                      choices=md.defined_reanalyses,
                      help="Reanalysis name (case sensitive)")
    
    cmd = prsr.parse_args()
    
    yr_s, yr_e = md.reanalysis_year_range[cmd.reanalysis]
    ny, nx = md.reanalysis_grid_dims_atm[cmd.reanalysis]
    ens_members = ["r1i1p1f1"]
    
    nt = yr_e - yr_s + 1
    
    # Load coordinates from the processed areacella data for
    # reanalyses, using the raw-data loading routines (which are
    # more generic than the processed-data loading routines and
    # can be adapted here for the reanalysis data):
    lon, lat, lon_bnds, lat_bnds = lrd._load_coordinate_arrays(
        Path(md.dir_out_nc_data_reanalyses, "areacella",
             f"areacella_{cmd.reanalysis}.nc"),
        [nf.nc_lon_1d_name, nf.nc_lat_1d_name,
         f"{nf.nc_lon_1d_name}_bnds",
         f"{nf.nc_lat_1d_name}_bnds"])
    
    # Load raw tas data (manually, as there is no generic
    # routine in src).
    # 
    # (1) Get netCDF file list (assumes common format as
    #     specified here and in metadata.py):
    nc_files = [
        Path(md.dir_raw_nc_data_reanalyses, cmd.reanalysis,
             md.reanalysis_nc_file_fmt["tas"].format(y))
        for y in range(yr_s, yr_e+1, 1)
    ]
    # (2) Load data. Need to also load latitude to check if the
    #     corresponding dimension needs to be reversed. This
    #     will already have been checked and corrected in the
    #     version of latitude loaded for areacella.
    print("Loading raw data...")
    lat_raw, tas = utils.nc_get_arrays(nc_files,
        [md.reanalysis_nc_coord_names[cmd.reanalysis][1]],
        [md.reanalysis_nc_names["tas"][cmd.reanalysis]])
    
    # Latitude must be increasing:
    if lat_raw[1] < lat_raw[0]:
        tas = tas[:,::-1,:]
        # Again don't actually need to correct lat or lat_raw,
        # since the former will have been checked and corrected
        # during computation of areacella, the nc file of which
        # it has been loaded from, and the latter here is only
        # needed for checking whether it *had* been done, so
        # that the latitude dimension of tas can also be
        # reversed as above.
    
    if nc_var_attrs["units"] == "K":
        if not all(tas[~np.isnan(tas)] > 100.0):
            # Safe bet (!) that we are in degrees_celsius
            tas = md.degree_C_to_K(tas)
    else:  # want degrees_C
        if all(tas[~np.isnan(tas)] > 100.0):
            # Safe bet (!) that we are in K
            tas = md.K_to_degree_C(tas)
            tas = md.K_to_degree_C(tas)
    
    # Data is monthly averaged; convert to annual mean:
    tas = diags.year_mean_2D(tas)
    
    # Re-shape into 4 dimensions, adding the ensemble member
    # dimension of size 1, for consistency with CMIP data and
    # other routines here [so that array now has shape
    # (nt, n_ens, ny, nx)]:
    tas = np.expand_dims(tas, axis=1)
    
    # ------------------------------------------------------- #
    
    print("Saving to NetCDF...")
    
    save_nc_kw = {
        "model_id": cmd.reanalysis,
        "member_ids": ens_members,
        "experiment_id": "reanalysis",
        "year_range": (yr_s, yr_e),
        "longitude": lon,
        "latitude": lat,
        "longitude_bnds": lon_bnds,
        "latitude_bnds": lat_bnds,
        "nc_global_attrs": {
            "comment": "Atmospheric reanalysis diagnostics "
                       + "for the analysis presented in Aylmer "
                       + "et al. 2024 [1]. This dataset "
                       + "contains one diagnostic for "
                       + "one reanalysis "
                       + f"({nf.nc_file_attrs_model_name}) "
                       + "[2,3].",
            "external_variables": "areacella",
            "title": "Atmospheric reanalysis diagnostics: "
                     + f"{cmd.reanalysis}: {nc_title_str}"
        }
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
        nc_field_type=tas.dtype,
        nc_field_attrs=nc_var_attrs,
        save_dir=Path(md.dir_out_nc_data_reanalyses,
                      nf.diag_name(**diag_kw)),
        **save_nc_kw)



if __name__ == "__main__":
    main()
