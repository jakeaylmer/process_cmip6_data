from argparse import ArgumentParser
from pathlib import Path
import numpy as np
import netCDF4 as nc

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
nc_var_name = "tas"

nc_var_attrs = {
    "cell_measures": "area: areacella",
    "cell_methods" : f"area: mean {nf.nc_time_name}: mean",
    "long_name"    : "Near-surface air temperature",
    "standard_name": "air_temperature",
    "units"        : nf.field_units["temperature"]
}

nc_title_str = "near-surface air temperature"


def save_monthly_one_member(field, diagnostic_id, nc_field_name,
        model_id, member_ids, experiment_id,
        year_range=(1950, 2014),
        longitude=None, latitude=None,
        longitude_bnds=None, latitude_bnds=None,
        nc_field_type=np.float64,
        nc_field_attrs={},
        nc_global_attrs={},
        nc_title_str="",
        save_dir=None, file_name=None
    ):
    """Save a 1, 2, or 3D field of yearly-mean data.
    
    Parameters
    ----------
    field : numpy array (2D, 3D, or 4D)
        Field data to save:
            4D : (nt, n_ens, ny, nx)
            3D : (nt, n_ens, ny)
            2D : (nt, n_ens)
    
    diagnostic_id : str (for automated directory structure)
    nc_field_name : str (netCDF variable name)
    model_id      : str (for metadata)
    experiment_id : str (for metadata)
    
    
    Optional parameters
    -------------------
    year_range : iterable of int (year_start, year_end)
        Start and end year (defines the length of the time
        axis, and the time coordinates of field along axis=0).
        End year is included.
    
    longitude, latitude: 2D (ny, nx) or 1D (ny,), (nx,) arrays
        Longitude and latitude coordinates in degrees_east and
        degrees_north, respectively. Can also be None if not
        required. Both must be provided if field is 4D;
        latitude must be provided if field is 3D.
    
    longitude_bnds, latitude_bnds:
        3D (ny, nx, 4) or 2D (ny, 2), (nx, 2) arrays
        Longitude and latitude cell bounds. Can also be None if
        not required.
    
    nc_field_type : default = np.float64
        Datatype for field.
    
    nc_field_attrs : dict of str, default = {}
        The keys are attributes set for the netCDF field.
    
    nc_global_attrs : dict of str, default = {}
        Global netCDF attributes set in addition to (overriding
        if applicable) the defaults (title, author, etc.).
    
    save_dir : str or pathlib.Path or None
        Absolute save directory path. If None, the default is
        used (set in metadata and netcdf modules).
    
    file_name : str or None
        File name to save to. If None, the default is used (set
        in metadata and netcdf modules). If the file already
        exists in save_dir, note that it is overwritten without
        warning.
    
    """
    
    save_dir, file_name = \
        nf.prepare_to_save(save_dir, file_name, diagnostic_id,
                           experiment_id, model_id)
    
    nc_file_attrs = nf.prepare_nc_global_attrs(model_id,
        experiment_id, nc_title_str, nc_global_attrs)
    
    with nc.Dataset(Path(save_dir, file_name), "w") as ncout:
        
        nf._set_attributes(ncout, nc_file_attrs)
        
        nf.create_time_variables(ncout,
                                nf.nc_time_units[experiment_id],
                                nf.nc_calendar[experiment_id])
        nf.set_time_data_monthly(ncout, year_range)
        nf.set_ripf_data(ncout, member_ids)
        
        # Assume field is (time, ens, y, x) for one member
        if np.ndim(longitude) == 2:
            xdim_name = nf.nc_x_dim_name
            ydim_name = nf.nc_y_dim_name
            nvert = 4
            lon_name = nf.nc_lon_2d_name
            lat_name = nf.nc_lat_2d_name
        else:
            xdim_name = nf.nc_lon_dim_name
            ydim_name = nf.nc_lat_dim_name
            nvert = 2
            lon_name = nf.nc_lon_1d_name
            lat_name = nf.nc_lat_1d_name
        
        ncout.createDimension(ydim_name, np.shape(field)[2])
        ncout.createDimension(xdim_name, np.shape(field)[3])
        ncout.createDimension(nf.nc_vert_dim_name, nvert)
        
        if nvert == 4:
            nf.set_lon_2d_data(ncout, longitude, longitude_bnds)
            nf.set_lat_2d_data(ncout, latitude, latitude_bnds)
        else:
            nf.set_lon_1d_data(ncout, longitude,longitude_bnds)
            nf.set_lat_1d_data(ncout, latitude, latitude_bnds)
        
        ncfield = ncout.createVariable(nc_field_name,
            nc_field_type, (nf.nc_time_dim_name,
                nf.nc_member_dim_name, ydim_name, xdim_name),
            fill_value=nf.nc_fill_value, zlib=True
        )
        ncfield[:,:,:,:] = field
        
        if (nvert == 4 and
                "coordinates" not in nc_field_attrs.keys()):
            nc_field_attrs["coordinates"] = \
                f"{lon_name} {lat_name}"
    
        nf._set_attributes(ncfield, nc_field_attrs)
    
    print(f"Saved: {str(Path(save_dir, file_name))}")



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
    n_ens = len(ens_members)
    
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
    
    # Re-shape into 4 dimensions, adding the ensemble member
    # dimension of size 1, for consistency with CMIP data and
    # other routines here [so that array now has shape
    # (nt, n_ens, ny, nx)]:
    tas = np.expand_dims(tas, axis=1)
    
    # ------------------------------------------------------- #
    
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
        "time_methods": nf.diag_nq_monthly
    }
    
    nc_var_kw = {
        "name"        : diag_name,
        "time_methods": nf.nc_var_nq_monthly
    }
    
    save_monthly_one_member(tas,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(hemi="", **nc_var_kw),
        nc_field_type=tas.dtype,
        nc_field_attrs=nc_var_attrs,
        save_dir=Path(md.dir_out_nc_data_reanalyses,
                      nf.diag_name(**diag_kw)),
        **save_nc_kw)



if __name__ == '__main__':
    main()
