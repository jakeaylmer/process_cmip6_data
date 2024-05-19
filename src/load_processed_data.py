"""Load processed data, usually intermediate data, so that
further (final) diagnostics may be computed. Also includes
routines to prepare sea ice concentration data ready for
ice edge latitude calculation (masking lakes, ensuring units are
fraction, etc).
"""

from pathlib import Path
import numpy as np

from .load_raw_data import _load_coordinate_arrays
from . import metadata as md
from .netcdf import (nc_lon_2d_name, nc_lat_2d_name,
                     nc_ref_lat_n_name, nc_ref_lat_s_name,
                     nc_ref_lat_single_name, nc_ripf_labels,
                     nc_file_path)
from . import utils

from ice_edge_latitude.utilities import regions as iel_regions


def areacello(model_id):
    """Load and return specified model ocean grid longitude,
    latitude, and cell area arrays.
    """
    return _load_coordinate_arrays(
        Path(md.dir_out_nc_data, "areacello",
             f"areacello_{model_id}.nc"),
        ["lon", "lat", "areacello"]
    )



def areacella(model_id):
    """Load and return specified model atmosphere grid
    longitude, latitude, and cell area arrays.
    """
    return _load_coordinate_arrays(
        Path(md.dir_out_nc_data, "areacella",
             f"areacella_{model_id}.nc"),
        ["lon", "lat", "areacella"])



def bndscello(model_id):
    """Load and return specified model ocean grid cell
    longitude and latitude bounds as arrays."""
    return _load_coordinate_arrays(
        Path(md.dir_out_nc_data, "areacello",
             f"areacello_{model_id}.nc"),
        ["lon_bnds", "lat_bnds"])



def bndscella(model_id):
    """Load and return specified model atmosphere grid cell
    longitude and latitude bounds as arrays."""
    return _load_coordinate_arrays(
        Path(md.dir_out_nc_data, "areacella",
             f"areacella_{model_id}.nc"),
        ["lon_bnds", "lat_bnds"])



def time_series_single_ref_lat(diagnostic_id, nc_field_name,
                               model_id, experiment_id):
    """Load a time series data that is based on a single (non-
    north vs. south) reference latitude data. Most commonly this
    is OHT derived from hfbasin where there is no concept of
    north and south diagnostics.
    
    
    Parameters
    ----------
    diagnostic_id : str, full diagnostic id.
    nc_field_name : str, name of variable in NetCDF file.
    model_id      : str, model id.
    experiment_id : str, model id.
    
    
    Returns
    -------
    ref_lat : 1D array of float
        Reference latitude (coordinate) data.
    
    r_vals, i_vals, p_vals, f_vals : 1D arrays of int
        Ensemble member axes labels.
    
    data : 3D array of float (nt, n_ens, n_lat)
        Data as a function of time, ensemble member, and
        latitude.
    
    """
    
    print(f"Loading {model_id} {experiment_id} "
        + f"{diagnostic_id} data")
    
    data_file = nc_file_path(diagnostic_id, experiment_id,
                                model_id)
    
    nc_coords = [nc_ref_lat_single_name]
    nc_ripf = [nc_ripf_labels[x] for x in "ripf"]
    
    ref_lat, r_vals, i_vals, p_vals, f_vals, data = \
        utils.nc_get_arrays([data_file], nc_coords + nc_ripf,
                            [nc_field_name])
    
    return ref_lat, r_vals, i_vals, p_vals, f_vals, data



def time_series_ref_lats(diagnostic_id, nc_field_name_n,
                         nc_field_name_s, model_id,
                         experiment_id):
    """Load a time series data that is based on reference
    latitude data separately for the northern (n) and southern
    (s) hemispheres.
    
    
    Parameters
    ----------
    diagnostic_id   : str, full diagnostic id.
    nc_field_name_n : str, name of variable (n) in NetCDF file.
    nc_field_name_s : str, name of variable (s) in NetCDF file.
    model_id        : str, model id.
    experiment_id   : str, model id.
    
    
    Returns
    -------
    ref_lat_n, ref_lat_s : 1D array of float
        Reference latitudes (n, s; coordinate) data.
    
    r_vals, i_vals, p_vals, f_vals : 1D arrays of int
        Ensemble member axes labels.
    
    data_n, data_s : 3D array of float (nt, n_ens, n_lat)
        Data as a function of time, ensemble member, and
        reference_latitude_* for n and s, respectively.
    
    """
    
    print(f"Loading {model_id} {experiment_id} "
        + f"{diagnostic_id} data")
    
    data_file = nc_file_path(diagnostic_id, experiment_id,
                                model_id)
    
    nc_coords = [nc_ref_lat_n_name, nc_ref_lat_s_name]
    nc_ripf = [nc_ripf_labels[x] for x in "ripf"]
    
    ref_lat_n, ref_lat_s, r_vals, i_vals, p_vals, f_vals, \
        data_n, data_s = utils.nc_get_arrays([data_file],
            nc_coords + nc_ripf,
            [nc_field_name_n, nc_field_name_s])
    
    return ref_lat_n, ref_lat_s, r_vals, i_vals, p_vals, \
        f_vals, data_n, data_s



def field_2D(diagnostic_id, nc_field_name, model_id,
             experiment_id):
    """Load a 2D processed data field.
    
    
    Parameters
    ----------
    diagnostic_id : str, full diagnostic id.
    nc_field_name : str, name of variable in NetCDF file.
    model_id      : str, model id.
    experiment_id : str, model id.
    
    
    Returns
    -------
    lon, lat : 1D or 2D array of float
        Longitude and latitude data.
    
    r_vals, i_vals, p_vals, f_vals : 1D arrays of int
        Ensemble member axes labels.
    
    field : 4D array of float (nt, n_ens, ny, nx)
        Data as a function of time, ensemble member, and
        space (y, x).
    
    """
    
    print(f"Loading {model_id} {experiment_id} "
        + f"{diagnostic_id} data")
    
    data_file = nc_file_path(diagnostic_id, experiment_id,
                                model_id)
    
    # Note that nf.nc_l??_1d_name happen to be set the same as
    # the 2d versions -- so this works for atmospheric fields as
    # well as the ocean fields, for now:
    nc_coords = [nc_lon_2d_name, nc_lat_2d_name]
    nc_ripf = [nc_ripf_labels[x] for x in "ripf"]
    
    lon, lat, r_vals, i_vals, p_vals, f_vals, field = \
        utils.nc_get_arrays([data_file], nc_coords + nc_ripf,
                            [nc_field_name])
    
    return lon, lat, r_vals, i_vals, p_vals, f_vals, field



def prepare_siconc_remapped(
        variable_id         = "siconc",
        member_id           = md.default_member_id,
        model_id            = md.default_model_id,
        experiment_id       = md.default_experiment_id,
        mask_reg_names      = ["lakes", "baltic_sea", "black_sea"],
        remapped_data_dir   = "/storage/basic/cpom/gb919150/CMIP6/_swap",
        remap_method        = "bil",
        remap_res           = "0.5",
        original_miss_value = md.default_original_missing_value,
        set_miss_value      = md.default_new_missing_value
    ):
    """Load remapped sea ice concentration and make necessary
    data transformations for sea ice extent, area, and ice-edge
    latitude calculations. Specifically, ensure that values are
    in fraction, not percentage, remove unphysical values, and
    mask out major lake bodies.
    
    
    Optional parameters
    -------------------
    variable_id : str, default = "siconc"
        Can also be "siconca"
    
    mask_reg_names : list of str
        Names of regions defined in ice_edge_latitude.utilities
        regions module, to apply masks to data.
    
    
    Returns
    -------
    siconc : array
    
    """
    
    print("Preparing remapped siconc ("
        + f"{model_id}, {experiment_id}, {member_id})")
    
    data_files_in = sorted(
        [str(x) for x in Path(remapped_data_dir).glob(
                f"remapped_{remap_method}_{remap_res}deg_"
                + f"{variable_id}*{model_id}*"
                + f"{experiment_id}*{member_id}*.nc"
            )
        ]
    )
    
    #if variable_id == "siconc":
    #    coords = md.lonlat_ocn_nc_names[model_id]
    #else:  # siconca
    #    coords = md.lonlat_atm_nc_names[model_id]
    # Pretty sure cdo remap* with global uniform grid replaces
    # lon/lat names with "lon" and "lat"?
    coords = ["lon", "lat"]
    
    lon, lat, siconc = utils.nc_get_arrays(data_files_in,
                                           coords,
                                           [variable_id])
    
    if original_miss_value != set_miss_value:
        siconc = np.where(siconc >= original_miss_value,
                          set_miss_value, siconc)
    
    # Need 2D coordinates for lake masking:
    if np.ndim(lon) == 1 and np.ndim(lat) == 1:
        lon_m, lat_m = np.meshgrid(lon, lat)
    else:
        lon_m = lon.copy()
        lat_m = lat.copy()
    
    lon_m = lon_m % 360.0  # need 0-360 range for lake masking
    
    if not (siconc < 10.0).all():
        # fair bet we are a percentage, not fraction [check is
        # with 10 rather than 1 in case there are "overshoot"
        # values like 1.02 (this can occur from interpolation)
        # even if the units truly are fraction]:
        siconc /= 100.0
    
    # Restrict range to [0.0, 1.0]:
    siconc = np.maximum(0.0, siconc)
    siconc = np.minimum(1.0, siconc)
    
    # Mask lakes
    for reg in mask_reg_names:
        
        reg_def = getattr(iel_regions, reg)
        
        for sub_reg in reg_def.keys():
            sub_reg_mask = np.where(
                (lon_m >= reg_def[sub_reg][0]) &
                (lon_m <= reg_def[sub_reg][1]) &
                (lat_m >= reg_def[sub_reg][2]) &
                (lat_m <= reg_def[sub_reg][3]),
                np.nan, 1.0)
            
            siconc *= sub_reg_mask[np.newaxis,:,:]
    
    return lon, lat, siconc
