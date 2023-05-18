"""Load raw data. Also includes routines to prepare sea ice
concentration data ready for sea ice extent calculation (masking
lakes, ensuring units are fraction, etc).
"""

from pathlib import Path
import numpy as np

from src import metadata as md

from my_python_utilities.data_tools import nc_tools as nct
from ice_edge_latitude.utilities import regions as iel_regions


def _load_coordinate_arrays(nc_file, nc_coord_names):
    """Generic function to load coordinate arrays from
    one netCDF file using nct.get_arrays() routine."""
    return nct.get_arrays([nc_file], nc_coord_names, [])



def areacello(model_id):
    """Load and return specified model ocean grid longitude,
    latitude, and cell area arrays.
    """
    return _load_coordinate_arrays(
        md.areacello_file(model_id),
        list(md.lonlat_ocn_nc_names[model_id]) + ["areacello"])



def areacella(model_id):
    """Load and return specified model atmosphere grid
    longitude, latitude, and cell area arrays.
    """
    return _load_coordinate_arrays(
        md.areacella_file(model_id),
        list(md.lonlat_atm_nc_names[model_id]) + ["areacella"])



def bndscello(model_id):
    """Load and return specified model ocean grid cell
    longitude and latitude bounds as arrays."""
    return _load_coordinate_arrays(
        md.lonlat_bnds_ocn_file(model_id),
        md.lonlat_bnds_ocn_nc_names[model_id])



def bndscella(model_id):
    """Load and return specified model atmosphere grid cell
    longitude and latitude bounds as arrays."""
    return _load_coordinate_arrays(
        md.lonlat_bnds_atm_file(model_id),
        md.lonlat_bnds_atm_nc_names[model_id])



def field_2D(variable_id,
        member_id           = md.default_member_id,
        model_id            = md.default_model_id,
        experiment_id       = md.default_experiment_id,
        original_miss_value = md.default_original_missing_value,
        set_miss_value      = md.default_new_missing_value
    ):
    """Load a 2D field (+ time -> 3D) as an array. Returns the
    longitude, latitude, and field arrays. Resets missing values
    as specified.
    
    Assumes raw netCDF files from CMIP archive have not been
    renamed, and are stored under the root directory specified
    in metadata, then subdirectories corresponding to model_id,
    experiment_id, then variable_id:
    
    ./MODEL_NAME/piControl/siconc/*.nc
    
    This is critical when there are multiple files for the same
    field, as the sorted-order is assumed to be datetime order
    (which works if the files are not renamed).
    
    """
    data_dir = Path(md.dir_raw_nc_data, model_id,
                    experiment_id, variable_id)
    
    data_files_in = sorted(
        [str(x) for x in Path(data_dir).glob(
            f"{variable_id}*{model_id}*"
            + f"{experiment_id}*{member_id}*.nc")])
    
    if md.variable_domain[variable_id] == "ocn":
        coords = md.lonlat_ocn_nc_names[model_id]
    else:
        coords = md.lonlat_atm_nc_names[model_id]
    
    lon, lat, fld = nct.get_arrays(data_files_in, coords,
                                   [variable_id])
    
    if original_miss_value != set_miss_value:
        fld = np.where(fld >= original_miss_value,
                       set_miss_value, fld)
    
    return lon, lat, fld



def field_3D(variable_id,
        member_id           = md.default_member_id,
        model_id            = md.default_model_id,
        experiment_id       = md.default_experiment_id,
        original_miss_value = md.default_original_missing_value,
        set_miss_value      = md.default_new_missing_value
    ):
    """Load a 3D field (+ time -> 4D) as an array. Returns the
    longitude, latitude, and field arrays. Depth coordinate is
    *not* returned. Resets missing values as specified.
    
    Assumes raw netCDF files from CMIP archive have not been
    renamed, and are stored under the root directory specified
    in metadata, then subdirectories corresponding to model_id,
    experiment_id, then variable_id:
    
    ./MODEL_NAME/piControl/siconc/*.nc
    
    This is critical when there are multiple files for the same
    field, as the sorted-order is assumed to be datetime order
    (which works if the files are not renamed).
    
    """
    data_dir = Path(md.dir_raw_nc_data, model_id,
                    experiment_id, variable_id)
    
    data_files_in = sorted(
        [str(x) for x in Path(data_dir).glob(
            f"{variable_id}*{model_id}*"
            + f"{experiment_id}*{member_id}*.nc")])
    
    if md.variable_domain[variable_id] == "ocn":
        coords = md.lonlat_ocn_nc_names[model_id]
    else:
        coords = md.lonlat_atm_nc_names[model_id]
    
    lon, lat, fld = nct.get_arrays(data_files_in, coords,
                                   [variable_id])
    
    if original_miss_value != set_miss_value:
        fld = np.where(fld >= original_miss_value,
                       set_miss_value, fld)
    
    return lon, lat, fld



def prepare_siconc(variable_id="siconc",
        member_id        = md.default_member_id,
        model_id         = md.default_model_id,
        experiment_id    = md.default_experiment_id,
        load_field_2d_kw = {},
        mask_reg_names   = ["lakes", "baltic_sea", "black_sea"]
    ):
    """Load sea ice concentration and make necessary data
    transformations for sea ice extent, area, and ice-edge
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
    
    
    load_field_2d_kw["member_id"] = member_id
    load_field_2d_kw["model_id"] = model_id
    load_field_2d_kw["experiment_id"] = experiment_id
    
    print(f"Preparing siconc ({load_field_2d_kw['model_id']}, "
        + f"{load_field_2d_kw['experiment_id']}, "
        + f"{load_field_2d_kw['member_id']})")
    
    lon, lat, siconc = field_2D(variable_id, **load_field_2d_kw)
    
    lon = lon % 360.0  # need 0-360 range for lake masking
    
    # Need 2D coordinates for lake masking:
    if np.ndim(lon) == 1 and np.ndim(lat) == 1:
        lon, lat = np.meshgrid(lon, lat)
    
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
                (lon >= reg_def[sub_reg][0]) &
                (lon <= reg_def[sub_reg][1]) &
                (lat >= reg_def[sub_reg][2]) &
                (lat <= reg_def[sub_reg][3]),
                np.nan, 1.0)
            
            siconc *= sub_reg_mask[np.newaxis,:,:]
    
    return siconc
