"""Settings and writing routines for NetCDF output data."""

from datetime import datetime as dt, timezone as tz
from pathlib import Path

from netCDF4 import Dataset, date2num
import numpy as np

from src import metadata as md


def nc_save_dir(diagnostic_id, experiment_id, model_id):
    """Absolute path to default directory for saving data.
    """
    return Path(md.dir_out_nc_data, diagnostic_id,
                experiment_id)


def nc_file_name(diagnostic_id, experiment_id, model_id):
    """Default file name for each """
    return f"{diagnostic_id}_{experiment_id}_{model_id}.nc"


def nc_file_path(diagnostic_id, experiment_id, model_id):
    """Absolute path to specified diagnostic, experiment,
    and model saved data.
    """
    return Path(
        nc_save_dir(diagnostic_id, experiment_id, model_id),
        nc_file_name(diagnostic_id, experiment_id, model_id)
    )
    


nc_file_attrs_model_name = "source"
nc_file_attrs_experiment_name = "experiment_id"
nc_file_attrs_member_name = "member_id"
nc_file_attrs_grid_name = "grid_label"

default_nc_file_attrs = {
    "author"     : "Jake R. Aylmer",
    "contact"    : "j.aylmer@reading.ac.uk",
    "coauthors"  : "David G. Ferreira, Daniel L. Feltham",
    "comment"    : "Climate model diagnostics derived from the CMIP6 "
                   + "archive for the analysis in the work, 'Modulation "
                   + "of future sea ice loss by ocean heat transport'. "
                   + "This dataset contains one diagnostic for one model "
                   + f"({nc_file_attrs_model_name}) and one experiment "
                   + f"({nc_file_attrs_experiment_name}).",
    "references" : "The associated work is not published; see the "
                   + "author's PhD thesis for further information: "
                   + "doi:10.48683/1926.00108418",
    "institution": "Centre for Polar Observation and Modelling (CPOM), "
                   + "Department of Meteorology, University of Reading, UK",
    "title"      : "CMIP6 diagnostics: {}, {}",
}

# Dimension names
# ------------------------- #
nc_time_dim_name           = "time"       # (coord. var.)
nc_member_dim_name         = "member"     # for the 4 ripf member axes (not constructed as coord. vars.)
nc_y_dim_name              = "y"          # for 2D lon/lat data (aux. coord. var.)
nc_x_dim_name              = "x"          # for 2D lon/lat data (aux. coord. var.)
nc_lon_dim_name            = "lon"        # for 1D lon data (coord. var.)
nc_lat_dim_name            = "lat"        # for 1D lat data (coord. var.)
nc_vert_dim_name           = "vert"       # vertices, for cell bounds
nc_tbnd_dim_name           = "bnd"        # to distinguish from vertices
nc_ref_lat_n_dim_name      = "ref_lat_n"  # for reference latitude (north; coord. var.)
nc_ref_lat_s_dim_name      = "ref_lat_s"  # for reference latitude (south; coord. var.)
nc_ref_lat_single_dim_name = "ref_lat"    # for datasets with single reference latitudes (coord. var.)

# Time data
# ------------------------- #
nc_time_name = nc_time_dim_name  # must match dim. for coord. var
nc_time_type = np.float64

nc_time_units = {
    "piControl" : "days since 0001-01-01",
    "historical": "days since 1850-01-01",
    "ssp126"    : "days since 1850-01-01",
    "ssp245"    : "days since 1850-01-01",
    "ssp370"    : "days since 1850-01-01",
    "ssp585"    : "days since 1850-01-01"
}

nc_calendar = {
    "piControl" : "proleptic_gregorian",
    "historical": "standard",  # "gregorian" is deprecated
    "ssp126"    : "standard",
    "ssp245"    : "standard",
    "ssp370"    : "standard",
    "ssp585"    : "standard"
}


# Coordinate data
# ----------------------------------------------- #
nc_lon_2d_name  = "lon"
nc_lon_2d_type  = np.float64
nc_lon_2d_attrs = {
    "long_name"    : "longitude",
    "standard_name": "longitude",
    "units"        : "degrees_east"
}

nc_lat_2d_name  = "lat"
nc_lat_2d_type  = np.float64
nc_lat_2d_attrs = {
    "long_name"    : "latitude",
    "standard_name": "latitude",
    "units"        : "degrees_north"
}

nc_lon_1d_name  = nc_lon_dim_name  # must match for coord. variable
nc_lon_1d_type  = np.float64
nc_lon_1d_attrs = nc_lon_2d_attrs.copy()

nc_lat_1d_name  = nc_lat_dim_name  # must match for coord. variable
nc_lat_1d_type  = np.float64
nc_lat_1d_attrs = nc_lat_2d_attrs.copy()

# -- reference latitudes (e.g., for area-
#    integrated quantities at selected
#    latitudes). These coordinate variables
#    are set to unlimited dimensions.
nc_ref_lat_n_name  = nc_ref_lat_n_dim_name
nc_ref_lat_s_name  = nc_ref_lat_s_dim_name
nc_ref_lat_n_type  = np.float64
nc_ref_lat_s_type  = np.float64
nc_ref_lat_n_attrs = {
    "long_name"    : "reference_latitude_north",
    "standard_name": "latitude",
    "units"        : "degrees_north"
}
nc_ref_lat_s_attrs = nc_ref_lat_n_attrs.copy()
nc_ref_lat_s_attrs["long_name"] = "reference_latitude_south"

# -- for datasets with single reference latitudes
#    (no north/south sets; e.g., hfbasin)
nc_ref_lat_single_name = nc_ref_lat_single_dim_name
nc_ref_lat_single_type = np.float64
nc_ref_lat_single_attrs = {
    "long_name"    : "reference_latitude",
    "standard_name": "latitude",
    "units"        : "degrees_north"
}

# Ensemble member label axes names:
# ----------------------------------------------- #
nc_ripf_type   = np.int32
nc_ripf_labels = {
    "r": "realisation_id",
    "i": "initialisation_id",
    "p": "physics_id",
    "f": "forcing_id"
}

# Standard units for different types of quantity
# (not set by any function here: they must be
# set in the main processing scripts and passed
# in the nc_variable_attrs dict):
# 
# ----------------------------------------------- #
field_units = {
    "cellarea"     : "m2",
    "heatflux"     : "W m-2",
    "heattransport": "PW",
    "seaicearea"   : "1e6 km2",
    "seaiceedge"   : "degrees_north",
    "seaiceextent" : "1e6 km2",
    "temperature"  : "K"
}

# For sea ice diagnostics:
nc_siconc_threshold_name = "siconc_threshold"


# Common diagnostic (for filenames) and variable (in the
# netCDF files) name qualifiers.
# 
# These should be used in each data generation script in
# conjunction with the nc_diag_name() and nc_var_name()
# functions defined below in this module.
# 
# They are also used in the api._diagnostic_definitions to
# load saved diagnostics from aliases/keywords.
# --------------------------------------------------------- #

# --- time methods (frequency/average)
diag_nq_monthly  = "_mon"
diag_nq_yearly   = "_yr"

nc_var_nq_monthly = ""
nc_var_nq_yearly  = ""


# --- space methods (mean, integral, etc.)
diag_nq_vertical_integral   = "_ver_int"
diag_nq_zonal_mean          = "_zm"
diag_nq_area_mean           = "_area_mean"
diag_nq_area_integral       = "_hor_int"

nc_var_nq_vertical_integral = "_ver_int"
nc_var_nq_zonal_mean        = "_zm"
nc_var_nq_area_mean         = ""
nc_var_nq_area_integral     = ""


# --- other methods (interpolation, etc.)
diag_nq_native               = "_gn"
diag_nq_native_interp        = "_gn_interp"
diag_nq_cell_center_approx   = "_cc_approx"

nc_var_nq_native             = ""
nc_var_nq_native_interp      = ""
nc_var_nq_cell_center_approx = ""

# ------ for ice edge latitude interpolation:
diag_nq_4deg_bil   = "_4deg_bil"
diag_nq_2deg_bil   = "_2deg_bil"
diag_nq_1deg_bil   = "_1deg_bil"
diag_nq_05deg_bil  = "_05deg_bil"   # 0.5 deg
diag_nq_025deg_bil = "_025deg_bil"  # 0.25 deg

nc_var_nq_4deg_bil   = ""
nc_var_nq_2deg_bil   = ""
nc_var_nq_1deg_bil   = ""
nc_var_nq_05deg_bil  = ""
nc_var_nq_025deg_bil = ""


# Other common NetCDF settings:
# ----------------------------------------------- #
nc_fill_value = md.default_new_missing_value


# =========================================================== #


def _format_nc_name(name="diag", time_methods="",
                    space_methods="", hemi="", other_methods="",
                    keep_case=False):
    """Set name from specified attributes in standardised order
    [to be used for the diagnostic (first) part of netCDF file
    names and variable names only]. All inputs must be strings
    and should, with the exception of name and hemi, be one of
    the name qualifiers defined in the netcdf module.
    """
    
    dname = name  # start with basic name (e.g., "sie")
    
    # Apply in the following order:
    for x in [space_methods, time_methods,
              other_methods, hemi
        ]:
        
        if x.startswith("_"):
            x = x[1:]
        if x.endswith("_"):
            x = x[:-1]
        
        x = x.replace("__", "_")
        
        if len(x) > 0:
            dname += f"_{x}"
    
    dname = dname.replace("__", "_")
    
    if not keep_case:
        dname = dname.lower()
    
    return dname



def diag_name(name="diag", time_methods="", space_methods="",
              hemi="", other_methods="", keep_case=False):
    """Get netCDF file diagnostic according to standardised
    format (overwrites hemi as no hemi information is needed for
    file names).
    """
    return _format_nc_name(name=name,
        time_methods=time_methods,
        space_methods=space_methods, hemi="",
        other_methods=other_methods, keep_case=keep_case)



def nc_var_name(name="diag", time_methods="", space_methods="",
                hemi="", other_methods="",keep_case=False):
    """Get netCDF variable name according to standardised
    format.
    """
    return _format_nc_name(name=name,
        time_methods=time_methods,
        space_methods=space_methods, hemi=hemi,
        other_methods=other_methods, keep_case=keep_case)



def _set_attributes(item, attrs):
    """Generic attribute setter for arbitrary item and
    dictionary of attributes. Applies each key in alphabetical
    order.
    """
    for k in sorted(attrs.keys()):
        setattr(item, k, attrs[k])



def create_time_variables(ncdat,
        units="days since 1850-01-01",
        calendar="gregorian"
    ):
    """Create time and time bounds variable on a given netCDF4
    Dataset instance opened in r+ or w mode. Does everything
    except set the data.
    """
    
    # Create dimensions for time and time bounds:
    ncdat.createDimension(nc_time_dim_name, None)
    ncdat.createDimension(nc_tbnd_dim_name, 2)
    
    # Create time variable:
    ncdat.createVariable(nc_time_name, nc_time_type,
                         (nc_time_dim_name,))
    
    # Create time bounds variable [should have same name (plus
    # "_bnds") and type as the ordinary time variable]. Note
    # also the CF convention is not to repeat other attributes
    # already in the ordinary time variable, including units
    # and calendar, as the bounds are themselves considered
    # metadata of the corresponding variable:
    ncdat.createVariable(f"{nc_time_name}_bnds", nc_time_type,
        (nc_time_dim_name, nc_tbnd_dim_name))
    
    nc_time_attrs = {"units": units, "calendar": calendar,
                     "bounds": f"{nc_time_name}_bnds"}
    
    _set_attributes(ncdat.variables[nc_time_name],
                    nc_time_attrs)



def set_time_data_yearly(ncdat, year_range):
    """Set time data for yearly averages accounting for time
    units and calendar in an existing netCDF4 Dataset instance.
    
    Requires a year_range = (year_start, year_end), and sets
    the time coordinates at the center of bounds y and y+1
    for y = year_start ... year_end inclusive.
    """
    
    cal_dates = []
    cal_dates_bnds = []
    for yr in range(year_range[0], year_range[1]+1):
        dt_0 = dt(yr, 1, 1, 0, 0, 0, 0)
        dt_1 = dt(yr+1, 1, 1, 0, 0, 0, 0)
        
        cal_dates.append(dt_0 + (dt_1 - dt_0)/2)
        cal_dates_bnds.append([dt_0, dt_1])
    
    d2n_kw = {"units": ncdat.variables[nc_time_name].units,
        "calendar": ncdat.variables[nc_time_name].calendar}
    
    ncdat.variables[nc_time_name][:] = \
        date2num(cal_dates, **d2n_kw)
    
    ncdat.variables[f"{nc_time_name}_bnds"][:,:] = \
        date2num(cal_dates_bnds, **d2n_kw)



def set_time_data_monthly(ncdat, year_range):
    """Set time data for yearly averages accounting for time
    units and calendar in an existing netCDF4 Dataset instance.
    
    Requires a year_range = (year_start, year_end), and sets
    the time coordinates at the center of bounds y and y+1
    for y = year_start ... year_end inclusive.
    """
    
    cal_dates = []
    cal_dates_bnds = []
    for yr in range(year_range[0], year_range[1]+1):
        for mn in range(1, 13):
            dt_0 = dt(yr, mn, 1, 0, 0, 0, 0)
            dt_1 = dt(yr+(mn==12), 1 if mn==12 else (mn+1),
                      1, 0, 0, 0, 0)
            
            cal_dates.append(dt_0 + (dt_1 - dt_0)/2)
            cal_dates_bnds.append([dt_0, dt_1])
    
    d2n_kw = {"units": ncdat.variables[nc_time_name].units,
        "calendar": ncdat.variables[nc_time_name].calendar}
    
    ncdat.variables[nc_time_name][:] = date2num(cal_dates,
                                                **d2n_kw)
    
    ncdat.variables[f"{nc_time_name}_bnds"][:,:] = date2num(
        cal_dates_bnds, **d2n_kw)



def _set_coord_2d_data(ncdat, coord_data, coord_bnds_data,
        coord_name=nc_lon_2d_name,
        coord_type=nc_lon_2d_type,
        coord_attrs=nc_lon_2d_attrs):
    """Set 2D (i.e., mesh grid) coordinate variables and
    associated cell bounds.
    
    
    Parameters
    ----------
    ncdat : netCDF4.Dataset instance in r+ or w mode
    
    coord_data : 2D array (ny, nx) of coordinate values
        (cell centres)
    
    coord_bnds_data: 3D array (ny, nx, 4) of coordinate
        bounds (cell vertices). Can also be None (for
        point data), in which case no bounds variable is
        created.
    
    
    Optional parameters
    -------------------
    Defaults to settings for longitude:
    
    coord_name: str, name of coordinate variable
    coord_type: coordinate variable data type
    coord_attrs: dict of str, coordinate variable attributes.
    
    """
    
    cvar_dims = (nc_y_dim_name, nc_x_dim_name)
    
    cvar = ncdat.createVariable(coord_name, coord_type,
                                cvar_dims)
    
    cvar[:,:] = coord_data
    
    if coord_bnds_data is not None:
        # Create bounds variable without metadata
        # (as is the CF convention):
        cvar_bnds = ncdat.createVariable(
            f"{coord_name}_bnds", coord_type,
            (nc_y_dim_name, nc_x_dim_name, nc_vert_dim_name)
        )
        cvar_bnds[:,:,:] = coord_bnds_data
        coord_attrs["bounds"] = f"{coord_name}_bnds"
    
    _set_attributes(cvar, coord_attrs)


def set_lon_2d_data(ncdat, lon_2d_data, lon_2d_bnds_data):
    """Set 2D (i.e., mesh-grid) longitude coordinates and
    associated cell bounds (wrapper for _set_coord_2d_data).
    """
    _set_coord_2d_data(ncdat, lon_2d_data, lon_2d_bnds_data,
        coord_name=nc_lon_2d_name, coord_type=nc_lon_2d_type,
        coord_attrs=nc_lon_2d_attrs)



def set_lat_2d_data(ncdat, lat_2d_data, lat_2d_bnds_data):
    """Set 2D (i.e., mesh-grid) latitude coordinates and
    associated cell bounds (wrapped for _set_coord_2d_data).
    """
    _set_coord_2d_data(ncdat, lat_2d_data, lat_2d_bnds_data,
        coord_name=nc_lat_2d_name, coord_type=nc_lat_2d_type,
        coord_attrs=nc_lat_2d_attrs)


def _set_coord_1d_data(ncdat, coord_data, coord_bnds_data,
        coord_dim_name=nc_lon_dim_name,
        coord_name=nc_lon_1d_name,
        coord_type=nc_lon_1d_type,
        coord_attrs=nc_lon_1d_attrs):
    """Set 1D coordinate variables and associated cell
    bounds, if applicable.
    
    
    Parameters
    ----------
    ncdat : netCDF4.Dataset instance in r+ or w mode
    
    coord_data : 1D array (n,) of coordinate values
        (cell centres)
    
    coord_bnds_data: 2D array (n, 2) of coordinate
        bounds (cell vertices). Can also be None (for
        point data), in which case no bounds variable is
        created.
    
    
    Optional parameters
    -------------------
    Defaults to settings for longitude:
    
    coord_dim_name: str, name of corresponding dimension
    coord_name: str, name of coordinate variable
    coord_type: coordinate variable data type
    coord_attrs: dict of str, coordinate variable attributes.
    
    """
    
    cvar_dims = (coord_dim_name,)
    
    cvar = ncdat.createVariable(coord_name, coord_type,
                                cvar_dims)
    
    cvar[:] = coord_data
    
    if coord_bnds_data is not None:
        # Create bounds variable without metadata
        # (as is the CF convention):
        cvar_bnds = ncdat.createVariable(
            f"{coord_name}_bnds", coord_type,
            (coord_dim_name, nc_vert_dim_name)
        )
        cvar_bnds[:,:] = coord_bnds_data
        coord_attrs["bounds"] = f"{coord_name}_bnds"
    
    _set_attributes(cvar, coord_attrs)



def set_lon_1d_data(ncdat, lon_1d_data, lon_1d_bnds_data):
    """Set 1D longitude coordinates and associated cell bounds
    (wrapper for _set_coord_1d_data).
    """
    _set_coord_1d_data(ncdat, lon_1d_data, lon_1d_bnds_data,
        coord_dim_name=nc_lon_dim_name,
        coord_name=nc_lon_1d_name, coord_type=nc_lon_1d_type,
        coord_attrs=nc_lon_1d_attrs)



def set_lat_1d_data(ncdat, lat_1d_data, lat_1d_bnds_data):
    """Set 1D longitude coordinates and associated cell bounds
    (wrapper for _set_coord_1d_data).
    """
    _set_coord_1d_data(ncdat, lat_1d_data, lat_1d_bnds_data,
        coord_dim_name=nc_lat_dim_name,
        coord_name=nc_lat_1d_name, coord_type=nc_lat_1d_type,
        coord_attrs=nc_lat_1d_attrs)



def set_ref_lat_data(ncdat, ref_lat_n_data, ref_lat_s_data,
                     ref_lat_n_bnds_data=None,
                     ref_lat_s_bnds_data=None):
    """Set 1D reference latitude coordinates and data."""
    
    for name, type, dim_name, data, data_bnds, attrs in zip(
            [nc_ref_lat_n_name, nc_ref_lat_s_name],
            [nc_ref_lat_n_type, nc_ref_lat_s_type],
            [nc_ref_lat_n_dim_name, nc_ref_lat_s_dim_name],
            [ref_lat_n_data, ref_lat_s_data],
            [ref_lat_n_bnds_data, ref_lat_s_bnds_data],
            [nc_ref_lat_n_attrs, nc_ref_lat_s_attrs]
        ):
        
        ncdat.createVariable(name, type, (dim_name,))
        _set_attributes(ncdat.variables[name], attrs)
        ncdat.variables[name][:] = data
        
        if data_bnds is not None:
            bnds = ncdat.createVariable(f"{name}_bnds", type,
                (dim_name, nc_vert_dim_name))
            bnds[:,:] = data_bnds
            setattr(ncdat.variables[name], "bounds",
                    f"{name}_bnds")



def _get_rx_val(ripf, i='r'):
    """Convert ripf string to integer value of specified index.
    
    Examples
    --------
    >>> get_e_val("r14i1p2f1", "r")
    14
    >>> get_e_val("r14i1p2f1", "p")
    2
    >>> type(get_e_val("r1i1p1f1"))
    int
    
    """
    e = "ripf"
    j = [ripf.index(x) for x in e]
    
    if i == "f":
        return int(ripf[j[3]+1:])
    else:
        return int(ripf[j[e.index(i)]+1:j[e.index(i)+1]])


def set_ripf_data(ncdat, ens_members):
    """Set ensemble member label (discrete axes) data.
    The data are determined from the 'ripf' labels directly.
    """
    
    ncdat.createDimension(nc_member_dim_name, None)
    
    for x in nc_ripf_labels.keys():
        ncdat.createVariable(nc_ripf_labels[x], np.int32,
                             (nc_member_dim_name,))
        ncdat.variables[nc_ripf_labels[x]][:] = \
            [_get_rx_val(m, x) for m in ens_members]



def set_nc_history_attr_timestamp(ncfa):
    """Set the history global attribute to the current
    time. Updates the 'history' key in the attributes
    dictionary, not the actual netCDF file.
    """
    ncfa["history"] = "created " + \
        dt.now(tz.utc).strftime("%H:%M UTC %d %b %Y")



def set_nc_source_attr(ncfa, model_id):
    """Set the source global attribute. Updates the dictionary
    key, not the actual netCDF file."""
    ncfa[nc_file_attrs_model_name] = model_id



def set_nc_experiment_attr(ncfa, experiment_id):
    """Set the experiment global attribute. Updates the
    dictionary key, not the actual netCDF file."""
    ncfa[nc_file_attrs_experiment_name] = experiment_id



def set_nc_member_attr(ncfa, member_id):
    """Set the member_id global attribute (only for saving
    fields of single members, without member dimension). Updates
    the dictionary key, not the actual netCDF file."""
    ncfa[nc_file_attrs_member_name] = member_id


def set_nc_grid_attr(ncfa, grid_label):
    """Set the grid_label global attribute (only for saving
    fields of single members). Updates the dictionary key, not
    the actual netCDF file."""
    ncfa[nc_file_attrs_grid_name] = grid_label


def set_nc_title_attr(ncfa, model_id, experiment_id):
    """Format the title global attribute with the model
    and experiment names. Updates the dictionary key,
    not the actual netCDF file.
    """
    ncfa["title"] = ncfa["title"].format(model_id,
                                         experiment_id)



def prepare_to_save(save_dir_in=None, file_name_in=None,
        diagnostic_id="", experiment_id="", model_id=""
    ):
    """Determine directory and file name (if defaults), create
    directory if needed, and return verified directory and
    file names.
    """
    
    if save_dir_in is None:  # set to default
        save_dir_out = nc_save_dir(diagnostic_id,
                                        experiment_id,
                                        model_id)
    else:
        save_dir_out = Path(save_dir_in)
    
    save_dir_out.mkdir(parents=True, exist_ok=True)
    
    if file_name_in is None:  # set to default
        file_name_out = nc_file_name(diagnostic_id,
                                     experiment_id,
                                     model_id)
    
    return save_dir_out, file_name_out


# =========================================================== #
# Saving routines
# =========================================================== #

def save_yearly(field, diagnostic_id, nc_field_name, model_id,
        member_ids, experiment_id,
        year_range=(1, 500),
        longitude=None, latitude=None,
        longitude_bnds=None, latitude_bnds=None,
        nc_field_type=np.float64,
        nc_field_attrs={},
        nc_global_attrs={},
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
    
    if np.ndim(field) < 2 or np.ndim(field) > 4:
        raise Exception("Not programmed for fields of shape "
                        + f"{np.shape(field)}")
    
    save_dir, file_name = \
        prepare_to_save(save_dir, file_name, diagnostic_id,
                        experiment_id, model_id)
    
    # Prepare the global attributes: ------------------------ #
    nc_file_attrs = default_nc_file_attrs.copy()
    set_nc_title_attr(nc_file_attrs, model_id, experiment_id)
    set_nc_history_attr_timestamp(nc_file_attrs)
    set_nc_source_attr(nc_file_attrs, model_id)
    set_nc_experiment_attr(nc_file_attrs, experiment_id)
    
    for extra_attr in nc_global_attrs.keys():
        nc_file_attrs[extra_attr] = nc_global_attrs[extra_attr]
    # ------------------------------------------------------- #
    
    with Dataset(Path(save_dir, file_name), "w") as ncout:
        
        _set_attributes(ncout, nc_file_attrs)
        
        create_time_variables(ncout,
                              nc_time_units[experiment_id],
                              nc_calendar[experiment_id])
        set_time_data_yearly(ncout, year_range)
        set_ripf_data(ncout, member_ids)
        
        if np.ndim(field) == 4:
            # Field is (time, member, y, x)
            
            if np.ndim(longitude) == 2:
                xdim_name = nc_x_dim_name
                ydim_name = nc_y_dim_name
                nvert = 4
                lon_name = nc_lon_2d_name
                lat_name = nc_lat_2d_name
            else:
                xdim_name = nc_lon_dim_name
                ydim_name = nc_lat_dim_name
                nvert = 2
                lon_name = nc_lon_1d_name
                lat_name = nc_lat_1d_name
            
            ncout.createDimension(ydim_name, np.shape(field)[2])
            ncout.createDimension(xdim_name, np.shape(field)[3])
            ncout.createDimension(nc_vert_dim_name, nvert)
            
            if nvert == 4:
                set_lon_2d_data(ncout, longitude,longitude_bnds)
                set_lat_2d_data(ncout, latitude, latitude_bnds)
            else:
                set_lon_1d_data(ncout, longitude,longitude_bnds)
                set_lat_1d_data(ncout, latitude, latitude_bnds)
            
            ncfield = ncout.createVariable(nc_field_name,
                nc_field_type, (nc_time_dim_name,
                    nc_member_dim_name, ydim_name,
                    xdim_name),
                fill_value=nc_fill_value, zlib=True
            )
            ncfield[:,:,:,:] = field
            
            if (nvert == 4 and
                    "coordinates" not in nc_field_attrs.keys()):
                nc_field_attrs["coordinates"] = \
                    f"{lon_name} {lat_name}"
        
        elif np.ndim(field) == 3:
            # Field is (time, member, y)
            # 
            # NOTE: for reference latitude, there is a separate
            # function (to account for integration from both
            # poles). This one just does proper latitude.
            
            ncout.createDimension(nc_lat_dim_name,
                                  np.shape(field)[2])
            set_lat_1D_data(ncout, latitude, latitude_bnds)
            
            ncfield = ncout.createVariable(nc_field_name,
                nc_field_type,
                (nc_time_dim_name, nc_member_dim_name,
                    nc_lat_dim_name),
                fill_value=nc_fill_value, zlib=True
            )
            ncfield[:,:,:] = field
        
        else:  # np.ndim(field) == 2 (time, member)
            ncfield = ncout.createVariable(
                nc_field_name, nc_field_type,
                (nc_time_dim_name, nc_member_dim_name),
                fill_value=nc_fill_value, zlib=True
            )
            ncfield[:,:] = field
        
        _set_attributes(ncfield, nc_field_attrs)
    
    print(f"Saved: {str(Path(save_dir, file_name))}")



def save_time_series_by_hemi(field_n, field_s,
        diagnostic_id, nc_field_name_n,
        nc_field_name_s, model_id,
        member_ids, experiment_id,
        year_range=(1, 500),
        nc_field_type=np.float64,
        nc_field_attrs_n={},
        nc_field_attrs_s={},
        scalar_coord_var_n_kw={"name": None},
        scalar_coord_var_s_kw={"name": None},
        nc_global_attrs={},
        monthly=False, save_dir=None, file_name=None
    ):
    """Save time series data, one for each hemisphere.
    
    
    Parameters
    ----------
    field : numpy array (nt, n_ens)
        Field data to save: 2D
    
    diagnostic_id   : str (for automated directory structure)
    nc_field_name_n : str (netCDF variable name, north)
    nc_field_name_s : str (netCDF variable name, south)
    model_id        : str (for metadata)
    experiment_id   : str (for metadata)
    
    
    Optional parameters
    -------------------
    year_range : iterable of int (year_start, year_end)
        Start and end year (defines the length of the time
        axis, and the time coordinates of field along axis=0).
        End year is included.
    
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
    
    for x in [field_n, field_s]:
        if np.ndim(x) != 2:
            raise Exception("Not programmed for fields of "
                            + f"shape {np.shape(x)}")
    
    save_dir, file_name = \
        prepare_to_save(save_dir, file_name, diagnostic_id,
                        experiment_id, model_id)
    
    # Prepare the global attributes: ------------------------ #
    nc_file_attrs = default_nc_file_attrs.copy()
    set_nc_title_attr(nc_file_attrs, model_id, experiment_id)
    set_nc_history_attr_timestamp(nc_file_attrs)
    set_nc_source_attr(nc_file_attrs, model_id)
    set_nc_experiment_attr(nc_file_attrs, experiment_id)
    
    for extra_attr in nc_global_attrs.keys():
        nc_file_attrs[extra_attr] = nc_global_attrs[extra_attr]
    # ------------------------------------------------------- #
    
    with Dataset(Path(save_dir, file_name), "w") as ncout:
        
        _set_attributes(ncout, nc_file_attrs)
        
        create_time_variables(ncout,
                              nc_time_units[experiment_id],
                              nc_calendar[experiment_id])
        if monthly:
            set_time_data_monthly(ncout, year_range)
        else:
            set_time_data_yearly(ncout, year_range)
        set_ripf_data(ncout, member_ids)
        
        # This script allows the creation of scalar coordinate
        # variables as it is used for sea ice extent:
        for scv_kw in [scalar_coord_var_n_kw,
                       scalar_coord_var_s_kw]:
            
            if scv_kw["name"] is not None:
                ncout.createVariable(scv_kw["name"],
                                     scv_kw["dtype"])
                ncout.variables[scv_kw["name"]][:] = \
                    scv_kw["value"]
                _set_attributes(ncout.variables[scv_kw["name"]],
                                scv_kw["attrs"])
        
        for name, data, attrs in zip(
                [nc_field_name_n, nc_field_name_s],
                [field_n, field_s],
                [nc_field_attrs_n, nc_field_attrs_s]):
            
            ncout.createVariable(name, nc_field_type,
                (nc_time_dim_name, nc_member_dim_name),
                fill_value=nc_fill_value, zlib=True
            )
            ncout.variables[name][:,:] = data
            _set_attributes(ncout.variables[name], attrs)
            
    
    print(f"Saved: {str(Path(save_dir, file_name))}")

    

def save_yearly_ref_lat(fint_n, fint_s, ref_lat_n, ref_lat_s,
        diagnostic_id, nc_field_name_n, nc_field_name_s,
        model_id, member_ids, experiment_id,
        year_range=(1, 500),
        nc_field_type_n=np.float64,
        nc_field_type_s=np.float64,
        ref_lat_n_bnds=None, ref_lat_s_bnds=None,
        unlimited_ref_lat_dim=True,
        nc_field_attrs_n={},
        nc_field_attrs_s={},
        nc_global_attrs={},
        save_dir=None, file_name=None
    ):
    """Save a yearly data for integrations/averages above a
    reference latitude (one variable for north, one for south).
    """
    
    for x in [fint_n, fint_s]:
        if np.ndim(x) != 3:
            raise Exception("Field should be 3 dimensional "
                            + "(nt, n_ens, n_ref_lat); instead,"
                            + " got array of shape "
                            + f"{np.shape(x)}.")
    
    save_dir, file_name = \
        prepare_to_save(save_dir, file_name, diagnostic_id,
                        experiment_id, model_id)
    
    # Prepare the global attributes: ------------------------ #
    nc_file_attrs = default_nc_file_attrs.copy()
    set_nc_title_attr(nc_file_attrs, model_id, experiment_id)
    set_nc_history_attr_timestamp(nc_file_attrs)
    set_nc_source_attr(nc_file_attrs, model_id)
    set_nc_experiment_attr(nc_file_attrs, experiment_id)
    
    for extra_attr in nc_global_attrs.keys():
        nc_file_attrs[extra_attr] = nc_global_attrs[extra_attr]
    # ------------------------------------------------------- #
    
    with Dataset(Path(save_dir, file_name), "w") as ncout:
        
        _set_attributes(ncout, nc_file_attrs)
        
        create_time_variables(ncout,
                              nc_time_units[experiment_id],
                              nc_calendar[experiment_id])
        set_time_data_yearly(ncout, year_range)
        set_ripf_data(ncout, member_ids)
        
        ncout.createDimension(nc_ref_lat_n_dim_name,
            None if unlimited_ref_lat_dim else len(ref_lat_n))
        ncout.createDimension(nc_ref_lat_s_dim_name,
            None if unlimited_ref_lat_dim else len(ref_lat_s))
        
        if (ref_lat_n_bnds is not None
                or ref_lat_s_bnds is not None
            ):
            ncout.createDimension(nc_vert_dim_name, 2)
        
        set_ref_lat_data(ncout, ref_lat_n, ref_lat_s,
                         ref_lat_n_bnds, ref_lat_s_bnds)
        
        for f_name, f_type, dim_name, data, attrs in zip(
                [nc_field_name_n, nc_field_name_s],
                [nc_field_type_n, nc_field_type_s],
                [nc_ref_lat_n_dim_name, nc_ref_lat_s_dim_name],
                [fint_n, fint_s],
                [nc_field_attrs_n, nc_field_attrs_s]
            ):
            
            ncout.createVariable(f_name,
                f_type, (nc_time_dim_name,
                    nc_member_dim_name, dim_name),
                fill_value=nc_fill_value, zlib=True
            )
            ncout.variables[f_name][:,:,:] = data
            _set_attributes(ncout.variables[f_name],
                            attrs)
    
    print(f"Saved: {str(Path(save_dir, file_name))}")


# =========================================================== #
# Saving specific fields
# =========================================================== #


def save_areacell(areacell, model_id, member_id, experiment_id,
        grid_label="gn", which="a",
        longitude=None, latitude=None,
        longitude_bnds=None, latitude_bnds=None,
        nc_field_type=np.float64,
        nc_field_attrs={},
        nc_global_attrs={},
        save_dir=None, file_name=None
    ):
    """Save areacella or areacello data in a standardised
    format.
    
    
    Parameters
    ----------
    areacell: array (ny, nx)
        Horizontal grid cell areas (areacella or areacello).
    
    model_id     : str (for metadata)
    member_id    : str (for metadata)
    experiment_id: str (for metadata)
    
    
    Optional parameters
    -------------------
    which : str, {'a' or 'o'}, default = 'a'
        Whether this is atmosphere (a) or ocean (o) data.
    
    longitude, latitude: 2D (ny, nx) or 1D (ny,), (nx,) arrays
        Longitude and latitude coordinates in degrees_east and
        degrees_north, respectively.
    
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
    
    if np.ndim(areacell) != 2:
        raise Exception("Grid cell areas should be 2D (got "
                        + f"{np.shape(areacell)})")
    
    which = which.lower()
    if which not in ["a", "o"]:
        raise Exception("Parameter \'which\' must be \'a\' or "
                        + "\'o\' (got \'" + which + "\')")
    
    if save_dir is None:
        save_dir_out = Path(md.dir_out_nc_data,
                            f"areacell{which}")
    else:
        save_dir_out = Path(save_dir)
    save_dir_out.mkdir(parents=True, exist_ok=True)
    
    if file_name is None:  # set to default
        file_name = f"areacell{which}_{model_id}.nc"
    
    
    # Prepare the global attributes: ------------------------ #
    nc_file_attrs = default_nc_file_attrs.copy()
    set_nc_title_attr(nc_file_attrs, model_id, experiment_id)
    set_nc_history_attr_timestamp(nc_file_attrs)
    set_nc_source_attr(nc_file_attrs, model_id)
    set_nc_member_attr(nc_file_attrs, member_id)
    set_nc_grid_attr(nc_file_attrs, grid_label)
    set_nc_experiment_attr(nc_file_attrs, experiment_id)
    
    for extra_attr in nc_global_attrs.keys():
        nc_file_attrs[extra_attr] = nc_global_attrs[extra_attr]
    # ------------------------------------------------------- #
    
    with Dataset(Path(save_dir_out,file_name), "w") as ncout:
        
        _set_attributes(ncout, nc_file_attrs)
        
        if np.ndim(longitude) == 2:
            xdim_name = nc_x_dim_name
            ydim_name = nc_y_dim_name
            nvert = 4
            lon_name = nc_lon_2d_name
            lat_name = nc_lat_2d_name
        else:
            xdim_name = nc_lon_dim_name
            ydim_name = nc_lat_dim_name
            nvert = 2
            lon_name = nc_lon_1d_name
            lat_name = nc_lat_1d_name
        
        ncout.createDimension(ydim_name, np.shape(areacell)[0])
        ncout.createDimension(xdim_name, np.shape(areacell)[1])
        ncout.createDimension(nc_vert_dim_name, nvert)
        
        if nvert == 4:
            set_lon_2d_data(ncout, longitude, longitude_bnds)
            set_lat_2d_data(ncout, latitude, latitude_bnds)
        else:
            set_lon_1d_data(ncout, longitude, longitude_bnds)
            set_lat_1d_data(ncout, latitude, latitude_bnds)
        
        nc_areacell = ncout.createVariable(f"areacell{which}",
            nc_field_type, (ydim_name, xdim_name),
            fill_value=nc_fill_value, zlib=True
        )
        nc_areacell[:,:] = areacell
        
        if nvert == 4 and "coordinates" not in nc_field_attrs.keys():
            nc_field_attrs["coordinates"] = \
                f"{lon_name} {lat_name}"
        
        _set_attributes(nc_areacell, nc_field_attrs)
    
    print(f"Saved: {str(Path(save_dir_out, file_name))}")



def save_yearly_ref_lat_single(fint, ref_lat,
        diagnostic_id, nc_field_name,
        model_id, member_ids, experiment_id,
        year_range=(1, 500),
        nc_field_type=np.float64,
        ref_lat_bnds=None,
        unlimited_ref_lat_dim=True,
        nc_field_attrs={},
        nc_global_attrs={},
        save_dir=None, file_name=None
    ):
    """Save a yearly data for integrations/averages above a
    reference latitude (one variable only, no north/south
    labelling). This is primarily for hfbasin diagnostics.
    """
    
    if np.ndim(fint) != 3:
        raise Exception("Field should be 3 dimensional "
                        + "(nt, n_ens, n_ref_lat); instead,"
                        + " got array of shape "
                        + f"{np.shape(fint)}.")
    
    save_dir, file_name = \
        prepare_to_save(save_dir, file_name, diagnostic_id,
                        experiment_id, model_id)
    
    # Prepare the global attributes: ------------------------ #
    nc_file_attrs = default_nc_file_attrs.copy()
    set_nc_title_attr(nc_file_attrs, model_id, experiment_id)
    set_nc_history_attr_timestamp(nc_file_attrs)
    set_nc_source_attr(nc_file_attrs, model_id)
    set_nc_experiment_attr(nc_file_attrs, experiment_id)
    
    for extra_attr in nc_global_attrs.keys():
        nc_file_attrs[extra_attr] = nc_global_attrs[extra_attr]
    # ------------------------------------------------------- #
    
    with Dataset(Path(save_dir, file_name), "w") as ncout:
        
        _set_attributes(ncout, nc_file_attrs)
        
        create_time_variables(ncout,
                              nc_time_units[experiment_id],
                              nc_calendar[experiment_id])
        set_time_data_yearly(ncout, year_range)
        set_ripf_data(ncout, member_ids)
        
        ncout.createDimension(nc_ref_lat_single_dim_name,
            None if unlimited_ref_lat_dim else len(ref_lat))
        
        if ref_lat_bnds is not None:
            ncout.createDimension(nc_vert_dim_name, 2)
        
        _set_coord_1d_data(ncout, ref_lat, ref_lat_bnds,
            coord_dim_name=nc_ref_lat_single_dim_name,
            coord_name=nc_ref_lat_single_name,
            coord_type=nc_ref_lat_single_type,
            coord_attrs=nc_ref_lat_single_attrs)
        
        ncout.createVariable(nc_field_name,
            nc_field_type, (nc_time_dim_name,
                nc_member_dim_name,
                nc_ref_lat_single_dim_name),
            fill_value=nc_fill_value, zlib=True
        )
        ncout.variables[nc_field_name][:,:,:] = fint
        
        _set_attributes(ncout.variables[nc_field_name],
                        nc_field_attrs)
    
    print(f"Saved: {str(Path(save_dir, file_name))}")



def save_ice_edge_latitude_per_longitude(iel_n, iel_s, lon,
        diagnostic_id, nc_field_name_n, nc_field_name_s,
        model_id, member_ids, experiment_id,
        year_range=(1, 500),
        nc_field_type=np.float64,
        lon_bnds=None, monthly=True,
        nc_field_attrs_n={},
        nc_field_attrs_s={},
        scalar_coord_var_n_kw={"name": None},
        scalar_coord_var_s_kw={"name": None},
        nc_global_attrs={},
        save_dir=None, file_name=None
    ):
    """Save ice edge latitude data as a function of time,
    ensemble member and longitude.
    """
    
    for x in [iel_n, iel_s]:
        if np.ndim(x) != 3:
            raise Exception("Field should be 3 dimensional "
                            + "(nt, n_ens, n_lon); instead,"
                            + " got array of shape "
                            + f"{np.shape(x)}.")
    
    save_dir, file_name = \
        prepare_to_save(save_dir, file_name, diagnostic_id,
                        experiment_id, model_id)
    
    # Prepare the global attributes: ------------------------ #
    nc_file_attrs = default_nc_file_attrs.copy()
    set_nc_title_attr(nc_file_attrs, model_id, experiment_id)
    set_nc_history_attr_timestamp(nc_file_attrs)
    set_nc_source_attr(nc_file_attrs, model_id)
    set_nc_experiment_attr(nc_file_attrs, experiment_id)
    
    for extra_attr in nc_global_attrs.keys():
        nc_file_attrs[extra_attr] = nc_global_attrs[extra_attr]
    # ------------------------------------------------------- #
    
    with Dataset(Path(save_dir, file_name), "w") as ncout:
        
        _set_attributes(ncout, nc_file_attrs)
        
        create_time_variables(ncout,
                              nc_time_units[experiment_id],
                              nc_calendar[experiment_id])
        if monthly:
            set_time_data_monthly(ncout, year_range)
        else:
            set_time_data_yearly(ncout, year_range)
        set_ripf_data(ncout, member_ids)
        
        # This script allows the creation of scalar coordinate
        # variables as it is used for sea ice extent:
        for scv_kw in [scalar_coord_var_n_kw,
                       scalar_coord_var_s_kw
            ]:
            if scv_kw["name"] is not None:
                ncout.createVariable(scv_kw["name"],
                                     scv_kw["dtype"])
                ncout.variables[scv_kw["name"]][:] =\
                    scv_kw["value"]
                _set_attributes(ncout.variables[scv_kw["name"]],
                                scv_kw["attrs"])
        
        ncout.createDimension(nc_lon_dim_name, len(lon))
        
        if lon_bnds is not None:
            ncout.createDimension(nc_vert_dim_name, 2)
        
        set_lon_1d_data(ncout, lon, lon_bnds)
        
        for name, data, attrs in zip(
                [nc_field_name_n, nc_field_name_s],
                [iel_n, iel_s],
                [nc_field_attrs_n, nc_field_attrs_s]
            ):
            
            ncout.createVariable(name, nc_field_type,
                (nc_time_dim_name, nc_member_dim_name,
                    nc_lon_dim_name),
                fill_value=nc_fill_value, zlib=True
            )
            ncout.variables[name][:,:,:] = data
            
            _set_attributes(ncout.variables[name], attrs)
    
    print(f"Saved: {str(Path(save_dir, file_name))}")
