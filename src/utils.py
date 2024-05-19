"""General python utilities."""

from datetime import datetime as dt
import warnings

import netCDF4 as nc
import numpy as np


def cftime_to_datetime(dt_cf):
    """Convert a list or array of cf datetimes to regular python
    datetimes.
    """
    return np.array([dt(x.year, x.month, x.day, x.hour,
        x.minute, x.second, x.microsecond) for x in dt_cf])



def ncdump(nc_fid, indent_size=4):
    """Rough implementation of NCAR ncdump utility in Python.
    Adapted from:
    
    schubert.atmos.colostate.edu/~cslocum/nf_example.html
    [Accessed May 2023; not accessible on 19 May 2024]
    
    Parameters
    ----------
    ncfid : str or pathlib.Path instance
        Absolute path to a netCDF file
    
    """
    
    indent = " "*indent_size
    
    def _print_ncattr(key):
        try:
            for ncattr in ncdat.variables[key].ncattrs():
                print(f"{indent*2}{var}:{ncattr} = " + repr(
                    ncdat.variables[key].getncattr(ncattr)))
        except:
            pass
    
    print("\n" + str(nc_fid) + "\n")
    
    with nc.Dataset(nc_fid, mode="r") as ncdat:
        
        print("global attributes:")
        for nc_attr in ncdat.ncattrs():
            print(indent + "%s:" % nc_attr,
                  repr(ncdat.getncattr(nc_attr)))
        
        nc_dims = list(ncdat.dimensions.keys())
        print("dimensions:")
        for dim in nc_dims:
            print(f"{indent}{dim} = {len(ncdat.dimensions[dim])}")
        
        nc_vars = [x for x in ncdat.variables.keys()]
                   #if x not in nc_dims]
        print("variables:")
        for var in nc_vars:
            print(f"{indent}{ncdat.variables[var].dtype} {var}"
                  + f"{ncdat.variables[var].dimensions}")
            _print_ncattr(var)
    
    print("")



def nc_get_time_units(ncfile):
    """Determine the units and calendar of a specified NetCDF
    file (input as string or pathlib.Path instance). Returns
    the units and calendar attribute. If the time variable
    cannot be identified, or if the attribute(s) does not exist,
    returns None in place (warnings are raised).
    """

    with nc.Dataset(ncfile, mode="r") as ncdat:

        # Determine the variable name of the time coordinate.
        # This is almost always either "time" or "t".
        # Raise warning if it doesn't exist and return None for
        # each attribute:
        if "time" in ncdat.variables:
            nc_t_name = "time"
        elif "t" in ncdat.variables:
            nc_t_name = "t"
        else:
            nc_t_name = None
            warnings.warn("Unable to identify time variable"
                          + f" in file {ncfile}")

        if nc_t_name is None:
            t_units = None
            t_calendar = None
        else:
            if hasattr(ncdat.variables[nc_t_name], "units"):
                t_units = ncdat.variables[nc_t_name].units
            else:
                t_units = None
                warnings.warn("variable \'{nc_t_name}\' has "
                              + f"no units in file {ncfile}")

            if hasattr(ncdat.variables[nc_t_name], "calendar"):
                t_calendar = ncdat.variables[nc_t_name].calendar
            else:
                t_calendar = None
                warnings.warn("variable \'{nc_t_name}\' has "
                              + "no calendar specified in file "
                              + ncfile)

    return t_units, t_calendar



def nc_get_arrays(ncfiles, coordinate_vars=[],
        diagnostic_vars=[],
        np_concatenate_kwargs={}):
    """
    Load an arbitrary set of coordinates and variables
    from a set of NetCDF data files as NumPy arrays.
    
    Variables specified in coordinate_vars are loaded from
    one file only and not concantenated, while those specified
    in diagnostic_vars are loaded from all files and are
    concatenated (along axis=0 by default).
    
    
    Parameters
    ----------
    ncfiles : list of str
        
        List of paths to each NetCDF file. Note that files
        are loaded and diagnostic variables (see below) are
        concatenated in the order of this list. There is no
        checking of time coordinates.
    
    coordinate_vars : list of str
        
        Names (matching variable name in the NetCDF metadata)
        of the coordinates to load from the ncfiles[0]. This
        list is for non-concatenated variables (latitude,
        pressure level, etc., although note that the time, if
        required, should be put as a diagnostic variable).
    
    diagnostic_vars : list of str
        
        Names (matching variable name in the NetCDF metadata)
        of the diagnostic variables to be loaded from all
        ncfiles and concatenated. These are concatenated in
        the order of ncfiles.
    
    
    Optional parameters
    -------------------
    np_concatenate_kwargs : dict, default = {}
        
        Keyword arguments passed to NumPy.concatenate(). Note
        that the default concatenation axis, 0, usually
        corresponds to time and thus, usually, does not need
        to be changed.
    
    
    Returns
    -------
    coordinate_var_1, ..., diagnostic_var_1, ...
    
    NumPy arrays of the requested coordinate variables, in
    the order specified by coordinate_vars, followed by the
    requested diagnostic variables, in the order specified by
    diagnostic_vars.
    
    
    Example
    -------
    Suppose the variable "t2m", as a function of time,
    latitude, and longitude, on a global 1 degree fixed grid,
    is stored across 5 NetCDF files containing monthly averages
    for one year each. This can be loaded using:
    
        >>> import numpy as np
        >>> import nc_tools as nct
        >>> 
        >>> ncfiles = ['./file%i.nc' % j for j in range(5)']
        >>>
        >>> lat, lon, time, t2m = nct.gather_data_arrays(
        ...     ncfiles, ['lat', 'lon'], ['time', 't2m']
        ... )
        >>> 
        >>> np.shape(lon)
        (360,)
        >>> np.shape(lat)
        (180,)
        >>> np.shape(time)
        (60,)
        >>> np.shape(t2m)
        (60, 180, 360)
    
    """
    
    n_cvars = len(coordinate_vars)
    n_dvars = len(diagnostic_vars)
    
    # Each array will be stored in a list (separately for
    # coordinates (cvar_arrays) and diagnostics (dvar_arrays):
    cvar_arrays = [None for j in range(n_cvars)]
    dvar_arrays = [None for j in range(n_dvars)]
    
    # Load the first dataset explicitly to get dimension
    # information, coordinates, and initialise arrays:
    with nc.Dataset(ncfiles[0], mode="r") as ncds_0:
        
        for k in range(n_cvars):
            cvar_arrays[k] = np.array(
                ncds_0.variables[coordinate_vars[k]])
        
        for k in range(n_dvars):
            dvar_arrays[k] = np.array(
                ncds_0.variables[diagnostic_vars[k]])
    
    # Now loop over remaining datasets and concatenate
    # the diagnostics:
    for j in range(1, len(ncfiles)):
        
        with nc.Dataset(ncfiles[j], mode="r") as ncds_j:
            
            for k in range(len(diagnostic_vars)):
                
                dvar_arrays[k] = np.concatenate(
                    (dvar_arrays[k],
                     np.array(ncds_j.variables[
                              diagnostic_vars[k]])
                    ), **np_concatenate_kwargs
                )
    
    return tuple(cvar_arrays + dvar_arrays)
