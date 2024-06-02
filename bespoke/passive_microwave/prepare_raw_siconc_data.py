# Prepare raw sea ice concentration data from passive microwave
# observations for processing by other scripts to calculate
# required diagnostics.
# 
# The idea is that this script is run first, which interprets
# the data exactly as obtained from the National Snow and Ice
# Data Center (NSIDC), and puts them into a consistent format.
# 
# In particular, monthly data is available with one netCDF file
# per month of data. This script loads in an arbitrary number of
# them and concatenates in time, adds the longitude/latitude
# coordinates, adds a proper time coordinate, and saves to a
# specified netCDF file/location. One complication, also dealt
# with here, is that the variable name corresponding to sea ice
# concentration in these files is the specific satellite sensor
# name and this differs per input file.
# 
# ============================================================ #
from argparse import ArgumentParser
from calendar import monthrange
from datetime import datetime as dt, timedelta
from pathlib import Path

import netCDF4 as nc
import numpy as np

from process_cmip6_data.src import (
    metadata as md, netcdf as nf, utils)


def get_datetime_from_nc_metadata(ncdat):
    """Get datetime coordinate from an input netCDF file (as a
    netCDF4 Dataset instance).
    
    Raw data is one time step (one month), but has a time
    coordinate which can be read. However, some files do not
    (e.g., missing data prior). In such cases, this function
    attempts to read the global attribute 'time_coverage_start'
    and work out the date/time from that. If that fails for some
    unforeseen reason, then this function returns None.
    """
    if "time" in ncdat.variables.keys():
        # Data contains a time variable; should be straight-
        # forward to get the date/time:
        ncdat_dt = nc.num2date(ncdat.variables["time"][0],
            units=ncdat.variables["time"].units,
            calendar=ncdat.variables["time"].calendar)
    
    elif hasattr(ncdat, "time_coverage_start"):
        
        # No time coordinate; get from the global attribute
        # called "time_coverage_start". This depends on the
        # attribute having a consistent format, so wrap this in
        # a try/except block and return None if it does not
        # work:
        try:
            ncdat_dt = dt.strptime(
                getattr(ncdat, "time_coverage_start"),
                "%Y-%m-%d %H:%M:%S")
        except Exception:
            ncdat_dt = None
    
    else:
        ncdat_dt = None
    
    return ncdat_dt



def get_datetime_from_file_name(fname, strfmt="%Y%m",
                                strsep="_", dtindex=-2):
    """Get the datetime for one input file name (fname, str),
    which should work if files are downloaded from the NSIDC and
    their names are not changed.
    
    
    Optional parameters
    -------------------
    strfmt : str, default = "%Y%m"
        Format string for interpreting datetime, using standard
        format codes (%Y for year etc.).
    
    strsep : str, default = "_"
        Character to split fname by. Should be the character
        surrounding the date string in the file name
     
    dtindex : int, default = -2
        Index of the resulting split string corresponding to the
        date part of the file name. By default, select the
        second counting from the end.
    
    
    Returns
    -------
    datetime.datetime
    
    """
    return dt.strptime(fname.split(strsep)[dtindex], strfmt)



def load_raw_data(input_files,
        dimensions=md.nsidc_grid_dims["n"],
        time_units=md.nsidc_nc_time_units,
        time_calendar="standard",
        new_missing_value=md.default_new_missing_value,
        old_flag_values=[251, 252, 253, 254],
        new_flag_values=[md.default_new_missing_value]*4):
    """Load raw monthly data from specified input file names,
    concatenate with respect to time, and fix y-coordinate
    order.
    
    
    Parameters
    ----------
    input_files : length n_files list of str or pathlib.Path
        Input netCDF files, full paths. It is assumed these are
        the raw NSIDC data files which each contain one time
        step, occassionally being empty for pre-1987 (i.e.,
        missing days for daily data). They are loaded in order
        after sorting, but are sorted with respect to time after
        concatenating anyway.
    
    
    Optional parameters
    -------------------
    dimensions : tuple of int, (ny, nx)
        Grid dimensions.
    
    time_units : str
        NetCDF time unit string for the output data.
    
    time_calendar : str
        NetCDF calendar string for the output data.
    
    old_missing_values : list of length n_miss_flags
        Missing value flag values in original data. Note these 
        are different between NSIDC-0051 (NASA Team) and NSIDC-
        0079 (Bootstrap).
    
    new_missing_values : list of length n_miss_flags of float
        New missing value flags to set corresponding to those in
        old_missing_values. Each list element can also be None,
        in which case the original value is left as-is.
    
    
    Returns
    -------
    time : array (n_files,)
        Time coordinates in units/calendar specified by inputs
        time_units and time_calendar.
    
    time_bnds : array (n_files, 2)
        Time coordinate bounds.
    
    siconc : array (n_files, ny, nx)
        Sea ice concentration data as loaded, with all missing
        value flags set to NaN and axis=1 reversed. Raw NSIDC
        data is already in units of fraction (not percentage).
    
    """
    
    # Sea ice concentration variable is different in each input
    # file because it is named by the specific sensor used to
    # collect the data. So in the loop below, need to work out
    # what it is. Do this by eliminating known variables that 
    # are presumably always present and correspond to the
    # coordinates or something irrelevant:
    known_var_names = ["time", "x", "y", "crs"]
    
    # Sorting the input files by name presumably puts then in
    # the correct time order, and does if raw NSIDC file names
    # are not modified, but this is checked at the end using the
    # tiem data anyway:
    input_files = sorted(list(input_files))
    n_files = len(input_files)
    
    # Prepare output arrays (this is the only reason to require
    # the dimensions up front; could potentially re-design this
    # so that the dimensions are worked out in place):
    date = np.zeros(n_files, dtype=dt)
    date_bnds = np.zeros((n_files, 2), dtype=dt)
    
    # Set siconc to missing initially and fill in data as input
    # files are iterated over. Also set the data type to double
    # precision (numpy/netCDF4 converts this from ubyte):
    siconc = new_missing_value*np.ones((n_files, *dimensions),
                                       dtype=np.float64)
    
    for n in range(n_files):
        
        with nc.Dataset(input_files[n], "r") as nc_n:
            
            # Get time coordinate as datetime, ideally from the
            # netCDF time variable or global attribute:
            dt_n = get_datetime_from_nc_metadata(nc_n)
            
            # ... failing that, get it from the file name:
            if dt_n is None:
                dt_n = get_datetime_from_file_name(
                    input_files[n],
                    strfmt="%Y%m%d" if daily else "%Y%m")
            
            # ... that could potentially still go wrong, but
            # no point handling it here as... well, I don't know
            # what the problem would be until it happens anyway.
            
            # Set datetime bounds (monthly averages assumed):
            date_bnds[n,0] = dt(dt_n.year, dt_n.month, 1)
            date_bnds[n,1] = date_bnds[n,0] + timedelta(
                days=monthrange(dt_n.year, dt_n.month)[1])
            
            # Set datetime coordinate at the center of bounds:
            date[n] = date_bnds[n,0] \
                + 0.5*(date_bnds[n,1] - date_bnds[n,0])
            
            # Get list of netCDF variable names except for those
            # in the "rule out being siconc" case above:
            var_names_n = [x for x in nc_n.variables.keys()
                           if x not in known_var_names]
            
            # It should be that var_names_n either contains no
            # variable names (missing data) or it contains one
            # variable name (which must correspond to sea ice
            # concentration). Otherwise, not sure what's
            # happening but skip anyway:
            if len(var_names_n) == 1:
                siconc[n,:,:] = np.array(
                    nc_n.variables[var_names_n[0]])[0,:,:]
    
    # In case input files were not in the correct time order,
    # sort with respect to the datetimes determined above:
    jsort = np.argsort(date)
    
    date = date[jsort]
    date_bnds = date_bnds[jsort,:]
    siconc = siconc[jsort,:,:]
    
    # Convert datetimes to time values based on specified output
    # units/calendar:
    time = nc.date2num(date, units=time_units,
                       calendar=time_calendar)
    
    time_bnds = nc.date2num(date_bnds, units=time_units,
                            calendar=time_calendar)
    
    # Final processing of siconc data: set missing values to
    # their new missing values:
    for j in range(len(old_flag_values)):
        if new_flag_values[j] is not None:
            siconc = np.where(siconc == old_flag_values[j],
                              new_flag_values[j], siconc)
    
    # And reverse y-direction (axis=1; this also is done for
    # auxilliary data, lon/lat/masks/area) so that it is the
    # correct way around:
    siconc = siconc[:,::-1,:]
    
    return time, time_bnds, siconc



def load_lon_lat(hemi="n",
        coordinate_var_names=["longitude", "latitude"]):
    """Load coordinates (longitude/latitude), and reverse the
    axis=0 dimension direction. These come directly from raw
    data (NSIDC dataset ID: NSIDC-0771), the file paths of which
    are set in metadata.nsidc_areacell_file.
    """
    with nc.Dataset(md.nsidc_lonlat_file[hemi], "r") as ncdat:
        lon = np.array(ncdat.variables[coordinate_var_names[0]])
        lat = np.array(ncdat.variables[coordinate_var_names[1]])
    
    lon = lon[::-1,:]
    lat = lat[::-1,:]
    
    return lon, lat



def main():
    
    # Arguments required are input netCDF files to prepare,
    # which should be the raw files exactly as downloaded from
    # the NSIDC (see data references), an output file path, a
    # string option to identify which dataset corresponds to
    # input data (to determine missing flag values; NSIDC-0051
    # or NSIDC-0079), a string to determine for which hemisphere
    # sea ice concentration is being prepared (to determine
    # which coordinates to load and add to the output; "n" or
    # "s"):
    prsr = ArgumentParser()
    prsr.add_argument("-i", "--infiles", type=str, nargs="*",
                      default=[], help="Input file paths")
    prsr.add_argument("-o", "--outfile", type=str,
                      default="~/output.nc",
                      help="Output file name (full path)")
    prsr.add_argument("--hemisphere", type=str, default="n",
                      choices="ns")
    prsr.add_argument("-d", "--datasetid", type=str,
                      default="NSIDC-0051_nasateam_v2")
    prsr.add_argument("--newmissvals", type=float, nargs="*",
                      default=[np.nan]*4)
    cmd = prsr.parse_args()
    
    # -------------------------------------------------------- #
    # Output netCDF attributes: most don't need to be set here
    # since this data is intermediate.
    # 
    # However, these time attributes (particularly the units and
    # calendar) are carried throughout the rest of the 
    # processing chain; i.e., sea ice edge final data will also
    # have these units calendar:
    nc_out_time_attrs = {
        "bounds"  : f"{nf.nc_time_name}_bnds",
        "calendar": "standard",
        "units"   : md.nsidc_nc_time_units}
    
    # Also, "coordinates" is needed for interpolation purposes:
    nc_out_siconc_attrs = {
        "coordinates": f"{nf.nc_lon_2d_name} {nf.nc_lat_2d_name}",
        "long_name": "sea ice concentration",
        "standard_name": "sea_ice_area_fraction",
        "units": "1"}
    
    # No other variable or global attributes need be set
    # -------------------------------------------------------- #
    
    ny, nx = md.nsidc_grid_dims[cmd.hemisphere]
    lon, lat = load_lon_lat(cmd.hemisphere)
    
    # For NSIDC-0051 or NSIDC-0081 (NASA Team), flag values
    # correspond to:
    #     251: pole_hole_mask
    #     252: unused
    #     253: coast
    #     254: land
    # 
    # For NSIDC-0079 (Bootstrap), the flag values are different
    # because the data type is different (short as opposed to
    # ubyte)! So here, flag values (also less detailed than
    # 0051/0081) are
    #    1100: missing
    #    1200: land
    
    if "NSIDC-0079" in cmd.datasetid:
        old_flag_values = [1100, 1200]
    else:
        old_flag_values = [251, 252, 253, 254]
    
    # Default is just to set them all to be missing (NaN):
    new_flag_values = list(cmd.newmissvals)
    
    time, time_bnds, siconc = load_raw_data(cmd.infiles,
        dimensions=(ny, nx),
        time_units=nc_out_time_attrs["units"],
        time_calendar=nc_out_time_attrs["calendar"],
        new_missing_value=md.default_new_missing_value,
        old_flag_values=old_flag_values,
        new_flag_values=new_flag_values)
    
    # Additional processing needed for the northern hemisphere:
    if cmd.hemisphere == "n":
        
        # Apply NSIDC-0622 "valid ice" mask. There is no
        # corresponding mask for the southern hemisphere. There
        # is a different mask for each month of the year, so
        # need the dates to apply these (without having to
        # assume Jan, Feb, ..., Dec, Jan, ... and exact years
        # input to script):
        date = nc.num2date(time,
            units=nc_out_time_attrs["units"],
            calendar=nc_out_time_attrs["calendar"])
        
        # Create a time-dependent mask:
        vmask = np.zeros((np.shape(siconc)))
        
        # Loop over months:
        for j in range(12):
            # Determine where to apply mask for month j+1
            # (j = 0 --> January, etc.):
            jt = np.array([x.month == j+1 for x in date])
            
            # File for mask for month j+1:
            mask_j = \
                md.nsidc_valid_ice_mask_n_nc_file_fmt.format(
                    j+1)
            
            with nc.Dataset(mask_j, "r") as ncdat_vmask:
                vmask[jt,:,:] = np.array(
                    ncdat_vmask.variables["valid_ice_flag"]
                    )[::-1,:]
                # also reversing the array direction as is done
                # done for the coordinates/sea ice concentration
        
        # The mask values are:
        # 
        #     0   ocean
        #     1   valid_ice
        #     2   coastline (should be implicitly masked)
        #     3   land (should be implicitly masked)
        #     4   lakes
        #
        # The vmask array should load fine as an integer such
        # that, e.g., vmask == 3 should correctly identify land,
        # but to be safe, just check against some arbitrary
        # threshold (here 0.1):
        # 
        for vval in [3, 4]:  # land and lakes
            siconc = np.where(abs(vmask - vval) < 0.1,
                              md.default_new_missing_value,
                              siconc)
        
        # valid_ice flag (1) note the reversed True/False
        # assignment: here, we *want* to include valid_ice and
        # exclude everywhere else. Logically, this should make
        # the other cases above redundant (?), but doesn't hurt
        # to cover all bases anyway:
        siconc = np.where(abs(vmask - 1) < 0.1, siconc,
                          md.default_new_missing_value)
    # -------------------------------------------------------- #
    
    # Save the prepared data:
    with nc.Dataset(cmd.outfile, "w") as nc_out:
        
        nc_out.createDimension(nf.nc_x_dim_name, nx)
        nc_out.createDimension(nf.nc_y_dim_name, ny)
        nc_out.createDimension(nf.nc_time_dim_name, None)
        nc_out.createDimension(nf.nc_tbnd_dim_name, 2)
        
        # Time variable:
        nc_out.createVariable(nf.nc_time_name, time.dtype,
                              (nf.nc_time_dim_name,))
        # Time bounds variable:
        nc_out.createVariable(nc_out_time_attrs["bounds"],
                              time_bnds.dtype,
                              (nf.nc_time_dim_name,
                               nf.nc_tbnd_dim_name))
        # Longitude:
        nc_out.createVariable(nf.nc_lon_2d_name, lon.dtype,
                              (nf.nc_y_dim_name,
                               nf.nc_x_dim_name))
        # Latitude:
        nc_out.createVariable(nf.nc_lat_2d_name, lat.dtype,
                              (nf.nc_y_dim_name,
                               nf.nc_x_dim_name))
        # Sea ice concentration (not a final field and won't
        # keep so don't need to compress here):
        nc_out.createVariable("siconc", siconc.dtype,
                              (nf.nc_time_dim_name,
                               nf.nc_y_dim_name,
                               nf.nc_x_dim_name),
                              zlib=False)
        
        nf._set_attributes(nc_out.variables[nf.nc_time_name],
                           nc_out_time_attrs)
        
        nf._set_attributes(nc_out.variables[nf.nc_lon_2d_name],
                           nf.nc_lon_2d_attrs)
        
        nf._set_attributes(nc_out.variables[nf.nc_lat_2d_name],
                           nf.nc_lat_2d_attrs)
        
        nf._set_attributes(nc_out.variables["siconc"],
                           nc_out_siconc_attrs)
        
        nc_out.variables[nf.nc_time_name][:] = time
        nc_out.variables[nc_out_time_attrs["bounds"]][:,:] = \
            time_bnds
        
        nc_out.variables[nf.nc_lon_2d_name][:,:] = lon
        nc_out.variables[nf.nc_lat_2d_name][:,:] = lat
        
        nc_out.variables["siconc"][:,:,:] = siconc
    
    print(f"Saved: {cmd.outfile}")


if __name__ == "__main__":
    main()
