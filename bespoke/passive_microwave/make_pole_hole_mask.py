# Create north pole hole masks for SSMR/SSMI/S passive microwave
# data, either native polar-stereographic grid or regridded data
# depending on inputs.
# ============================================================ #

from argparse import ArgumentParser
import calendar
from datetime import datetime as dt
from pathlib import Path

import netCDF4 as nc
import numpy as np

from process_cmip6_data.src import metadata as md, netcdf as nf


# Actual/known latitudes of pole hole boundaries for each
# instrument, in degrees_north, taken from the documentation
# (doi:10.5067/MPYG15WAA4WX):
pole_hole_lat_bound = {"SMMR": 84.50, "SSMI": 87.20,
                       "SSMIS": 89.18}

# Start and end date/times, as strings to appear as-is in the
# output netCDF global attributes, corresponding to
# applicability of mask for each instrument:
pmask_time_valid = {
    "SMMR" : ("1978-11-01 00:00:00", "1987-08-21 00:00:00"),
    "SSMI" : ("1987-08-21 00:00:00", "2008-01-01 00:00:00"),
    "SSMIS": ("2008-01-01 00:00:00",
              dt.now().strftime("%Y-%m-%d 00:00:00"))}
# End date of SSMIS is really "present", but prefer to keep
# consistent format so just set to whatever date (not time) that
# the script is run.


def main():
    
    # Arguments required are input netCDF file (x1) to determine
    # mask from (should be output from "prepare_raw_siconc_data"
    # python script), and several string options: one to
    # identify which mask (SMMR, SSMI, or SSMIS), one to
    # identify which dataset corresponds to input netCDF data
    # (for data references; NSIDC-0051 or NSIDC-0079), and a
    # long and a short description string for interpolation
    # method (only used when interpolated data is detected,
    # based on dimensions of coordinate arrays in input data):
    prsr = ArgumentParser()
    prsr.add_argument("-i", "--infile", type=str, default="",
                      help="Sea ice concentration file to "
                           + "determine mask from.")
    prsr.add_argument("--whichmask", type=str, default="SMMR",
                      choices=list(pole_hole_lat_bound.keys()))
    prsr.add_argument("-d", "--datasetid", type=str,
                      default="NSIDC-0051",
                      choices=md.defined_nsidc_datasets)
    prsr.add_argument("--interplongdescription", type=str,
                      default="This mask is computed after "
                              + "interpolating raw data to a "
                              + "0.5 degree regular lon-lat "
                              + "grid using bilinear "
                              + "interpolation.")
    prsr.add_argument("--interpshortdescription", type=str,
                      default="0.5 degree regular lon-lat")
    cmd = prsr.parse_args()
    # -------------------------------------------------------- #
    
    # Load data (need coordinates and siconc, first time step
    # only):
    with nc.Dataset(cmd.infile, "r") as ncdat:
        lon = np.array(ncdat.variables["lon"])
        lat = np.array(ncdat.variables["lat"])
        # Save datetime of first time step for output netCDF
        # global attributes:
        date_used = nc.num2date(ncdat.variables["time"][0],
            units=ncdat.variables["time"].units,
            calendar=ncdat.variables["time"].calendar)
        siconc = np.array(ncdat.variables["siconc"])[0]
    
    # Determine whether input data is on the native grid or has
    # been interpolated based on the number of dimensions of the
    # coordinate arrays:
    if np.ndim(lon) == 1 and np.ndim(lat) == 1:
        # Must be re-gridded (regular lon-lat). Save a string
        # for output netCDF "grid_label" attribute and generate
        # 2D equivalent arrays for calculation below (when
        # saving, the original 1D lon/lat arrays are saved):
        grid_label = "gr"
        lon_use, lat_use = np.meshgrid(lon, lat)
    elif np.ndim(lon) == 2 and np.ndim(lat) == 2:
        # Must be raw data (i.e., native grid).
        lon_use = lon
        lat_use = lat
        grid_label = "gn"
    else:
        raise Exception("Unsupported coordinates of shape "
                        + repr(np.shape(lon)) + " (lon) or "
                        + repr(np.shape(lat)) + " (lat).")
    
    lat_bound = pole_hole_lat_bound[cmd.whichmask]
    
    if grid_label == "gn":
        # Here, find where the sea ice concentration data is
        # missing and latitude exceeds the pole-hole boundary.
        # Note that this is necessarily approximate: 
        pmask = np.where(
            (np.isnan(siconc)) & (lat_use > lat_bound),
            md.default_new_missing_value, 1.0)
    else:
        # Here (gr), can just use the known latitude boundaries
        # of the pole hole to determine where to set missing,
        # noting that latitudes are at the center of grid cells
        # of uniform width. This is therefore exact, although
        # there is still a possibility that data originally
        # missing due to the pole hole gets remapped to a
        # slightly southward grid point and thus is not captured
        # by this mask:
        dlat = abs(lat[1] - lat[0])
        # Set to missing wherever the upper latitude limit
        # exceeds (i.e., grid cell overlaps) the actual pole-
        # hole boundary:
        pmask = np.where(lat_use + 0.5*dlat > lat_bound,
                         md.default_new_missing_value, 1.0)
    
    ny, nx = np.shape(pmask)
    
    # ------------------------------------------------------ #
    # NetCDF attributes:
    # ------------------------------------------------------ #
    
    nc_global_attrs = nf.default_nc_file_attrs.copy()
    
    nf._set_nc_history_attr_timestamp(nc_global_attrs)
    nc_global_attrs["grid_label"] = grid_label
    nc_global_attrs["time_coverage_start"] = \
        pmask_time_valid[cmd.whichmask][0]
    nc_global_attrs["time_coverage_end"] = \
        pmask_time_valid[cmd.whichmask][1]
    nc_global_attrs["source_id"] = cmd.datasetid
    nc_global_attrs["source"] = (md.nsidc_source[cmd.datasetid]
                                 + " [2]")
    
    nc_global_attrs["title"] = (
        f"{cmd.whichmask} north pole-hole mask ("
        + ("native NSIDC 25 km north polar stereographic grid"
           if grid_label == "gn" else
           f"interpolated {cmd.interpshortdescription} grid")
        + ")")
    
    nc_global_attrs["comment"] = (
        "Observations of sea ice concentration, diagnostics "
        + "derived from raw data obtained from the National "
        + "Snow and Ice Data Center (NSIDC), for the analysis "
        + "presented in Aylmer et al. 2024 [1]. This dataset "
        + "contains a north \'pole-hole\' mask, which is used "
        + "when processing the raw data to calculate other "
        + "diagnostics and determined from NSIDC dataset ID "
        + f"{cmd.datasetid[-4:]} [2], monthly data for year "
        + f"{date_used.year} and month {date_used.month} "
        + f"({calendar.month_name[date_used.month]}). ")
    
    if grid_label == "gr":
        nc_global_attrs["comment"] += cmd.interplongdescription
    else:
        nc_global_attrs["comment"] += ("This mask is computed "
            + "on the native NSIDC 25 km north polar "
            + "stereographic grid.")
    
    if cmd.whichmask == "SMMR":
        nc_global_attrs["comment"] += (
            " \'time_coverage_end\' is based on an inspection "
            + "of daily concentration data, which shows that "
            + "the pole hole corresponding to this instrument "
            + "occurs until 20 August 1987, after which it "
            + "switches to that corresponding to the SSMI pole "
            + "hole.")
    
    elif cmd.whichmask == "SSMI":
        nc_global_attrs["comment"] += (
            " \'time_coverage_start\' is based on an "
            + "inspection of daily concentration data, which "
            + "shows that the (larger) pole hole corresponding "
            + "to the SMMR instrument occurs until 20 August "
            + "1987, after which it switches to that "
            + "corresponding to the SSMI pole hole (i.e., this "
            + "mask).")
    
    elif cmd.whichmask == "SSMIS":
        nc_global_attrs["comment"] += (
            " \'time_coverage_end\' is just included for "
            + "consistency with other pole mask datasets; this "
            + "mask is valid beyond the stated end time.")
    
    nc_global_attrs["references"] = ("[1] " + nf.paper_reference
        + nf.nc_list_delimiter + "[2] "
        + md.nsidc_data_reference[cmd.datasetid])
    
    nc_pmask_attrs = {
        "comment": f"pmask = "
                   + repr(md.default_new_missing_value)
                   + " where masked, 1.0 otherwise",
        "units"  : "1",
        "long_name": f"{cmd.whichmask} north pole-hole mask"}
    
    if grid_label == "gn":
        nc_pmask_attrs["coordinates"] = \
            f"{nf.nc_lon_2d_name} {nf.nc_lat_2d_name}"
    # "coordinates" attribute not needed if regridded (because
    # lon and lat are coordinate variables in this case)
    
    # Output NetCDF file name and directory:
    nc_out_dir = Path(md.dir_out_nc_data_nsidc, "auxilliary")
    nc_out_dir.mkdir(parents=True, exist_ok=True)
    nc_out_file = Path(nc_out_dir,
        f"pmask_{grid_label}_{cmd.whichmask.lower()}.nc")
    
    # ------------------------------------------------------ #
    
    with nc.Dataset(nc_out_file, "w") as nc_out:
        
        nf._set_attributes(nc_out, nc_global_attrs)
        
        # Dimensions and shapes of coordinates depend on whether
        # data is on the native grid or interpolated:
        if grid_label == "gn":
            
            nc_out.createDimension(nf.nc_x_dim_name, nx)
            nc_out.createDimension(nf.nc_y_dim_name, ny)
            
            # Longitude:
            nc_out.createVariable(nf.nc_lon_2d_name, lon.dtype,
                                  (nf.nc_y_dim_name,
                                   nf.nc_x_dim_name))
            nf._set_attributes(
                nc_out.variables[nf.nc_lon_2d_name],
                nf.nc_lon_2d_attrs)
            
            nc_out.variables[nf.nc_lon_2d_name][:,:] = lon
                
            # Latitude:
            nc_out.createVariable(nf.nc_lat_2d_name, lat.dtype,
                                  (nf.nc_y_dim_name,
                                   nf.nc_x_dim_name))
            nf._set_attributes(
                nc_out.variables[nf.nc_lat_2d_name],
                nf.nc_lat_2d_attrs)
            
            nc_out.variables[nf.nc_lat_2d_name][:,:] = lat
            
            # Pole mask variable (set data outside if block
            # as it's the same in both cases):
            nc_out.createVariable("pmask", pmask.dtype,
                                 (nf.nc_y_dim_name,
                                  nf.nc_x_dim_name))
        else:
            # Interpolated grid (1D coordinate variables)
            nc_out.createDimension(nf.nc_lon_dim_name, nx)
            nc_out.createDimension(nf.nc_lat_dim_name, ny)
            
            # Longitude:
            nc_out.createVariable(nf.nc_lon_1d_name, lon.dtype,
                                  (nf.nc_lon_dim_name,))
            nf._set_attributes(
                nc_out.variables[nf.nc_lon_1d_name],
                nf.nc_lon_1d_attrs)
            
            nc_out.variables[nf.nc_lon_1d_name][:] = lon
                
            # Latitude:
            nc_out.createVariable(nf.nc_lat_1d_name, lat.dtype,
                                  (nf.nc_lat_dim_name,))
            nf._set_attributes(
                nc_out.variables[nf.nc_lat_1d_name],
                nf.nc_lat_1d_attrs)
            
            nc_out.variables[nf.nc_lat_1d_name][:] = lat
            
            # Pole mask variable (set data outside of if block
            # as it is the same in both cases):
            nc_out.createVariable("pmask", pmask.dtype,
                                 (nf.nc_lat_dim_name,
                                  nf.nc_lon_dim_name))
        
        nf._set_attributes(nc_out.variables["pmask"],
                           nc_pmask_attrs)
        
        nc_out.variables["pmask"][:,:] = pmask
    
    print(f"Saved: {str(nc_out_file)}")


if __name__ == "__main__":
    main()
