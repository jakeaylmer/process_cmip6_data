# Calculate observations of sea ice-edge latitude, as a function
# of longitude and zonal mean, from passive microwave
# observations of sea ice concentration interpolated to a
# regular lon-lat grid.
# ============================================================ #

from argparse import ArgumentParser
from datetime import datetime as dt
from pathlib import Path

import netCDF4 as nc
import numpy as np

from process_cmip6_data.src import (
    metadata as md, netcdf as nf, utils)

from ice_edge_latitude.diagnostics import ice_edge_latitude as iel
from ice_edge_latitude.utilities import regions as iel_regions


# ------------------------------------------------------------ #

# This script needs the north pole-hole masks calculated from
# a separate script. The paths to these pole-hole masks data:
_pmask_path = Path(md.dir_out_nc_data_nsidc, "auxilliary")

nc_file_pole_hole_mask = {
    "SMMR" : Path(_pmask_path, "pmask_gr_smmr.nc"),
    "SSMI" : Path(_pmask_path, "pmask_gr_ssmi.nc"),
    "SSMIS": Path(_pmask_path, "pmask_gr_ssmis.nc")}

# datetime (start, end) for when each mask applies (input
# concentration data is monthly, so these are not precise nor
# correspond to the values in the pmask metadata -- e.g., note
# the overlap of SMMR and SSMI in August 1987 here):
pmask_valid_range = {
    "SMMR" : (dt(1978, 10, 1), dt(1987,  8, 31)),
    "SSMI" : (dt(1987,  8, 1), dt(2008, 1, 1)),
    "SSMIS": (dt(2008, 1, 1), dt(2023, 12, 31))}

# ------------------------------------------------------------ #

# This metadata for netCDF output copied from CMIP6 sea ice
# calculation script:
# 
# Basic diagnostic names (details such as time averaging, zonal
# mean, etc., and the corresponding nf variable name, are
# added automatically as specified in the nf module):
diag_name = "iel"

nc_long_name = (
    "{} hemisphere{} sea ice-edge latitude{}")

# "coordinates" used for a scalar coordinate variable
# corresponding to the threshold concentration for ice edge:
nc_var_attrs = {
    "cell_methods": f"{nf.nc_time_name}: mean "
                    + f"{nf.nc_lon_1d_name}: mean",
    "coordinates": f"{nf.nc_siconc_threshold_name} "
                   + f"{nf.nc_lon_1d_name}",
    "standard_name": "latitude",
    "units": nf.field_units["seaiceedge"]
}

nc_var_zm_attrs = {
    "cell_methods": f"{nf.nc_time_name}: mean",
    "coordinates": nf.nc_siconc_threshold_name,
    "standard_name": "latitude",
    "units"      : nf.field_units["seaiceedge"]
}


def load_pole_hole_mask(mask="SMMR"):
    """Load north pole-hole mask data for specified mask
    satellite instrument period, mask = "SMMR", "SSMI", or
    "SSMIS", regridded only. These datasets are computed in a
    separate script.
    """
    with nc.Dataset(nc_file_pole_hole_mask[mask]) as ncdat:
        pmask = np.array(ncdat.variables["pmask"])
    return pmask


def main():
    
    # Arguments required are input netCDF files to calculate
    # from (should be the regridded output from
    # "prepare_raw_siconc_data" python script, a string option
    # to identify which dataset corresponds to input data (for
    # data references and filename; NSIDC-0051 or NSIDC-0079), a
    # string to determined which hemisphere is being computed
    # ("n" or "s"), and a string description of the
    # interpolation method/resolution to be added to the output
    # netCDF global attributes:
    prsr = ArgumentParser()
    prsr.add_argument("-i", "--infiles", type=str, nargs="*",
                      default=[], help="Input file paths")
    prsr.add_argument("--hemisphere", type=str, default="n",
                      choices="ns")
    prsr.add_argument("-d", "--datasetid", type=str,
                      default="NSIDC-0051")
    prsr.add_argument("--threshold", type=float,
                      default=md.default_sie_threshold)
    prsr.add_argument("--interpdescription", type=str,
                      default="Diagnostics computed after "
                              + "interpolating data to a 0.5 "
                              + "degree regular lon-lat grid "
                              + "using bilinear interpolation.")
    cmd = prsr.parse_args()
    # -------------------------------------------------------- #
    
    # Prepare siconc here as in the src.load_raw_data function
    # prepare_siconc(). Load input data files:
    lon, lat, time_mon, time_bnds_mon, siconc = \
        utils.nc_get_arrays(sorted(list(cmd.infiles)),
        ["lon", "lat"], ["time", "time_bnds", "siconc"])
    
    # Inherit time variable properties/attributes:
    t_units, t_calendar = utils.nc_get_time_units(cmd.infiles[0])
    
    nc_out_time_attrs = {
        "bounds"  : f"{nf.nc_time_name}_bnds",
        "calendar": t_calendar,
        "units"   : t_units}
    
    # Will need the date to determine the time/year range input,
    # set the 'valid ice' and pole-hole masks if processing the
    # northern hemisphere, and to set missing data for the 1987-
    # 1988 missing data period:
    date = nc.num2date(time_mon, units=t_units,
                       calendar=t_calendar)
    
    nt, ny, nx = np.shape(siconc)
    yr_s = date[0].year
    yr_e = date[-1].year
    
    # Need longitude bounds for the saved data (ice edge as a
    # function of longitude; technically, the value at each
    # longitude represents an average ice edge latitude over the
    # longitude 'grid cell' width):
    lon_res = abs(lon[1] - lon[0])
    lon_bnds = np.array([[lon[k] - lon_res/2.0,
                          lon[k] + lon_res/2.0]
                        for k in range(nx)])
    
    # Assume units are fraction (not percentage); this should be
    # dealt with in the prepare_raw_siconc_data.py script:
    siconc = np.maximum(0.0, siconc)
    siconc = np.minimum(1.0, siconc)
    
    # Additional processing needed for the northern hemisphere:
    if cmd.hemisphere == "n":

        # Need 2D coordinate arrays for lake masking (but lon
        # needs to be saved for the ice edge as a function of
        # longitude, so don't overwrite the original 1D arrays):
        lon_mg, lat_mg = np.meshgrid(lon, lat)
        
        # Mask out lakes. Determine the mask (1 if allowed,
        # nan if masked). This is adapted from the
        # src.load_raw_data prepare_siconc script:
        lon_mg = lon_mg % 360.0
        lake_mask = np.ones((ny, nx), dtype=np.float64)
        
        for reg in ["lakes", "baltic_sea", "black_sea"]:
            
            reg_def = getattr(iel_regions, reg)
            
            for sub_reg in reg_def.keys():
                lake_mask *= np.where(
                    (lon_mg >= reg_def[sub_reg][0]) &
                    (lon_mg <= reg_def[sub_reg][1]) &
                    (lat_mg >= reg_def[sub_reg][2]) &
                    (lat_mg <= reg_def[sub_reg][3]),
                    np.nan, 1.0)
        
        # Apply the mask:
        siconc *= lake_mask[np.newaxis,:,:]
        
        # Deal with pole holes. These vary in size throughout
        # the time period as instruments change. Auxilliary nc
        # files are created with masks for these different
        # periods.
        # 
        # The pole hole area is set to 1 (100% concentration)
        # which is standard practice. For ice edge, this likely
        # has minimal impact anyway:
        for p in pmask_valid_range.keys():
            
            # Time indices corresponding to this mask (if any;
            # by default, the whole time period 1980-2021 is
            # loaded so each mask is required over its entire
            # valid range):
            jt = (date >= pmask_valid_range[p][0]) \
               & (date <= pmask_valid_range[p][1])
            
            if jt.any():
                
                pmask = load_pole_hole_mask(p)
                
                # Set pole mask to 100% concentration:
                siconc[jt,:,:] = np.where(
                    np.isnan(pmask[np.newaxis,:,:]),
                    1.0, siconc[jt,:,:])
    
    else:
        # Southern hemisphere doesn't require any extra
        # processing as such, but the ice_edge_latitude code
        # requires positive increasing latitudes, so sort that
        # here (i.e., put in degrees_south and reverse direction
        # of arrays):
        lat = -lat[::-1]
        siconc = siconc[:,::-1,:]
    
    # Options to get_ice_edge function as keyword arguments.
    # 
    # Here, because we load data for one hemisphere only, need
    # to use the "get_ice_edge()" function [rather than the
    # 'api' function "ice_edge()" of the ice_edge_latitude
    # package].
    # 
    # This requires an explicit indication of the values to set
    # for ice-free meridians (not to be confused with longitudes
    # that are set to NaN/missing because of the land constraint
    # effect -- see external package documentation). These are
    # determined here using the 'hidden' _get_max_lat() function
    # which is the default behaviour of "ice_edge()", is also
    # used for the CMIP6 data, and sets such ice-free meridians
    # to the maximum latitude in the input data (i.e., north
    # pole for northern hemisphere and Antarctic coastline for
    # the southern hemisphere). That only works for missing
    # values set to NaN, which is the case here again by
    # default:
    ie_kw = {"threshold": cmd.threshold,
             "ice_free_meridian_values":
                 iel._get_max_lat(lat, siconc[0,:,:])}
    
    # Pre-allocate and explicitly loop over time -- again this
    # is necessary for the aforementioned reasons (input data is
    # one hemisphere at a time) and the get_ice_edge() function
    # takes one time step only:
    iel_mon = np.zeros((nt, nx))
    
    for k in range(nt):
        iel_mon[k,:] = iel.get_ice_edge(lat, siconc[k,:,:],
                                        **ie_kw)
        # Routine is slightly time consuming, so print progress
        # percentage:
        print("Calculating sea ice-edge latitude "
              + f"({100.0*(k+1)/nt:.0f}%)", end="\r")
    print("")  # cancels carriage return for next print()
    
    # The zonal mean ice edge latitude as a function of time:
    iel_zm_mon = np.nanmean(iel_mon, axis=1)
    
    if cmd.hemisphere == "s":
        # Put back in degrees_north
        iel_mon *= -1.0
        iel_zm_mon *= -1.0
    
    # Set December 1987 and January 1988 to missing.    
    # 
    # There is a period of missing daily data from 03 December
    # 1987 to 13 January 1988, so monthly mean extents for both
    # months cannot be computed. For the Bootstrap data, 
    # NSIDC-0079, the raw monthly concentration fields are
    # already set to missing and so the above calculation gives
    # the ice-free meridian value (missing) for all data anyway.
    # For the NASA Team data, NSIDC-0051, raw monthly
    # concentration fields are not missing (presumably are the
    # average of available days for each month), so the
    # calculation above gives a misleading value:
    jtmiss = (date >= dt(1987,12,1)) & (date <= dt(1988,1,31))
    iel_mon[jtmiss,:] = md.default_new_missing_value
    iel_zm_mon[jtmiss] = md.default_new_missing_value
    
    # Annual means:
    iel_yr = np.nanmean(np.reshape(iel_mon, (nt//12, 12, nx)),
                        axis=1)
    iel_zm_yr = np.nanmean(np.reshape(iel_zm_mon, (nt//12, 12)),
                           axis=1)
    
    # ------------------------------------------------------ #
    # NetCDF attributes:
    # ------------------------------------------------------ #
    
    print("Saving to NetCDF...")
    
    hemi_long = ("Northern" if cmd.hemisphere == "n"
                 else "Southern")
    
    nc_global_attrs = nf.default_nc_file_attrs.copy()
    nf._set_nc_history_attr_timestamp(nc_global_attrs)
    nc_global_attrs["grid_label"] = "gr"
    nc_global_attrs["hemisphere"] = hemi_long
    nc_global_attrs["source_id"] = cmd.datasetid
    nc_global_attrs["source"] = (md.nsidc_source[cmd.datasetid]
                                 + " [2]")
    
    # Different title per output as includes "" or
    # "zonal_mean " and "monthly" or "annual-mean"; format when
    # writing to nc:
    nc_title = hemi_long + " hemisphere {}sea ice-edge latitude"
    nc_title += " {}-mean time series from passive microwave "
    nc_title += "data ("
    nc_title += ("Bootstrap" if "0079" in cmd.datasetid
                 else "NASA Team") + ")"
    # .format({"", "zonal-mean "}, {"monthly", "annual"}
    
    nc_global_attrs["comment"] = (
        "Observations of sea ice concentration, diagnostics "
        + "derived from raw data obtained from the National "
        + "Snow and Ice Data Center (NSIDC), for the analysis "
        + "presented in Aylmer et al. 2024 [1]. Diagnostics "
        + "are computed from monthly concentration data "
        + "obtained from the NSIDC dataset ID "
        + f"{cmd.datasetid[-4:]} [2].")
    
    nc_global_attrs["references"] = ("[1] " + nf.paper_reference
        + nf.nc_list_delimiter + "[2] "
        + md.nsidc_data_reference[cmd.datasetid])
    
    if cmd.hemisphere == "n":
        nc_global_attrs["comment"] += (
            " Additional National Ice Center (NIC) northern "
            + "hemisphere \'valid ice\' masks, NSIDC dataset "
            + f"ID 0622 [3], applied to raw concentration data "
            + "before further processing.")
        
        nc_global_attrs["references"] += (
            nf.nc_list_delimiter + f"[3] "
            + md.nsidc_data_reference["NSIDC-0622"])
    
    nc_global_attrs["comment"] += " " + cmd.interpdescription
    
    if cmd.hemisphere == "n":
        nc_global_attrs["comment"] += (
            " Region of satellite north \'pole hole\', which "
            + "varies in size with instrument source (SMMR, "
            + "SSM/I, and SSMIS), i.e., time period, assumed "
            + "to have a concentration 1 (100%).")
    
    out_dir = md.dir_out_nc_data_nsidc
    
    # For simplicity, hard-code the "other_methods":
    diag_kw = {
        "name": diag_name,
        "time_methods": nf.diag_nq_monthly,
        "other_methods": nf.diag_nq_05deg_bil}
    
    nc_var_kw = {
        "name": diag_name,
        "time_methods": nf.nc_var_nq_monthly,
        "other_methods": nf.nc_var_nq_05deg_bil}
    
    # Need hemisphere information in the file names; then, this
    # follows an analogue structure for the CMIP6 outputs, with
    # "passive_microwave" taking the spot of <experiment_id>
    # and the NSIDC dataset ID taking the spot of <model_id>:
    diag_nq_bespoke = (f"_{cmd.hemisphere}h_passive_microwave"
                       + f"_{cmd.datasetid}.nc")
    
    # Ice-edge latitude (time, longitude), monthly mean:
    dir_iel_mon = Path(out_dir, nf.diag_name(**diag_kw))
    file_iel_mon = Path(dir_iel_mon,
        nf.diag_name(**diag_kw) + diag_nq_bespoke)
    nc_var_name_iel_mon = nf.nc_var_name(**nc_var_kw)
    
    # Ice-edge latitude (time, longitude) yearly mean:
    diag_kw["time_methods"] = nf.diag_nq_yearly
    nc_var_kw["time_methods"] = nf.nc_var_nq_yearly
    dir_iel_yr = Path(out_dir, nf.diag_name(**diag_kw))
    file_iel_yr = Path(dir_iel_yr,
        nf.diag_name(**diag_kw) + diag_nq_bespoke)
    nc_var_name_iel_yr = nf.nc_var_name(**nc_var_kw)
    
    # Ice-edge latitude (time), zonal mean, monthly mean:
    diag_kw["time_methods"] = nf.diag_nq_monthly
    nc_var_kw["time_methods"] = nf.nc_var_nq_monthly
    diag_kw["space_methods"] = nf.diag_nq_zonal_mean
    nc_var_kw["space_methods"] = nf.nc_var_nq_zonal_mean
    dir_iel_zm_mon = Path(out_dir, nf.diag_name(**diag_kw))
    file_iel_zm_mon = Path(dir_iel_zm_mon,
        nf.diag_name(**diag_kw) + diag_nq_bespoke)
    nc_var_name_iel_zm_mon = nf.nc_var_name(**nc_var_kw)
    
    # Ice-edge latitude (time), zonal mean, yearly mean:
    diag_kw["time_methods"] = nf.diag_nq_yearly
    nc_var_kw["time_methods"] = nf.nc_var_nq_yearly
    dir_iel_zm_yr = Path(out_dir, nf.diag_name(**diag_kw))
    file_iel_zm_yr = Path(dir_iel_zm_yr,
        nf.diag_name(**diag_kw) + diag_nq_bespoke)
    nc_var_name_iel_zm_yr = nf.nc_var_name(**nc_var_kw)
    
    for x in [dir_iel_mon, dir_iel_yr, dir_iel_zm_mon,
              dir_iel_zm_yr]:
        x.mkdir(parents=True, exist_ok=True)
    
    with (nc.Dataset(file_iel_mon, "w") as nc_m,
          nc.Dataset(file_iel_yr, "w") as nc_y,
          nc.Dataset(file_iel_zm_mon, "w") as nc_zm,
          nc.Dataset(file_iel_zm_yr, "w") as nc_zy):
        
        # Set global attributes:
        nf._set_attributes(nc_m,
            {"title": nc_title.format("", "monthly"),
             **nc_global_attrs})
        
        nf._set_attributes(nc_y,
            {"title": nc_title.format("", "annual"),
             **nc_global_attrs})
        
        nf._set_attributes(nc_zm,
            {"title": nc_title.format("zonal-mean ", "monthly"),
             **nc_global_attrs})
        
        nf._set_attributes(nc_zy,
            {"title": nc_title.format("zonal-mean ", "annual"),
             **nc_global_attrs})
        
        # Create/set common variables/attributes:
        for x in [nc_m, nc_y, nc_zm, nc_zy]:
            
            # Dimensions (time and time_bounds):
            x.createDimension(nf.nc_time_dim_name, None)
            x.createDimension(nf.nc_tbnd_dim_name, 2)
            
            # Time variable:
            x.createVariable(nf.nc_time_name, time_mon.dtype,
                             (nf.nc_time_dim_name,))
            
            # Time bounds variable:
            x.createVariable(nc_out_time_attrs["bounds"],
                             time_bnds_mon.dtype,
                             (nf.nc_time_dim_name,
                              nf.nc_tbnd_dim_name))
            
            # Attributes for time (no attributes for time bounds
            # variable as per CF conventions):
            nf._set_attributes(x.variables[nf.nc_time_name],
                               nc_out_time_attrs)
        
        # Set longitude dimension/variable for non-zonal-means:
        for x in [nc_m, nc_y]:
            x.createDimension(nf.nc_lon_dim_name, nx)
            x.createDimension(nf.nc_vert_dim_name, 2)
            x.createVariable(nf.nc_lon_1d_name, lon.dtype,
                             (nf.nc_lon_dim_name,))
            
            x.createVariable(nf.nc_lon_1d_name + "_bnds",
                             lon_bnds.dtype,
                             (nf.nc_lon_dim_name,
                              nf.nc_vert_dim_name))
            
            nf._set_attributes(x.variables[nf.nc_lon_1d_name],
                {"bounds": nf.nc_lon_1d_name + "_bnds",
                 **nf.nc_lon_1d_attrs})
            
            x.variables[nf.nc_lon_1d_name][:] = lon
            x.variables[nf.nc_lon_1d_name+"_bnds"][:] = lon_bnds
        
        # Set siconc threshold scalar coordinate variable:
        for x in [nc_m, nc_y, nc_zm, nc_zy]:
            x.createVariable(nf.nc_siconc_threshold_name,
                             np.float64)
            x.variables[nf.nc_siconc_threshold_name][:] = \
                cmd.threshold
            nf._set_attributes(
                x.variables[nf.nc_siconc_threshold_name],
                {"standard_name": "sea_ice_area_fraction",
                 "units": "1"})
        
        # Set monthly time data:
        for x in [nc_m, nc_zm]:
            x.variables[nf.nc_time_name][:] = time_mon
            x.variables[nc_out_time_attrs["bounds"]][:,:] = \
                time_bnds_mon
        
        # Set yearly time data:
        for x in [nc_y, nc_zy]:
            nf.set_time_data_yearly(x, (yr_s, yr_e))
        
        # Create ice-edge latitude variables:
        nc_iel_mon = nc_m.createVariable(nc_var_name_iel_mon,
            iel_mon.dtype,
            (nf.nc_time_dim_name, nf.nc_lon_dim_name))
        nc_iel_mon[:,:] = iel_mon
        nf._set_attributes(nc_iel_mon,
            {"long_name": nc_long_name.format(hemi_long,"",""),
             **nc_var_attrs})
        
        nc_iel_yr = nc_y.createVariable(nc_var_name_iel_yr,
            iel_yr.dtype,
            (nf.nc_time_dim_name, nf.nc_lon_dim_name))
        nc_iel_yr[:,:] = iel_yr
        nf._set_attributes(nc_iel_yr,
            {"long_name": nc_long_name.format(
                hemi_long, "", ", annual mean"),
             **nc_var_attrs})
        
        nc_iel_zm_mon = nc_zm.createVariable(
            nc_var_name_iel_zm_mon,
            iel_zm_mon.dtype, (nf.nc_time_dim_name,))
        nc_iel_zm_mon[:] = iel_zm_mon
        nf._set_attributes(nc_iel_zm_mon,
            {"long_name": nc_long_name.format(
                hemi_long, " zonal-mean", ""),
             **nc_var_zm_attrs})
        
        nc_iel_zm_yr = nc_zy.createVariable(
            nc_var_name_iel_zm_yr,
            iel_zm_yr.dtype, (nf.nc_time_dim_name,))
        nc_iel_zm_yr[:] = iel_zm_yr
        nf._set_attributes(nc_iel_zm_yr,
            {"long_name": nc_long_name.format(
                hemi_long, " zonal-mean", ", annual mean"),
             **nc_var_zm_attrs})
    
    for x in [file_iel_mon, file_iel_yr, file_iel_zm_mon,
              file_iel_zm_yr]:
        print(f"Saved: {x}")


if __name__ == "__main__":
    main()
