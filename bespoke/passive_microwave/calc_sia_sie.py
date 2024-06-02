# Calculate observations of sea ice area and extent from passive
# microwave observations of sea ice concentration.
# ============================================================ #

from argparse import ArgumentParser
from datetime import datetime as dt
from pathlib import Path

import netCDF4 as nc
import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    metadata as md,
    netcdf as nf,
    utils)

from ice_edge_latitude.utilities import regions as iel_regions


# ------------------------------------------------------------ #

# This script needs the north pole-hole masks calculated from
# a separate script. The paths to these pole-hole masks data:
_pmask_path = Path(md.dir_out_nc_data_nsidc, "auxilliary")

nc_file_pole_hole_mask = {
    "SMMR" : {"gn": Path(_pmask_path, "pmask_gn_smmr.nc"),
              "gr": Path(_pmask_path, "pmask_gr_smmr.nc")},
    "SSMI" : {"gn": Path(_pmask_path, "pmask_gn_ssmi.nc"),
              "gr": Path(_pmask_path, "pmask_gr_ssmi.nc")},
    "SSMIS": {"gn": Path(_pmask_path, "pmask_gn_ssmis.nc"),
              "gr": Path(_pmask_path, "pmask_gr_ssmis.nc")}
}

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
# mean, etc., and the corresponding netCDF variable name, are
# added automatically as specified in the nf module):
diag_area_name   = "sia"
diag_extent_name = "sie"

nc_extent_long_name = (
    "{}ern hemisphere sea ice extent (area of "
    f"ocean where siconc >= {nf.nc_siconc_threshold_name})"
    + "{}")

nc_area_long_name = "{}ern hemisphere sea ice area{}"

# Extent has "coordinates" used for a scalar coordinate variable
# corresponding to the threshold concentration for extent:
nc_var_extent_attrs = {
    "standard_name": "sea_ice_extent",
    "cell_methods" : f"{nf.nc_time_name}: mean",
    "coordinates"  : nf.nc_siconc_threshold_name,
    "units"        : nf.field_units["seaiceextent"]
}

nc_var_area_attrs = {
    "standard_name": "sea_ice_area",
    "cell_methods" : f"{nf.nc_time_name}: mean",
    "units"        : nf.field_units["seaicearea"]
}


def load_native_ps_cell_area(hemi="n", nc_var_name="cell_area"):
    """Load cell areas on the native grid. These come directly
    from raw data (NSIDC dataset ID: NSIDC-0771), the file paths
    of which are set in metadata.nsidc_areacell_file.
    """
    with nc.Dataset(md.nsidc_areacell_file[hemi]) as ncdat:
        cell_area = np.array(ncdat.variables[nc_var_name])
    return cell_area


def load_pole_hole_mask(mask="SMMR", grid="gn"):
    """Load north pole-hole mask data for specified mask
    satellite instrument period, mask = "SMMR", "SSMI", or
    "SSMIS", on the native (grid="gn") or regridded (grid="gr")
    grid. These datasets are computed in a separate script.
    """
    with nc.Dataset(nc_file_pole_hole_mask[mask][grid]) as ncdat:
        pmask = np.array(ncdat.variables["pmask"])
    return pmask


def main():
    
    # Arguments required are input netCDF files to calculate
    # from (should be the native output from
    # "prepare_raw_siconc_data" python script or those that have
    # been regridded using the bash script, a string option to
    # identify which dataset corresponds to input data (for data
    # references and filename; NSIDC-0051 or NSIDC-0079), a
    # string to determined which hemisphere is being computed
    # ("n" or "s"), and a string description of the
    # interpolation method/resolution to be added to the output
    # netCDF global attributes, if interpolated data is passed:
    prsr = ArgumentParser()
    prsr.add_argument("-i", "--infiles", type=str, nargs="*",
                      default=[], help="Input file paths")
    prsr.add_argument("--hemisphere", type=str, default="n",
                      choices="ns")
    prsr.add_argument("-d", "--datasetid", type=str,
                      default="NSIDC-0051",
                      choices=["NSIDC-0051", "NSIDC-0079"])
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
    
    nt = len(time_mon)
    yr_s = date[0].year
    yr_e = date[-1].year
    
    # Assume units are fraction (not percentage); this should be
    # dealt with in the prepare_raw_siconc_data.py script:
    siconc = np.maximum(0.0, siconc)
    siconc = np.minimum(1.0, siconc)
    
    # Determine whether input data is on the native grid or has
    # been interpolated based on the number of dimensions of the
    # coordinate arrays:
    if np.ndim(lon) == 1 and np.ndim(lat) == 1:
        # Must be re-gridded (regular lon-lat). Save a string
        # for output netCDF "grid_label" attribute, generate
        # 2D equivalent coordinate arrays for calculation below,
        # and work out grid cell areas (note by default raw data
        # is not interpolated to a global grid):
        lon_bnds, lat_bnds = \
            diags.estimate_lon_lat_bnds_1D(lon, lat)
        areacello = diags.grid_cell_area(lon_bnds, lat_bnds)
        lon, lat = np.meshgrid(lon, lat)
        grid_label = "gr"
        
    elif np.ndim(lon) == 2 and np.ndim(lat) == 2:
        # Data is on the native grid; load corresponding native
        # grid cell areas:
        areacello = load_native_ps_cell_area(cmd.hemisphere)
        grid_label = "gn"
    
    # Additional processing needed for the northern hemisphere:
    if cmd.hemisphere == "n":
        
        # Mask out lakes. Determine the mask (1 if allowed,
        # nan if masked). Again this is copied from the
        # src.load_raw_data prepare_siconc script:
        lon = lon % 360.0
        lake_mask = np.ones(np.shape(areacello),
                            dtype=np.float64)
        
        for reg in ["lakes", "baltic_sea", "black_sea"]:
            
            reg_def = getattr(iel_regions, reg)
            
            for sub_reg in reg_def.keys():
                lake_mask *= np.where(
                    (lon >= reg_def[sub_reg][0]) &
                    (lon <= reg_def[sub_reg][1]) &
                    (lat >= reg_def[sub_reg][2]) &
                    (lat <= reg_def[sub_reg][3]),
                    np.nan, 1.0)
        
        # Apply the mask:
        siconc *= lake_mask[np.newaxis,:,:]
        
        # Deal with pole holes. These vary in size throughout
        # the time period as instruments change. Auxilliary nc
        # files are created with masks for these different
        # periods.
        # 
        # The pole hole area is set to 1 (100% concentration)
        # which is standard practice. For extent this is fine
        # but caution should be taken with sea ice area,
        # particularly in summer:
        for p in pmask_valid_range.keys():
            
            # Time indices corresponding to this mask (if any;
            # by default, the whole time period 1980-2021 is
            # loaded so each mask is required over its entire
            # valid range):
            jt = (date >= pmask_valid_range[p][0]) \
               & (date <= pmask_valid_range[p][1])
            
            if jt.any():
                
                pmask = load_pole_hole_mask(p, grid_label)
                
                # Set pole-hole area to 100% concentration:
                siconc[jt,:,:] = np.where(
                    np.isnan(pmask[np.newaxis,:,:]),
                    1.0, siconc[jt,:,:])
    
    # Calculate monthly sea ice area and extent in 1e6 km2:
    sia_mask = siconc
    sie_mask = np.where(siconc >= cmd.threshold, 1.0, 0.0)
    
    sia_mon = np.nansum(sia_mask*areacello[np.newaxis,:],
                        axis=(1,2))/1.0E12
    
    sie_mon = np.nansum(sie_mask*areacello[np.newaxis,:],
                        axis=(1,2))/1.0E12
    
    # Set December 1987 and January 1988 to missing.    
    # 
    # There is a period of missing daily data from 03 December
    # 1987 to 13 January 1988, so monthly mean extents for both
    # months cannot be computed. For the Bootstrap data, 
    # NSIDC-0079, the raw monthly concentration fields are
    # already set to missing but the above calculation then
    # gives zero for these months. For the NASA Team data,
    # NSIDC-0051, raw monthly concentration fields are not
    # missing (presumably are the average of available days for
    # each month), so the calculation above gives a misleading
    # value:
    jtmiss = (date >= dt(1987,12,1)) & (date <= dt(1988,1,31))
    sia_mon[jtmiss] = md.default_new_missing_value
    sie_mon[jtmiss] = md.default_new_missing_value
    
    # Yearly-mean sea ice area:
    sia_mon_rs = np.reshape(sia_mon, (nt//12, 12))
    sia_yr = np.nanmean(sia_mon_rs, axis=1)
    
    # Yearly-mean sea ice extent:
    sie_mon_rs = np.reshape(sie_mon, (nt//12, 12))
    sie_yr = np.nanmean(sie_mon_rs, axis=1)
    
    
    # ------------------------------------------------------ #
    # NetCDF attributes:
    # ------------------------------------------------------ #
    
    hemi_long = ("Northern" if cmd.hemisphere == "n"
                 else "Southern")
    
    nc_global_attrs = nf.default_nc_file_attrs.copy()
    nf._set_nc_history_attr_timestamp(nc_global_attrs)
    nc_global_attrs["grid_label"] = grid_label
    nc_global_attrs["hemisphere"] = hemi_long
    nc_global_attrs["source_id"] = cmd.datasetid
    nc_global_attrs["source"] = (md.nsidc_source[cmd.datasetid]
                                 + " [2]")
    
    # Different title per output as includes "area" or "extent";
    # format when writing to nc:
    nc_title = hemi_long + " hemisphere sea ice {} {}-mean "
    nc_title += "time series from passive microwave data ("
    nc_title += ("Bootstrap" if "0079" in cmd.datasetid
                 else "NASA Team") + ")"
    # .format({"area", "extent"}, {"monthly", "annual"})
    
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
    
    ref = 3
    
    if cmd.hemisphere == "n":
        nc_global_attrs["comment"] += (
            " Additional National Ice Center (NIC) northern "
            + "hemisphere \'valid ice\' masks, NSIDC dataset "
            + f"ID 0622 [{ref}], applied to raw concentration "
            + "data before further processing.")
        
        nc_global_attrs["references"] += (
            nf.nc_list_delimiter + f"[{ref}] "
            + md.nsidc_data_reference["NSIDC-0622"])
        
        ref += 1
    
    if grid_label == "gn":
        nc_global_attrs["comment"] += (
            " Grid cell areas for computation of diagnostics "
            + f"on the native grid obtained from [{ref}].")
        
        nc_global_attrs["references"] += (
            nf.nc_list_delimiter + f"[{ref}] "
            + md.nsidc_data_reference["NSIDC-0771"])
        
    else:
        nc_global_attrs["comment"] += \
            " " + cmd.interpdescription
    
    if cmd.hemisphere == "n":
        nc_global_attrs["comment"] += (
            " Region of satellite north \'pole hole\', which "
            + "varies in size with instrument source (SMMR, "
            + "SSM/I, and SSMIS), i.e., time period, assumed "
            + "to have a concentration 1 (100%).")
    
    out_dir = md.dir_out_nc_data_nsidc
    
    diag_kw = {"time_methods": nf.diag_nq_monthly}
    nc_var_kw = {"time_methods": nf.nc_var_nq_monthly}
    
    if grid_label == "gn":
        diag_kw["other_methods"] = nf.diag_nq_native
        nc_var_kw["other_methods"] = nf.nc_var_nq_native
    else:
        diag_kw["other_methods"] = nf.diag_nq_regridded
        nc_var_kw["other_methods"] = nf.nc_var_nq_regridded
    
    # Need hemisphere information in the file names; then, this
    # follows an analogue structure for the CMIP6 outputs, with
    # "passive_microwave" taking the spot of <experiment_id> and
    # the NSIDC dataset ID taking the spot of <model_id>:
    diag_nq_bespoke = (f"_{cmd.hemisphere}h_passive_microwave_"
                       + f"{cmd.datasetid}.nc")
    
    # Sea ice area, monthly:
    dir_sia_mon = Path(out_dir, nf.diag_name(diag_area_name,
                                             **diag_kw))
    file_sia_mon = Path(dir_sia_mon,
        nf.diag_name(diag_area_name, **diag_kw)
        + diag_nq_bespoke)
    
    nc_var_name_sia_mon = nf.nc_var_name(diag_area_name,
                                         **nc_var_kw)
    
    # Sea ice extent, monthly:
    dir_sie_mon = Path(out_dir, nf.diag_name(diag_extent_name,
                                             **diag_kw))
    file_sie_mon = Path(dir_sie_mon,
        nf.diag_name(diag_extent_name, **diag_kw)
        + diag_nq_bespoke)
    
    nc_var_name_sie_mon = nf.nc_var_name(diag_extent_name,
                                         **nc_var_kw)
    
    # Sea ice area, yearly:
    diag_kw["time_methods"] = nf.diag_nq_yearly
    nc_var_kw["time_methods"] = nf.nc_var_nq_yearly
    
    dir_sia_yr = Path(out_dir, nf.diag_name(diag_area_name,
                                             **diag_kw))
    file_sia_yr = Path(dir_sia_yr,
        nf.diag_name(diag_area_name, **diag_kw)
        + diag_nq_bespoke)
    
    nc_var_name_sia_yr = nf.nc_var_name(diag_area_name,
                                        **nc_var_kw)
    
    # Sea ice extent, yearly:
    dir_sie_yr = Path(out_dir, nf.diag_name(diag_extent_name,
                                             **diag_kw))
    file_sie_yr = Path(dir_sie_yr,
        nf.diag_name(diag_extent_name, **diag_kw)
        + diag_nq_bespoke)
    
    nc_var_name_sie_yr = nf.nc_var_name(diag_extent_name,
                                        **nc_var_kw)
    
    for x in [dir_sia_mon, dir_sia_yr, dir_sie_mon, dir_sie_yr]:
        x.mkdir(parents=True, exist_ok=True)
    
    with (nc.Dataset(file_sia_mon, "w") as nc_am,
          nc.Dataset(file_sia_yr, "w") as nc_ay,
          nc.Dataset(file_sie_mon, "w") as nc_em,
          nc.Dataset(file_sie_yr, "w") as nc_ey):
        
        # Set global attributes:
        nf._set_attributes(nc_am,
            {"title": nc_title.format("area", "monthly"),
             **nc_global_attrs})
        
        nf._set_attributes(nc_ay,
            {"title": nc_title.format("area", "annual"),
             **nc_global_attrs})
        
        nf._set_attributes(nc_em,
            {"title": nc_title.format("extent", "monthly"),
             **nc_global_attrs})
        
        nf._set_attributes(nc_ey,
            {"title": nc_title.format("extent", "annual"),
             **nc_global_attrs})
        
        # Create/set common variables/attributes:
        for x in [nc_am, nc_ay, nc_em, nc_ey]:
            
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
        
        
        # Create threshold scalar coordinate variable for extent
        # variables:
        for x in [nc_em, nc_ey]:
            x.createVariable(nf.nc_siconc_threshold_name,
                             np.float64)
            x.variables[nf.nc_siconc_threshold_name][:] = \
                cmd.threshold
            nf._set_attributes(
                x.variables[nf.nc_siconc_threshold_name],
                {"standard_name": "sea_ice_area_fraction",
                 "units": "1"})
        
        # Set monthly time data:
        for x in [nc_am, nc_em]:
            x.variables[nf.nc_time_name][:] = time_mon
            x.variables[nc_out_time_attrs["bounds"]][:,:] = \
                time_bnds_mon
        
        # Set yearly time data:
        for x in [nc_ay, nc_ey]:
            nf.set_time_data_yearly(x, (yr_s, yr_e))
        
        # Create sea ice area/extent variables:
        nc_sia_mon = nc_am.createVariable(nc_var_name_sia_mon,
            sia_mon.dtype, (nf.nc_time_dim_name,))
        nc_sia_mon[:] = sia_mon
        nf._set_attributes(nc_sia_mon, nc_var_area_attrs)
        
        nc_sia_yr = nc_ay.createVariable(nc_var_name_sia_yr,
            sia_yr.dtype, (nf.nc_time_dim_name,))
        nc_sia_yr[:] = sia_yr
        nf._set_attributes(nc_sia_yr, nc_var_area_attrs)
        
        nc_sie_mon = nc_em.createVariable(nc_var_name_sie_mon,
            sie_mon.dtype, (nf.nc_time_dim_name,))
        nc_sie_mon[:] = sie_mon
        nf._set_attributes(nc_sie_mon, nc_var_extent_attrs)
        
        nc_sie_yr = nc_ey.createVariable(nc_var_name_sie_yr,
            sie_yr.dtype, (nf.nc_time_dim_name,))
        nc_sie_yr[:] = sie_yr
        nf._set_attributes(nc_sie_yr, nc_var_extent_attrs)
    
    for x in [file_sia_mon, file_sia_yr, file_sie_mon,
              file_sie_yr]:
        print(f"Saved: {x}")


if __name__ == "__main__":
    main()
