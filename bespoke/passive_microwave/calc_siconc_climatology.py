# Create fixed annual climatology sea ice concentration maps
# from passive microwave observations.
# ============================================================ #

from argparse import ArgumentParser
from datetime import datetime as dt
from pathlib import Path

import netCDF4 as nc
import numpy as np

from process_cmip6_data.src import (
    metadata as md, netcdf as nf, utils)


def main():
    
    # Arguments required are input netCDF files to determine
    # climatology from (should be the native output from
    # "prepare_raw_siconc_data" python script -- this script
    # does not support, though could easily be adapted to,
    # producing climatologies from remapped data), a string
    # option to identify which dataset corresponds to input
    # data (for data references and filename; NSIDC-0051 or
    # NSIDC-0079), and a string to determined which hemisphere
    # is being computed ("n" or "s"):
    prsr = ArgumentParser()
    prsr.add_argument("-i", "--infiles", type=str, nargs="*",
                      default=[], help="Input file paths")
    prsr.add_argument("--hemisphere", type=str, default="n",
                      choices="ns")
    prsr.add_argument("-d", "--datasetid", type=str,
                      default="NSIDC-0051",
                      choices=["NSIDC-0051", "NSIDC-0079"])
    cmd = prsr.parse_args()
    # -------------------------------------------------------- #
    
    # Load data (need coordinates and siconc):
    lon, lat, time_bnds, siconc = utils.nc_get_arrays(
        sorted(list(cmd.infiles)), ["lon", "lat"],
        ["time_bnds", "siconc"])
    
    # Get datetime bounds from time_bnds to save the start and
    # end points of the climatology as global netCDF attributes:
    t_units, t_calendar = utils.nc_get_time_units(cmd.infiles[0])
    date_bnds = nc.num2date(time_bnds, units=t_units,
                            calendar=t_calendar)
    
    ny, nx = np.shape(siconc)[1:]
    
    clim_year_start = date_bnds[0,0].year
    clim_year_end = date_bnds[-1,0].year
    
    # Set December 1987 and January 1988 to missing.    
    # 
    # There is a period of missing daily data from 03 December
    # 1987 to 13 January 1988. For the Bootstrap data, 
    # NSIDC-0079, the raw monthly concentration fields are
    # already set to missing, while for the NASA Team data,
    # NSIDC-0051, raw monthly concentration fields are not
    # missing (presumably are the average of available days for
    # each month). Either way, identify if and which time
    # indices correspond to these months and set to missing:
    jtmiss = (date_bnds[:,0] >= dt(1987,12,1)) & (
              date_bnds[:,1] <= dt(1988,2,1))
    siconc[jtmiss,:,:] = md.default_new_missing_value
    
    # Compute climatology (average over all time):
    siconc_clim = np.nanmean(siconc, axis=0)
    
    
    # ------------------------------------------------------ #
    # NetCDF attributes:
    # ------------------------------------------------------ #
    
    diag_kw = {"time_methods": "climatology",
               "space_methods": nf.diag_nq_native}
    nc_var_kw = {"space_methods": nf.nc_var_nq_native}
    
    # Need hemisphere and time-range information in the file
    # names; then, this follows an analogue structure for the
    # CMIP6 outputs, with "passive_microwave" taking the spot of
    # <experiment_id> and the NSIDC dataset ID taking the spot
    # of <model_id>:
    diag_nq_bespoke = (f"_{cmd.hemisphere}h_{clim_year_start}-"
                       + f"{clim_year_end}_passive_microwave_"
                       + f"{cmd.datasetid}.nc")
    
    nc_out_dir = Path(md.dir_out_nc_data_nsidc,
                   nf.diag_name("siconc", **diag_kw))
    nc_out_dir.mkdir(parents=True, exist_ok=True)
    
    nc_out_file = Path(nc_out_dir,
        nf.diag_name("siconc", **diag_kw) + diag_nq_bespoke)
    
    nc_var_name = nf.nc_var_name("siconc", **nc_var_kw)
    
    hemi_long = ("Northern" if cmd.hemisphere == "n" else
                 "Southern")
    
    nc_global_attrs = nf.default_nc_file_attrs.copy()
    nf._set_nc_history_attr_timestamp(nc_global_attrs)
    nc_global_attrs["grid_label"] = "gn"
    nc_global_attrs["time_climatology_start"] = \
        date_bnds[0,0].strftime("%Y-%m-%d %H:%M:%S")
    nc_global_attrs["time_climatology_end"] = \
        date_bnds[-1,1].strftime("%Y-%m-%d %H:%M:%S")
    nc_global_attrs["hemisphere"] = hemi_long
    nc_global_attrs["source_id"] = cmd.datasetid
    nc_global_attrs["source"] = (md.nsidc_source[cmd.datasetid]
                                 + " [2]")
    
    nc_global_attrs["title"] = (hemi_long + " hemisphere sea "
        + "ice concentration annual-mean climatology from "
        + "passive microwave data ("
        + ("Bootstrap" if "0079" in cmd.datasetid
           else "NASA Team")
        + f"), {clim_year_start}-{clim_year_end}")
    
    nc_global_attrs["comment"] = (
        "Observations of sea ice concentration, diagnostics "
        + "derived from raw data obtained from the National "
        + "Snow and Ice Data Center (NSIDC), for the analysis "
        + "presented in Aylmer et al. 2024 [1]. This dataset "
        + f"contains {hemi_long.lower()} hemisphere sea ice "
        + "concentration, annual-mean climatology for years "
        + f"{clim_year_start}-{clim_year_end} inclusive, "
        + "obtained from monthly data from the NSIDC dataset "
        + f"ID {cmd.datasetid[-4:]} [2]. Grid coordinates "
        + "obtained from [3].")
    
    nc_global_attrs["references"] = ("[1] " + nf.paper_reference
        + nf.nc_list_delimiter + "[2] "
        + md.nsidc_data_reference[cmd.datasetid]
        + nf.nc_list_delimiter + "[3] "
        + md.nsidc_data_reference["NSIDC-0771"])
    
    if cmd.hemisphere == "n":
        nc_global_attrs["comment"] += (
            " Additional National Ice Center (NIC) northern "
            + "hemisphere \'valid ice\' mask, NSIDC dataset ID "
            + "0622 [4], has been applied to the this data "
            + "before computation of the climatology.")
        
        nc_global_attrs["references"] += (nf.nc_list_delimiter
            + "[4] " + md.nsidc_data_reference["NSIDC-0622"])
    
    nc_out_siconc_attrs = {
        "coordinates": f"{nf.nc_lon_2d_name} {nf.nc_lat_2d_name}",
        "long_name": "sea ice concentration",
        "standard_name": "sea_ice_area_fraction",
        "units": "1"}
    
    with nc.Dataset(nc_out_file, "w") as nc_out:
        
        nf._set_attributes(nc_out, nc_global_attrs)
        
        nc_out.createDimension(nf.nc_x_dim_name, nx)
        nc_out.createDimension(nf.nc_y_dim_name, ny)
        
        # Longitude:
        nc_out.createVariable(nf.nc_lon_2d_name, lon.dtype,
                              (nf.nc_y_dim_name,
                               nf.nc_x_dim_name))
        # Latitude:
        nc_out.createVariable(nf.nc_lat_2d_name, lat.dtype,
                              (nf.nc_y_dim_name,
                               nf.nc_x_dim_name))
        # Sea ice concentration:
        nc_out.createVariable(nc_var_name, siconc_clim.dtype,
                              (nf.nc_y_dim_name,
                               nf.nc_x_dim_name),
                              zlib=True)
        
        nf._set_attributes(nc_out.variables[nf.nc_lon_2d_name],
                           nf.nc_lon_2d_attrs)
        
        nf._set_attributes(nc_out.variables[nf.nc_lat_2d_name],
                           nf.nc_lat_2d_attrs)
        
        nf._set_attributes(nc_out.variables[nc_var_name],
                           nc_out_siconc_attrs)
        
        nc_out.variables[nf.nc_lon_2d_name][:,:] = lon
        nc_out.variables[nf.nc_lat_2d_name][:,:] = lat
        nc_out.variables[nc_var_name][:,:] = siconc_clim
    
    print(f"Saved: {str(Path(nc_out_file))}")


if __name__ == "__main__":
    main()
