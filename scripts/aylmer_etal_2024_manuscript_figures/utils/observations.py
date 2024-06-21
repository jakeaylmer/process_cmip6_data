"""Load and process observation data for plots."""

from pathlib import Path

import netCDF4 as nc
import numpy as np
from tabulate import tabulate

from .ebm import area_factor as calc_area_factor
from .script_tools import get_path


# ------------------------------------------------------------ #
# Metadata for loading passive microwave derived observations of
# sea ice edge change
# ------------------------------------------------------------ #
defined_passive_microwave_datasets = ["NSIDC-0051","NSIDC-0079"]
_path_passive_microwave = \
    get_path("path_nsidc_processed_data.txt")

# ------------------------------------------------------------ #


# ------------------------------------------------------------ #
# Metadata for loading reanalysis data for surface air
# temperature area averages:
# ------------------------------------------------------------ #
defined_reanalyses = ["CFSR", "CFSv2", "CFSR/CFSv2", "ERA5",
                      "JRA-55", "MERRA-2"]

_path_reanalyses = \
    get_path("path_atmospheric_reanalyses_processed_data.txt")

# --- directories and file formats depend on diagnostic type
#     (keys of the following two dictionaries):
reanalysis_tas_diagnostic = {
    "gn"       : "tas_area_mean_yr_gn",
    "gn_interp": "tas_area_mean_yr_gn_interp",
    "cc_approx": "tas_area_mean_yr_cc_approx"}

reanalysis_tas_file_fmt = {
    "gn"       : "tas_area_mean_yr_gn_reanalysis_{}.nc",
    "gn_interp": "tas_area_mean_yr_gn_interp_reanalysis_{}.nc",
    "cc_approx": "tas_area_mean_yr_cc_approx_reanalysis_{}.nc"}
# ------------------------------------------------------------ #


# ------------------------------------------------------------ #
# Metadata for loading ECCO v4r5 OHT data
# ------------------------------------------------------------ #
# Note this data is raw/unprocessed, obtained directly from
# doi:10.5281/zenodo.7869067
# 
ecco_label = "ECCO v4r5 (RC2)"

_path_ecco_oht = get_path("path_ecco_oht.txt")

# The path above points to "MHT.jld2", which is actually a
# netCDF file despite the extention. It contains only the OHT,
# no time or coordinates:
# 
# >>> ncdump -h MHT.jld2
# netcdf MHT {
# dimensions:
#       phony_dim_0 = 336 ;
#       phony_dim_1 = 179 ;
# variables:
#       double single_stored_object(phony_dim_0, phony_dim_1) ;
# }
# 
# The above is 12*28 = 336 monthly data, starting at 1992 and
# ending in 2019 (= 28 years) and 179 latitudes as in the
# originally published data (doi:10.7910/DVN/AVVGYX) which are
# -89 to 89 N in steps of 1 degree. Units of data are PW.
yr_range_ecco        = (1992, 2019)
ecco_latitude        = np.arange(-89.0, 90.0, 1.0)
nc_var_name_ecco_oht = "single_stored_object"
ecco_oht_unit_factor = 1000.0

# ------------------------------------------------------------ #


def _get_hemi(ref_lats):
    """Determine the hemisphere ("n" or "s") that is being
    evaluated based on input reference latitudes.
    """
    return ("s" if (ref_lats[0] < 0.0 and ref_lats[1] < 0.0)
            else "n")



def get_passive_microwave_delta_iel(hemi="n",
        time_period_1=(1980, 2000),
        time_period_2=(2001, 2021), verbose=True):
    """Load passive microwave sea ice-edge latitude data and
    calculate the change between two time period averages.
    
    The input data filenames are hardcoded, except for the top-
    level directory containing all passive-microwave related
    diagnostics, assuming they have been saved in the format as
    specified in the processing scripts. It is assumed that the
    datasets exist.
    
    
    Optional parameters
    -------------------
    hemi : str "n" (default) or "s"
        Which hemisphere to analyse.
    
    time_period_1, time_period_2 : 2-tuple of int,
        defaults: (1980, 2000) and (2001, 2021) respectively.
        Year ranges to average, inclusive.
    
    
    Returns
    -------
    delta_iel : 2-tuple of float
        Lower and upper bounds on the estimate of the zonal-mean
        sea-ice edge latitude change between the two time-
        averaging periods. These are defined to be the minimum
        and maximum possible changes that can be computed by
        combining the two input datasets NSIDC-0051 and NSIDC-
        0079 (different algorithms for converting the satellite
        brightness temperature into sea ice concentration,
        "NASA Team" and "Bootstrap", respectively). Units are
        degrees_north if hemi == "n", otherwise, degrees_south.
    
    """
    
    # Load datasets for the zonal-mean, annual mean sea ice-
    # edge latitude, assuming they are saved in the format and
    # directory structure as specified in the processing code:
    nc_data_0051 = Path(_path_passive_microwave,
                        "iel_zm_yr_05deg_bil",
                        f"iel_zm_yr_05deg_bil_{hemi}h_passive_"
                        + f"microwave_NSIDC-0051.nc")
    
    nc_data_0079 = Path(_path_passive_microwave,
                        "iel_zm_yr_05deg_bil",
                        f"iel_zm_yr_05deg_bil_{hemi}h_passive_"
                        + f"microwave_NSIDC-0079.nc")
    
    if verbose:
        print(f"Loading: {str(nc_data_0051)}")
    
    with nc.Dataset(nc_data_0051, "r") as ncdat:
        
        t_units = ncdat.variables["time"].units
        t_calendar = ncdat.variables["time"].calendar
        date = nc.num2date(ncdat.variables["time"],
                           units=t_units, calendar=t_calendar)
        
        yr = np.array([x.year for x in date])
        
        iel_0051 = np.array(ncdat.variables["iel_zm"])
    
    if verbose:
        print(f"Loading: {str(nc_data_0079)}")
    
    with nc.Dataset(nc_data_0079, "r") as ncdat:
        iel_0079 = np.array(ncdat.variables["iel_zm"])
        # assume time/units etc. are as in nc_data_0051
    
    if hemi == "s":  # convert to degrees_south
        iel_0051 *= -1.0
        iel_0079 *= -1.0
    
    # Indices to average in time:
    jt_1 = (yr >= time_period_1[0]) & (yr <= time_period_1[1])
    jt_2 = (yr >= time_period_2[0]) & (yr <= time_period_2[1])
    
    # Averages over the two periods, for each dataset:
    iel_0051_1 = np.nanmean(iel_0051[jt_1])
    iel_0051_2 = np.nanmean(iel_0051[jt_2])
    
    iel_0079_1 = np.nanmean(iel_0079[jt_1])
    iel_0079_2 = np.nanmean(iel_0079[jt_2])
    
    # Compute delta_iel for each combination of the above.
    # 
    # Note that, with defaults, d_iel are very similar (within
    # a few percent) for each dataset considered independently
    # (i.e., the first and third entries in the list below),
    # giving no meaningful information about uncertainty because
    # the two time series are roughly offset by a constant value
    # but have similar trends. Hence, the minimum and maximum
    # that are reported as the range do come from the 'mixed'
    # differences (i.e., the second and fourth entries in the
    # list below):
    d_iel = [iel_0051_2 - iel_0051_1,
             iel_0051_2 - iel_0079_1,
             iel_0079_2 - iel_0079_1,
             iel_0079_2 - iel_0051_1]
    
    delta_iel = (min(d_iel), max(d_iel))
    
    if verbose:
        print(f"delta_iel_obs = [{delta_iel[0]:.3f}, "
              + f"{delta_iel[1]:.3f}] degrees_"
              + f"{'nor' if hemi=='n' else 'sou'}th "
              + f"(range = {delta_iel[1]-delta_iel[0]:.3f} "
              + f"degrees_{'nor' if hemi=='n' else 'sou'}th)")
    
    return delta_iel



def get_passive_microwave_siconc_climatology(hemi="n",
    time_period=(1980, 2000), dataset_id="NSIDC-0079"):
    """Load passive microwave sea ice concentration climatology
    data for a specified hemisphere (hemi: "n" or "s"), time-
    averaging period [time_period: 2-tuple, default
    (1980, 2000)], and dataset (dataset_id: str, default =
    "NSIDC-0079").
    
    The input data filenames are hardcoded, except for the top-
    level directory containing all passive-microwave related
    diagnostics, assuming they have been saved in the format as
    specified in the processing scripts. It is assumed that
    the datasets exist.
    
    
    Returns
    -------
    lon, lat, siconc
        Arrays of shape (ny, nx) corresponding to the longitude
        coordinates, latitude coordinates, and sea ice
        concentration data, respectively.
    
    """
    
    # Load data for the multi-year, annual climatology of sea
    # ice concentration, assuming they are saved in the format
    # and directory structure as specified in the processing code:
    nc_data = Path(_path_passive_microwave,
                   "siconc_gn_climatology",
                   f"siconc_gn_climatology_{hemi}h_"
                   + f"{time_period[0]}-{time_period[1]}"
                   + f"_passive_microwave_{dataset_id}.nc")
    
    with nc.Dataset(nc_data, "r") as ncdat:
        siconc = np.array(ncdat.variables["siconc"])
        lon = np.array(ncdat.variables["lon"])
        lat = np.array(ncdat.variables["lat"])
    
    return lon, lat, siconc



def get_reanalysis_delta_tas(reanalyses, ref_lats=(65.0, 90.0),
        global_mean=False, time_period_1=(1980, 2000),
        time_period_2=(2001, 2021), diagnostic_type="gn_interp",
        verbosity=1):
    """Load data for one or more reanalyses and compute their
    change(s) in near-surface air temperature, delta_tas,
    between specified time average periods, where tas is
    averaged between specified reference latitudes or is the
    global mean. Uses the saved polar-cap area average data
    to determine these. It is assumed that the relevant datasets
    exist.
    
    
    Parameters
    ----------
    reanalyses: str or list of str
        The names of the reanalysis or reanalyses to load data
        and compute delta_tas for.
    
    
    Optional parameters
    -------------------
    ref_lats: 2-tuple of float, default = (65.0, 90.0)
        Reference latitudes in degrees_north for computing the
        average of tas between. Unless "diagnostic_type" == "gn"
        these are restricted to 50.0-90.0 in magnitude.
    
    global_mean: bool, default = False
        If True, ignores ref_lats and diagnostic_type and
        computes the change in global mean tas.
    
    time_period_1, time_period_2 : 2-tuple of int,
        defaults: (1980, 2000) and (2001, 2021) respectively.
        Year ranges to average, inclusive.
    
    diagnostic_type: str, default = "gn_interp"
        Refers to the method of calculating tas polar-cap area
        averages which are loaded to compute these diagnostics
        (see data processing scripts documentation).
    
    verbosity: int, default = 1
        How much information to print to the console, with 0
        being no output, 1 being key steps in the routine, and 2
        being full details.
    
    
    Returns
    -------
    delta_tas: float or array of float
        Change in tas with specified parameters for each
        reanalysis corresponding to input (if one provided as a
        single string, returns a single value, otherwise returns
        an array of floats corresponding to each reanalysis).
    
    """
    
    if global_mean:
        ref_lats = (-100.0, 100.0)
        # ^ to ensure getting the lowest possible ref_lat_n;
        # second tuple value is arbitrary as it is unused.
        hemi = "n"
        diagnostic_type = "gn"  # necessary!
    else:
        # Determine which hemisphere to analyse from the input
        # reference latitudes:
        hemi = _get_hemi(ref_lats)
    
    if type(reanalyses) == str:
        reanalyses = [reanalyses]
        unlist = True  # flag to convert back to float
    else:
        unlist = False
    
    n_reanalyses = len(reanalyses)
    
    diagnostic = reanalysis_tas_diagnostic[diagnostic_type]
    file_fmt = reanalysis_tas_file_fmt[diagnostic_type]
    
    if verbosity > 0:
        print(f"Loading {diagnostic} data for reanalys"
              + ("e" if n_reanalyses > 1 else "i") + "s: "
              + ", ".join(reanalyses))
    
    delta_tas = np.zeros(n_reanalyses)
    
    if not global_mean and ref_lats[1] < 85.0:
        # Work out weights for converting loaded polar-cap
        # average data to averages between reference latitudes.
        # 
        # First, data are converted to area integrals by
        # multiplying by the total area of the corresponding
        # polar caps, which are
        # 
        #     2*pi*(a**2)*[1 - sin(phi)]
        # 
        # for reference latitude phi and Earth radius a. Then,
        # the difference is taken and divided by the area
        # between reference latitudes to get the new average.
        # This amounts to a weighted difference (sum) with
        # weighting factors given by the [1 - sin(phi)] terms.
        # 
        # Work out the weights, also absorbing the
        # minus from differencing into the second weight w1:
        w0 = 1.0 - np.sin(abs(ref_lats[0])*np.pi/180.0)
        w1 = np.sin(abs(ref_lats[1])*np.pi/180.0) - 1.0
        w = w0 + w1
    else:
        # Shouldn't be needed, but define these just in case
        # of any missed cases...
        w0 = 1.0
        w1 = 0.0
        w = 1.0
    
    def _get_data_time_only(nc_file_i):
        """Loads the data from an netCDF file and sub-selects
        the correct latitude, converting from polar-cap to
        latitude-band average if required, returning data with
        only the time dimension. See main loop for usage.
        """
        if verbosity > 1:
            print(f"Loading: {nc_file_i}")
        
        with nc.Dataset(nc_file_i, "r") as ncdat_i:
            
            date_i = nc.num2date(ncdat_i.variables["time"],
                units=ncdat_i.variables["time"].units,
                calendar=ncdat_i.variables["time"].calendar)
            
            yr_i = np.array([x.year for x in date_i])
            
            ref_lat_i = np.array(
                ncdat_i.variables[f"ref_lat_{hemi}"])
            
            tas_i = np.array(ncdat_i.variables[f"tas_{hemi}"])
            
            # remove length 1 ensemble member dimension, which
            # was only added for consistency with the CMIP data
            # processing/saving routines:
            tas_i = tas_i[:,0,:]  # now shape (time, ref_lat)
            
            # Index corresponding to specified first
            # (equatorward) reference latitude (or, the latitude
            # closest to it):
            jl0_i = np.argmin(abs(ref_lat_i - ref_lats[0]))
            
            if not global_mean and abs(ref_lats[1]) <= 85.0:
                
                # Convert to latitude-band average (rather than
                # polar cap or global mean).
                # 
                # Only do this if we are *not* computing global
                # mean *and* the second reference latitude is
                # not the pole (more precisely, 85.0, because
                # data was never saved more poleward than that).
                
                # Index corresponding to second (more poleward)
                # reference latitude (or, the closest to it):
                jl1_i = np.argmin(abs(ref_lat_i - ref_lats[1]))
                
                # Use weights computed above (w0, w1) to get
                # tas averaged in the ref_lats[0] to ref_lats[1]
                # latitude band:
                tas_i = (w0*tas_i[:,jl0_i]+w1*tas_i[:,jl1_i])/w
            
            else:
                # Otherwise, we just want the polar-cap averages
                # or global mean as loaded:
                tas_i = tas_i[:,jl0_i]
        
        # [ End with nc.Dataset(...) ]
        return yr_i, tas_i
    
    
    for r in range(n_reanalyses):
        
        if reanalyses[r] in ["CFSR", "CFSv2", "CFSR/CFSv2"]:
            
            # Need to combine CFSR and CFSv2 data into one
            # time series [this is why the _get_data_time_only
            # function evaluates the latitude first; if
            # diagnostic_type is gn (native), then the latitude
            # dimensions are different for CFSR and CFSv2]:
            
            yr = []
            tas = []
            
            for rx in ["CFSR", "CFSv2"]:
                yr_rx, tas_rx = _get_data_time_only(
                    Path(_path_reanalyses, diagnostic,
                         file_fmt.format(rx)))
                
                yr.append(yr_rx.copy())
                tas.append(tas_rx.copy())
            
            yr = np.concatenate(yr)
            tas = np.concatenate(tas)
            
        else:
            # Otherwise, data is one file:
            yr, tas = _get_data_time_only(
                Path(_path_reanalyses, diagnostic,
                     file_fmt.format(reanalyses[r])))
        
        # [ End if reanalyses[r] in ["CFSR", ... ] ]
        
        # Time-axis indices for each time-average period:
        j1 = (yr >= time_period_1[0]) & (yr <= time_period_1[1])
        j2 = (yr >= time_period_2[0]) & (yr <= time_period_2[1])
        
        delta_tas[r] = np.nanmean(tas[j2]) - np.nanmean(tas[j1])
    
    # [ End for r in range(n_reanalyses) ]
    
    if verbosity > 0:
        print(tabulate([[f"{x:.2f} K" for x in delta_tas]],
              headers=reanalyses, tablefmt="plain"))
        print(f"(mean = {np.mean(delta_tas):.3f} K, range = "
              + f"{np.max(delta_tas)-np.min(delta_tas):.3f} K)")
    
    if unlist:
        delta_tas = delta_tas[0]
    
    return delta_tas



def get_ecco_delta_oht(ref_lats=(65.0, 90.0),
        time_period_1=(1980, 2000),
        time_period_2=(2001, 2021),
        convergence=False,
        n_std_clim_extrapolate=2,
        return_full_data=False,
        verbosity=1):
    """Estimate the change in ocean heat transport (convergence)
    OHT(C) from ECCO data [1]. For the default time periods,
    this involves a minimal-assumption extrapolation
    (persistence of climatology) of the data for the first time-
    average period; see Methods section of [2] for details.
    
    
    Optional parameters
    -------------------
    ref_lats: 2-tuple of float, default = (65.0, 90.0)
        Reference latitudes in degrees_north for computing the
        change in OHTC between. If ref_lats[1] >= 85.0, then the
        change in poleward OHT at ref_lats[0] is computed
        instead.
    
    time_period_1, time_period_2 : 2-tuple of int,
        defaults: (1980, 2000) and (2001, 2021) respectively.
        Year ranges to average, inclusive.
    
    convergence: bool, default = False
        Only used if ref_lats[1] < 85.0. If False (default),
        computes the change in difference of poleward OHTs at
        each reference latitude in (T)W. If True, computes an
        approximate OHT convergence change between reference
        latitudes in W m-2. It is approximate because it only
        divides the former calculation by the area of the globe
        between reference latitudes, not the area of the ocean,
        and is intended just to convert to a more intuitive
        scale.
    
    n_std_clim_extrapolate : int, default = 2
        The extrapolation of M missing years is done by assuming
        they exhibit a climatology that matches the nearest
        available M years of data. That climatology is defined
        based on a mean plus or minus n_std_clim_extrapolate
        standard deviations.
    
    return_full_data : bool, default = False
        If False (default), returns the delta_oht central
        estimate and uncertainty range interval estimate. If
        True, returns additional intermediate diagnostics,
        including the full time series of raw data.
    
    verbosity : int, default = 1
        How much information to print to the console, with 0
        being no output, 1 being key steps in the routine, and 2
        being full details.
    
    
    Returns
    -------
    If return_full_data == False (default), returns:
        
        doht_central : float
            The central estimate of change in OHT(C)
        
        doht_err : array of float
            The surrounding uncertainty interval arising from
            extrapolation of data assuming persistence of
            climatology.
    
    If return_full_data == True, instead returns the following:
        
        yr_actual : array of int
            The years of the full raw, annually averaged, data.
        
        oht : array of float
            The full raw, annually averaged, OHT(C) data.
        
        j_partial_tp1 : length-2 list of int
            The time indices of yr_actual corresponding to the
            partial average of time period 1.
        
        oht_partial_tp1 : float
            The partial average of OHT(C) for time period 1
            
        j_partial_tp2, oht_partial_tp2
            As above for the second time average period.
        
        yr_miss : 2-tuple of int
            Start and end year of the period that has been
            extrapolated using the climatology of years
            yr_extrap_clim.
        
        yr_extrap_clim : 2-tuple of int
            Start and end year of the period that has been used
            to define a climatology for extrapolation of missing
            years (yr_miss).
        
        oht_miss_mean : float
            Mean of raw OHT(C) data over the climatology years
            yr_extrap_clim, assumed to also be that over the
            missing years yr_miss.
        
        oht_miss_sd : float
            As above, but the standard deviation.
        
        oht_extrap_range_avg1 : 2-tuple of float
            The minimum and maximum mean OH(C) values estimated
            for the first time averaging period (time_period_1). 
    
    
    References
    ----------
    [1] Forget, G., 2023: ECCO v4 standard analysis sample
        (v4r5-rc2) [Data set], Zenodo,
        doi:10.5281/zenodo.7869067
    
    [2] Aylmer, J. R., D. Ferreira, and D. L. Feltham, 2024:
        Impact of ocean heat transport on sea ice captured by a
        simple energy balance model, Commun. Earth. Environ.,
        accepted in principle May 2024.
    
    """
    
    n_yr_avg_1 = time_period_1[1] - time_period_1[0] + 1
    n_yr_avg_2 = time_period_2[1] - time_period_2[0] + 1
    
    if verbosity > 0:
        print(f"Loading {str(_path_ecco_oht)}")
    
    # Load raw OHT data (full time series):
    with nc.Dataset(_path_ecco_oht, "r") as ncdat:
        oht = ecco_oht_unit_factor*np.array(
            ncdat.variables[nc_var_name_ecco_oht])
    
    # Determine where to evaluate for specified latitudes:
    jl0 = np.argmin(abs(ecco_latitude - ref_lats[0]))
    jl1 = np.argmin(abs(ecco_latitude - ref_lats[1]))
    
    # Need to know which hemisphere we are evaluating (in case
    # we need to compute a difference):
    hemi = _get_hemi(ref_lats)
    lat_units = f"degrees_{'nor' if hemi == 'n' else 'sou'}th"
    area_factor = calc_area_factor(ref_lats)  # for OHTC
    oht_units = "W m-2" if convergence else "TW"
    oht_factor = 1.0E12 if convergence else 1.0
    
    if abs(ref_lats[1]) < 85.0:
        
        if convergence:
            
            if verbosity > 0:
                print(f"Evaluating {ecco_label} OHTC over "
                      + "latitudes "
                      + f"{abs(ecco_latitude[jl0]):.1f}-"
                      + f"{abs(ecco_latitude[jl1]):.1f} "
                      + lat_units)
            
            # If hemi == "s", it's the same because we want
            # poleward OHT so the minus from different
            # order of latitudes is reversed by switching
            # northward OHT to southward OHT:
            oht = (oht[:,jl1] - oht[:,jl0]) / area_factor
        
        else:  # not convergence
            
            if verbosity > 0:
                print(f"Evaluating {ecco_label} poleward "
                      + f"OHT at {abs(ecco_latitude[jl1]):.1f} "
                      + f"{lat_units} minus that at "
                      + f"{abs(ecco_latitude[jl0]):.1f} "
                      + lat_units)
            
            # If hemi == "s", it's the same because we want
            # poleward OHT so the minus from different
            # order of latitudes is reversed by switching
            # northward OHT to southward OHT:
            oht = oht[:,jl1] - oht[:,jl0]
    
    else:
        
        # Just evaluate at lower (more equatorward) reference
        # latitude:
        oht = oht[:,jl0]
        
        if hemi == "s":  # make poleward
            oht *= -1.0
        
        if convergence:
            oht /= area_factor
            
            if verbosity > 0:
                print(f"Evaluating {ecco_label} OHTC over"
                      + f"{abs(ecco_latitude[jl0]):.1f} "
                      + f"{lat_units} to the "
                      + f"{'nor' if hemi == 'n' else 'sou'}th "
                      + "pole")
        else:
            if verbosity > 0:
                print(f"Evaluating {ecco_label} poleward "
                      + f"OHT at {abs(ecco_latitude[jl0]):.1f} "
                      + lat_units)
    
    # Now determine time averages and difference. First, the
    # time (year) axis for the actual data:
    yr_actual = np.arange(yr_range_ecco[0],
                            yr_range_ecco[1]+1, 1)
    n_yr_actual = len(yr_actual)
    
    if len(oht) != n_yr_actual*12:
        raise Exception(f"{ecco_label} data does not match "
                        + f"specified time span "
                        + f"{yr_range_ecco[0]}-"
                        + f"{yr_range_ecco[1]} (data is "
                        + f"actually {len(oht)} months ("
                        + f"{len(oht)//12} years) long).")
    
    # Second, the time (year) axis for the period that we would
    # like (i.e., spanning the specified time averaging periods;
    # with defaults, this would be 1980..2021 inclusive):
    yr_required = np.arange(time_period_1[0],time_period_2[1],1)
    
    # Convert monthly data to annual means:
    oht = np.nanmean(np.reshape(oht, (n_yr_actual, 12)), axis=1)
    
    # Indices of *required* years matching input time average
    # periods (with defaults, this would be the indices
    # corresponding to 1980, 2000, 2001, and 2021 from the
    # array yr_required which spans 1980-2021 inclusive):
    j_rq_t1_0 = np.argmin(abs(yr_required - time_period_1[0]))
    j_rq_t1_1 = np.argmin(abs(yr_required - time_period_1[1]))
    j_rq_t2_0 = np.argmin(abs(yr_required - time_period_2[0]))
    j_rq_t2_1 = np.argmin(abs(yr_required - time_period_2[1]))
    
    # Indices of *actual* years nearest to the input time
    # average periods (with defaults, and assuming ECCO data is
    # 1992-2019, this would be the indices corresponding to
    # 1992, 2000, 2001, and 2019 from the array yr_actual which
    # spans 1992-2019 inclusive):
    j_t1_0 = np.argmin(abs(yr_actual - time_period_1[0]))
    j_t1_1 = np.argmin(abs(yr_actual - time_period_1[1]))
    j_t2_0 = np.argmin(abs(yr_actual - time_period_2[0]))
    j_t2_1 = np.argmin(abs(yr_actual - time_period_2[1]))
    
    # Determine the averages from the data matching the time
    # periods requested as closely as possible (i.e., the
    # 'partial averages', without any extrapolation). With
    # defaults, these would be the averages of the data over the
    # periods 1992-2000 and 2001-2019, inclusive, respectively:
    oht_partial_tp1 = np.nanmean(oht[j_t1_0:j_t1_1+1])
    oht_partial_tp2 = np.nanmean(oht[j_t2_0:j_t2_1+1])
    
    # From here, the extrapolation part is mostly 'hard-coded'
    # for the defaults.
    # 
    # Specifically, assume that the climatologies of the missing
    # 12 years (1980-1991) are the same as the climatologies of
    # the nearest 12 years of data (1992-2003). Neglect the
    # missing 2 years at the end of the second time average
    # period (i.e., approximate the average over years 2001-
    # 2021 as being the same as that over years 2001-2019).
    # 
    n_yr_missing_t1 = max(0, yr_actual[0] - yr_required[0])
    
    # Climatology (mean and standard deviation) of missing
    # years:
    oht_miss_mean = np.nanmean(oht[:n_yr_missing_t1+1])
    oht_miss_sd = np.nanstd(oht[:n_yr_missing_t1+1])
    
    if verbosity > 1:
        print(f"Climatology of missing years {yr_required[0]}"
              + f"-{yr_required[n_yr_missing_t1-1]} assumed "
              + "to match that of "
              + f"{yr_actual[0]}-"
              + f"{yr_actual[n_yr_missing_t1-1]}")
        print(f"Mean = {oht_factor*oht_miss_mean:.3f} "
              + f"{oht_units}, {n_std_clim_extrapolate}x std. "
              + f"dev. = {oht_factor*oht_miss_sd:.3f} "
              + oht_units)
    
    # Assuming this climatological range of possible OHTs for
    # the missing years, calculate the minimum and maximum
    # possible averages for the whole of the first period
    # 1980-2000.
    # 
    # In other words, for 1980-1991 (missing data), assume they
    # all took the maximum of climatology. Combine that with the
    # actual data for the remaining years 1992-2000 and get the
    # average for 1980-2000. This amounts to a weighted average
    # of the average for the known years 1992-2000 and the
    # extrapolated missing years 1980-1991:
    oht_extrap_max_avg1 = (
        n_yr_missing_t1*(oht_miss_mean
                         + n_std_clim_extrapolate*oht_miss_sd)
        + (n_yr_avg_1 - n_yr_missing_t1)*oht_partial_tp1
    ) / n_yr_avg_1
    
    # And similarly for the minimum:
    oht_extrap_min_avg1 = (
        n_yr_missing_t1*(oht_miss_mean
                         - n_std_clim_extrapolate*oht_miss_sd)
        + (n_yr_avg_1 - n_yr_missing_t1)*oht_partial_tp1
    ) / n_yr_avg_1
    
    # The central estimate of the OHT(C) change is based on the
    # partial averages, i.e., data only, no explicit
    # extrapolation:
    doht_central = oht_factor*(oht_partial_tp2 - oht_partial_tp1)
    
    # The uncertainty interval is the extrapolated minimum and
    # maximum possible differences based on the assumption of
    # persistent climatologies above:
    doht_err = oht_factor*np.array(sorted([
        oht_partial_tp2 - oht_extrap_min_avg1,
        oht_partial_tp2 - oht_extrap_max_avg1]))
    
    if verbosity > 0:
        print(f"delta_oht{'c' if convergence else ''}"
              + f"_obs = {doht_central:.3f} {oht_units}, "
              + f"uncertainty interval "
              + f"[{doht_err[0]:.3f}, "
              + f"{doht_err[1]:.3f}] {oht_units} "
              + f"(range = {doht_err[1]-doht_err[0]:.3f} "
              + f"{oht_units})")
    
    if return_full_data:
        
        # Return diagnostics for illustration of method
        
        j_partial_tp1 = [j_t1_0, j_t1_1]
        j_partial_tp2 = [j_t2_0, j_t2_1]
        
        yr_miss = (yr_actual[0], yr_actual[n_yr_missing_t1-1])
        yr_extrap_clim = (yr_required[0], yr_actual[0]-1)
        
        oht_extrap_range_avg1 = (oht_extrap_min_avg1,
                                 oht_extrap_max_avg1)
        
        return (yr_actual, oht, j_partial_tp1,
                oht_partial_tp1, j_partial_tp2, oht_partial_tp2,
                yr_miss, yr_extrap_clim, oht_miss_mean,
                oht_miss_sd, oht_extrap_range_avg1)
    else:
        return doht_central, doht_err
