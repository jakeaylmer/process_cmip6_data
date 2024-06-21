"""General routines used by plotting scripts: loading and
preparing data and parsing command-line arguments.
"""

from argparse import ArgumentParser
from datetime import datetime as dt, timezone as tz

import numpy as np
from tabulate import tabulate

from . import maths, model_lists

from process_cmip6_data import load_data
from process_cmip6_data.src.metadata import _get_path


# ------------------------------------------------------------ #
# Data loading functions/wrappers
# ============================================================ #

def get_path(path_file):
    """Reads the first line of a file called path_file located
    in the "paths" package directory and returns a corresponding
    pathlib.Path instance. If the file is not found, returns a
    default Path (the current working directory). If the
    directory written in the text file is not absolute (e.g.,
    starting /) then it is interpreted (or created) relative to
    the current working directory.
    
    This is using the function defined in processing code
    src/metadata.py module, i.e., points to the same "paths"
    directory. In the context of making the manuscript figures,
    it is used by utils/observations.py and utils/plot_tools.py
    modules to find observational data and know where to save
    figures to, respectively. So, those locations can and should
    be set in the corresponding paths/path_*.txt files,
    specifically:
        
        paths/path_atmospheric_reanalyses_processed_data.txt
        paths/path_cmip6_processed_data.txt
        paths/path_ecco_oht.txt
        paths/path_manuscript_figures.txt
        paths/path_nsidc_processed_data.txt
    
    """
    return _get_path(path_file)



# Constants to multiple different diagnostics by (north, south),
# used by loading routines below:
scale_factors = {
    "aht" : (1000.0, -1000.0),
    "dhdt": (1000.0, 1000.0),
    "iel" : (1.0, -1.0),
    "oht" : (1000.0, -1000.0)
}
has_scale_factor = list(scale_factors.keys())


def load_data_multi_model(models, diagnostics,
        lat_eval=(65.0, -60.0), n_year_averages=21,
        time_period_1=(1980, 2000), time_period_2=(2001, 2021),
        datetime_lims=(dt(1850, 1, 1), dt(2100, 12, 31)),
        full_time_series=False,
        experiment_ids=["historical", "ssp370"], verbosity=1):
    """Load data for multiple models, for specified diagnostics,
    experiments, for both hemispherse, computing differences
    between two time-averaging periods or fixed year averages,
    or returning the full time series data, as specified.
    
    This utilises the API module load_data from the processing
    code.
    
    
    Parameters
    ----------
    models : list of str (length n_models > 0)
        Names of models to load data for.
    
    diagnostics : list of str (length n_diags > 0)
        Names of diagnostics to load. These should be "keywords"
        as defined by the API module load_data. Typical choices
        are, e.g., "aht", "oht", "iel", "f_sw_surf", "f_olr",
        "dhdt", "f_up", "f_down", among others.
    
    
    Optional parameters
    -------------------
    lat_eval : 2-tuple of float, default = (65.0, -60.0)
        North and south latitudes to evaluate diagnostics at,
        where applicable.
    
    n_year_averages : int, default = 21
        Number of years making up the contiguous anomaly
        averages of data that is returned if experiment_ids
        (see below) contains "piControl".
    
    time_period_1, time_period_2 : 2-tuple of float
        Start and end years to average over (inclusive) for the
        first and second time-averaging periods used to compute
        differences of specified diagnostics where the argument
        experiment_ids (see below) does not contain "piControl".
        Defaults: (1980, 2000) and (2001, 2021) respectively.
        Only used if experiment_ids (see below) does not
        contain "piControl".
    
    datetime_lims : 2-tuple of datetime.datetime
        Time span of data to load. Only used when the argument
        full_time_series = True (see below).
    
    full_time_series : bool, default = False
        Whether to return the full time series data rather than
        time differences or averages.
    
    experiment_ids : list of str, length > 0
                     default = ["historical", "ssp370"]
        Names of experiments to load and concatenate
    
    verbosity : int, default = 1
        Level of detail to print to the console (0 = none).
    
    
    Returns
    -------
    save_dict_n, save_dict_s : dict of list of array
        Dictionaries where each key is an entry in the list
        diagnostics and the values are lists of length n_models
        containing data such that, for instance assuming "aht"
        is in diagnostics:
        
        save_dict_n["aht"][0]
        
        is an array of length equal to one of:
            
            a) the number of ensemble members for the first
               model in the input list models, for combined
               experiments, or;
            b) the number of n_year_averages year averages in
               the piControl simulation for that model, or;
            c) the number of time steps in the data for that
               model, for combined experiments, for specified
               datetime_lims,
        
        depending on input options. One dictionary is returned
        for each hemisphere (north and south respectively).
        
        In the case of option (c), which occurs with input
        full_time_series = True, each dictionary also contains
        the entry "date", containing the datetime coordinates.
    
    """
    
    n_models = len(models)
    n_diags = len(diagnostics)
    
    if verbosity > 0:
        print(f"Loading {n_diags} diagnostics for {n_models} "
              + f"models, experiment"
              + f"{'s' if len(experiment_ids)>1 else ''}: "
              + " + ".join(experiment_ids))
    
    # If a single latitude provided, put into a tuple as this
    # is required by the load_data function:
    if type(lat_eval) in [float, np.float32, np.float64]:
        lat_eval = (lat_eval, -lat_eval)
    
    # Keyword arguments passed to the load_data function:
    load_kw = {"lat_eval": lat_eval, "verbose": verbosity > 1}
    
    # Determine which load function to use:
    if len(experiment_ids) == 1:
        load_func = load_data.one_experiment
        load_kw["experiment_id"] = experiment_ids[0]
    else:
        load_func = load_data.two_experiments
        load_kw["experiment_ids"] = experiment_ids
    
    # Only presence of "piControl" in experiment_ids switches
    # behaviour of this function from time differencing between
    # time_period_1 and time_period_2 to contiguous
    # n_year_averages averages:
    anomalies = "piControl" in experiment_ids
    
    # Get datetime limits for each time averaging period (data
    # is at most monthly frequency so these are sufficient):
    dt1_0 = dt(time_period_1[0], 1, 1)
    dt1_1 = dt(time_period_1[1], 12, 31)
    dt2_0 = dt(time_period_2[0], 1, 1)
    dt2_1 = dt(time_period_2[1], 12, 31)
    
    # Set up the dictionaries that are returned:
    save_dict_n = {diag: list() for diag in diagnostics}
    save_dict_s = {diag: list() for diag in diagnostics}
    
    # Switch behaviour of this function is specified to return
    # full data (regardless of experiment):
    if full_time_series:
        load_kw["datetime_lims"] = datetime_lims
        save_dict_n["date"] = list()
        save_dict_s["date"] = list()
    
    for m in range(n_models):
        
        if verbosity > 0:
            print(f"Loading {m+1} of {n_models}: {models[m]}")
        
        # Load the data:
        dates_m, _, _, _, diags_n_m, diags_s_m, _ = \
            load_func(diagnostics, models[m], **load_kw)
        
        if full_time_series:
            save_dict_n["date"].append(dates_m.copy())
            save_dict_s["date"].append(dates_m.copy())
        
        for k in range(n_diags):
            
            # Convert units if required:
            if diagnostics[k] in has_scale_factor:
                diags_n_m[k] *= scale_factors[diagnostics[k]][0]
                diags_s_m[k] *= scale_factors[diagnostics[k]][1]
            
            if full_time_series:
                
                # Just save data as loaded:
                save_dict_n[diagnostics[k]].append(
                    diags_n_m[k].copy())
                save_dict_s[diagnostics[k]].append(
                    diags_s_m[k].copy())
            
            else:
                if anomalies:
                    
                    # Compute non-overlapping n_year_averages.
                    # The number of such averages:
                    n_avgs_mk = len(dates_m[k]) // n_year_averages
                    
                    # If the number of time steps is not an
                    # exact multiple of the averaging length, we
                    # discard the remainder at the end of the
                    # time series:
                    n_discard = len(dates_m[k]) % n_year_averages
                    
                    if n_discard > 0:
                        diags_n_m[k] = diags_n_m[k][:-n_discard]
                        diags_s_m[k] = diags_s_m[k][:-n_discard]
                    
                    # Remove ensemble member axis (should be
                    # only one member in this case, for
                    # piControl, anyway):
                    diags_n_m[k] = diags_n_m[k][:,0]
                    diags_s_m[k] = diags_s_m[k][:,0]
                    
                    # Convert to anomalies:
                    diags_n_m[k] -= np.nanmean(diags_n_m[k])
                    diags_s_m[k] -= np.nanmean(diags_s_m[k])
                    
                    # Compute the averages and append to the
                    # list for this diagnostic:
                    save_dict_n[diagnostics[k]].append(
                        np.nanmean(
                            np.reshape(diags_n_m[k],
                                (n_avgs_mk, n_year_averages)),
                            axis=1))
                    
                    save_dict_s[diagnostics[k]].append(
                        np.nanmean(
                            np.reshape(diags_s_m[k],
                                (n_avgs_mk, n_year_averages)),
                            axis=1))
                
                else:
                    
                    # Simply compute differences between the two
                    # time averaging periods; time averaging
                    # indices for each period:
                    j1 = (dates_m[k] >= dt1_0) \
                       & (dates_m[k] <= dt1_1)
                    j2 = (dates_m[k] >= dt2_0) \
                       & (dates_m[k] <= dt2_1)
                    
                    # Append the differences for each ensemble
                    # member to the results dictionary, list for
                    # this diagnostic:
                    save_dict_n[diagnostics[k]].append(
                        np.nanmean(diags_n_m[k][j2], axis=0)
                        - np.nanmean(diags_n_m[k][j1], axis=0))
                    
                    save_dict_s[diagnostics[k]].append(
                        np.nanmean(diags_s_m[k][j2], axis=0)
                        - np.nanmean(diags_s_m[k][j1], axis=0))
                # [ End if anomalies ]
            
            # [ End if full_time_series ]
        # [ End for k in range(n_diags) ]
    # [ End for m in range(n_models) ]
    
    return save_dict_n, save_dict_s



def load_data_multi_model_one_hemi_two_ref_lats(models,
        diagnostics, hemi, ref_lats=(-60.0, -70.0),
        fix_integral_diagnostics=["aht", "dhdt", "oht"],
        fix_averaged_diagnostics=["f_down", "f_olr",
                                  "f_sw_surf", "f_up", "tas"],
        **kwargs):
    """Load model diagnostics, differenced between two time
    averaging periods, for a single hemisphere. It also re-
    calculates for a second reference latitude that isn't the
    north or south pole using the existing polar-cap data.
    
    This utilises the API module load_data from the processing
    code.
    
    
    Parameters
    ----------
    models : list of str (length n_models > 0)
        Names of models to load data for.
    
    diagnostics : list of str (length n_diags > 0)
        Names of diagnostics to load. These should be "keywords"
        as defined by the API module load_data. Typical choices
        are, e.g., "aht", "oht", "iel", "f_sw_surf", "f_olr",
        "dhdt", "f_up", "f_down", among others.
    
    hemi : str, "n" or "s"
        Which hemisphere to analyse (north or south).
    
    
    Optional parameters
    -------------------
    ref_lats : length-2 iterable of float,
               default = (-60.0, -70.0)
        Reference latitudes, (lower_latitude, upper_latitude)
        in degrees north (i.e., should be negative if
        hemi == "s").
    
    fix_integral_diagnostics : list of str
                               default = ["aht", "dhdt", "oht"]
        Diagnostics that, if included in diagnostics, should be
        re-calculated for two reference latitudes, assuming they
        are area-integrated quantities (e.g., OHT, AHT, dh/dt).
    
    fix_averaged_diagnostics : list of str
        Diagnostics that, if included in diagnostics, should be
        re-calculated for two reference latitudes, assuming they
        are area-averaged quantities (e.g., vertical heat flux).
        
        Note that if ref_lats[1] is the pole (greater than
        +/- 85), this part of the code is ignored and this
        routine is essentially a single-hemisphere wrapper for
        script_tools.load_data_multi_model().
    
    Additional keyword arguments are passed to
    script_tools.load_data_multi_model().
    
    
    Returns
    -------
    data : dict of list of array
        Dictionary of data as output from script_tools function
        load_data_multi_model but for one hemisphere only,
        including any corrections for two reference latitudes
        as required.
    
    """
    
    if hemi not in "ns":
        print("WARNING: assuming you meant southern hemisphere"
              + " [choose hemi = \"n\" OR \"s\"; utils.script_"
              + "tools.load_data_multi_model_one_hemi_two_ref_"
              + "lats()]")
        hemi = "s"
    
    n_models = len(models)
    
    # For load_data_multi_model(), the lat_eval argument expects
    # a northern and a southern reference latitude. So here,
    # where we only want one hemisphere, need to put in a dummy
    # reference latitude for the other, unwanted hemisphere:
    if hemi == "n":
        # Put a dummy reference latitude for the southern
        # hemisphere of 60.0 degrees_south:
        load_lat_eval_0 = (abs(ref_lats[0]), -60.0)
    else:
        # Put a dummy reference latitude for the northern
        # hemisphere of 65.0 degrees_north:
        load_lat_eval_0 = (65.0, -abs(ref_lats[0]))
    
    # Load data for the first (lower) reference latitude:
    data = load_data_multi_model(models, diagnostics,
                                 lat_eval=load_lat_eval_0,
                                 **kwargs)[int(hemi=="s")]
    
    # Account for the second reference latitude, if it is not
    # the north or south pole. This means, for integrated
    # diagnostics such as heat transport, subtracting that
    # evaluated at reference latitude 1 from that evaluated at
    # reference latitude 0, and for averages such as vertical
    # heat fluxes and temperature, computing a weighted sum of
    # the existing diagnostics.
    # 
    # Thus, here the data for the second reference latitude is
    # loaded and used to modify the data loaded above, so that
    # the rest of the code can run using the array called "data"
    # regardless of whether this is needed.
    # 
    # The initial check is with +/- 85.0 as that is the maximum
    # latitude that averages/integrals were computed from:
    # 
    if abs(ref_lats[1]) <= 85.0:
        
        # Modify the loading keyword arguments to load from the
        # average between reference latitude 1 (i.e., the non-
        # pole latitude) and the pole:
        if hemi == "n":
            load_lat_eval_1 = (abs(ref_lats[1]), -60.0)
        else:
            load_lat_eval_1 = (65.0, -abs(ref_lats[1]))
        
        data_lat1 = load_data_multi_model(models, diagnostics,
            lat_eval=load_lat_eval_1, **kwargs)[int(hemi=="s")]
        
        # Correct integral diagnostics: simply subtract that
        # evaluated at the second (higher) reference latitude
        # from that evaluated at the first (lower) reference 
        # latitude:
        for d in fix_integral_diagnostics:
            if d in diagnostics:
                for m in range(n_models):
                    data[d][m] -= data_lat1[d][m]
        
        # Correct averaged diagnostics. First convert them to
        # area integrals by multiplying by the total area of the
        # corresponding polar caps, which are
        # 
        #     2*pi*(a**2)*[1 - sin(phi)]
        # 
        # for reference latitude phi and Earth radius a. Then,
        # take the difference and divide by the area between
        # reference latitudes to get the new average. This
        # amounts to a weighted difference (sum) with weighting
        # factors given by the [1 - sin(phi)] terms.
        # 
        # Work out the weights in advance, also absorbing the
        # minus from differencing into the second weight w1:
        w0 = 1.0 - np.sin(abs(ref_lats[0])*np.pi/180.0)
        w1 = np.sin(abs(ref_lats[1])*np.pi/180.0) - 1.0
        w = w0 + w1
        
        for d in fix_averaged_diagnostics:
            if d in diagnostics:
                for m in range(n_models):
                    data[d][m] = \
                        (w0*data[d][m] + w1*data_lat1[d][m])/w
    
    return data



def get_data_for_lagged_correlation_plots(models,
        diagnostics_combinations, year_span=(1956, 2045),
        fixed_yr_avg_1=(1980, 2000),
        fixed_yr_avg_2=(2001, 2021), ref_lats_n=(65.0, 90.0),
        ref_lats_s=(-60.0, -90.0), future_experiment="ssp370",
        as_ensemble=True):
    """Load data and calculate lagged correlations between time
    changes in specified diagnostics, for specified
    models, and lags, for both hemispheres. 
    
    Here, lagged correlations are computed between delta X(t)
    and delta Y, where the Y variable is fixed, always computed
    as the difference between the same two time periods, and the
    X(t) variable is lagged by t years relative to those time
    periods.
    
    
    Parameters
    ----------
    models : list of str (length n_models > 0)
        Names of models to load data for.
    
    diagnostics_combinations : list of length n_combos
                               of 2-tuples of str
        Names of combinations of diagnostics (X(t), Y) to
        compute lagged correlations between. These should be
        "keywords" as defined by the API module load_data.
        Typical choices are, e.g., "aht", "oht", "iel",
        "f_sw_surf", "f_olr", "dhdt", "f_up", "f_down",
        among others.
    
    
    Optional parameters
    -------------------
    year_span : 2-tuple of int, default = (1956, 2045)
        Total year range to compute lagged correlations over
        (determines the overall number of lags, n_lags, and
        should generally be centred on the fixed averaging
        periods; see below).
    
    fixed_yr_avg_1, fixed_yr_avg_2 : 2-tuple of int
        The fixed time averaging periods for the differences in
        the Y variables (second of each tuple in
        diagnostics_combinations).
    
    ref_lats_n, ref_lats_s : 2-tuple of float,
                             default = (65.0, 90.0)
                                   and (-60.0, -90.0)
                             respectively
        Reference latitudes, (lower_latitude, upper_latitude)
        in degrees north for north and southern hemispheres
        respectively.
    
    future_experiment : str, default = "ssp370"
        Which future experiment to combine with historical.
    
    as_ensemble : bool, default = True
        If True, combines all models and their ensemble members
        and computes the lagged correlations in combination. If
        False, calculates separately for each input model (it is
        assumed that each input model has more than one ensemble
        member).
    
    
    Returns
    -------
    t_lag : array of int, shape (n_lags,)
        The time lags in years relative to the fixed time-
        averaging periods (corresponds to axis=-1 in r_n and
        r_s; see below). Negative lags imply that the first
        diagnostic in a given combination leads the second.
    
    r_n, r_s : arrays of shape (n_combos, n_lags)
               if as_ensemble = True, otherwise arrays of shape
               (n_combos, n_models, n_lags)
        The correlations as a function of lag (second axis) for
        each diagnostic combination (first axis). If as_ensemble
        is False, i.e., correlations are computed for individual
        models, an extra axis is inserted (axis=1) corresponding
        to the model.
    
    n_ens : int or array of int, shape (n_models,)
        The number of ensemble members entering the lagged
        correlations; either the total number of members summed
        across all models if as_ensemble = True, otherwise
        an array containing the number of members for each model
        in input list models.
    
    """
    
    # Get unique set of diagnostics from the input diagnostics
    # combinations (X, Y):
    diagnostics = sorted(list(set(
        [x[0] for x in diagnostics_combinations]
        + [x[1] for x in diagnostics_combinations])))
    
    n_models = len(models)
    n_diags = len(diagnostics)
    n_combos = len(diagnostics_combinations)
    
    # First get the delta_X (as a function of lag) and delta_Y
    # (difference between the fixed time periods), saving data
    # to dictionaries:
    dX_n = {k: list() for k in diagnostics}
    dX_s = {k: list() for k in diagnostics}
    dY_n = {k: list() for k in diagnostics}
    dY_s = {k: list() for k in diagnostics}
    
    # Number of years in each averaging interval:
    n_yr_avg = fixed_yr_avg_1[1] - fixed_yr_avg_1[0] + 1
    
    # Shifted time axes, as years, for the start/end of each 
    # lagged averaging period, used in the loop below when
    # computing differences:
    yr_avg_0 = np.arange(year_span[0],
                         year_span[1] - 2*n_yr_avg + 2)
    yr_avg_1 = np.arange(year_span[0] + n_yr_avg - 1,
                         year_span[1] - n_yr_avg + 1)
    yr_avg_2 = np.arange(year_span[0] + n_yr_avg,
                         year_span[1] - n_yr_avg + 2)
    yr_avg_3 = np.arange(year_span[0] + 2*n_yr_avg - 1,
                         year_span[1] + 1)
    
    # Negative t_lag => dX(t) leads dY:
    t_lag = yr_avg_0 - fixed_yr_avg_1[0]
    n_lags = len(yr_avg_0)
    
    # Load data:
    load_kw = {
        "experiment_ids": ("historical", future_experiment),
        "datetime_lims": (dt(year_span[0], 1, 1),
                          dt(year_span[1], 12, 31)),
        "full_time_series": True, "verbosity": 1}
    
    data_n = load_data_multi_model_one_hemi_two_ref_lats(models,
        diagnostics, "n", ref_lats=ref_lats_n, **load_kw)
    
    data_s = load_data_multi_model_one_hemi_two_ref_lats(models,
        diagnostics, "s", ref_lats=ref_lats_s, **load_kw)
    
    # Loop over each model and diagnostic to compute delta_X and
    # delta_Y. Afterwards, compute the lagged correlations:
    for m in range(n_models):
        for d in diagnostics:
            
            # Determine years as time axis for this model/
            # diagnostic (it should generally be the same,
            # anyway):
            years_dm = np.array([dm.year for dm in
                                 data_n["date"][m][0]])
            
            # Indices for the fixed time interval, for delta_Y:
            t_avg_fx_01 = (years_dm >= fixed_yr_avg_1[0]) \
                        & (years_dm <= fixed_yr_avg_1[1])
            t_avg_fx_23 = (years_dm >= fixed_yr_avg_2[0]) \
                        & (years_dm <= fixed_yr_avg_2[1])
            
            # Indices for lagged time intervals, for delta_X:
            t_avg_g_01 = [(years_dm >= yr_avg_0[g])
                          & (years_dm <= yr_avg_1[g])
                          for g in range(n_lags)]
            
            t_avg_g_23 = [(years_dm >= yr_avg_2[g])
                          & (years_dm <= yr_avg_3[g])
                          for g in range(n_lags)]
            
            # For differences as a function of lag (lag w.r.t.
            # fixed time period, shifting both averaging
            # intervals), need arrays of shape (n_lags, number
            # of ens. members):
            n_ens_m = len(data_n[d][m][0,:])
            
            # Differences over the fixed time period:
            dY_n[d].append(
                np.nanmean(data_n[d][m][t_avg_fx_23], axis=0)
                - np.nanmean(data_n[d][m][t_avg_fx_01], axis=0))
            
            dY_s[d].append(
                np.nanmean(data_s[d][m][t_avg_fx_23], axis=0)
                - np.nanmean(data_s[d][m][t_avg_fx_01], axis=0))
            
            # Differences over the variable/lagged time periods:
            dX_n_dm = np.zeros((n_ens_m, n_lags))
            dX_s_dm = np.zeros((n_ens_m, n_lags))
            
            for g in range(n_lags):
                
                dX_n_dm[:,g] = \
                    np.nanmean(data_n[d][m][t_avg_g_23[g]],
                               axis=0) \
                    - np.nanmean(data_n[d][m][t_avg_g_01[g]],
                                 axis=0)
                
                dX_s_dm[:,g] = \
                    np.nanmean(data_s[d][m][t_avg_g_23[g]],
                               axis=0) \
                    - np.nanmean(data_s[d][m][t_avg_g_01[g]],
                                 axis=0)
            
            dX_n[d].append(dX_n_dm)
            dX_s[d].append(dX_s_dm)
            
        # [ End for d in diagnostics ]
    # [ End for m in range(n_models) ]
    
    if as_ensemble:
        
        # Total number of ensemble members, all models:
        n_ens = sum([len(dX_n[diagnostics[0]][m])
                     for m in range(n_models)])
        
        # Combine prepared diagnostics above for all models
        # (concatenate along ensemble member axes):
        for diag in diagnostics:
            dX_n[diag] = np.concatenate(dX_n[diag], axis=0)
            dX_s[diag] = np.concatenate(dX_s[diag], axis=0)
            dY_n[diag] = np.concatenate(dY_n[diag], axis=0)
            dY_s[diag] = np.concatenate(dY_s[diag], axis=0)
        
        r_n = np.zeros((n_combos, n_lags))
        r_s = np.zeros((n_combos, n_lags))
        
        for j in range(n_combos):
            r_n[j,:] = maths.cross_correlation(
                dX_n[diagnostics_combinations[j][0]],
                dY_n[diagnostics_combinations[j][1]])
            r_s[j,:] = maths.cross_correlation(
                dX_s[diagnostics_combinations[j][0]],
                dY_s[diagnostics_combinations[j][1]])
    else:
        
        # Total number of ensemble members per model:
        n_ens = np.array([len(dX_n[diagnostics[0]][m])
                          for m in range(n_models)])
        
        r_n = np.zeros((n_models, n_combos, n_lags))
        r_s = np.zeros((n_models, n_combos, n_lags))
        
        for m in range(n_models):
            for j in range(n_combos):
                r_n[m,j,:] = maths.cross_correlation(
                    dX_n[diagnostics_combinations[j][0]][m],
                    dY_n[diagnostics_combinations[j][1]][m])
                r_s[m,j,:] = maths.cross_correlation(
                    dX_s[diagnostics_combinations[j][0]][m],
                    dY_s[diagnostics_combinations[j][1]][m])
    
    return t_lag, r_n, r_s, n_ens


# ------------------------------------------------------------ #
# Command-line argument parsing
# 
# Provide an ArgumentParser and add groups of related arguments
# via the add_*_cmd_args functions. Thus, sometimes arguments
# are available for a script that doesn't use it, or is
# interpreted slightly differently but this should be obvious
# from context
# ============================================================ #

def argument_parser(description=""):
    """Provides an argparse.ArgumentParser() instance to define
    command line arguments in. Takes an optional description
    (str) which shows when using --help.
    """
    return ArgumentParser(description=description)



def add_cmip_selection_cmd_args(prsr, default_future=False,
                                both_hemispheres=False):
    """Add common, general command-line arguments to an argparse
    ArgumentParser instance (prsr) for selection and analysis
    of CMIP6 data.
    
    
    Optional parameters
    -------------------
    default_future : bool, default = False
        Set defaults for future analysis if True or historical
        analysis if False (currently only affects year selection
        arguments).
    
    both_hemispheres : bool, default = False
        Adds for analysing both hemispheres if True,
        specifically, a separate set of reference latitudes for
        each. If False (default), only includes one reference
        latitude argument.
    
    """
    prsr.add_argument("--yravg1", type=int, nargs=2,
                      default=(1980, 2000),
                      help="First average (y1, y2) inclusive")
    
    prsr.add_argument("--yravg2", type=int, nargs=2,
                      default=(2030, 2050) if default_future
                              else (2001, 2021),
                      help="Second average (y1, y2) inclusive")
    
    if both_hemispheres:
        for h, default in zip("ns", [(65., 90.),(-60., -90.)]):
            prsr.add_argument(f"--reflats{h}", type=float,
                              nargs=2, default=default,
                              help="Reference latitudes (phi_0,"
                                   + f" phi_1), {h}h")
    else:
        prsr.add_argument("-l", "--reflats", type=float,
                          nargs=2, default=(65., 90.),
                          help="Reference latitudes (phi_0, "
                               + "phi_1")
    
    prsr.add_argument("-e", "--exp", type=str, default="ssp370",
                      choices=["ssp370", "ssp585"],
                      help="For future/extension to historical")
    
    prsr.add_argument("--exclude", type=str, nargs="*",
                      default=[],
                      help="Models to exclude from analysis")



def add_ebm_fitting_cmd_args(prsr):
    """Add common command-line arguments to an argparse
    ArgumentParser instance (prsr) relevant to EBM fitting.
    """
    prsr.add_argument("--piparams", type=str, nargs="*",
        default=["beta", "bup", "S"],
        choices=["bc", "beta", "bup", "h", "S"],
        help="Which, if any, piControl parameters to use for "
             + "EBM fitting")
    
    prsr.add_argument("--exclude-from-ebm", type=str, nargs="*",
                      default=[],
                      help="Models to exclude from EBM fitting")



def add_plotting_cmd_args(prsr, n_panels=2,
        text_labels=["ebm", "cmip6", "obs"], n_cbars=2,
        default_cbar_cats=[(-60., 60., 20.), (0., 3.5, 0.5)],
        default_cbar_ticks=[(-60., 60., 20.), (0., 3., 1.)]):
    """Add common command-line arguments to an argparse
    ArgumentParser instance (prsr) relevant to producing
    and saving figures.
    
    
    Optional parameters
    -------------------
    n_panels : int, default = 2
        Total number of panels in the figure. This determines
        the number of arguments to add for setting axes limits
        and tick labels. Panels here are referred to as a, b,
        ... (n_panels letters) regardless of what they are
        actually labelled as on the figures.
    
    text_labels : list of str, default = ["ebm", "cmip6", "obs"]
        An argument is added called "--{p}-dr-{x}-text" (taking
        2 float arguments), for each panel p (i.e., p =a, b, ...
        up to n_panels) and text label x in text_labels, used
        for offsetting the positions of annotations labelled as
        x in the script code.
    
    n_cbars : int, default = 2
        Total number of color bars in the figure. Like n_panels,
        this determines the number of arguments to add for
        setting colorbar axes limits and tick labels.
    
    default_cbar_cats : list of length n_cbars of 3-tuple of
                        float
        Default color bar category boundaries for each colorbar
        (lowest category lower bound, highest category upper
        bound, step size). Defaults are set for manuscript Fig.
        1.
    
    default_cbar_ticks : list of length n_cbars of 3-tuple of
                         float
        As in default_cbar_cats but for the displayed tick
        labels.
    
    
    Note
    ----
    The following arguments
        
        --{p}-xlim (xmin, xmax)
        --{p}-xticks (start, stop inclusive, step)
        --{p}-ylim (ymin, ymax)
        --{p}-yticks (start, stop inclusive, step)
    
    are added for each panel of n_panels p = a, b, ...
    The defaults are set to zero (in which case the
    corresponding function (plot_tools.set_axis_ticks_from_cmd)
    does not set limits at all, instead using matplotlib's
    default axis limit/ticks).
    
    """
    
    prsr.add_argument("--savefig", action="store_true")
    prsr.add_argument("--savefigname", type=str, default="",
                      help="File name (without extension) for "
                           + "saved figure(s)")
    prsr.add_argument("--savefigtitle", type=str, default="",
                      help="Title for file metadata if saving")
    
    pl = "abcdefghijkl"
    
    for p in range(n_panels):
        # For setting axis limits/ticks on each panel (p) and
        # axis (x and y):
        for x in "xy":
            prsr.add_argument(f"--{pl[p]}-{x}lim", type=float,
                              nargs=2, default=(0.0, 0.0))
            prsr.add_argument(f"--{pl[p]}-{x}ticks", type=float,
                              nargs=3, default=(0.0, 0.0, 0.0))
        
        # For offsetting the positions of the text labels
        # (one for each panel):
        for tl in text_labels:
            prsr.add_argument(f"--{pl[p]}-dr-{tl.lower()}-text",
                              type=float, nargs=2,
                              default=(0.0, 0.0))
    
    cl = pl
    
    for cj in range(n_cbars):
        # For setting the color bar categories and ticks:
        prsr.add_argument(f"--{cl[cj]}-cbar-cats", type=float,
            nargs=3,
            default=default_cbar_cats[cj%len(default_cbar_cats)])
        
        prsr.add_argument(f"--{cl[cj]}-cbar-ticks", type=float,
            nargs=3,
            default=default_cbar_ticks[cj%len(default_cbar_ticks)])



def get_args(prsr, models=None, piControl=False,
        suppress_options=[], ad_hoc_options={}):
    """Prints selected parsed command-line arguments where
    available. Requires the Namespace output from
    argparse.ArgumentParser().parse_args(), cmd, and optionally
    model list models. Returns model list again, possibly
    modified if excluding any models using the
    --exclude MODEL1 MODEL2 [...] option (see
    add_cmip6_selection_args).
    
    
    Optional parameters
    -------------------
    piControl : bool, default = False
        Changes the way time selection options are interpreted
        if piControl is being analysed rather than historical/
        future simulations.
    
    ad_hoc_options : dict, default = {}
        Manually add options from cmd (keys) to print values of.
        Options for formatting are limited; the values 2-tuples
        of string for the option label and format string for
        the value, e.g., with:
        
        dict = {"my_cmd_arg_1": ("Arg 1", "{:.2f}"),
                "my_cmd_arg_2": ("Arg 2", "{}")}
        
        prints:
        
        Arg 1  < arg1 value formatted .2f >
        Arg 2  < arg2 value formatted using repr() >
    
    
    Returns
    -------
    model_list : list of str or None
        List of model names or None.
    
    """
    
    cmd = prsr.parse_args()
    
    # Collect options to print in a dictionary. Always print the
    # time/date of running the script:
    info = {"Time and date":
            dt.now(tz.utc).strftime("%H:%M UTC, %a %d %b %Y")}
    
    sup = suppress_options  # alias
    
    # Print info about time averaging/selection. The standard
    # arguments are "yravg1" and "yravg2", tuples of int for
    # year ranges, each denoting an averagine period. However,
    # if analysis is on piControl data, then the difference of
    # values in "yravg1" is interpreted as the n-year contiguous
    # averaging interval applied across all piControl
    # simulations.
    if piControl:
        # Using n-year averages, not difference of two averaging
        # periods and this is set from "yravg1":
        if hasattr(cmd, "yravg1") and "yravg1" not in sup:
            info["n-year averages"] = \
                cmd.yravg1[1] - cmd.yravg1[0] + 1
    else:
        # Interpret "yravg1" and "yravg2" as fixed averaging
        # periods:
        for j in "12":
            opt = f"yravg{j}"
            if hasattr(cmd, opt) and opt not in sup:
                yrs = getattr(cmd, opt)
                info[f"Average period {j}"] = (
                    f"{yrs[0]}-{yrs[1]} ({yrs[1]-yrs[0]+1}-"
                    + "year average)")
    
    # Print info about reference latitudes. The arguments
    # available are either "reflatsn" and "reflatss" if both
    # hemispheres are being analysed in one script, or just
    # "reflats" if one hemisphere is being analysed, in which
    # case the hemisphere is assumed to be southern if both
    # values are negative (in degrees_north):
    if hasattr(cmd, "reflats") and "reflats" not in sup:
        if cmd.reflats[0] < 0.0 and cmd.reflats[1] < 0.0:
            info["Analysing"] = "Southern Ocean"
            info["Ref. latitudes"] = (
                f"{-cmd.reflats[0]:.0f}-"
                + f"{-cmd.reflats[1]:.0f}S")
        else:
            info["Analysing"] = "Arctic"
            info["Ref. latitudes"] = (
                f"{cmd.reflats[0]:.0f}-"
                + f"{cmd.reflats[1]:.0f}N")
    
    elif (hasattr(cmd, "reflatsn") and hasattr(cmd, "reflatss")
          and "reflatsn" not in sup and "reflatss" not in sup):
        
        info["Analysing"] = "Both hemispheres"
        info["Ref. latitudes"] = (f"{cmd.reflatsn[0]:.0f}-"
            + f"{cmd.reflatsn[1]:.0f}N, "
            + f"{abs(cmd.reflatss[0]):.0f}-"
            + f"{abs(cmd.reflatss[1]):.0f}S")
    
    # Print info about model selection. If a models list is not
    # provided (models = None), determine it if possible from
    # the experiment ("exp") input arg, otherwise leave as None
    # for the next part:
    if models is None:
        if piControl:
            models = model_lists.by_experiment["piControl"]
            if "exp" not in sup:
                info["Experiment"] = "piControl"
        elif hasattr(cmd, "exp"):
            models = model_lists.by_experiment[
                f"historical+{cmd.exp}"]
            if "exp" not in sup:
                info["Experiment"] = f"historical + {cmd.exp}"
    
    # If we have a model list (possibly after generating it
    # above -- the below cannot be combined as and "else" to the
    # previous block), remove any explicitly excluded models:
    if models is not None:
        if hasattr(cmd, "exclude"):
            for x in cmd.exclude:
                if x in models:
                    kx = models.index(x)
                    del models[kx]
        
        if "exclude" not in sup and len(cmd.exclude) > 0:
            info["Excluded models"] = ", ".join(cmd.exclude)
        
        if "n_models" not in sup:
            info["Number of models"] = f"{len(models)}"
    
    # Print EBM related options:
    if (hasattr(cmd, "exclude_from_ebm")
        and len(cmd.exclude_from_ebm) > 0
        and "exclude_from_ebm" not in sup):
        
        info["Excluded from EBM"] = \
            ", ".join(cmd.exclude_from_ebm)
    
    if (hasattr(cmd, "piparams")
        and len(cmd.piparams) > 0
        and "piparams" not in sup):
        
        info["piControl parameters"] = ", ".join(cmd.piparams)
    
    # Print save figure options:
    if hasattr(cmd, "savefig") and "savefig" not in sup:
        if cmd.savefig:
            if hasattr(cmd, "savefigname"):
                if cmd.savefigname == "":
                    info["Saving"] = "<no file name provided>"
                else:
                    info["Saving"] = f"{cmd.savefigname}.*"
        else:
            info["Saving"] = "Not saving"
    
    
    # Print any ad-hoc options:
    for k in ad_hoc_options.keys():
        if hasattr(cmd, k):
            info[ad_hoc_options[k][0]] = \
                ad_hoc_options[k][1].format(getattr(cmd, k))
    
    info_table = tabulate(
        [[k, info[k]] for k in sorted(list(info.keys()))],
        headers="firstrow", tablefmt="plain",
        maxcolwidths=[None, 80])
    
    print("\n" + prsr.prog + "\n" + "-"*len(prsr.prog) + "\n"
          + prsr.description + "\n"
          + ("\n"*(len(prsr.description)>0))
          + info_table, end="\n\n")
    
    return cmd, models
