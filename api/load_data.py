"""Load processed diagnostic data."""

from datetime import datetime as dt
from netCDF4 import num2date
import numpy as np

from . import model_diagnostics as mdiags
from ..src.netcdf import (nc_file_path, nc_ripf_labels,
                        nc_time_name)

from my_python_utilities.data_tools.datetime_tools import (
    cftime_to_datetime)
from my_python_utilities.data_tools.nc_tools import (
    get_arrays, get_nc_time_units)


def _get_j_t_data(data):
    """Determine which index of data list construct corresponds
    to time. TODO: make more robust?
    """
    return (4 if len(data) == 7 else 6)


def _time_slice(data,
        datetime_lims=(dt(1, 1, 1), dt(3000, 1, 1)),
        verbose=True):
    """Take time slice of data from specified datetime lims,
    where data is the list
        
        [ripf, datetime, [coords], var_n, var_s]
    
    structure used by the main api functions.
    """
    
    j_t_data = _get_j_t_data(data)
    
    # Time slicing:
    t_want = (data[j_t_data] >= datetime_lims[0]
             ) & (data[j_t_data] <= datetime_lims[1])
    
    # Only need to slice the datetime data and variables:
    data[j_t_data] = data[j_t_data][t_want]  # datetime
    data[j_t_data+1] = data[j_t_data+1][t_want]  # var_n
    data[j_t_data+2] = data[j_t_data+2][t_want]  # var_s
    
    if verbose:
        print("    datetime range: "
            + f"{data[j_t_data][0].strftime('%d-%m-%Y')}"
            + " -- "
            + f"{data[j_t_data][-1].strftime('%d-%m-%Y')}")
    
    return data



def _load_diag_data(data_path, coord_names, var_names,
        lat_eval=(65.0, -65.0),
        id_names=[nc_ripf_labels[x] for x in "ripf"],
        verbose=True):
    """Load diagnostic data, including ripf values,
    convert time data to datetime, and evaluate reference
    latitudes where applicable.
    
    
    Parameters
    ----------
    data_path : str or pathlib.Path instance
        Full path to data to load.
    
    coord_names : list of str
        Spatial coordinate variables (list should be empty or
        length 2, north and south coordinate names).
    
    var_names : list of str, length 2
        NetCDF variable names, north and south.
    
    
    Optional parameters
    -------------------
    lat_eval : tuple (float, float), default = (65.0, -65.0)
        North and south latitudes to evaluate, degrees_north.
        Can also be None, in which case no latitude evaluation
        is done (and the coordinates are returned).
    
    id_names : list of str
        NetCDF variable names for the ripf member axes.
    
    verbose : bool, default = True
        Print information to console.
    
    
    Returns
    -------
    data_ret : list (length 7) of arrays:
        
        data_ret[0-3] : ripf values (n_ens,)
        
    Then, if lat_eval is None:
        data_ret[4] : latitude coordinates north (n_lat_n)
        data_ret[5] : latitude coordinates south (n_lat_s)
    
    Then (or continue from here otherwise):
        data_ret[6 (4)] : datatime values (n_t,)
        data_ret[7 (5)] : variable, north (n_t, n_ens[, n_lat])
        data_ret[8 (6)] : variable, south (n_t, n_ens[, n_lat])
    
    lat_eval_n, lat_eval_s : float
        Actual latitudes of evaluation (or NaN if not
        applicable).
    
    """
    
    t_units, calendar = get_nc_time_units(data_path)
    
    # Load data arrays, including the ensemble member id numbers
    # (r, i, p, f), then coordinate variables (if applicable)
    # and then time and the actual variable names of diagnostic
    # d. Keep these altogther in a list (mutable):
    data = list(get_arrays([data_path], id_names + coord_names,
                [nc_time_name] + var_names))
    
    # Do not use function for this (we are manipulating
    # the data construct here -- that function is for the "final
    # form" version after such manipulation):
    j_t_data = len(id_names) + len(coord_names)
    
    data[j_t_data] = cftime_to_datetime(
        num2date(data[j_t_data], units=t_units,
                 calendar=calendar)
    )
    
    # Latitude evaluation:
    if len(coord_names) == 2 and lat_eval is not None:
        # In this case, the data indices len(id_names) + 1, 2
        # correspond to northern and southern variable
        # coordinates, so we have to evaluate but do not need to
        # include the coordinates themselves in the returned
        # data list.
        # 
        # Work out indices where to evaluate coords:
        cn = np.argmin(abs(data[len(id_names)] - lat_eval[0]))
        cs = np.argmin(abs(data[len(id_names)+1] - lat_eval[1]))
        
        lat_eval_n = data[len(id_names)][cn]
        lat_eval_s = data[len(id_names)+1][cs]
        
        if verbose:
            print(f"    evaluating at lat_n = "
                + f"{lat_eval_n:.2f}, lat_s = "
                + f"{lat_eval_s:.2f}")
        
        data[j_t_data+1] = data[j_t_data+1][:,:,cn]
        data[j_t_data+2] = data[j_t_data+2][:,:,cs]
    
    else:
        lat_eval_n = np.nan
        lat_eval_s = np.nan
    
    # Construct a new list to return:
    data_ret = [data[j].copy() for j in range(len(id_names))]
    
    if lat_eval is None and len(coord_names) == 2:
        data_ret += [data[len(id_names)+j].copy()
                     for j in range(2)]
    # Else skip any spatial coordinate data
    
    # Next is datetime and then the variables:
    data_ret += [data[len(id_names)+len(coord_names)+j].copy()
                 for j in range(3)]
    
    return data_ret, lat_eval_n, lat_eval_s
    



def _diagnostics_intersection_by_ensemble_members(data_list,
                                                  verbose=True):
    """Return sub-setted diagnostic list containing only
    ensemble members that exist for each one.
    
    
    Parameters
    ----------
    data_list : list of length n_diag of list of array with
        structure:
        
        data_list = [
            
            [ r(n_ens,), i(n_ens,), p(n_ens,), f(n_ens,),
              time(n_t), var_n(n_t, n_ens, (n_lat),),
              var_s(n_t, n_ens, [n_lat,])
            ]
            
            for d in range(n_diag)
        ]
        
        where n_ens is the number of ensemble members, n_t is
        the number of time steps, and n_lat is the length of the
        latitude dimension (usually only present if OHT or AHT
        profile). The first four arrays (indices = 0-3) are the
        ensemble member labels (ripf), the next (index = 4) is
        time (not modified by this routine), and the next (5 and
        6) are the north and south diagnostic data.
        
        If n_diag == 1, the input is not modified.
    
    
    Optional parameter
    ------------------
    verbose : bool, default = True
        Print information to console.
    
    """
    
    n_diag = len(data_list)
    
    # Create a boolean array for each variable. Here, we now
    # know the lengths of each array for each diagnostic.
    # Starting with all of these False, it doesn't matter which
    # array we reference the others to, since there's no risk of
    # missing an ensemble member (each one needs to be in every
    # array for every diagnostic). Unless, of course, there is
    # only one diagnostic in the list, in which case this
    # function should not be in use, but just set everything to
    # True and the returned data unmodified:
    if n_diag == 1:
        ens_want = [np.ones(len(data_list[d][0]), dtype=bool)
                    for d in range(n_diag)]
    else:
        ens_want = [np.zeros(len(data_list[d][0]), dtype=bool)
                    for d in range(n_diag)]
    
    # Loop over each ensemble member of the first diagnostic:
    for m0 in range(len(ens_want[0])):
        
        # Check remaining diagnostics:
        for d in range(1, n_diag):
            
            # Check realisation (r):
            if data_list[0][0][m0] in data_list[d][0]:
                
                md = np.nonzero(
                    data_list[d][0] == data_list[0][0][m0]
                )[0][0]
                
                # Now check corresponding i, p, and f values:
                if (data_list[0][1][m0] == data_list[d][1][md]
                    and data_list[0][2][m0] == data_list[d][2][md]
                    and data_list[0][3][m0] == data_list[d][3][md]
                    ):
                    
                    ens_want[0][m0] = True
                    ens_want[d][md] = True
    
    # Now we know which ensemble members to take for all
    # diagnostics. First need to figure out the time data index
    # (in case there are latitude coordinate present, which do
    # not need to be touched):
    j_t_data = _get_j_t_data(data_list[0])
    
    for d in range(n_diag):
        
        # Start with ripf values (always indices 0-3):
        for j in range(4):
            data_list[d][j] = data_list[d][j][ens_want[d]]
        
        # Just variables remaining (don't need to touch
        # datetime data or latitude coordinates if present):
        for j in range(j_t_data+1, len(data_list[0])):
            data_list[d][j] = data_list[d][j][:,ens_want[d]]
    
    return data_list



def _format_data_to_return(data, data_paths, lat_eval_n,
                           lat_eval_s,  verbose=True):
    """Convert processed data construct into 4 or 6 return
    values (dependent on whether spatial coordinates are
    present):
    
    datetime_list, ripf_list, [coords_n_list, coords_s_list],
    diags_n_list, diags_s_list
    
    """
    
    n_diag = len(data)
    
    # Needed to determine if there are coordinates: it will be
    # either 7 (no coordinates, data has been evaluated at
    # reference latitudes) or 9 (data has coordinates for the
    # northern and southern hemispheres after the ripf values
    # and before time, and there has been no evaluation at
    # reference latitudes):
    dat_size = len(data[0])
    j_t_data = _get_j_t_data(data[0])
    
    # Ensemble members are now the same per diagnostic:
    ripf = [f"r{data[0][0][m0]}i{data[0][1][m0]}"
            + f"p{data[0][2][m0]}f{data[0][3][m0]}"
            for m0 in range(len(data[0][0]))]
    
    if verbose:
        print(f"  returning {len(ripf)} ensemble members")
    # TODO: make more robust / not hardcoded?
    if dat_size == 9:
        # Data contains latitude coordinates:
        coords_n_ret = [data[d][4] for d in range(n_diag)]
        coords_s_ret = [data[d][5] for d in range(n_diag)]
    # Else there are no coordinates to return
    
    # All datetimes are the same range but in the unlikely event
    # that they are different (could be the case if monthly data
    # is loaded for one diagnostic and yearly for another),
    # save and return all datetime data for each diagnostic:
    date_ret = [data[d][j_t_data] for d in range(n_diag)]
    
    diags_ret_n = [data[d][j_t_data+1] for d in range(n_diag)]
    diags_ret_s = [data[d][j_t_data+2] for d in range(n_diag)]
    
    if dat_size == 9:
        # Return full latitudes, rather than lat_evals:
        return date_ret, ripf, coords_n_ret, coords_s_ret, \
               diags_ret_n, diags_ret_s, data_paths
    else:
        # Return latitude evals, rather than full latitudes
        # (but in the same place, so it's consistent):
        return date_ret, ripf, lat_eval_n, lat_eval_s, \
               diags_ret_n, diags_ret_s, data_paths



def one_experiment(diagnostics, model_id,
        experiment_id="piControl",
        lat_eval=(65.0, -65.0),
        datetime_lims=(dt(1, 1, 1), dt(3000, 12, 31)),
        select_diags=0, verbose=True
    ):
    """Load multiple diagnostics for a specified model_id,
    experiment_id, slicing in time and evaluating at reference
    latitudes if applicable.
    
    
    Parameters
    ----------
    diagnostics : list of str of length n_diag >= 1
        Names of diagnostics to load. These can either be full
        diagnostic names, their aliases, or a keyword which,
        with the model_id, is used to identify the required
        diagnostic to load. The list can contain a mixture of
        such identifiers. See api.model_diagnostics for
        explanation.
    
    model_id : str, model_id
    
    
    Optional parameters
    -------------------
    experiment_id : str, experiment_id (default = "piControl")
    
    lat_eval : tuple (float, float)
        Northern and southern hemisphere reference latitude to
        evaluate at (if required). Both are interpreted as
        degrees_north. Default is (65.0, -65.0). Can also be
        None, in which case no evaluation is done.
    
    datetime_lims : tuple (datetime.datetime, datetime.datetime)
        Start and end of datetime range for time slicing.
        Default is an extensive range that (usually) results in
        all data being returned (no time slicing).
    
    select_diags : int or list of int, default = 0
        If a diagnostic specified in parameter diagnostics is
        a keyword, this can be used to switch between multiple
        choices if a given model_id has them for the specified
        keyword (see api.model_diagnostics -- 0 is the default,
        highest-priority for each model_id and keyword).
        
        If an int, the same index is used for all diagnostics.
        Alternatively, can be a list of up to length n_diag so
        that a different index is used per diagnostic.
    
    verbose : bool, default = True
        Print information to console.
    
    
    Returns
    -------
    date_ret : list of length n_diag or array (nt,)
        Datetime arrays for each diagnostic.
    
    ripf : list of length n_ens of str
        Ensemble member labels (e.g., "r1i1p1f1")
    
    diags_ret_n, diags_ret_s : list of length n_diag of arrays
                               (nt, n_ens, [n_lat])
        Data arrays (north, south) for each diagnostic.
    
    data_paths : list of length n_diag of pathlib.Path
        Full data paths that have been loaded for each
        diagnostic.
    
    lat_eval_n, lat_eval_s : array(n_diag,) of float
        Actual latitudes data have been evaluated at (north,
        south). Set to NaN if no evaluation for the
        corresponding diagnostic.
    
    
    Notes
    -----
    [1] There is no exception handling in case of data not being
        found. A FileNotFoundError is raised if the data does
        not exist or a KeyError if the specified keyword/alias/
        diagnostic is not defined or available for model_id.
        Such exceptions are raised by api.model_diagnostics.
        This is by design, as it may be desirable to "try"
        loading a diagnostic anyway and handle these cases by
        non-fatal means (e.g., in quality control plotting
        scripts).
    
    """
    
    if verbose:
        print(f"Loading {experiment_id} data for {model_id}")
    
    n_diag = len(diagnostics)
    
    diagnostics, coord_names, var_names = \
        mdiags.get_load_info(diagnostics, model_id,
                             experiment_id, select_diags)
    
    # Create a list to append all diagnostics. The lengths/
    # shapes of the arrays are not known before loading:
    data = []
    data_paths = []
    
    # Actual latitudes evaluated at (set to NaN for diagnostics
    # which do not have a latitude evaluation):
    lat_eval_n = np.full(n_diag, np.nan, dtype=np.float64)
    lat_eval_s = np.full(n_diag, np.nan, dtype=np.float64)
    
    for d in range(n_diag):
        
        if verbose:
            print(f"  ({d+1}/{n_diag}) {diagnostics[d]}")
        
        # Load data arrays, including the ensemble member
        # id numbers (r, i, p, f), then coordinate variables
        # and then time and the actual variable names of
        # diagnostic d. Keep these altogther in the tuple
        # data_e1_d (data for experiment 1, diagnostic d):
        data_paths.append(nc_file_path(diagnostics[d],
            experiment_id, model_id))
        
        data_d, lat_eval_n[d], lat_eval_s[d] = _load_diag_data(
            data_paths[d], coord_names[d], var_names[d],
            lat_eval=lat_eval, verbose=verbose)
        
        data_d = _time_slice(data_d, verbose=verbose,
                             datetime_lims=datetime_lims)
        
        # Append to the main data list:
        data.append(data_d)
    
    # All diagnostics are loaded and sliced in time and space,
    # but still need to sub-set the intersection of ensemble
    # members per diagnostic:
    data = _diagnostics_intersection_by_ensemble_members(data,
        verbose=verbose)
    
    return _format_data_to_return(data, data_paths, lat_eval_n,
                                  lat_eval_s, verbose=verbose)



def two_experiments(diagnostics, model_id,
        experiment_ids=("historical", "ssp370"),
        lat_eval=(65.0, -65.0),
        datetime_lims=(dt(1850, 1, 1), dt(2100, 12, 31)),
        select_diags=0, verbose=True
    ):
    """Load multiple diagnostics for a specified model_id for
    two experiments, concatenating in time, slicing to the
    specified time range, and evaluating at reference
    latitudes if applicable.
    
    
    Parameters
    ----------
    diagnostics : list of str of length n_diag >= 1
        Names of diagnostics to load. These can either be full
        diagnostic names, their aliases, or a keyword which,
        with the model_id, is used to identify the required
        diagnostic to load. The list can contain a mixture of
        such identifiers. See api.model_diagnostics for
        explanation.
    
    model_id : str, model_id
    
    
    Optional parameters
    -------------------
    experiment_ids : tuple (str, str)
        Experiment ids, default = ("historical", "ssp370"). Note
        there is no checking that two experiments make sense to
        combine or their order [e.g., ("ssp370", "piControl")
        will not raise any exceptions].
    
    lat_eval : tuple (float, float)
        Northern and southern hemisphere reference latitude to
        evaluate at (if required). Both are interpreted as
        degrees_north. Default is (65.0, -65.0). Can also be
        None, in which case no evaluation is done.
    
    datetime_lims : tuple (datetime.datetime, datetime.datetime)
        Start and end of datetime range for time slicing.
        Default is an extensive range that (usually) results in
        all data being returned (no time slicing).
    
    select_diags : int or list of int, default = 0
        If a diagnostic specified in parameter diagnostics is
        a keyword, this can be used to switch between multiple
        choices if a given model_id has them for the specified
        keyword (see api.model_diagnostics -- 0 is the default,
        highest-priority for each model_id and keyword).
        
        If an int, the same index is used for all diagnostics.
        Alternatively, can be a list of up to length n_diag so
        that a different index is used per diagnostic.
    
    verbose : bool, default = True
        Print information to console.
    
    
    Returns
    -------
    date_ret : list of length n_diag or array (nt,)
        Datetime arrays for each diagnostic.
    
    ripf : list of length n_ens of str
        Ensemble member labels (e.g., "r1i1p1f1")
    
    diags_ret_n, diags_ret_s : list of length n_diag of arrays
                               (nt, n_ens, [n_lat])
        Data arrays (north, south) for each diagnostic.
    
    data_paths : list of length n_diag of [pathlib.Path]*2
        Full data paths that have been loaded for each
        diagnostic, and each experiment.
    
    lat_eval_n, lat_eval_s : array(n_diag, 2) of float
        Actual latitudes data have been evaluated at (north,
        south), for each experiment. Set to NaN if no evaluation
        for the corresponding diagnostic.
    
    
    Notes
    -----
    [1] There is no exception handling in case of data not being
        found. A FileNotFoundError is raised if the data does
        not exist or a KeyError if the specified keyword/alias/
        diagnostic is not defined or available for model_id.
        Such exceptions are raised by api.model_diagnostics.
        This is by design, as it may be desirable to "try"
        loading a diagnostic anyway and handle these cases by
        non-fatal means (e.g., in quality control plotting
        scripts).
    
    """
    
    if verbose:
        print(f"Loading {' + '.join(experiment_ids)} data for "
            + f"model {model_id}")
    
    n_diag = len(diagnostics)
    
    # Currently in api.model_diagnostics, the routine called
    # get_load_info does not utilise the experiment_id argument
    # (it is included in an attempt at "future proofing" the
    # code in case it becomes necessary to have a different list
    # of aliases per model AND per experiment -- not currently
    # the case). So here, for now, just give the first
    # experiment:
    diagnostics, coord_names, var_names = \
        mdiags.get_load_info(diagnostics, model_id,
                             experiment_ids[0], select_diags)
    
    # Create a list to append all diagnostics (for both
    # experiments concatenated along time). The lengths/shapes
    # of the arrays are not known before loading:
    data = []
    data_paths = []
    
    # Actual latitudes evaluated at (set to NaN for diagnostics
    # which do not have a latitude evaluation, and determine
    # for each experiment):
    lat_eval_n = np.full((n_diag, 2), np.nan, dtype=np.float64)
    lat_eval_s = np.full((n_diag, 2), np.nan, dtype=np.float64)
    
    for d in range(n_diag):
        
        if verbose:
            print(f"  ({d+1}/{n_diag}) {diagnostics[d]}")
        
        # Load data arrays, including the ensemble member
        # id numbers (r, i, p, f), then coordinate variables
        # and then time and the actual variable names of
        # diagnostic d. Keep these altogther in the tuple
        # data_x1_d (data for experiment 1, diagnostic d):
        data_paths.append([
            nc_file_path(diagnostics[d], exp, model_id)
            for exp in experiment_ids])
        
        data_x1_d, lat_eval_n[d,0], lat_eval_s[d,0] = \
            _load_diag_data(data_paths[d][0], coord_names[d],
                            var_names[d], lat_eval=lat_eval,
                            verbose=verbose)
        
        # and similarly for experiment 2 (in general, a
        # different number of ensemble members):
        data_x2_d, lat_eval_n[d,1], lat_eval_s[d,1] = \
            _load_diag_data(data_paths[d][1], coord_names[d],
                            var_names[d], lat_eval=lat_eval,
                            verbose=verbose)
        
        # Determine and match the ensemble members loaded for
        # each experiment. Create lists to append matches to:
        ens_want_x1_d = []
        ens_want_x2_d = []
        
        # Now loop over ensemble members of experiment 1 (m1)
        # and check whether the same ensemble member id exists
        # in the second experiment (it shouldn't matter if
        # experiment 1 and 2 are switched here):
        for m1 in range(len(data_x1_d[0])):
            
            # first check realisation (r, index 0).
            # In data_x1_d[0][m1], the [0] index accesses
            # the realisation_id of experiment 1 data, then
            # [m1] index accesses the m1-th ensemble member:
            if data_x1_d[0][m1] in data_x2_d[0]:
                
                # found a matching realisation_id (we assume
                # there's only 1!): find its index in the
                # second experiment array, then cross check
                # the other id numbers:
                m2 = np.nonzero(
                    data_x2_d[0] == data_x1_d[0][m1])[0][0]
                
                # data indices 1-3 correspond to i, p, f:
                if (data_x1_d[1][m1] == data_x2_d[1][m2]
                    and data_x1_d[2][m1] == data_x2_d[2][m2]
                    and data_x1_d[3][m1] == data_x2_d[3][m2]
                    ):
                    ens_want_x1_d.append(m1)
                    ens_want_x2_d.append(m2)
        
        # Convert lists of indices into arrays (to obtain
        # boolean indices when sorted):
        ens_want_x1_d = np.array(ens_want_x1_d)
        ens_want_x2_d = np.array(ens_want_x2_d)
        
        # Sort by realisation index of data 1 (again, should
        # not matter if experiments are switched):
        asort = np.argsort(data_x1_d[0][ens_want_x1_d])
        ens_want_x1_d = ens_want_x1_d[asort]
        ens_want_x2_d = ens_want_x2_d[asort]
        
        # Create the combined-experiment data array. First set
        # the intersected ensemble members (doesn't matter which
        # one):
        data_d = [data_x1_d[j][ens_want_x1_d] for j in range(4)]
        
        # Add coordinates, if they exist (assume these are the
        # same for both experiments -- they are, at least, the
        # same size, otherwise we would not get this far):
        j_t_data = _get_j_t_data(data_x1_d)
        if j_t_data == 6:
            # north and south latitudes, respectively
            data_d += [data_x1_d[4], data_x1_d[5]]
        
        # Next is datetime (which requires concatenation along
        # time but no ensemble member selection):
        data_d += [np.concatenate((data_x1_d[j_t_data],
                                   data_x2_d[j_t_data]),
                                  axis=0)]
        
        # Next and finally are variables north and south (which
        # require concatenation along time and ensemble member
        # selection):
        data_d += [np.concatenate(
            (data_x1_d[j_t_data+j][:,ens_want_x1_d],
             data_x2_d[j_t_data+j][:,ens_want_x2_d]), axis=0)
            for j in [1,2]]
        
        # Get time slice (could probably be done earlier, but
        # seems safer to do it last):
        data_d = _time_slice(data_d, verbose=verbose,
                             datetime_lims=datetime_lims)
        
        # Append processed data for this diagnostic, then loop
        # will repeat for remaining diagnostics:
        data.append(data_d)
    
    # All diagnostics are loaded and sliced in time and space,
    # but still need to sub-set the intersection of ensemble
    # members per diagnostic:
    data = _diagnostics_intersection_by_ensemble_members(data,
        verbose=verbose)
    
    return _format_data_to_return(data, data_paths, lat_eval_n,
                                  lat_eval_s, verbose=verbose)

