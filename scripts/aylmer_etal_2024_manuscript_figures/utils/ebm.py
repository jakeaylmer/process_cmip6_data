"""Calculation of energy-balance model (EBM) parameters and
estimates of sensitivities dphi_i/dOHT, etc. See Aylmer et al.
(2024) for details:
    
[1] Aylmer, J. R., D. Ferreira, and D. L. Feltham, 2024: Impact
    of ocean heat transport on sea ice capture by a simple
    energy balance model, Commun. Earth Environ., accepted in
    principle May 2024.
"""

import numpy as np

from . import maths, model_lists
from .script_tools import (
    load_data_multi_model_one_hemi_two_ref_lats as ldmmohtrl)

default_ref_lats = (65.0, 90.0)  # for north
earth_radius_m   = 6.371E6

def default_calc_param_func(x, y):
    """Must take two arrays x and y, and return a regression
    slope (y on x) and its error (floats). Default is to use
    orthogonal distance regression.
    """
    return maths.orthogonal_distance_regression(x, y)[2:]

# EBM parameter values: either fixed values (value, error), the
# tuple (None, None) to calculate from input historical/future
# data, or the string "piControl" to calculate as fixed values
# from pre-industrial control run data:
default_ebm_params = {
    "bc"  : (None, None),  # i.e., calculate from input data
    "beta": "piControl",   # i.e., calculate from piControl
    "bup" : "piControl",   # i.e., calculate from piControl
    "h"   : (None, None),  # i.e., calculate from input data
    "S"   : "piControl"}   # i.e., calculate from piControl


# Diagnostics for each EBM parameter for regression: [x, y].
# These are "keywords" defined by the "load data API":
ebm_param_diagnostics = {
    "S"   : ["iel", "f_sw_surf"],
    "bc"  : ["oht", "aht"],
    "beta": ["f_down", "f_olr"],
    "bup" : ["tas", "f_up"],
    "h"   : ["oht", "dhdt"]}

# Units of parameters, for reporting values:
ebm_param_units = {
    "S"   : "W m-2 {}-1",  # format string for "N" or "S"
    "bc"  : "",
    "beta": "",
    "bup" : "W m-2 K-1",
    "h"   : ""}

# Default format strings (precision) for value reporting:
default_fmt_str = {
    "S"   : ".2f",
    "bc"  : ".2f",
    "beta": ".3f",
    "bup" : ".2f",
    "h"   : ".2f"}

# Format string for correlation coefficients:
default_fmt_str_corr = ".2f"


def _print_param_value(param, value, corr=None,
        description_before="", description_after="",
        hemi="n", val_fmt=None, override_units=None):
    """Print parameter (or some other) value and error. Sub-
    routine for other functions in module.
    
    
    Parameters
    ----------
    param : str
        Parameter name.
    
    value : tuple of float (value, error)
    
    
    Optional parameters
    -------------------
    corr : float or None (default)
        Correlation value of data underlying parameter.
    
    description_before, description_after: str, default = ""
        Text(s) to display before and/or after parameter value.
    
    hemi : str "n" (default) or "s"
        Which hemisphere is being analysed.
    
    val_fmt : str or None (default)
        Format string for parameter value/error. If None, uses
        the defaults defined in module.
    
    override_units : str or None (default)
        Units of parameter value. If None, uses the defaults
        defined in module.
    
    """
    
    # Determine units string:
    if override_units is None:
        if param == "S":
            units = ebm_param_units["S"].format(hemi.upper())
        else:
            units = ebm_param_units[param]
    else:
        units = override_units
    
    # For formatting values:
    fmt = default_fmt_str[param] if val_fmt is None else val_fmt
    
    x = (f"{param}: {description_before} "
        + ("(" if units != "" else "")
        + f"{value[0]:{fmt}} +/- {abs(value[1]):{fmt}}"
        + (") " if units != "" else " ") + units)
    
    if corr is not None:
        x += f"; r = {corr:{default_fmt_str_corr}}"
    
    x += f" {description_after}"
    
    print(" ".join(x.split()))
    # [ End function _print_param_value(); nothing returned ]



def area_factor(ref_lats, radius=earth_radius_m):
    """Surface area of Earth between two reference latitudes in
    units of [radius]**2 (with default: m2).
    """
    return 2.0*np.pi*(radius**2)*abs(
        np.sin(abs(ref_lats[1])*np.pi/180.0)
        - np.sin(abs(ref_lats[0])*np.pi/180.0))



def radiation_coefficient(ebm_params):
    """Compute the 'radiation coefficient', and its error, from
    EBM parameter values. In the EBM equation in [1], this is
    called R and appears in the first term on the right-hand
    side, multiplying delta_tas:
    
        R = bup*beta / (1 + beta)
    
    """
    R = ebm_params["bup"][0]*ebm_params["beta"][0]
    R /= 1.0 + ebm_params["beta"][0]
    
    # Assume bup and beta are never zero:
    R_err = R*np.sqrt(
        (ebm_params["bup"][1]/ebm_params["bup"][0])**2 +
        (ebm_params["beta"][1]/ebm_params["beta"][0])**2)
    
    return R, R_err



def ocean_correction_factor(ebm_params, A):
    """Compute the 'ocean correction factor', and its error,
    from EBM parameter values and area factor. In the EBM
    equation in [1], this is called C and appears in the second
    term on the right-hand side, multiplying delta_OHT:
    
        C = (1/A)*(1 + bc/(1+beta) - h)
    
    Here, input A (area factor) should include any units change
    for OHT.
    """
    
    C = (1.0
         + ebm_params["bc"][0]/(1.0 + ebm_params["beta"][0])
         - ebm_params["h"][0]) / A
    
    # Here, can't neglect that bc and/or h could be zero. If
    # that's the case, then the EBM isn't going to work anyway,
    # but best to handle the divide-by-zero here (by simply
    # ignoring it) so we can at least see the result through):
    # 
    C_err_root = 0.0  # bit inside the square root (to add to)
    
    for x in ["bc", "h"]:
        if ebm_params[x][0] != 0.0:
            C_err_root += (ebm_params[x][1]/ebm_params[x][0])**2
    
    # Error from (1 + beta) part of bc term:
    C_err_root += (ebm_params["beta"][1] /
                   (1.0 + ebm_params["beta"][0]))**2
    
    C_err = C*np.sqrt(C_err_root)
    
    return C, C_err



def get_piControl_parameter(param, ref_lats=default_ref_lats,
        hemi="n", verbosity=1, n_year_averages=21,
        regression_function=maths.orthogonal_distance_regression):
    """Derive pre-industrial control value of EBM parameter
    param (str).
    
    
    Optional parameters
    -------------------
    ref_lats : tuple of float (ref_lat_1, ref_lat_2)
        Reference latitude(s) in degrees_north (with
        ref_lat_2 >= 85 effectively being one reference
        latitude, ref_lat_1).
    
    hemi : str "n" (default) or "s"
        Which hemisphere is being analysed.
    
    verbosity: int, default = 1
        How much information on progress to print to console.
        With 1, prints parameter value only.
    
    n_year_averages : int, default = 21
        Number of years making up contiguous averages taken from
        pre-industrial control runs.
    
    regression_function : function
        Function of two input arrays that returns a regression
        slope estimate of the second input on the first and the
        error of this estimate.
    
    
    Returns
    -------
    pval : float
        Parameter value derived from piControl
    
    pval_err : float
        Error/uncertainty on estimate pval
    
    corr : float
        Correlation of input datasets underlying derivation of
        pval and pval_err.
    
    """
    
    # Load required pre-industrial control data:
    pi_data = ldmmohtrl(model_lists.by_experiment["piControl"],
        ebm_param_diagnostics[param], hemi=hemi,
        ref_lats=ref_lats, n_year_averages=n_year_averages,
        experiment_ids=["piControl"],
        verbosity=0 if verbosity < 2 else verbosity)
    
    pval, pval_error = regression_function(
        np.concatenate(pi_data[ebm_param_diagnostics[param][0]]),
        np.concatenate(pi_data[ebm_param_diagnostics[param][1]]))
    
    corr = maths.correlation_coefficient(
        np.concatenate(pi_data[ebm_param_diagnostics[param][0]]),
        np.concatenate(pi_data[ebm_param_diagnostics[param][1]]))
    
    if verbosity > 0:
        _print_param_value(param, (pval, pval_error), corr=corr,
            description_before="piControl-derived", hemi=hemi)
    
    return pval, pval_error, corr



def fit(data_changes,
        ebm_params=default_ebm_params,
        calc_param_func=default_calc_param_func,
        ref_lats=default_ref_lats, hemi="n",
        piparams_n_year_averages=21):
    """Determines EBM estimate of the slope of ice edge latitude
    against OHT and against surface temperature. Also computes
    the actual slopes from the data and EBM parameters. This is
    the main function that should be called externally.
    
    
    Parameters
    ----------
    data_changes: dict of list of arrays
        The CMIP6 data for the historical/future changes in ice
        edge latitude, OHT, etc., as loaded by the script_tools
        module; e.g., data_changes["iel"][m][:] is the change
        in sea ice edge latitude over some time period for model
        m (all its ensemble members).
    
    
    Optional parameters
    -------------------
    ebm_params : dict of tuple of float
        Values of EBM parameters (keys) "S", "beta", "bup",
        "bc", and "h", (value, error) in the right units
        assumed. Alternatively, entries can be (None, None) to
        calculate them from input data, or "piControl" to
        calculate them from pre-industrial control simulation
        data. Default is to calculate them 'on the fly' from
        piControl data (S, beta, bup) and input data (bc, h).
    
    calc_param_func : function
        Function of two input arrays that returns a regression
        slope estimate of the second input on the first and the
        error of this estimate. This is used to calculate
        parameter values and data sensitivities (e.g., ice-edge
        latitude against OHT).
    
    ref_lats : tuple of float
        Reference latitude(s) in degrees_north (with
        ref_lat_2 >= 85 effectively being one reference
        latitude, ref_lat_1).
    
    hemi : str "n" (default) or "s"
        Which hemisphere is being analysed.
    
    piparams_n_year_averages : int, default = 21
        Number of years making up contiguous averages taken from
        pre-industrial control runs where those are being used
        to calculate EBM parameter values.
    
    
    Returns
    -------
    ebm_params : dict {"param": (value, error)}
        The EBM parameter values.
    
    ebm_results : dict {"result": (value, error)}
        The EBM fitting results; specifically, the keys:
            
            "diel/dtas_data"
            "corr_tas_iel"
            "diel/dtas_ebm"
            "diel/doht_data"
            "corr_oht_iel"
            "diel/doht_ebm"
    """
    
    # Keyword arguments passed to piControl calculation routine,
    # common for all parameters, in case this is needed in the
    # loop below:
    pi_kw = {"ref_lats": ref_lats, "hemi": hemi,
             "verbosity": 1,
             "n_year_averages": piparams_n_year_averages,
             "regression_function": calc_param_func}
    
    # Loop over required EBM parameters, which are the keys of
    # default_ebm_params, and determine the values as
    # appropriate/required:
    for p in default_ebm_params.keys():
        
        if (type(ebm_params[p]) == str
                and ebm_params[p] == "piControl"):
            
            # Calculate piControl values now.
            # 
            # NOTE: this uses the piControl model set, which is
            # not necessarily the exact same model set as input
            # data.
            ebm_params[p] = \
                get_piControl_parameter(p, **pi_kw)[:2]
        
        elif type(ebm_params[p]) in [tuple, list, np.ndarray]:
            
            if (ebm_params[p][0] is None
                    and ebm_params[p][1] is None):
                
                # Calculate from input historical or future
                # change data:
                ebm_params[p] = calc_param_func(
                    np.concatenate(data_changes
                        [ebm_param_diagnostics[p][0]]),
                    np.concatenate(data_changes
                        [ebm_param_diagnostics[p][1]]))
                
                corr = maths.correlation_coefficient(
                    np.concatenate(data_changes
                        [ebm_param_diagnostics[p][0]]),
                    np.concatenate(data_changes
                        [ebm_param_diagnostics[p][1]]))
                
                _print_param_value(p, ebm_params[p], corr=corr,
                    description_before="historical/future data",
                    hemi=hemi)
                
            else:
                # Do nothing; input values provided:
                _print_param_value(p, ebm_params[p],
                    description_before="input value", hemi=hemi)
        
        else:
            
            raise TypeError("Invalid input type for "
                            + "ebm_params[\"" + f"{p}" + "\"]: "
                            + f"{type(ebm_params[p])}")
        
        # [ End if type(ebm_params) == str ]
    # [ End for p in default_ebm_params.keys() ]
    
    
    # Use the same EBM parameter routines for the slope of
    # delta_OHT against delta_tas. This is needed for the EBM
    # estimate of delta_iel against delta_tas.
    # 
    # (At some point, we wondered if the PI control data could
    # determine this; it can't, so now this is always calculated
    # in this way from the input historical/future data.)
    # 
    # NOTE: OHT is either the OHT evaluated at ref_lats[0] if
    # ref_lats[1] is the pole, otherwise it is the difference
    # between OHT at those two latitudes (i.e., convergence
    # multiplied by area). Same goes for AHT and dhdt.
    # 
    doht_dtas = calc_param_func(
        np.concatenate(data_changes["tas"]),
        np.concatenate(data_changes["oht"]))
    
    r_doht_dtas = maths.correlation_coefficient(
        np.concatenate(data_changes["tas"]),
        np.concatenate(data_changes["oht"]))
    
    _print_param_value("doht/dtas", doht_dtas, corr=r_doht_dtas,
        description_before="historical/future data",
        description_after="", hemi=hemi, val_fmt=".2f",
        override_units="TW K-1")
    
    # Now do it the other way around (slope of delta_tas against
    # delta_OHT), for getting the EBM estimate of delta_iel
    # against delta_oht. If using ODR regression (default), this
    # is just one divided by doht_dtas, anyway:
    dtas_doht = calc_param_func(
        np.concatenate(data_changes["oht"]),
        np.concatenate(data_changes["tas"]))
    
    r_dtas_doht = maths.correlation_coefficient(
        np.concatenate(data_changes["oht"]),
        np.concatenate(data_changes["tas"]))
    
    _print_param_value("dtas/doht",
        100.0*np.array(dtas_doht), corr=r_dtas_doht,
        description_before="historical/future data",
        description_after="", hemi=hemi, val_fmt=".3f",
        override_units="K per 100 TW")
    
    # A is the area of the spherical cap between the lower
    # reference latitude, ref_lats[0], and the pole or higher
    # reference latitude, ref_lats[1], if it is less than 85
    # degrees_north or degrees_south. It also includes factors
    # to account for units of OHT (TW):
    A = area_factor(ref_lats) / 1.0E12
    
    # Compute the coefficients C and R from the EBM equation:
    R = radiation_coefficient(ebm_params)
    C = ocean_correction_factor(ebm_params, A)
    
    # Now calculate the EBM estimate of the slope delta_iel
    # against delta_tas:
    ebm_diel_dtas = (R[0]-C[0]*doht_dtas[0])/ebm_params["S"][0]
    
    ebm_diel_dtas_err = ebm_diel_dtas*(
        np.sqrt(
            (0.0 if R[0]==0.0 else R[1]/R[0])**2 +
            (0.0 if ebm_params["S"][0]==0.0
             else ebm_params["S"][1]/ebm_params["S"][0])**2
        )
        + np.sqrt(
            (0.0 if C[0]==0.0 else C[1]/C[0])**2 +
            (0.0 if ebm_params["S"][0]==0.0
             else ebm_params["S"][1]/ebm_params["S"][0])**2 +
            (0.0 if doht_dtas[0]==0.0
             else doht_dtas[1]/doht_dtas[0])**2
        )
    )
    
    # Repeat for the EBM estimate of the slope diel/doht
    # (i.e., dividing the EBM equation through by dOHT):
    ebm_diel_doht = (R[0]*dtas_doht[0]-C[0])/ebm_params["S"][0]
    
    ebm_diel_doht_err = ebm_diel_doht*(
        np.sqrt(
            (0.0 if R[0]==0.0 else R[1]/R[0])**2 +
            (0.0 if ebm_params["S"]==0.0
             else ebm_params["S"][1]/ebm_params["S"][0])**2 +
            (0.0 if dtas_doht[0]==0.0
             else dtas_doht[1]/dtas_doht[0])**2
        )
        + np.sqrt(
            (0.0 if C[0]==0.0 else C[1]/C[0])**2 +
            (0.0 if ebm_params["S"]==0.0
             else ebm_params["S"][1]/ebm_params["S"][0])**2
        )
    )
    
    # Combine results into a dictionary to return:
    ebm_results = {
        "diel/dtas_data": calc_param_func(
            np.concatenate(data_changes["tas"]),
            np.concatenate(data_changes["iel"])),
        "corr_tas_iel": maths.correlation_coefficient(
            np.concatenate(data_changes["tas"]),
            np.concatenate(data_changes["iel"])),
        "diel/dtas_ebm": (ebm_diel_dtas, ebm_diel_dtas_err),
        "diel/doht_data": calc_param_func(
            np.concatenate(data_changes["oht"]),
            np.concatenate(data_changes["iel"])),
        "corr_oht_iel": maths.correlation_coefficient(
            np.concatenate(data_changes["oht"]),
            np.concatenate(data_changes["iel"])),
        "diel/doht_ebm": (ebm_diel_doht, ebm_diel_doht_err)
    }
    
    # Print EBM fitting results/summary:
    _print_param_value("diel/dtas",
        ebm_results["diel/dtas_data"],
        corr=ebm_results["corr_tas_iel"],
        description_before="historical/future data",
        description_after="", hemi=hemi, val_fmt=".3f",
        override_units=f"{hemi.upper()} K-1")
    
    _print_param_value("diel/dtas",
        ebm_results["diel/dtas_ebm"],
        corr=None, description_before="EBM fit",
        description_after="", hemi=hemi, val_fmt=".3f",
        override_units=f"{hemi.upper()} K-1")
    
    _print_param_value("diel/doht",
        100.0*np.array(ebm_results["diel/doht_data"]),
        corr=ebm_results["corr_oht_iel"],
        description_before="historical/future data",
        description_after="", hemi=hemi, val_fmt=".3f",
        override_units=f"{hemi.upper()} per 100 TW")
    
    _print_param_value("diel/doht",
        100.0*np.array(ebm_results["diel/doht_ebm"]),
        corr=None, description_before="EBM fit",
        description_after="", hemi=hemi, val_fmt=".3f",
        override_units=f"{hemi.upper()} per 100 TW")
    
    return ebm_params, ebm_results
