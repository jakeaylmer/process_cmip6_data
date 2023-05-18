"""Provide data in the form of dictionaries listing what
diagnostics are available for each model, using keywords and
aliases as two levels of abstraction for the true diagnostic
names. These are defined in _diagnostic_definitions.

Keywords are used along with model_id to automatically determine
which diagnostic to load. For example, keyword = "dhdt" and
model_id = "HadGEM3-GC31-LL" will load heat content tendency
data for that model, but there are multiple possible diagnostics
corresponding to dhdt. This code will determine which one to
load for that model. More than one diagnostic can therefore have
the same keyword.

Aliases are just short names for the true diagnostic names, and
are unique per diagnostic. For example, heat content tendency
computed from the residual of hfds and OHT (itself computed from
hfx and hfy) might be called:

opottemptend_ver_int_from_hfds_hfx_hfy_hor_int_yr_cc_approx

which is fully descriptive but awkward to use, so it has the
alias (short name):

dhdt_pot_residual_hfx

This data can be loaded by specifying the full name, the alias,
or just "dhdt" in this particular example if that is the data
available for HadGEM3-GC31-LL.

"""

from api import _diagnostic_definitions as ddef
from src import metadata as md


print_defined_diagnostics = ddef.print_defined_diagnostics


# ============================================================ #
# Dictionary structure:
# ---------------------
# 
#     "keyword" : { "model_id": [alias_1, ..., alias_n] }
# 
# where
# 
#     keyword  : short "code-name" for basic diagnostic type
#     model_id : full name of model/source_id
#     alias_n  : alias for full diagnostic name; if multiple are
#                available, they should be listed in order of
#                priority starting at highest (default). The
#                get_diagnostic() routine has an option "select"
#                (0 by default) which indexes the alias list.
# 
# Note that keywords and aliases are defined in the
# _diagnostic_definitions module.
# 
# ============================================================ #
diagnostic_per_model = {
    "iel"      : {
        "CESM2-FV2"    : [f"iel_{x}deg_bil" for x in ["1", "05", "025", "2", "4"]],
        "IPSL-CM6A-LR" : [f"iel_{x}deg_bil" for x in ["1", "05", "025", "2", "4"]]
    },
    "oht"      : {
        "CanESM5-CanOE": ["oht_residual_pot", "oht_from_hfbasin"],
        "CESM2-FV2"    : ["oht_residual_pot"],
        "CNRM-ESM2-1"  : ["oht_from_hfx_hfy", "oht_from_hfbasin"],
        "IPSL-CM6A-LR" : ["oht_residual_con", "oht_from_hfbasin", "oht_from_hfx_hfy"]
    },
    "dhdt"    : {
        "CanESM5-CanOE": ["dhdt_pot_direct", "dhdt_pot_residual_hfbasin"],
        "CESM2-FV2"    : ["dhdt_pot_direct", "dhdt_pot_residual_hfbasin"],
        "CNRM-ESM2-1"  : ["dhdt_pot_residual_hfx", "dhdt_pot_residual_hfbasin"],
        "IPSL-CM6A-LR" : ["dhdt_con_direct", "dhdt_con_residual_hfbasin", "dhdt_con_residual_hfx"]
    },
    "aht"      : dict.fromkeys(md.defined_models, ["aht_from_net_flux_gn_interp"]),
    "f_down"   : dict.fromkeys(md.defined_models, ["f_down_gn_interp"]),
    "f_olr"    : dict.fromkeys(md.defined_models, ["f_olr_gn_interp"]),
    "f_sw_surf": dict.fromkeys(md.defined_models, ["f_sw_surf_gn_interp"]),
    "f_sw_toa" : dict.fromkeys(md.defined_models, ["f_sw_toa_gn_interp"]),
    "f_up"     : dict.fromkeys(md.defined_models, ["f_up_gn_interp"]),
    "hfds"     : dict.fromkeys(md.defined_models, ["hfds"]),
    "hfds_wm-2": dict.fromkeys(md.defined_models, ["hfds_wm-2"]),
    "sia"      : dict.fromkeys(md.defined_models, ["sia"]),
    "sia_mon"  : dict.fromkeys(md.defined_models, ["sia_mon"]),
    "sie"      : dict.fromkeys(md.defined_models, ["sie"]),
    "sie_mon"  : dict.fromkeys(md.defined_models, ["sie_mon"]),
    "tas"      : dict.fromkeys(md.defined_models, ["tas_gn_interp"])
}

diagnostic_per_model["dhdt_wm-2"] = {}
for model in diagnostic_per_model["dhdt"].keys():
    diagnostic_per_model["dhdt_wm-2"][model] = [
        x + "_wm-2" for x in diagnostic_per_model["dhdt"][model]
    ]

diagnostic_per_model["iel_mon"] = {}
for model in diagnostic_per_model["iel"].keys():
    diagnostic_per_model["iel_mon"][model] = [
        x + "_mon" for x in diagnostic_per_model["iel"][model]
    ]


# ============================================================ #
# Quick validation on import (saves such errors getting lost
# too far down the call stack...):
# ------------------------------------------------------------ #
for keyword in diagnostic_per_model.keys():
    
    if keyword not in ddef.defined_keywords:
        ddef.print_defined_diagnostics()
        raise ValueError(f"Check {__name__}: keyword " + "\""
            + keyword + "\" not defined")
    
    for model_id in diagnostic_per_model[keyword].keys():
        
        if model_id not in md.defined_models:
            ddef.print_defined_diagnostics()
            raise ValueError(f"Check {__name__}: model_id "
                + "\"" + model_id + "\" not defined (found in "
                + "keyword:\"" + keyword + "\")")
        
        for alias in diagnostic_per_model[keyword][model_id]:
            
            if alias not in ddef.defined_aliases:
                ddef.print_defined_diagnostics()
                raise ValueError(f"Check {__name__}: alias "
                    + "\"" + alias + "\" not defined (found in "
                    + "keyword:\"" + keyword + "\", model:\""
                    + model_id + "\")")
# ============================================================ #


def get_valid_aliases_from_keyword(keyword):
    """List of all valid alises for a specified keyword."""
    return sorted([a for a in ddef.alias_keywords.keys()
                  if ddef.alias_keywords[a] == keyword])


# Provide main API functions: give me a keyword ("oht", "iel",
# etc.) and a model_id ("CESM2", "HadGEM3-GC31-LL", etc.) and I
# will tell you what diagnostic to load, what the coordinate
# variables are, what the nc variable names are, or all at
# once, with the option of specifying a non-default diagnostic
# in the case of multiple choices for the specified model_id
# via a "alias_select" index parameter
# ------------------------------------------------------------ #
def get_diagnostic_from_keyword(keyword, model_id,
                                alias_select=0):
    """Full diagnostic name from input keyword and model_id
    (optionally select non-default diagnostic via alias_select).
    """
    return ddef.alias_diags[
        diagnostic_per_model[keyword][model_id][
            alias_select % len(
                diagnostic_per_model[keyword][model_id])
        ]
    ]


def get_coords_from_keyword(keyword, model_id, alias_select=0):
    """Diagnostic coordinates from a specified keyword and
    model_id (optionally select non-default diagnostic via
    alias_select).
    """
    return ddef.alias_coords[
        diagnostic_per_model[keyword][model_id][
            alias_select % len(
                diagnostic_per_model[keyword][model_id])
        ]
    ]


def get_nc_vars_from_keyword(keyword, model_id, alias_select=0):
    """Diagnostic variable names (north, south) from a specified
    keyword and model_id (optionally select non-default
    diagnostic via alias_select).
    """
    return ddef.alias_vars[
        diagnostic_per_model[keyword][model_id][
            alias_select % len(
                diagnostic_per_model[keyword][model_id])
        ]
    ]


# Also defined for getting from aliases:
def get_diagnostic_from_alias(alias):
    """Full diagnostic name from input alias."""
    return ddef.alias_diags[alias]


def get_coords_from_alias(alias):
    """Diagnostic coordinates (north, south) from a specified
    alias."""
    return ddef.alias_coords[alias]


def get_nc_vars_from_alias(alias):
    """Diagnostic variable names (north, south) from a specified
    alias.
    """
    return ddef.alias_vars[alias]



# All three at once (as is usually needed, thus this is the main
# API function). Here diag_list can contain keywords, aliases,
# or full diagnostic names, and is converted to the latter:
def get_load_info(diag_list, model_id,
        experiment_id="piControl", select=[0]):
    """Load one or more diagnostics for a specified model_id and
    experiment_id.
    
    
    Parameters
    ----------
    diag_list : list (length n_diag) of str
        Diagnostics to load. This can be a mixture of keywords,
        aliases, and/or full diagnostic names.
        
        For each entry, if a keyword is detected, the model_id
        is used to determine the full diagnostic name (see
        select parameter). If instead an alias is detected,
        this directly corresponds to a full diagnostic name
        (independent of model_id). Otherwise, it is assumed that
        a full diagnostic name is provided and nothing is done
        to it.
    
    model_id: str, model_id
    
    
    Optional parameters
    -------------------
    experiment_id: str, experiment_id, default = "piControl"
    
    select: int or list of int of length n_diag, default = [0]
        In the case of a keyword, some models may have multiple
        corresponding diagnostics available. The default is to
        select the first (highest priority) one. This parameter
        can be used to select a different index.
        
        Can be an individual choice for multiple diagnostics. If
        select values run out (length of this list less than
        n_diag) then the first value is used for remaining
        diagnostics.
    
    
    Returns
    -------
    diag_list_ret : list of length n_diag of str
        Full diagnostic names corresponding to input diag_list.
    
    coord_names_list : list of length n_diag of lists [str, ...]
        Spatial coordinate variable names for the corresponding
        diagnostics. Always either an empty list (if there are
        no applicable coordinates) or length-2 list.
    
    var_names_list : list of length n_diag of lists [str, ...]
        NetCDF variable names for each corresponding diagnostic.
        Always length-2 lists corresponding to [north, south]
        versions of the diagnostic.
    
    """
    
    diag_list_ret = diag_list.copy()
    
    if type(select) not in [list, tuple]:
        select = [select]
    
    for d in range(len(diag_list)):
        
        if diag_list_ret[d] in ddef.defined_keywords:
            diag_list_ret[d] = get_diagnostic_from_keyword(
                diag_list_ret[d], model_id,
                select[d%len(select)])
        
        elif diag_list_ret[d] in ddef.defined_aliases:
            diag_list_ret[d] = get_diagnostic_from_alias(
                diag_list_ret[d])
        
        # else assumed to be a full diagnostic name
    
    coord_names_list = [ddef.diag_coords[x]
                        for x in diag_list_ret]
    var_names_list = [ddef.diag_vars[x] for x in diag_list_ret]
    
    return diag_list_ret, coord_names_list, var_names_list
