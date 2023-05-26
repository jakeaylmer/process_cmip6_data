"""Define diagnostic keywords and aliases here. See
model_diagnostics module for explanation.
"""

from src import metadata as md
from src import netcdf as nf

from tabulate import tabulate


def _validate_input(x, y="diag"):
    """Validate input for diagnostic or netcdf variable name."""
    return (x if x == "" else getattr(nf, f"{y}_nq_{x}"))



def _gdn(name, tm="", sm="", om=""):
    """Wrapper function for netcdf.get_diag_name() -- method
    strings do not need to include 'diag_nq_'.
    """
    return nf.diag_name(name=name,
        time_methods=_validate_input(tm, "diag"),
        space_methods=_validate_input(sm, "diag"),
        other_methods=_validate_input(om, "diag"))



def _gvn(name, hemi="ns", tm="", sm="", om=""):
    """Wrapper function for netcdf.get_nc_var_name() -- method
    strings do not need to include 'nc_var_nq_'.
    """
    return [nf.nc_var_name(name=name, hemi=x,
        time_methods=_validate_input(tm, "nc_var"),
        space_methods=_validate_input(sm, "nc_var"),
        other_methods=_validate_input(om, "nc_var"))
        for x in hemi]


# With aliases as keys:
alias_diags    = {}  # Diagnostic (directory) names
alias_coords   = {}  # Corresponding nc coordinate vars.
alias_vars     = {}  # Corresponding nc variable names
alias_keywords = {}  # Corresponding keyword

# As above but with the diagnostic names as keys:
diag_coords   = {}
diag_vars     = {}
diag_keywords = {}


def _define_diagnostic(keyword, alias,
        diag_name="", nc_var_name="", hemi="ns",
        time_methods="", space_methods="", other_methods="",
        coord_vars=[nf.nc_ref_lat_n_name, nf.nc_ref_lat_s_name]
    ):
    """Append diagnostic alias and its actual (directory),
    name, coordinate and variable names (north and south), and
    corresponding keyword, to appropriate dictionaries.
    """
    
    if nc_var_name == "":
        nc_var_name = diag_name
    
    kw = {"tm": time_methods, "sm": space_methods,
          "om": other_methods}
    
    alias_diags[alias]    = _gdn(diag_name, **kw)
    alias_vars[alias]     = _gvn(nc_var_name, hemi=hemi, **kw)
    alias_coords[alias]   = coord_vars
    alias_keywords[alias] = keyword
    
    diag_coords[alias_diags[alias]]   = coord_vars
    diag_vars[alias_diags[alias]]     = alias_vars[alias]
    diag_keywords[alias_diags[alias]] = keyword



# Atmospheric diagnostics #
# ----------------------- #
for method, method_alias_part in zip(
        ["native", "native_interp", "cell_center_approx"],
        ["gn", "gn_interp", "cc"]):
    
    _define_diagnostic("aht",
        f"aht_from_net_flux_{method_alias_part}",
        diag_name="aht_from_net_flux", nc_var_name="aht",
        time_methods="yearly", other_methods=method)

    for x in ["down", "olr", "sw_surf", "sw_toa", "up"]:
        _define_diagnostic(f"f_{x}",
            f"f_{x}_{method_alias_part}",
            diag_name=f"f_{x}", time_methods="yearly",
            space_methods="area_mean",
            other_methods=method)

    _define_diagnostic("tas",
        f"tas_{method_alias_part}", diag_name="tas",
        time_methods="yearly", space_methods="area_mean",
        other_methods=method)

# Ocean diagnostics       #
# ----------------------- #
_define_diagnostic("hfds", "hfds",
    diag_name="hfds", time_methods="yearly",
    space_methods="area_integral",
    other_methods="cell_center_approx")

_define_diagnostic("hfds_wm-2", "hfds_wm-2",
    diag_name="hfds", time_methods="yearly",
    space_methods="area_mean",
    other_methods="cell_center_approx")

_define_diagnostic("oht", "oht_from_hfbasin",
    diag_name="oht_from_hfbasin", nc_var_name="oht",
    hemi=[""]*2, time_methods="yearly",
    other_methods="native_interp",
    coord_vars=[nf.nc_ref_lat_single_name]*2)

_define_diagnostic("oht", "oht_residual_pot",
    diag_name="oht_from_hfds_opottemptend" + nf.diag_nq_vertical_integral,
    nc_var_name="oht", time_methods="yearly",
    other_methods="cell_center_approx")

_define_diagnostic("oht", "oht_residual_con",
    diag_name="oht_from_hfds_ocontemptend" + nf.diag_nq_vertical_integral,
    nc_var_name="oht", time_methods="yearly",
    other_methods="cell_center_approx")

_define_diagnostic("oht", "oht_from_hfx_hfy",
    diag_name="oht_from_hfx_hfy", nc_var_name="oht",
    time_methods="yearly", other_methods="cell_center_approx")

for x in ["pot", "con"]:
    _define_diagnostic("dhdt", f"dhdt_{x}_direct",
        diag_name=f"o{x}temptend" + nf.diag_nq_vertical_integral,
        nc_var_name=f"o{x}temptend" + nf.nc_var_nq_vertical_integral,
        time_methods="yearly", space_methods="area_integral",
        other_methods="cell_center_approx")
    
    _define_diagnostic("dhdt", f"dhdt_{x}_residual_hfx",
        diag_name=f"o{x}temptend" + nf.diag_nq_vertical_integral + "_from_hfds_hfx_hfy",
        nc_var_name=f"o{x}temptend" + nf.nc_var_nq_vertical_integral,
        time_methods="yearly", space_methods="area_integral",
        other_methods="cell_center_approx")
    
    _define_diagnostic("dhdt", f"dhdt_{x}_residual_hfbasin",
        diag_name=f"o{x}temptend" + nf.diag_nq_vertical_integral + "_from_hfds_hfbasin",
        nc_var_name=f"o{x}temptend" + nf.nc_var_nq_vertical_integral,
        time_methods="yearly", space_methods="area_integral",
        other_methods="cell_center_approx")
    
    # Same as above 3, but area mean (wm-2) versions:
    _define_diagnostic("dhdt_wm-2", f"dhdt_{x}_direct_wm-2",
        diag_name=f"o{x}temptend" + nf.diag_nq_vertical_integral,
        nc_var_name=f"o{x}temptend" + nf.nc_var_nq_vertical_integral,
        time_methods="yearly", space_methods="area_mean",
        other_methods="cell_center_approx")
    
    _define_diagnostic("dhdt_wm-2", f"dhdt_{x}_residual_hfx_wm-2",
        diag_name=f"o{x}temptend" + nf.diag_nq_vertical_integral + "_from_hfds_hfx_hfy",
        nc_var_name=f"o{x}temptend" + nf.nc_var_nq_vertical_integral,
        time_methods="yearly", space_methods="area_mean",
        other_methods="cell_center_approx")
    
    _define_diagnostic("dhdt_wm-2", f"dhdt_{x}_residual_hfbasin_wm-2",
        diag_name=f"o{x}temptend" + nf.diag_nq_vertical_integral + "_from_hfds_hfbasin",
        nc_var_name=f"o{x}temptend" + nf.nc_var_nq_vertical_integral,
        time_methods="yearly", space_methods="area_mean",
        other_methods="cell_center_approx")
    


# Sea ice diagnostics     #
# ----------------------- #
for res in ["4deg", "2deg", "1deg", "05deg", "025deg"]:
    for method in ["bil", "dis"]:
        _define_diagnostic("iel", f"iel_{res}_{method}",
            diag_name="iel" + nf.diag_nq_zonal_mean,
            time_methods="yearly", coord_vars=[],
            other_methods=f"{res}_{method}")
        
        _define_diagnostic("iel_mon", f"iel_{res}_{method}_mon",
            diag_name="iel" + nf.diag_nq_zonal_mean,
            time_methods="monthly", coord_vars=[],
            other_methods=f"{res}_{method}")

for x in "ae":
    _define_diagnostic(f"si{x}", f"si{x}",
        diag_name=f"si{x}", time_methods="yearly",
        coord_vars=[])
    
    _define_diagnostic(f"si{x}_mon", f"si{x}_mon",
        diag_name=f"si{x}", time_methods="monthly",
        coord_vars=[])

# ----------------------------------------------------------- #


defined_keywords    = sorted(list(set(alias_keywords.values())))
defined_aliases     = sorted(list(alias_keywords.keys()))
defined_diagnostics = [alias_diags[k] for k in defined_aliases]


def print_defined_diagnostics():
    
    headers = ["Keyword", "Alias", "Diagnostic name",
               "Coordinates", "Variables"]
    
    rows = []
    
    # Group aliases per keyword (keep track of current keyword,
    # which is initially not known):
    current_keyword = "something_random_definitely_not_keyword"
    
    for keyword in defined_keywords:
        
        # Get ordered list of aliases for current keyword:
        aliases_k = sorted([a for a in alias_keywords.keys()
                   if alias_keywords[a] == keyword])
        
        kw_entry = keyword
        
        for a in aliases_k:
            rows.append([kw_entry, a, alias_diags[a],
                ", ".join(alias_coords[a]),
                ", ".join(alias_vars[a])
            ])
            kw_entry = "|"
    
    print("\n" + tabulate(rows, headers=headers) + "\n")
