"""Settings for the appearance of plots and customising
matplotlib rcParams."""

from pathlib import Path
import matplotlib as mpl

# Figure dimensions in millimetres. In plotting scripts,
# generally change the height as required but not the width (but
# a default height is needed here). This is so that sizes of
# fonts and other plot elements are always consistent when
# scaled to the same display (or print) width:
fig_width_single_mm = 90.0
fig_width_double_mm = 180.0
fig_height_mm = 101.25
# Here, 180.0 mm is the width required by the journal for our 
# manuscript and 101.25 mm sets the aspect ratio to be 16:9
# (for double-width figures)

# Matplotlib requires figure dimensions in inches:
mm_per_in = 25.40
fig_width_double_in = fig_width_double_mm / mm_per_in
fig_width_single_in = fig_width_single_mm / mm_per_in

# Default: double width:
fig_width_in = fig_width_double_in
fig_height_in = fig_height_mm / mm_per_in

# Resolution (dots per inch) for (embedded) raster graphics:
fig_raster_dpi = 300

# Set default sizes:
lw_0 = 0.75  # default line width
fs_0 = 12    # default font size
dfs  = 1     # default font-size increment
ms_0 = 3.5   # default marker size

panel_label_fmt = "{}"
panel_labels = "abcdefghij"

# ============================================================ #
# Create system for setting zorder of plot elements
# consistently
# ------------------------------------------------------------ #

_zorder = {
    "errorbar"       : -1,
    "regression_line": -2,
    "scatter_exclude": -3,
    "scatter"        : -22,
    "legend"         : -23,
    "ebm_line"       : -24,
    "gridlines"      : -25,
    "regression_fill": -26,
    "ebm_fill"       : -27,
    "obs_rectangle"  : -28,
    "rasterization"  : -29
}

def get_zorder(plot_element, zorder_dict=_zorder,
               default_zorder=0):
    """Get the z-order of a plot element from zorder_dict (or
    default value if plot_element is not in that dict).
    """
    if plot_element in zorder_dict.keys():
        return zorder_dict[plot_element]
    else:
        return default_zorder
# ------------------------------------------------------------ #


# ============================================================ #
# Default appearance/keyword arguments for specific plot
# elements
# ------------------------------------------------------------ #
default_annotate_kw = {}

default_cmap_delta_oht = "RdBu_r"
default_cmap_delta_tas = "inferno"

default_ebm_annotate_kw={"color": "tab:green", "ha": "center",
    "va": "center"}

default_cmip6_annotate_kw = {"color": "tab:gray"}

default_ebm_fill_kw = {"facecolor": [0.831, 0.925, 0.831],
    "edgecolor": "none", "zorder": get_zorder("ebm_fill")}

default_ebm_line_kw = {"color": "tab:green",
    "zorder": get_zorder("ebm_line")}

default_errorbar_kw = {"ecolor": "k", "capsize": 3,
    "zorder": get_zorder("errorbar"), "elinewidth": lw_0}

default_legend_kw = {"labelspacing": 0.275,
    "handletextpad": 0.035, "markerscale": 1.25,
    "fontsize": fs_0-3*dfs}

default_obs_annotate_kw = {"fontsize": fs_0-4*dfs,
    "color": "k", "fontstyle": "italic",
    "bbox": {"facecolor": "white", "edgecolor": "none",
             "pad": 0}}

default_obs_rectangle_kw = {"fill": True, "facecolor": [0.9]*3,
    "edgecolor": "none", "zorder": get_zorder("obs_rectangle")}

default_panel_label_title_kw = {"fontsize": fs_0,
    "fontweight": "bold", "loc": "left"}

default_panel_title_kw = {"fontsize": fs_0,
    "fontweight": "normal", "loc": "center"}

default_regression_line_kw = {"color": "tab:red",
    "zorder": get_zorder("regression_line")}

default_regression_fill_kw = {"facecolor": "tab:red",
    "edgecolor": "none", "alpha": 0.25,
    "zorder": get_zorder("regression_fill")}

default_regression_legend_kw = {"handlelength": 0,
    "handletextpad": 0, "facecolor": "white", "framealpha": 1,
    "loc": "upper left", "fontsize": fs_0 - 3*dfs}

default_scatter_kw = {"joinstyle": "miter",
    "zorder": get_zorder("scatter")}

default_scatter_exclude_kw = {"facecolors": "none",
    "edgecolors": "k", "joinstyle": "miter",
    "zorder": get_zorder("scatter_exclude")}

default_zero_gridline_kw = {
    "linewidth": lw_0, "color": [0.85]*3,
    "linestyle": "-", "zorder": get_zorder("gridlines")}
# ------------------------------------------------------------ #


# ============================================================ #
# Customize model markers/colours
# ------------------------------------------------------------ #

# Any "excluded" models are excluded from EBM fitting but
# still plotted -- this is a list (or iterable) of markers used
# to identify such models. This is used in the paper
# supplementary figure 6:
marker_exclude = "so^dX"

model_markers = {
    "AWI-CM-1-1-MR": "o",
    "CESM2"        : "o",
    "CESM2-FV2"    : "s",
    "CESM2-WACCM"  : "x",
    "CNRM-CM6-1"   : "o",
    "CNRM-CM6-1-HR": "s",
    "CNRM-ESM2-1"  : "x",
    "CanESM5"      : "o",
    "CanESM5-CanOE": "s",
    "GFDL-ESM4"    : "o",
    "GISS-E2-2-G"  : "x",
    "IPSL-CM6A-LR" : "o",
    "MIROC6"       : "x",
    "MPI-ESM1-2-HR": "o",
    "MPI-ESM1-2-LR": "s",
    "MRI-ESM2-0"   : "^",
    "NorESM2-LM"   : "o",
    "NorESM2-MM"   : "s",
    "UKESM1-0-LL"  : "o",
    "UKESM1-1-LL"  : "s"}


model_colors = {
    "AWI-CM-1-1-MR": "tab:orange",
    "CESM2"        : "tab:blue",
    "CESM2-FV2"    : "tab:blue",
    "CESM2-WACCM"  : "tab:blue",
    "CNRM-CM6-1"   : "tab:pink",
    "CNRM-CM6-1-HR": "tab:pink",
    "CNRM-ESM2-1"  : "tab:pink",
    "CanESM5"      : "tab:olive",
    "CanESM5-CanOE": "tab:olive",
    "GFDL-ESM4"    : "tab:brown",
    "GISS-E2-2-G"  : "tab:orange",
    "IPSL-CM6A-LR" : "tab:cyan",
    "MIROC6"       : "tab:brown",
    "MPI-ESM1-2-HR": "tab:purple",
    "MPI-ESM1-2-LR": "tab:purple",
    "MRI-ESM2-0"   : "tab:orange",
    "NorESM2-LM"   : "tab:red",
    "NorESM2-MM"   : "tab:red",
    "UKESM1-0-LL"  : "tab:green",
    "UKESM1-1-LL"  : "tab:green"}


def set_mpl_rcParams(font="Nimbus Sans"):
    """Customise matplotlib rcParams. This function must be
    called before any plotting calls.
    
    See matplotlib.org/stable/users/explain/customizing.html
    for explanation of these rcParams.
    """
    
    # -------------------------------------------------------- #
    # Set fonts
    # ======================================================== #
    # Sans-serif fonts is the default family for text and math
    # text already (which is what we want), so just need to set
    # the font for the sans-serif family.
    # 
    # rcParam "font.sans-serif" is a list of fonts of which the
    # first available is selected and used. So prepend the
    # desired font to that list:
    mpl.rcParams["font.sans-serif"].insert(0, font)
    
    # Update default font size:
    mpl.rcParams["font.size"] = fs_0
    
    # Set math text to use custom settings, which means that
    # each style (bold, italic, etc.) is specified individually.
    # The default values for those are fine (point to the
    # appropriate sans-serif fonts variants) so just need to
    # set that mathtext fontset to custom:
    mpl.rcParams["mathtext.fontset"] = "custom"
    
    # This is just to supress warnings if no cursive fonts are
    # available (they are not used in any plots anyway):
    mpl.rcParams["mathtext.cal"] = "sans:italic"
    # -------------------------------------------------------- #
    
    # Set various plot defaults. Line plots [ax.plot()]:
    mpl.rcParams["lines.linewidth"] = lw_0
    mpl.rcParams["lines.markeredgewidth"] = lw_0
    mpl.rcParams["lines.markersize"] = ms_0
    mpl.rcParams["lines.dash_joinstyle"] = "miter"
    mpl.rcParams["lines.dash_capstyle"] = "projecting"
    mpl.rcParams["lines.solid_joinstyle"] = "miter"
    mpl.rcParams["lines.solid_capstyle"] = "projecting"
    
    # Other types of plotting:
    mpl.rcParams["errorbar.capsize"] = \
        default_errorbar_kw["capsize"]
    mpl.rcParams["scatter.edgecolors"] = "none"
    
    mpl.rcParams["axes.titlesize"] = fs_0
    mpl.rcParams["axes.titlepad"] = 8
    mpl.rcParams["axes.titleweight"] = "bold"
    mpl.rcParams["axes.titlelocation"] = "left"
    mpl.rcParams["axes.labelsize"] = fs_0
    mpl.rcParams["axes.labelpad"] = 8
    mpl.rcParams["axes.spines.right"] = False
    mpl.rcParams["axes.spines.top"] = False
    mpl.rcParams["axes.axisbelow"] = True
    mpl.rcParams["axes.linewidth"] = lw_0
    
    # Grid off by default, but set anyway, e.g., for zero lines:
    mpl.rcParams["axes.grid"] = False
    mpl.rcParams["grid.color"] = [0.95]*3
    mpl.rcParams["grid.linewidth"] = lw_0
    mpl.rcParams["grid.linestyle"] = "-"
    
    for x in "xy":
        mpl.rcParams[f"{x}tick.direction"] = "in"
        mpl.rcParams[f"{x}tick.labelsize"] = fs_0
        
        mpl.rcParams[f"{x}tick.major.width"] = lw_0
        mpl.rcParams[f"{x}tick.major.size"] = 3
        mpl.rcParams[f"{x}tick.major.pad"] = 4
        
        mpl.rcParams[f"{x}tick.minor.visible"] = False
        mpl.rcParams[f"{x}tick.minor.width"] = lw_0
        mpl.rcParams[f"{x}tick.minor.size"] = 2.5
    
    mpl.rcParams["xtick.top"] = False
    mpl.rcParams["xtick.bottom"] = True
    mpl.rcParams["xtick.minor.top"] = False
    mpl.rcParams["xtick.minor.bottom"] = False
    mpl.rcParams["xtick.major.top"] = False
    mpl.rcParams["xtick.major.bottom"] = True
    
    mpl.rcParams["ytick.left"] = True
    mpl.rcParams["ytick.right"] = False
    mpl.rcParams["ytick.minor.left"] = False
    mpl.rcParams["ytick.minor.right"] = False
    mpl.rcParams["ytick.major.left"] = True
    mpl.rcParams["ytick.major.right"] = False
    
    mpl.rcParams["xaxis.labellocation"] = "right"
    mpl.rcParams["yaxis.labellocation"] = "top"
    
    mpl.rcParams["legend.fancybox"]  = False
    mpl.rcParams["legend.fontsize"]  = \
        default_legend_kw["fontsize"]
    mpl.rcParams["legend.facecolor"] = "none"
    mpl.rcParams["legend.edgecolor"] = "none"
    
    # Figure setup (figsize must be given in inches):
    mpl.rcParams["figure.figsize"] = (fig_width_in,
                                      fig_height_in)
    mpl.rcParams["savefig.dpi"]    = fig_raster_dpi
    
    # Miscelleneous:
    mpl.rcParams["axes.autolimit_mode"] = "round_numbers"
    mpl.rcParams["patch.linewidth"] = lw_0
    mpl.rcParams["svg.fonttype"] = "none"  # no text-to-paths
