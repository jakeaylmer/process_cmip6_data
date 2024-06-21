"""Common plotting routines."""

from datetime import datetime as dt, timezone as tz
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from . import maths, plot_style, script_tools


# ============================================================ #
# General
# ------------------------------------------------------------ #

def _add_default_kw(kw_in, kw_default):
    """Add default keyword arguments (kw_default, dict) to input
    keyword arguments (kw_in, dict) if not present. This is to
    allow defaults to be overridden in various routines below.
    Returns a new (copy, possibly modified) version of kw_in.
    """
    
    kw_out = kw_in.copy()
    
    for k in kw_default.keys():
        if k not in kw_in.keys():
            kw_out[k] = kw_default[k]
    
    return kw_out



def distribute_subplots(axs, l=0.15, r=0.05, b=0.15, t=0.10,
        s_hor=0.05, s_ver=0.05, vertically=False):
    """Uniformly distribute the subplots (array of matplotlib
    Axes instances) on a figure within specified margins and
    subplot spacing.
    
    
    Parameters
    ----------
    axs : array of matplotlib.axes.Axes instances [e.g., as
          returned by plt.subplots()]
        The subplots/axes to position, <= 2 dimensions.
    
    
    Optional parameters
    -------------------
    l, r, b, t: float, between 0 and 1
        Left, right, bottom, and top margins respectively, in
        figure coordinates.
    
    s_hor : float, default = 0.05
        Horizontal spacing between columns of subplots, in
        figure coordinates.
    
    s_ver : float, default = 0.05
        Vertical spacing between rows of subplots, in figure
        coordinates.
    
    vertically : bool, default = False
        This parameter is only used in cases where a 1D array
        of Axes is passed as axs. In this case the distribution
        direction is ambiguous (it could be a row or a column of
        subplots) so must be specified here: by default, the
        subplots are distributed horizontally.
    """
    
    # Input needs to be a 2D array. If a single Axes instance is
    # input, or a 1D array of Axes, expand the dimensions to
    # make it a 2D array, and set a flag to remove the expanded
    # dimensions at the end:
    if np.ndim(axs) == 0:
        axs = np.array([[axs]])
        squeeze = True
    elif np.ndim(axs) == 1:
        if vertically:
            axs = axs[:, np.newaxis]
        else:
            axs = axs[np.newaxis, :]
        squeeze = True
    else:
        squeeze = False
    
    ny, nx = np.shape(axs)
    
    # Compute the equal widths and heights of each Axes:
    axw = (1 - l - r - s_hor*(nx-1)) / nx
    axh = (1 - t - b - s_ver*(ny-1)) / ny
    
    for j in range(ny):
        for i in range(nx):
            axs[j,i].set_position([l + i*(s_hor + axw),
                                   b + (ny-j-1)*(s_ver + axh),
                                   axw, axh], which="both")
    
    if squeeze:
        axs = axs.squeeze()



def set_axis_ticks(ax, tick_start, tick_end, tick_step,
        lims_actual=None, which="x"):
    """Set an Axes (ax) major ticks manually on a specified axis
    (which = "x" or "y") from tick_start to tick_end (inclusive)
    in steps of tick_step.
    
    
    Optional parameters
    -------------------
    lims_actual: tuple of float (xmin, xmax) or None
        The actual axes limits required, if different from
        (tick_start, tick_end).
    
    which : str "x" or "y", default = "x"
        Which axis to apply ticks to.
    """
    
    ticks = np.arange(tick_start, tick_end+tick_step/2,
                      tick_step)
    
    if which == "x":
        ax.set_xticks(ticks)
        if lims_actual is None:
            ax.set_xlim((ticks[0], ticks[-1]))
        else:
            ax.set_xlim(lims_actual)
    else:
        ax.set_yticks(ticks)
        if lims_actual is None:
            ax.set_ylim((ticks[0], ticks[-1]))
        else:
            ax.set_ylim(lims_actual)



def set_axis_ticks_from_cmd(axs, cmd,
                            panel_labels=plot_style.panel_labels):
    """Set axes ticks and limits from command-line arguments.
    
    
    Parameters
    ----------
    axs : matplotlib.axes.Axes instance, or array of such
        The axes(s) to set ticks to.
    
    cmd : something with attributes:
              "p_xlim", "p_ylim" : 2-tuple of float
              "p_xticks", "p_yticks": 3-tuple of float
              and so on replacing "p" with panel labels for each
              subplot. Usually a return from
              argparse.ArgumentParser().parse_args().
        
        The x/y tick limits (min, max) and x/y tick
        (start, stop, step) values. Limits/ticks are not set if
        all values are 0.0 for a given panel.
        
        Subplots for each panel label are applied in order to
        axs.flatten().
    
    
    Optional parameters
    -------------------
    panel_labels : iterable of str
        Panel labels (p above).
    
    """
    
    if isinstance(axs, mpl.axes.Axes):
        axs_iter = np.array([axs])
    else:
        axs_iter = axs
    
    # Set main panel axes limits and ticks from the command-line
    # specified values for each panel j / axis x/y, if provided:
    for j in range(len(axs_iter.flatten())):
        for x in "xy":
            
            cmd_axlim = getattr(cmd,
                f"{panel_labels[j]}_{x}lim")
            cmd_axticks = getattr(cmd,
                f"{panel_labels[j]}_{x}ticks")
            
            if not (all([lim==0.0 for lim in cmd_axlim])
                    and all([t==0.0 for t in cmd_axticks])):
                
                set_axis_ticks(axs_iter.flatten()[j],
                    *cmd_axticks, lims_actual=cmd_axlim,which=x)



def shared_yaxis_label(axs, ylabel, left_frac=0.375,
                       ylabel_kw={}):
    """Add a common y-axis label to a vertical stack of plots,
    centred vertically. It is added as a figure (not axes)
    element.
    
    
    Parameters
    ----------
    axs : 1D array of matplotlib.axes.Axes instances
    ylabel : str
    
    
    Optional parameters
    -------------------
    left_frac : float, default = 0.375
        Where to position label horizontally, as a fraction of
        the space between the left-hand edge of the figure and
        the left-hand side of axs[0] position in figure
        coordinates.
    
    ylabel_kw : dict, default = {}
        Additional keyword arguments passed to fig.text().
    
    """
    
    y0 = axs[-1].get_position().y0
    y1 = axs[0].get_position().y1
    
    kw = _add_default_kw(ylabel_kw, {"ha": "right",
                                     "va": "center",
                                     "rotation": 90})
    
    axs[0].get_figure().text(left_frac*axs[0].get_position().x0,
                             y0 + 0.5*(y1 - y0), ylabel, **kw)



def add_zero_gridlines(axs, which="both", line_kw={}):
    """Draw gridlines on one or more subplots (axs, an array of
    Axes instances or a single instance) at x = 0 (which = "x"
    or "both") and/or y = 0 (which = "y" or "both").
    
    Note: lines are only drawn if visible on the plot. This also
    means that axis limit autoscaling is switched off by this
    routine, so that it should be called last or after data is
    plotted or fixed axis limits are set.
    
    
    Other parameters
    ----------------
    line_kw = dict, default = {}
        Keyword arguments passed to ax.plot() and specifying the
        gridline appearance, overriding defaults in
        plot_style.default_zero_gridline_kw where provided.
    """
    
    kw = _add_default_kw(line_kw,
                         plot_style.default_zero_gridline_kw)
    
    if type(axs) != np.ndarray:
        axs = np.array([axs])
        squeeze = True
    else:
        squeeze = False
    
    for ax in axs.flatten():
        
        # Get current axis limits:
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        
        if which in ["x", "both"]:
            ax.set_xlim(xlim)  # turns off auto-scaling
            if min(xlim) < 0 and max(xlim) > 0:
                ax.plot((0,0), ylim, **kw)
        
        # Repeat for y axis:
        if which in ["y", "both"]:
            ax.set_ylim(ylim)
            if min(ylim) < 0 and max(ylim) > 0:
                ax.plot(xlim, (0,0), **kw)
    
    if squeeze:
        axs = axs[0]



def add_subplot_panel_labels(axs,
        fmt=plot_style.panel_label_fmt,
        panel_labels=plot_style.panel_labels, title_kw={}):
    """Add subplot panel labels. Input axs must be an array of
    Axes instances (panel labels are not added to single panel
    figures).
    """
    for j in range(min(len(panel_labels), len(axs.flatten()))):
        axs.flatten()[j].set_title(fmt.format(panel_labels[j]),
            **_add_default_kw(title_kw,
                plot_style.default_panel_label_title_kw))



def add_subplot_panel_titles(axs, titles=[""], title_kw={}):
    """Add subplot panel titles (descriptive).
    """
    for j in range(min(len(titles), len(axs.flatten()))):
        axs.flatten()[j].set_title(titles[j],
            **_add_default_kw(title_kw,
                plot_style.default_panel_title_kw))



def add_text_label(ax, label, xy0, slope=0.0,
                   dr_label_text=(0.0, 0.0), annotate_kw={}):
    """Add a text label (label; str) to a plot (ax) at location
    xy0 (2-tuple of float, in data coordinates).
    
    
    Optional parameters
    -------------------
    slope : float, default = 0.0
        Rotate text so that it is parallel to a data line of
        this slope.
    
    dr_label_text : 2-tuple of float, default = (0.0, 0.0)
        Offset dr = (dx, dy), in data units, for the positioning
        of the annotation label, so that the actual position is
        (xy0[0] + dr_label_text[0], xy0[1] + dr_label-text[1]).
    
    annotate_kw : dict, default = {}
        Keyword arguments passed to ax.annotate(), overriding
        defaults in plot_style.default_annotate_kw where
        provided.
    """
    
    kw = _add_default_kw(annotate_kw,
                         plot_style.default_annotate_kw)
    
    if "rotation" not in kw.keys():
        # Calculate rotation from input slope. The rotation
        # argument of ax.annotate is in display coordinates, so
        # for example if slope is 1.0, this would be translated
        # to a rotation angle of +45 degrees just using
        # tan^-1(slope), but if the plot axes limits are not
        # equal aspect, then this angle is wrong. Therefore,
        # need to account for the aspect ratio of the plot axes:
        ax_ar = ax.get_xlim()[1] - ax.get_xlim()[0]
        ax_ar /= ax.get_ylim()[1] - ax.get_ylim()[0]
        kw["rotation"] = 180.0*np.arctan(slope*ax_ar)/np.pi
    
    ax.annotate(label, (xy0[0] + dr_label_text[0],
                        xy0[1] + dr_label_text[1]), **kw)


# ============================================================ #
# Specific
# ------------------------------------------------------------ #

def add_ebm_line(ax, ebm_slope, xlims, xy0,
                 ebm_label="EBM", dr_label_text=(0.0, 0.0),
                 line_kw={}, fill_kw={}, annotate_kw={}):
    """Add the EBM estimate (line with error shading) onto a
    plot (ax).
    
    
    Parameters
    ----------
    ebm_slope : tuple of float
        First element is the value of the EBM estimate of slope
        and the second element is the error (this is as output
        from the routines in the .ebm module).
    
    xlims : array of float
        The values [xmin, xmax] between which to draw the EBM
        line/shading (usually the min/max of the x-axis data).
    
    xy0 : 2-tuple of float
        The coordinates (in data units) to draw the EBM line/
        shading through (usually the ensemble mean of the data
        in the plot).
    
    
    Optional parameters
    -------------------
    ebm_label : str, default = "EBM"
        Text for the annotation label.
    
    dr_label_text : tuple of float, default = (0.0, 0.0)
        Offset dr = (dx, dy), in data units, for the positioning
        of the annotation label, which is set next to the
        maximum value on the EBM line at x = xlims[1].
    
    line_kw : dict, default = {}
        Keyword arguments passed to ax.plot() specifying the
        appearance of the line, overriding defaults in
        plot_style.default_ebm_line_kw where provided.
    
    fill_kw : dict, default = {}
        Keyword arguments passed to ax.fill_between() specifying
        the appearance of the shading, overriding defaults in
        plot_style.default_ebm_fill_kw where provided.
    
    annotate_kw : dict, default = {}
        Keyword arguments passed to ax.annotate() specifying the
        appearance of the annotation label, overriding defaults
        in plot_style.default_ebm_annotate_kw where provided.
    """
    
    ax.plot(xlims, xy0[1] + ebm_slope[0]*(xlims-xy0[0]),
            **_add_default_kw(line_kw,
                              plot_style.default_ebm_line_kw))
    
    ax.fill_between(xlims,
        xy0[1] + (ebm_slope[0] - ebm_slope[1])*(xlims - xy0[0]),
        xy0[1] + (ebm_slope[0] + ebm_slope[1])*(xlims - xy0[0]),
        **_add_default_kw(fill_kw,
                          plot_style.default_ebm_fill_kw))
    
    add_text_label(ax, ebm_label,
        (xlims[-1], xy0[1] + ebm_slope[0]*(xlims[-1] - xy0[0])),
        slope=0.0, dr_label_text=dr_label_text,
        annotate_kw=_add_default_kw(annotate_kw,
            plot_style.default_ebm_annotate_kw))



def add_regression_line(ax, xdata, ydata, add_legend=True,
                        slope_symbol="$m$", slope_units="",
                        decimal_places=(2, 2),
                        line_kw={}, fill_kw={}, legend_kw={}):
    """Add an orthogonal distance regression line to a subplot
    (ax) from input xdata and ydata (1D arrays). Labels with the
    value of the slope, the standard error, and correlation
    coefficient.
    
    
    Optional parameters
    -------------------
    add_legend : bool, default = True
        Whether to add a legend call to ax.
    
    slope_symbol : str, default = "$m$"
        For the text label/legend entry of the regression line,
        the symbol/label for the slope value.
    
    slope_units : str, default = ""
        For the text label/legend entry of the regression line,
        the units of the slope to be included after the value.
    
    decimal_places : 2-tuple of int
        Number of decimal places to report to the slope value
        and error, and correlation coefficient, respectively.
    
    line_kw, fill_kw, legend_kw : dict, default = {}
        Keyword arguments passed to ax.plot(),
        ax.fill_between(), and ax.legend(), respectively.
    
    """
    
    xplot = np.array([np.min(xdata), np.max(xdata)])
    
    x0 = np.nanmean(xdata)
    y0 = np.nanmean(ydata)
    
    m, m_err = maths.orthogonal_distance_regression(xdata,
                                                    ydata)[2:]
    
    if add_legend:
        
        # Construct label text for legend entry
        
        msign = u"\u2212" if m < 0 else ""
        mval = abs(m)
        
        r = maths.correlation_coefficient(xdata, ydata)
        
        rval = abs(r)
        rsign = u"\u2212" if r < 0 else "+"
        
        lbracket = "" if slope_units == "" else "("
        rbracket = "" if slope_units == "" else ")"
        
        label = slope_symbol
        
        if msign == "":
            label += "$=" + f"{lbracket}"
            label += f"{mval:.{decimal_places[0]}f}"
            label += r"\pm"
        else:
            label += f"$={lbracket}${msign}"
            label += f"{mval:.{decimal_places[0]}f}"
            label += r" $\pm\ "
        
        label += f"{m_err:.{decimal_places[0]}f}{rbracket}"
        label += "$ "
        label += slope_units
        label += "\n" + r"$r=$"
        label += f"{rsign}{rval:.{decimal_places[1]}f}"
        
        label = label.replace("$$", "")
        
    else:
        label = None
    
    l_kw = _add_default_kw(line_kw,
        plot_style.default_regression_line_kw)
    
    f_kw = _add_default_kw(fill_kw,
        plot_style.default_regression_fill_kw)
    
    ax.plot(xplot, y0 + m*(xplot - x0), label=label, **l_kw)
    ax.fill_between(xplot, y0 + (m - m_err)*(xplot - x0),
        y0 + (m + m_err)*(xplot - x0), **f_kw)
    
    if add_legend:
        leg = ax.legend(**_add_default_kw(legend_kw,
            plot_style.default_regression_legend_kw))
        leg.set_zorder(plot_style.get_zorder("legend"))



def add_range_as_error_bar(ax, error_range, offset, which="y",
                           errorbar_kw={}):
    """Add an x (or y) error bar to a plot (ax) from an error/
    uncertainty range [error_range; tuple of float (min, max) in
    data coordinates] at a specified y (or x) offset (float).
    
    
    Optional parameters
    -------------------
    which : str, "x" or (default) "y"
        Whether this is an x-data (i.e., horizontal) or y-data
        (i.e., vertical) error bar.
    
    errorbar_kw : dict, default = {}
        Keyword arguments passed to ax.errorbar(), overriding
        defaults in plot_style.default_errorbar_kw where
        provided. The argument "fmt" is also overridden to be
        "none" here.
    """
    
    eb_kw = _add_default_kw(errorbar_kw,
                            plot_style.default_errorbar_kw)
    
    if "zorder" in eb_kw:
        zorder = eb_kw["zorder"]
    else:
        zorder = plot_style.get_zorder("errorbar")
    
    eb_kw["fmt"] = "none"
    
    # ax.errorbar() requires a central value with symmetrical
    # +/- error as an array (2, N) of -/+ errors for each N
    # (here 1) data points, so determine these from input
    # uncertainty range:
    err = np.array([[0.5*abs(error_range[1] - error_range[0])]
                    for j in range(2)])
    
    avg = np.mean(error_range)
    
    # ax.errorbar() returns plotline, caplines, barlinecols. The
    # first is not needed since here we only draw the error bars
    # and nothing related to any underlying data (fmt="none").
    # The second two are needed (the error bar cap lines and the
    # barlines themselves), to set the zorders after creation;
    # it seems these cannot be set in the call to ax.errorbar():
    _, caplines, barlinecols = ax.errorbar(
        avg if which=="x" else offset,
        offset if which=="x" else avg,
        yerr=err if which=="y" else None,
        xerr=err if which=="x" else None,
        **eb_kw)
    
    for cap in caplines:
        cap.set_zorder(zorder)
    
    for col in barlinecols:
        col.set_zorder(zorder)



def add_points_along_error_bar(ax, values, offset, which="x",
        labels=[None], markers="os^d", facecolors=["k"],
        edgecolors=["none"], errorbar_kw={},
        scatter_kw={}, add_legend=True, legend_kw={}):
    """Add to a subplots (ax) a set of points (values) along a
    horizontal (which="x") or vertical (which="y") error-bar
    line.
    
    
    Parameters
    ----------
    values : iterable, length nval, of float
        Data points.
    
    offset : float
        Either x (if which = "y") or y (if which = "x") value to
        plot points along.
    
    
    Optional parameters
    -------------------
    which : str, "x" or "y"
        Whether to plot points along a horizontal ("x") or
        vertical ("y") bar.
    
    labels : list of str or None
        Labels for "label" argument of ax.scatter()
        corresponding to each value in values. If len(labels) <
        nval, indexing wraps back to the beginning of the list.
    
    markers : str or list of str
        Markers for ax.scatter(), selected in turn for each
        value in values. If len(markers) < nval, indexing wraps
        back to the beginning of the list.
    
    facecolors : list of matplotlib color identifiers
        Facecolors for ax.scatter(), selected in turn for each
        value in values. If len(facecolors) < nval, indexing
        wraps back to the beginning of the list.
    
    edgecolors : list of matplotlib color identifiers
        As in facecolors but for the marker edgecolors.
    
    errorbar_kw : dict, default = {}
        Keyword arguments passed to ax.errorbar(), overriding
        defaults in plot_style.default_errorbar_kw where
        provided. These arguments are passed via the
        plot_tools.add_range_as_error_bar() function, in which
        the argument "fmt" is also overridden to be "none".
    
    scattter_kw : dict, default = {}
        Keyword arguments passed to ax.scatter(), overriding
        defaults in plot_style.default_scatter_kw where
        provided. Here, the "zorder" is also overridden to be
        the zorder of the errorbar (in errorbar_kw) + 1 so that
        symbols appear on top.
    
    add_legend : bool, default = True
        Whether to add a legend (in which case labels should
        contain strings, not None).
    
    legend_kw : dict, default = {}
        Keyword arguments passed to ax.legend(), overriding
        defaults in plot_style.default_legend_kw provided.
    """
    
    eb_kw = _add_default_kw(errorbar_kw,
                            plot_style.default_errorbar_kw)
    
    sc_kw = scatter_kw.copy()
    
    # Override some properties:
    eb_kw["capsize"] = 0
    if "zorder" not in eb_kw:
        eb_kw["zorder"] = plot_style.get_zorder("errorbar")
    sc_kw["zorder"] = eb_kw["zorder"] + 1
    
    add_range_as_error_bar(ax, (min(values), max(values)),
                           offset, which=which,
                           errorbar_kw=eb_kw)
    
    nval = len(values)
    
    xdat = values if which == "x" else [offset]*nval
    ydat = [offset]*nval if which == "x" else values
        
    for j in range(nval):
        ax.scatter(xdat[j], ydat[j], marker=markers[j],
                   label=labels[j%len(labels)],
                   edgecolors=edgecolors[j%len(edgecolors)],
                   facecolors=facecolors[j%len(facecolors)],
                   **sc_kw)
    
    if add_legend:
        leg = ax.legend(**_add_default_kw(legend_kw,
                  plot_style.default_legend_kw))
        leg.set_zorder(plot_style.get_zorder("legend"))



# ============================================================ #
# Figure/suplot layouts
# ------------------------------------------------------------ #
# 
# The figure layouts here are straightforward grids of equal-
# sized subplots/panels.
# 
# Each fig_layout_*() function takes optional arguments "ax_l",
# "ax_r", "ax_t", "ax_b". These are 0 < floats < 1 specifying
# the left, right, top, and bottom margins on the figure canvas
# between the axes area. For example, ax_l=0.1 means 10% of the
# figure canvas width is the distance between the left-hand
# edge of the figure and the left-hand axis spine of the left-
# most column of axes.
# 
# Additionally, some take options "s_hor" and/or "s_ver", which
# indicate the distance between columns and rows of subplots
# respectively in the same fractional units.
# 
# Some shared defaults are here:
_ax_l  = 0.075
_ax_r  = 0.025
_s_hor = 0.090


def fig_layout_2_panel_2_cbar(ax_l=_ax_l, ax_r=_ax_r,
        ax_t=0.065, ax_b=0.30, s_hor=_s_hor):
    """Figure layout for two side-by-side panels each with a
    colour bar, with a space for a text legend on the left-hand
    colour bar. Used by Figs. 1, 3, S5, and S6. Returns fig,
    axs, cbar_axs.
    """
    
    fig, axs = plt.subplots(ncols=2,
        figsize=(plot_style.fig_width_in,
                 plot_style.fig_height_in))
    
    distribute_subplots(axs, l=ax_l, r=ax_r, t=ax_t, b=ax_b,
                        s_hor=s_hor)
    
    # Additional measures are required for the colorbar axes:
    # the location of the bottom edge and height, both in
    # figure units, here expressed as a fraction of the main
    # axes bottom margin:
    cbar_loc_b = 0.4*ax_b
    cbar_height = 0.14*ax_b
    
    # Also make space for the legend between the left colorbar.
    # First, the cbar widths *not* accounting for this (i.e.,
    # same as width of each axes above):
    cbar_width = 0.5*(1 - ax_l - ax_r - s_hor)
    
    # The right-hand colorbar has that width and sits directly
    # underneath the right-hand subplot. The left-hand colorbar
    # has a different left position (here set to the right-hand
    # axes margin, for symmetry) and a reduced overall width,
    # here expressed as a fraction of the unmodified width:
    cbar_0_left = ax_r
    cbar_0_width_frac = 0.875
    
    # Add the colorbar axes at the right locations on the figure
    # and collect into an array:
    cbar_axs = np.array([
        fig.add_axes([cbar_0_left, cbar_loc_b,
                      cbar_0_width_frac*cbar_width,
                      cbar_height]),
        fig.add_axes([ax_l + cbar_width + s_hor, cbar_loc_b,
                      cbar_width, cbar_height])])
    
    return fig, axs, cbar_axs



def fig_layout_2_panel(ax_l=_ax_l, ax_r=_ax_r, ax_t=0.095,
                       ax_b=0.155, s_hor=_s_hor,
                       fig_height_fraction_of_default=0.8):
    """Figure layout for two side-by-side panels. Used by Figs.
    4 and S1. The figure height is adjustable via the parameter
    fig_height_fraction_of_default.
    """
    
    fig, axs = plt.subplots(ncols=2,
        figsize=(plot_style.fig_width_in,
                 fig_height_fraction_of_default\
                 *plot_style.fig_height_in))
    
    distribute_subplots(axs, l=ax_l, r=ax_r, t=ax_t, b=ax_b,
                        s_hor=s_hor)
    
    return fig, axs



def fig_layout_2_panel_vertical(ax_l=0.17, ax_r=0.015,
                                ax_t=0.07, ax_b=0.065,
                                s_ver=0.13):
    """Figure layout for two stacked panels. Used by Fig. 8."""
    
    fig, axs = plt.subplots(nrows=2,
        figsize=(plot_style.fig_width_single_in,
                 1.125*plot_style.fig_width_single_in))
    
    distribute_subplots(axs, vertically=True,
        l=ax_l, r=ax_r, t=ax_t, b=ax_b, s_ver=s_ver)
    
    return fig, axs



def fig_layout_4_panel(ax_l=0.10, ax_r=_ax_r, ax_t=0.06,
                       ax_b=0.085, s_hor=0.10, s_ver=0.13):
    """Figure layout for a 2x2 grid of subplots. Used by Figs.
    6 and 7.
    """
    
    fig, axs = plt.subplots(ncols=2, nrows=2,
        figsize=(plot_style.fig_width_in,
                 0.75*plot_style.fig_width_in))
    
    distribute_subplots(axs, l=ax_l, r=ax_r, t=ax_t, b=ax_b,
                        s_hor=s_hor, s_ver=s_ver)
    
    return fig, axs



def fig_layout_6_panel(ax_l=0.09, ax_r=0.022, ax_t=0.04,
                       ax_b=0.06, s_hor=_s_hor, s_ver=0.08):
    """Figure layout for a 3 rows by 2 column grid of subplots.
    Used by Fig. 5.
    """
    
    fig, axs = plt.subplots(ncols=2, nrows=3,
        figsize=(plot_style.fig_width_in,
                 1.125*plot_style.fig_width_in))
    
    distribute_subplots(axs, l=ax_l, r=ax_r, t=ax_t, b=ax_b,
                        s_hor=s_hor, s_ver=s_ver)
    
    return fig, axs
    


# ============================================================ #
# Saving figures
# ------------------------------------------------------------ #

# Metadata is added to SVG output only. Only certain keys are
# allowed (by the backend SVG writer) and an error is raised if
# one is passed that is invalid. They are case sensitive. In
# the scripts, generally only the "Title" need be passed in to
# the save_figure function below. Other common values are
# defined in this dictionary, except for the "Date" which is
# added at the time of saving. I would add the coauthors, but 
# there seems to be a backend bug as only one "Contributor" is
# saved!
_default_svg_metadata = {
    "Contributor": ["Jake R. Aylmer"],
    "Language": "English"}


def save_figure(fig, file_name="", set_raster_level=False,
                file_fmts=[".svg", ".png"],
                fig_metadata={},
                save_dir=script_tools.get_path(
                    "path_manuscript_figures.txt")):
    """Save a figure (fig, matplotlib Figure instance) to a
    specified file(s) in one or more formats.
    
    
    Optional parameters
    -------------------
    file_name : str, default = ""
        File name excluding extension. If "" or None, the figure
        canvas window title is used.
    
    set_raster_level : bool, default = False
        Set to True if figure contains rasterized elements and
        saving to a vector format. In this case the savefig dpi
        is set so that those elements are rasterized.
    
    file_fmts : list of str, default = [".png", ".svg"]
        File extensions for image formats to save.
    
    fig_metadata : dict, default = {}
        Metadata added to SVG output only. Only specific keys
        are allowed; see
        
        https://matplotlib.org/stable/api/backend_svg_api.html
        #matplotlib.backends.backend_svg.FigureCanvasSVG.
        print_svg
        
        for details. Generally, only "Title" (case sensitive)
        should be provided here.
    
    save_dir : str or pathlib.Path
        Directory to save figures.
    """
    
    # Known raster and vector forms to check during loop below:
    raster_fmts = [".gif", ".jpg", ".png"]
    vector_fmts = [".eps", ".pdf", ".ps", ".svg"]
    
    # Make save directory if it doesn"t exist:
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    
    # Set fig file names to the window title if not provided:
    if file_name is None or len(file_name) == 0:
        file_name = fig.canvas.manager.get_window_title()
    
    file_name = file_name.replace(" ", "_")
    
    for fmt in file_fmts:
        
        kw = {}
        
        # Set DPI if saving as a raster format or with raster
        # elements embedded in a vector format:
        if (fmt in raster_fmts
                or (fmt in vector_fmts and set_raster_level)):
            kw["dpi"] = plot_style.fig_raster_dpi
        
        # Add metadata (SVG only):
        if "svg" in fmt:
            kw["metadata"] = _add_default_kw(fig_metadata,
                {"Date": dt.now(tz.utc).strftime(
                         "%H:%M UTC %d %b %Y"),
                 **_default_svg_metadata})
        
        fig.savefig(Path(save_dir, file_name + fmt), **kw)
        print(f"Saved: {str(Path(save_dir, file_name + fmt))}")



def finish_fig(fig, savefig=False, **kwargs):
    """Final steps for use in scripts that generate figures;
    either saves the figure or displays it interactively. Takes
    the Figure instance (fig), an optional parameter savefig
    (bool: whether to save), and additional keyword arguments
    are passed to plot_tools.save_figure().
    """
    if savefig:
        save_figure(fig, **kwargs)
    else:
        fig.show()
