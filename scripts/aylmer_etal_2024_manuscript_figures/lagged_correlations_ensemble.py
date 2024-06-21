"""Script to plot lagged correlations between various diagnostic
combinations, for CMIP6 historical simulations, all models
combined (supplementary information Fig. S2).
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from utils import maths, plot_style, plot_tools, script_tools



def main():
    
    # ======================================================== #
    # Parse command-line arguments and print script parameters
    # -------------------------------------------------------- #
    prsr = script_tools.argument_parser(
        "Lagged correlations (multiple models combined)")
    script_tools.add_cmip_selection_cmd_args(prsr,
        both_hemispheres=True)
    
    # Even though there are two panels, set n_panels=1 because
    # the axes limits are always the same in both cases (also
    # the y-axis ticks/limits are fixed and ignored from cmd):
    script_tools.add_plotting_cmd_args(prsr, n_panels=1, 
        n_cbars=0, text_labels=[])
    prsr.add_argument("--yearspan", type=int, nargs=2,
                      default=(1956, 2045))
    
    cmd, models = script_tools.get_args(prsr,
        ad_hoc_options={"yearspan": ("Time span", "{}")})
    
    n_models = len(models)
    n_yr_avg_1 = cmd.yravg1[1] - cmd.yravg1[0] + 1
    n_yr_avg_2 = cmd.yravg2[1] - cmd.yravg2[0] + 1
    
    
    # ======================================================== #
    # Some hard-coded settings:
    # -------------------------------------------------------- #
    # Combinations (X(t), Y), colors, whether to plot in
    # (n, s) hemi, label y offset for each hemisphere,
    # respectively; for computing lagged correlations:
    combo_settings = [
        [("oht", "iel"), "tab:blue", (True, True), (0, 0.025)],
        [("oht", "tas"), "tab:green", (False, True), (0, -0.15)],
        [("tas", "iel"), "tab:cyan", (True, True), (0, 0.15)],
        [("aht", "iel"), "tab:red", (True, True), (0, 0)],
        [("oht", "aht"), "tab:purple", (True, True), (0.05, 0.2)]
    ]
    
    # Text labels for legend entries:
    diag_labels = {
        "iel" : r"$\Delta\phi_\mathrm{i}$",
        "tas" : r"$\Delta{}T$",
        "oht" : r"$\Delta\mathrm{OHT}$",
        "aht" : r"$\Delta\mathrm{AHT}$"}
    
    
    # ======================================================== #
    # Load data
    # -------------------------------------------------------- #
    combos = [x[0] for x in combo_settings]
    colors = [x[1] for x in combo_settings]
    combos_hemi = [x[2] for x in combo_settings]
    labels_dy = [x[3] for x in combo_settings]
    
    # Load data:
    t_lag, r_n, r_s, n_ens = \
        script_tools.get_data_for_lagged_correlation_plots(
            models, combos, year_span=cmd.yearspan,
            fixed_yr_avg_1=cmd.yravg1,fixed_yr_avg_2=cmd.yravg2,
            ref_lats_n=cmd.reflatsn, ref_lats_s=cmd.reflatss,
            future_experiment=cmd.exp, as_ensemble=True)
    
    # For the correlation "insignificance" region (not exactly
    # as it should do some sort of weighting by number of
    # ensemble members; this is just for a rough indication):
    r_crit = maths.correlation_critical_value(n_models, 0.95)
    
    
    # ======================================================== #
    # Create figure
    # -------------------------------------------------------- #
    # Create a vertically-stacked 2-panel figure with the same
    # width as other two-column figures (note: this figure is
    # for supplementary information):
    fig, axs = plt.subplots(nrows=2,
        figsize=(plot_style.fig_width_in,
                 plot_style.fig_height_in))
    
    plot_tools.distribute_subplots(axs, vertically=True,
        l=0.25, r=0.30, b=0.025, t=0.065, s_ver=0.085)
    
    for ax in axs:
        # The horizontal axes are shared; set the horizontal
        # axis spine to be centered on y = 0:
        ax.spines["bottom"].set_position(("data", 0.0))
        
        # Set horizontal ticks on both sides (above and below):
        ax.tick_params(axis="x", which="major",
                       direction="inout", pad=2,
                       size=mpl.rcParams["xtick.major.size"]*2)
        
        # Set the axes spines to be above the data:
        ax.set_axisbelow(False)
        ax.spines["bottom"].set(zorder=1E10)
        ax.spines["left"].set(zorder=1E10)
    
        # Plot a rightward filled triangle to make the
        # horizontal axes spines have arrows:
        ax.plot(1, 0, ">k", transform=ax.get_yaxis_transform(),
                clip_on=False)
        
        # Set the axes limits/ticks which are always the same
        # for both panels (use cmd.a_xlim and cmd.a_xticks for
        # both and ignore cmd.a_ylim/cmd.a_yticks):
        plot_tools.set_axis_ticks_from_cmd(ax, cmd)
        plot_tools.set_axis_ticks(ax, -1, 1, 0.5,
            lims_actual=[-1, 1], which="y")
        
        # Correlation "insignificant" region:
        ax.fill_between(ax.get_xlim(), -r_crit, r_crit,
            facecolor=[0.9]*3, edgecolor="none", zorder=-100)
    
    # Create/setup the legends on each axes. Here everything is
    # drawn in data units (i.e., x is in units of years and y is
    # in units of correlation) and drawn as if on the axes but
    # with clipping set off so the legend ultimately appears
    # outside of the plot area. This means that changing time
    # spans etc. will probably mess this up (then again, in that
    # case the y-offsets specified in combo_settings which are 
    # hardcoded will also be messed up anyway).
    # 
    # Set some metrics:
    legend_left_width = 7.0  # years
    legend_right_width = 5.5  # years
    legend_ax_margin = 2.0  # years; gap between axes and legend
    legend_title_y0 = 1.0  # texts "X" and "Y" (on y-axis scale)
    legend_rule_y0 = 0.82  # rule height (on y-axis scale)
    
    # Need to use "top" vertical alignment for texts because mpl
    # uses the bounding box of the text to determine alignment,
    # and so using "center" would cause X/Y labels to be mis-
    # aligned depending on the character heights in each label
    # (e.g., \delta{}OHT would be higher than \delta\phi_{i}
    # with va="center"). Therefore, need to include a constant
    # offset to account for the use of "top" alignment, so that
    # y=0 in the ad-hoc combo_settings[-1] y-label positioning
    # offsets still sits roughly centred on the end point of
    # the corresponding curve. This works as long as none of the
    # text labels have superscripts:
    legend_dy_font = 0.06
    
    legend_left_col_x0 = axs[0].get_xlim()[1] + legend_ax_margin
    legend_right_col_x0 = legend_left_col_x0 + legend_left_width
    
    legend_text_kw = {"annotation_clip": False,
        "fontsize": mpl.rcParams["legend.fontsize"],
        "ha": "left", "va": "top"}
    
    for j in range(2):
        axs[j].annotate("$X$", (legend_left_col_x0,
                                legend_title_y0),
                        **legend_text_kw)
        axs[j].annotate("$Y$", (legend_right_col_x0,
                                legend_title_y0),
                        **legend_text_kw)
    
    rule_line_kw = {"color": "k", "clip_on": False,
        "linewidth": mpl.rcParams["axes.linewidth"]}
    
    for ax in axs:
        ax.plot([legend_left_col_x0,
            legend_left_col_x0 + legend_left_width +
            legend_right_width],
            [legend_rule_y0]*2, **rule_line_kw)
    
    for j in range(len(combos)):
        
        for h, data in zip([0, 1], [r_n[j], r_s[j]]):
            
            # Only plot if specified for this hemisphere:
            if combos_hemi[j][h]:
                
                # Lagged-correlation lines:
                axs[h].plot(t_lag, data, color=colors[j],
                            zorder=100)
                
                # Scatter points and lines at lags of maximum
                # (or minimum if negative) correlation:
                max_corr = np.max(data)
                min_corr = np.min(data)
                
                if min_corr < 0 and abs(min_corr) > max_corr:
                    xplot = t_lag[np.argmin(data)]
                    yplot = np.min(data)
                else:
                    xplot = t_lag[np.argmax(data)]
                    yplot = np.max(data)
                
                # Add drop line for point of maximum/minimum
                # correlation:
                axs[h].plot([xplot]*2, [0.0, yplot],
                            color=colors[j], zorder=90,
                            linestyle=(0, (5,5)))
                
                # Add marker for point of maximum/minimum
                # correlation:
                axs[h].scatter(xplot, yplot,
                    facecolor=colors[j], zorder=125)
                
                # Add legend X entry for this combination:
                axs[h].annotate(diag_labels[combos[j][0]],
                    (legend_left_col_x0,
                     data[-1] + labels_dy[j][h]+legend_dy_font),
                    color=colors[j], **legend_text_kw)
                
                # Add legend Y entry for this combination:
                axs[h].annotate(diag_labels[combos[j][1]],
                    (legend_right_col_x0,
                     data[-1] + labels_dy[j][h]+legend_dy_font),
                    color=colors[j], **legend_text_kw)
    
    # Add t_lag = 0 grid lines (must do this after plotting
    # data/setting any manual axes limits):
    plot_tools.add_zero_gridlines(axs, which="x")
    
    # x-axis label on top panel (a), above axis:
    axs[0].annotate(r"$t$ (years)",
                    (axs[0].get_xlim()[1], 0.1),
                    ha="right", va="bottom",
                    fontsize=mpl.rcParams["axes.labelsize"]-2)
    
    # Shared y-axis label:
    plot_tools.shared_yaxis_label(axs,
        r"Correlation between $X(t)$ and $Y$",
        left_frac=0.77)
    
    # Lead/lag labels (at hard-coded locations):
    kw = {"va": "center", "clip_on": False,
          "fontsize": mpl.rcParams["axes.labelsize"],
          "arrowprops": {"arrowstyle": "->",
                         "color": "tab:gray"},
          "color": "tab:gray"}
    
    axs[1].annotate("$X$ leads $Y$",
        xy=(-20, axs[1].get_ylim()[0]),
        xytext=(-1, axs[1].get_ylim()[0]), ha="right", **kw)
    axs[1].annotate("$Y$ leads $X$",
        xy=(20, axs[1].get_ylim()[0]),
        xytext=(1, axs[1].get_ylim()[0]), ha="left", **kw)
    
    # Set subplot panel labels (a/b) and titles:
    plot_tools.add_subplot_panel_labels(axs)
    plot_tools.add_subplot_panel_titles(axs,
        ["Arctic Ocean", "Southern Ocean"])
    
    plot_tools.finish_fig(fig, savefig=cmd.savefig,
        file_name=cmd.savefigname,
        fig_metadata={"Title": cmd.savefigtitle})



if __name__ == "__main__":
    main()
