"""Script to plot lagged correlations between various diagnostic
combinations, for CMIP6 historical simulations, for specified
models individually (supplementary information Fig. S3).
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
        "Lagged correlations (5 individual models)")
    script_tools.add_cmip_selection_cmd_args(prsr,
        both_hemispheres=True)
    
    # Even though there are ten panels, set n_panels=1 because
    # the axes limits are always the same in each case (also
    # the y-axis ticks/limits are fixed and ignored):
    script_tools.add_plotting_cmd_args(prsr, n_panels=1, 
        n_cbars=0, text_labels=[])
    prsr.add_argument("--yearspan", type=int, nargs=2,
                      default=(1956, 2045))
    prsr.add_argument("--models", type=str, nargs="*",
        default=["CanESM5", "IPSL-CM6A-LR", "MPI-ESM1-2-LR",
                 "MPI-ESM1-2-HR", "UKESM1-0-LL"])
    
    cmd, _ = script_tools.get_args(prsr, models=None,
        ad_hoc_options={"yearspan": ("Time span", "{}")})
    
    n_models = len(cmd.models)
    n_yr_avg_1 = cmd.yravg1[1] - cmd.yravg1[0] + 1
    n_yr_avg_2 = cmd.yravg2[1] - cmd.yravg2[0] + 1
    
    
    # ======================================================== #
    # Some hard-coded settings:
    # -------------------------------------------------------- #
    # Combinations (X(t), Y), colors, whether to plot in
    # (n, s) hemi; for computing lagged correlations:
    combo_settings = [
        [("oht", "iel"), "tab:blue", (True, True)],
        [("oht", "tas"), "tab:green", (False, True)],
        [("tas", "iel"), "tab:cyan", (True, True)],
        [("aht", "iel"), "tab:red", (True, True)],
        [("oht", "aht"), "tab:purple", (True, True)]
    ]
    
    # Text labels for legend entries:
    combos_label = [
        r"$r\left[\Delta\mathrm{OHT}(t),\Delta\phi_\mathrm{i}\right]$",
        r"$r\left[\Delta\mathrm{OHT}(t),\Delta{}T\right]$",
        r"$r\left[\Delta{}T(t),\Delta\phi_\mathrm{i}\right]$",
        r"$r\left[\Delta\mathrm{AHT}(t),\Delta\phi_\mathrm{i}\right]$",
        r"$r\left[\Delta\mathrm{OHT}(t),\Delta\mathrm{AHT}\right]$"
    ]
    
    
    # ======================================================== #
    # Load data
    # -------------------------------------------------------- #
    
    combos = [x[0] for x in combo_settings]
    colors = [x[1] for x in combo_settings]
    combos_hemi = [x[2] for x in combo_settings]
    
    # Load data:
    t_lag, r_n, r_s, n_ens = \
        script_tools.get_data_for_lagged_correlation_plots(
            cmd.models, combos, year_span=cmd.yearspan,
            fixed_yr_avg_1=cmd.yravg1,
            fixed_yr_avg_2=cmd.yravg2, ref_lats_n=cmd.reflatsn,
            ref_lats_s=cmd.reflatss, future_experiment=cmd.exp,
            as_ensemble=False)
    
    # For the correlation "insignificance" regions:
    r_crit = [maths.correlation_critical_value(n_ens[m], 0.95)
              for m in range(n_models)]
    
    
    # ======================================================== #
    # Create figure
    # -------------------------------------------------------- #
    
    fig_width = 207.965*(plot_style.fig_width_double_mm/147.084)
    fig_width /= plot_style.mm_per_in
    fig_height = (9.0/16.0)*fig_width
    
    fig, axs = plt.subplots(nrows=2, ncols=n_models,
        figsize=(fig_width, fig_height))
    
    plot_tools.distribute_subplots(axs, l=0.085, r=0.008,
        b=0.05, t=0.08, s_ver=0.075, s_hor=0.02)
    
    for ax in axs.flatten():
        # The horizontal axes are shared; set the horizontal
        # axis spine to be centered on y = 0:
        ax.spines["bottom"].set_position(("data", 0.0))
        
        # Set horizontal ticks on both sides (above and below):
        ax.tick_params(axis="x", which="major",
                       direction="inout", pad=2,
                       labelbottom=False,
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
    
    # Turn y-axis ticks off for all but the first column:
    for ax in axs[:,1:].flatten():
        ax.tick_params(labelleft=False)
    
    # Correlation "insignificant" regions:
    for col in range(n_models):
        for row in range(2):
            axs[row,col].fill_between(axs[row,col].get_xlim(),
                -r_crit[col], r_crit[col], facecolor=[0.9]*3,
                edgecolor="none", zorder=-100)
    
    for j in range(len(combos)):
        for m in range(n_models):
            for h, data in zip([0, 1], [r_n[m,j,:], r_s[m,j,:]]):
                
                # Only plot if specified for this hemisphere:
                if combos_hemi[j][h]:
                    
                    # Lagged-correlation lines:
                    axs[h,m].plot(t_lag, data, color=colors[j],
                                  label=combos_label[j],
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
                    axs[h,m].plot([xplot]*2, [0.0, yplot],
                                  color=colors[j], zorder=90,
                                  linestyle=(0, (5,5)))
                    
                    # Add marker for point of maximum/minimum
                    # correlation:
                    axs[h,m].scatter(xplot, yplot,
                                     facecolor=colors[j],
                                     zorder=125)
    
    # Add t_lag = 0 grid lines (must do this after plotting
    # data/setting any manual axes limits):
    plot_tools.add_zero_gridlines(axs, which="x")
    
    # Shared y-axis label:
    plot_tools.shared_yaxis_label(axs[:,0],
        r"Correlation between $X(t)$ and $Y$",
        left_frac=0.52)
    
    for j, txt in zip([0, 1], ["Arctic", "Southern Ocean"]):
        x0 = axs[j,0].get_position().x0
        y0 = axs[j,0].get_position().y0
        y1 = axs[j,0].get_position().y1
        fig.text(0.26*x0, y0 + 0.5*(y1-y0), txt,
            fontweight="bold", ha="right", va="center",
            rotation=90,
            fontsize=mpl.rcParams["axes.labelsize"])
    
    # Set the legend: x0 is as above from j = 1 in loop
    x1 = axs[-1,-1].get_position().x1
    y0 = axs[-1,-1].get_position().y0
    axs[1,0].legend(ncol=5, loc="center",
        bbox_to_anchor=(x0 + 0.5*(x1 - x0), 0.45*y0),
        bbox_transform=fig.transFigure,
        columnspacing=1.5,
        fontsize=mpl.rcParams["axes.labelsize"])
    
    # Set subplot panel labels (a, b, ...) and titles:
    plot_tools.add_subplot_panel_labels(axs)
    plot_tools.add_subplot_panel_titles(axs[0,:],
        [cmd.models[m] + "\n" + r"$n={}$".format(n_ens[m])
         for m in range(n_models)],
        title_kw={})
    
    plot_tools.finish_fig(fig, savefig=cmd.savefig,
        file_name=cmd.savefigname,
        fig_metadata={"Title": cmd.savefigtitle})



if __name__ == "__main__":
    main()
