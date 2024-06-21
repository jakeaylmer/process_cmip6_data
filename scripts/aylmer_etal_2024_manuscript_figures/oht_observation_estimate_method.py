"""Script to plot ocean heat transport into the Arctic and into
the Southern Ocean from ECCO data, demonstrating the method for
estimate the change in poleward OHT over the historical period
(manuscript Fig. 8).
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from utils import observations, plot_tools, script_tools



def main():
    
    # ======================================================== #
    # Parse command-line arguments and print script parameters
    # -------------------------------------------------------- #
    prsr = script_tools.argument_parser(
        "Plot ECCO data for demonstration of method to "
        "calculate observed delta_oht")
    script_tools.add_cmip_selection_cmd_args(prsr,
        both_hemispheres=True)
    script_tools.add_plotting_cmd_args(prsr, n_cbars=0,
        text_labels=[])
    prsr.add_argument("--nstddev", type=int, default=2)
    
    cmd, _ = script_tools.get_args(prsr, models=None,
        suppress_options=["exp"],
        ad_hoc_options={"nstddev": ("n_sigma for clim", "{}")})
    
    n_yr_avg_1 = cmd.yravg1[1] - cmd.yravg1[0] + 1
    n_yr_avg_2 = cmd.yravg2[1] - cmd.yravg2[0] + 1
    
    
    # ======================================================= #
    # Load data whilst creating figure 
    # ------------------------------------------------------- #
    # So that the method used to calculate delta_oht in the
    # other plotting scripts is consistent with the plot here,
    # the same function in the utils.observations module is
    # used, which has a switch "return_full_data" to return
    # everything required here rather than just the estimate of
    # delta_oht:
    load_kw = {"time_period_1": cmd.yravg1,
          "time_period_2": cmd.yravg2,
          "convergence"  : False,
          "n_std_clim_extrapolate": cmd.nstddev,
          "return_full_data": True,
          "verbosity": 2}
    
    fig, axs = plot_tools.fig_layout_2_panel_vertical(ax_b=0.17)
    
    plot_tools.set_axis_ticks_from_cmd(axs, cmd)
    
    titles = ["ECCO estimate of OHT at "]*2
    
    for j in range(2):
        
        # Load data:
        rlats = getattr(cmd, f"reflats{'ns'[j]}")
        
        (years_data, oht, j_partial_avg1, oht_partial_avg1,
            j_partial_avg2, oht_partial_avg2, years_clim,
            years_extrap_clim, oht_extrap_clim, oht_extrap_std,
            oht_extrap_range_avg1) = \
                observations.get_ecco_delta_oht(ref_lats=rlats,
                                                **load_kw)
        
        # Plot the raw OHT data:
        axs[j].plot(years_data, oht, color="k",
                    label="Raw data")
        
        # Plot thick horizontal lines to indicate the "partial"
        # averages over the two desired periods for the years
        # that are available in the data:
        line_kw = {"color": "tab:red",
            "linewidth": 2*mpl.rcParams["lines.linewidth"],
            "label": "Partial average 1 ("
                + f"{years_data[j_partial_avg1[0]]}" + u"\u2013"
                + f"{years_data[j_partial_avg1[1]]})"}
        
        axs[j].plot(years_data[j_partial_avg1],
                    [oht_partial_avg1]*2, **line_kw)
        
        line_kw["color"] = "tab:blue"
        line_kw["label"] = ("Partial average 2 ("
            + f"{years_data[j_partial_avg2[0]]}" + u"\u2013"
            + f"{years_data[j_partial_avg2[1]]})")
        
        axs[j].plot(years_data[j_partial_avg2],
                   [oht_partial_avg2]*2, **line_kw)
        
        # Shade the "unknown" area that is assumed to be the
        # climatology of the nearest available period of the same
        # length:
        axs[j].fill_between(years_extrap_clim,
            [oht_extrap_clim - cmd.nstddev*oht_extrap_std]*2,
            [oht_extrap_clim + cmd.nstddev*oht_extrap_std]*2,
            label=f"{years_clim[0]}" + u"\u2013"
                  + f"{years_clim[1]}" + " climatology",
            facecolor="tab:red", edgecolor="none", alpha=0.25)
        
        # Plot dashed lines to indicate the range of possible mean
        # OHTs for the first time period assumed based on
        # climatology for the nearest period of available data of
        # equal length to the number of missing years in average 1:
        axs[j].plot(cmd.yravg1, [oht_extrap_range_avg1[0]]*2,
                    color="tab:red", linestyle=(0, (5,5)),
                    label="Uncertainty range\n(average 1)")
        
        axs[j].plot(cmd.yravg1, [oht_extrap_range_avg1[1]]*2,
                    color="tab:red", linestyle=(0, (5,5)))
    
        titles[j] += f"{abs(rlats[j]):.0f}{u'\u00b0'}{'NS'[j]}"
    
    # Axes labels and titles:
    plot_tools.shared_yaxis_label(axs, "Poleward OHT (TW)")
    plot_tools.add_subplot_panel_labels(axs)
    plot_tools.add_subplot_panel_titles(axs, titles)
    
    # Legend below lower panel (b):
    axs[1].legend(fontsize=mpl.rcParams["legend.fontsize"]-1,
        ncol=2, handlelength=1.5, labelspacing=0.3,
        columnspacing=1.6, loc="upper left",
        bbox_to_anchor=(-0.2, -0.13))
    
    plot_tools.finish_fig(fig, savefig=cmd.savefig,
        file_name=cmd.savefigname,
        fig_metadata={"Title": cmd.savefigtitle})



if __name__ == "__main__":
    main()
