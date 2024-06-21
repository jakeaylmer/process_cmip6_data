"""Script to plot derivation of EBM parameters bc and h for both
hemispheres from historical or future simulation data
(manuscript Figs. 6 and 7 respectively).
"""

import matplotlib as mpl
import numpy as np

from utils import plot_style, plot_tools, script_tools

# Alias for long function name:
from utils.script_tools import (
    load_data_multi_model_one_hemi_two_ref_lats as ldmmohtrl)



def main():
    
    # ======================================================== #
    # Parse command-line arguments and print script parameters
    # -------------------------------------------------------- #
    prsr = script_tools.argument_parser(
        "Plots that show variable EBM parameters bc and h")
    script_tools.add_cmip_selection_cmd_args(prsr,
        both_hemispheres=True)
    script_tools.add_plotting_cmd_args(prsr, n_panels=4,
        n_cbars=0, text_labels=[])
    
    cmd, models = script_tools.get_args(prsr)
    
    n_models = len(models)
    n_yr_avg_1 = cmd.yravg1[1] - cmd.yravg1[0] + 1
    n_yr_avg_2 = cmd.yravg2[1] - cmd.yravg2[0] + 1
    
    
    # ======================================================== #
    # Load CMIP6 model data
    # -------------------------------------------------------- #
    load_kw = {"models"        : models,
               "time_period_1" : cmd.yravg1,
               "time_period_2" : cmd.yravg2,
               "experiment_ids": ["historical", cmd.exp],
               "diagnostics"   : ["aht", "dhdt", "oht"],
               "verbosity"     : 1}
    
    data_n = ldmmohtrl(hemi="n", ref_lats=cmd.reflatsn,
                       **load_kw)
    
    data_s = ldmmohtrl(hemi="s", ref_lats=cmd.reflatss,
                       **load_kw)
    
    
    # ======================================================= #
    # Create figure
    # ------------------------------------------------------- #
    
    # This code does Fig. 6 and 7, but need to slightly tweak
    # some elements between them. Set a flag to check whether
    # we are currently likely to be running Fig. 7 (in practice
    # I obviously know, but in the spirit of keeping the code
    # "generic" and all...):
    probably_fig_7 = cmd.yravg2[1] > 2021
    
    # Create and set the figure/axes layout:
    fig, axs = plot_tools.fig_layout_4_panel(
        ax_l=0.1 if probably_fig_7 else 0.085,
        s_hor=0.1 if probably_fig_7 else 0.115)
    
    plot_tools.set_axis_ticks_from_cmd(axs, cmd)
    
    # Plot data:
    for m in range(n_models):
        
        mcol = plot_style.model_colors[models[m]]
        mmrk = plot_style.model_markers[models[m]]
        
        fc = mcol   if mmrk == "x" else "none"
        ec = None if mmrk == "x" else mcol
        
        kw = {"s":(1.3*mpl.rcParams["lines.markersize"])**2,
              "marker": mmrk, "facecolors": fc,
              "edgecolors": ec,
              **plot_style.default_scatter_kw}
        
        # There is apparently a glitch where matplotlib (or more
        # likely the SVG backend) outputs some symbols with
        # rounded corners despite joinstyle='miter'. This is
        # alleviated (somehow) by drawing symbols one by one,
        # hence the otherwise unnecessary loop over ensemble
        # members e below:
        for row, y in zip([0, 1], ["aht", "dhdt"]):
            for col, data in zip([0, 1], [data_n, data_s]):
                for e in range(len(data["oht"][m])):
                    axs[row,col].scatter(data["oht"][m][e],
                        data[y][m][e], **kw)
    
    # Add y = 0 and x = 0 grid lines (must do this after
    # plotting data/setting any manual axes limits):
    plot_tools.add_zero_gridlines(axs, which="both")
    
    # Add regression lines with values/correlations printed as
    # the legend in each panel:
    line_kw = {"color": "k"}
    fill_kw = {"facecolor": "k", "alpha": 0.25}
    
    legend_kw = [
        {"bbox_to_anchor": (1.0, 1.0), "loc": "upper right",
         "facecolor": "none" if probably_fig_7 else "white"},
        {"bbox_to_anchor": (-0.015, 1.025),
         "facecolor": "none" if probably_fig_7 else "white"}]
    
    for col, data in zip([0, 1], [data_n, data_s]):
        plot_tools.add_regression_line(axs[0,col],
            np.concatenate(data["oht"]),
            np.concatenate(data["aht"]),
            slope_symbol=r"$b_\mathrm{c}$",
            line_kw=line_kw, fill_kw=fill_kw,
            legend_kw=legend_kw[0])
        
        plot_tools.add_regression_line(axs[1,col],
            np.concatenate(data["oht"]),
            np.concatenate(data["dhdt"]),
            slope_symbol=r"$h$",
            line_kw=line_kw, fill_kw=fill_kw,
            legend_kw=legend_kw[1])
    
    # Set subplot panel labels (a/b) and titles:
    plot_tools.add_subplot_panel_labels(axs)
    plot_tools.add_subplot_panel_titles(axs[0,:],
        ["Northern hemisphere", "Southern hemisphere"])
    
    # Axis labels:
    doht_label = r"$\Delta$OHT (TW)"
    daht_label = r"$\Delta$AHT (TW)"
    ddhdt_label = r"$\Delta\partial_tH$ (TW)" 
    
    # For future figure (Fig. 7), need to adjust the position of
    # the ylabel on panel c to align with that of panel a.
    # 
    # Similarly, for both figures need to adjust the position of
    # the ylabel on panel d to align with that of panel a (but
    # by different amounts):
    dlabelpadc = 7 if probably_fig_7 else 0
    dlabelpadd = 5 if probably_fig_7 else 7
    
    for col in range(2):
        axs[0,col].set_xlabel(doht_label, labelpad=2)
        axs[0,col].set_ylabel(daht_label, labelpad=0)
        axs[1,col].set_xlabel(doht_label, labelpad=2)
        axs[1,col].set_ylabel(ddhdt_label,
                              labelpad=0 + dlabelpadc*(col==0)
                              + dlabelpadd*(col==1))
    
    plot_tools.finish_fig(fig, savefig=cmd.savefig,
        file_name=cmd.savefigname,
        fig_metadata={"Title": cmd.savefigtitle})



if __name__ == "__main__":
    main()
