"""Script to plot derivation of EBM parameters S, Bup, and beta
for both hemispheres from pre-industrial control simulation
data (manuscript Fig. 5).
"""

import numpy as np

from utils import plot_tools, plot_style, script_tools

# Alias for long function name:
from utils.script_tools import (
    load_data_multi_model_one_hemi_two_ref_lats as ldmmohtrl)



def main():
    
    # ======================================================== #
    # Parse command-line arguments and print script parameters
    # -------------------------------------------------------- #
    prsr = script_tools.argument_parser(
        "Plots that show piControl EBM parameters S, beta, "
        + "and bup")
    script_tools.add_cmip_selection_cmd_args(prsr,
        both_hemispheres=True)
    script_tools.add_plotting_cmd_args(prsr, n_panels=6,
        n_cbars=0, text_labels=[])
     
    cmd, models = script_tools.get_args(prsr, piControl=True)
    
    n_models = len(models)
    n_yr_avg = cmd.yravg1[1] - cmd.yravg1[0] + 1
    
    
    # ======================================================== #
    # Load CMIP6 model data
    # -------------------------------------------------------- #
    diagnostics = ["f_down", "f_olr", "f_up", "f_sw_surf",
                   "iel", "tas"]
    
    load_kw = {"models"         : models,
               "n_year_averages": n_yr_avg,
               "experiment_ids" : ["piControl"],
               "diagnostics"    : diagnostics,
               "verbosity"      : 1}
    
    data_n = ldmmohtrl(hemi="n", ref_lats=cmd.reflatsn,
                       **load_kw)
    
    data_s = ldmmohtrl(hemi="s", ref_lats=cmd.reflatss,
                       **load_kw)
    
    # Only ever need concatenated arrays here (i.e., do not
    # need models separated):
    for d in diagnostics:
        data_n[d] = np.concatenate(data_n[d])
        data_s[d] = np.concatenate(data_s[d])
    
    
    # ======================================================= #
    # Create figure
    # ------------------------------------------------------- #
    
    fig, axs = plot_tools.fig_layout_6_panel()
    plot_tools.set_axis_ticks_from_cmd(axs, cmd)
    
    kw = {"marker": "x", "color": "k",
        "alpha": 0.35, "clip_on": False,
        "zorder": plot_style.get_zorder("scatter")}
    
    # Plot data:
    axs[0,0].scatter(data_n["tas"], data_n["f_up"], **kw)
    axs[0,1].scatter(data_s["tas"], data_s["f_up"], **kw)
    axs[1,0].scatter(data_n["f_down"], data_n["f_olr"], **kw)
    axs[1,1].scatter(data_s["f_down"], data_s["f_olr"], **kw)
    axs[2,0].scatter(data_n["iel"], data_n["f_sw_surf"], **kw)
    axs[2,1].scatter(data_s["iel"], data_s["f_sw_surf"], **kw)
    
    # Add y = 0 and x = 0 grid lines (must do this after
    # plotting data/setting any manual axes limits):
    plot_tools.add_zero_gridlines(axs, which="both")
    
    # Add regression lines with values/correlations printed as
    # the legend in each panel:
    legend_kw = {"bbox_to_anchor": (-0.015, 1.025),
                 "facecolor": "white"}
    
    for j, data in zip([0,1], [data_n, data_s]):
        plot_tools.add_regression_line(axs[0,j], data["tas"],
            data["f_up"], slope_symbol=r"$B_\mathrm{up}$",
            slope_units=r"W m$^{-2}$", legend_kw=legend_kw)
        plot_tools.add_regression_line(axs[1,j], data["f_down"],
            data["f_olr"], slope_symbol=r"$\beta$",
            decimal_places=(3,2), legend_kw=legend_kw)
        plot_tools.add_regression_line(axs[2,j], data["iel"],
            data["f_sw_surf"], slope_symbol=r"$S$",
            slope_units=r"W m$^{-2}$ $\degree$"
                        + "NS"[j] + r"$^{-1}$",
            legend_kw=legend_kw)
    
    # Set subplot panel labels (a/b) and titles:
    plot_tools.add_subplot_panel_labels(axs)
    plot_tools.add_subplot_panel_titles(axs[0,:],
        ["Northern hemisphere", "Southern hemisphere"])
    
    # Axis labels:
    dtas_label = r"$\Delta{}T$ (K)"
    dfup_label = r"$\Delta{}F_\mathrm{up}$ (W m$^{-2}$)"
    dfdn_label = r"$\Delta{}F_\mathrm{down}$ (W m$^{-2}$)"
    dfolr_label = r"$\Delta{}F_\mathrm{OLR}$ (W m$^{-2}$)"
    dfsw_label = r"$\Delta{}F_\mathrm{sw}$ (W m$^{-2}$)"
    dphi_label = r"$\Delta\phi_\mathrm{i}$ ($\degree$"
    # dphi_label is not incomplete; add N or S below:
    
    for j in range(2):
        axs[0,j].set_xlabel(dtas_label, labelpad=2)
        axs[0,j].set_ylabel(dfup_label, labelpad=1)
        axs[1,j].set_xlabel(dfdn_label, labelpad=2)
        axs[1,j].set_ylabel(dfolr_label, labelpad=1+7*(j==0))
        axs[2,j].set_xlabel(dphi_label+ "{})".format("NS"[j]),
                            labelpad=2)
        axs[2,j].set_ylabel(dfsw_label, labelpad=1+7*(j==0))
    
    plot_tools.finish_fig(fig, savefig=cmd.savefig,
        file_name=cmd.savefigname,
        fig_metadata={"Title": cmd.savefigtitle})



if __name__ == "__main__":
    main()
