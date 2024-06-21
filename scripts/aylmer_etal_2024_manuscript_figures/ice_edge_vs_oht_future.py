"""Script to plot change in sea ice edge against change in
poleward ocean heat transport in a near-future projection in
CMIP6 (manuscript Fig. 4).
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from utils import ebm, plot_style, plot_tools, script_tools

# Alias for long function name:
from utils.script_tools import (
    load_data_multi_model_one_hemi_two_ref_lats as ldmmohtrl)



def main():
    
    # ======================================================== #
    # Parse command-line arguments and print script parameters
    # -------------------------------------------------------- #
    prsr = script_tools.argument_parser(
        "Plot delta_iel against delta_oht for near-future "
        + "projection with observations/reanalyses")
    script_tools.add_cmip_selection_cmd_args(prsr,
        default_future=True, both_hemispheres=True)
    script_tools.add_ebm_fitting_cmd_args(prsr)
    script_tools.add_plotting_cmd_args(prsr, n_cbars=0,
        text_labels=["ebm"])
    
    cmd, models = script_tools.get_args(prsr)
    
    n_models = len(models)
    n_yr_avg_1 = cmd.yravg1[1] - cmd.yravg1[0] + 1
    n_yr_avg_2 = cmd.yravg2[1] - cmd.yravg2[0] + 1
    
    # Boolean array, whether model m is excluded from EBM
    # fitting (if True) or not:
    exclude_from_ebm = \
        np.array([x in cmd.exclude_from_ebm for x in models],
                       dtype=bool)
    
    # ======================================================== #
    # Load CMIP6 model data
    # -------------------------------------------------------- #
    
    diagnostics = ["aht", "dhdt", "f_down", "f_olr",
                   "f_sw_surf", "f_up", "iel", "oht", "tas"]
    
    kw = {"time_period_1": cmd.yravg1,
          "time_period_2": cmd.yravg2,
          "experiment_ids": ["historical", cmd.exp],
          "verbosity": 1}
    
    data_n = ldmmohtrl(models, diagnostics, "n",
        ref_lats=cmd.reflatsn, **kw)
    data_s = ldmmohtrl(models, diagnostics, "s",
        ref_lats=cmd.reflatss, **kw)
    
    # Create a copy of the data needed for the plot in case
    # we have to remove any data for models to exclude when
    # EBM fitting, next step:
    oht_n_plot = data_n["oht"].copy()
    iel_n_plot = data_n["iel"].copy()
    oht_s_plot = data_s["oht"].copy()
    iel_s_plot = data_s["iel"].copy()
    
    # Remove models from the data array if excluding any from
    # EBM fitting. Start at the end (reverse alphabetical order
    # by model name) so that deleting can be done on the fly:
    for x in sorted(cmd.exclude_from_ebm)[::-1]:
        j = models.index(x)
        for d in data_n.keys():  # i.e., every diagnostic
            del data_n[d][j]
            del data_s[d][j]
    
    # Dictionary of EBM parameters
    # 
    # "name": (value, error) [OR] (None, None) [OR] "piControl"
    # 
    # corresponding to options:
    #     1) providing value directly
    #     2) calculate from input historical/future data
    #     3) calculate from pre-industrial control data
    # respectively. Then they get over-written as appropriate
    # during the ebm.fit() routine:
    # 
    ebm_params = {"bc"  : (None, None),
                  "beta": (None, None),
                  "bup" : (None, None),
                  "h"   : (None, None),
                  "S"   : (None, None)}
    
    for p in ebm_params.keys():
        if p in cmd.piparams:
            ebm_params[p] = "piControl"
    
    ebm_params_n, ebm_results_n = ebm.fit(data_n,
        ebm_params=ebm_params.copy(), ref_lats=cmd.reflatsn,
        hemi="n", piparams_n_year_averages=n_yr_avg_1)
    
    ebm_params_s, ebm_results_s = ebm.fit(data_s,
        ebm_params=ebm_params.copy(), ref_lats=cmd.reflatss,
        hemi="s", piparams_n_year_averages=n_yr_avg_1)
    
    # Convert OHT to OHTC for plotting (NB: this is approximate
    # as it does not account for land, but can be considered a
    # scaling for units/to put in a more intuitive form for what
    # is being represented, i.e., the difference in OHT across
    # the two reference latitudes. It doesn't affect the
    # validity of the underlying calculations, it's just the
    # plot/axis scale which is approximate):
    
    # Labels and flag for each hemisphere (figure panel);
    # by default set to OHT case, then in loop below change
    # to OHTC case if required:
    oht_xlabel = [r"Poleward OHT change, $\Delta$OHT (TW)"]*2
    oht_convergence = [False, False]
    
    for h, oht_data, ebm_results in zip("ns",
        [oht_n_plot, oht_s_plot],
        [ebm_results_n, ebm_results_s]):
        
        reflats = getattr(cmd, f"reflats{h}")
        
        if abs(reflats[1]) < 85.0:
            
            # Approximate area of ocean between reference
            # latitudes:
            A_approx = ebm.area_factor(reflats)
            
            # Adjust/scale relevant data and EBM results:
            for m in range(len(oht_data)):
                oht_data[m] *= 1.0E12 / A_approx  # TW --> W m-2
            ebm_results["diel/doht_ebm"] = (
                ebm_results["diel/doht_ebm"][0]*A_approx*1.E-12,
                ebm_results["diel/doht_ebm"][1]*A_approx*1.E-12)
            
            oht_xlabel[h=="s"] = r"OHT convergence change, "
            oht_xlabel[h=="s"] += r"$\Delta$OHTC (W m$^{-2}$)"
            oht_convergence[h=="s"] = True
        
    # ======================================================= #
    # Create figure
    # ------------------------------------------------------- #
    
    fig, axs = plot_tools.fig_layout_2_panel(ax_r=0.24)
    plot_tools.set_axis_ticks_from_cmd(axs, cmd)
    
    # Plot data:
    for m in range(n_models):
        
        if exclude_from_ebm[m]:
            
            marker_excl = plot_style.marker_exclude[
                cmd.exclude_from_ebm.index(models[m])
                % len(plot_style.marker_exclude)]
            
            for xdata, ydata, j in zip(
                    [oht_n_plot, oht_s_plot],
                    [iel_n_plot, iel_s_plot], (0, 1)):
                
                axs[j].scatter(xdata[m], ydata[m],
                    s=(1.3*mpl.rcParams["lines.markersize"])**2,
                    marker=marker_excl, facecolors="none",
                    edgecolors="k",
                    label=models[m] if j==1 else None,
                    **plot_style.default_scatter_exclude_kw)
        else:
            
            mcol = plot_style.model_colors[models[m]]
            mmrk = plot_style.model_markers[models[m]]
            
            fc = mcol   if mmrk == "x" else "none"
            ec = None if mmrk == "x" else mcol
            
            kw = {"s":(1.3*mpl.rcParams["lines.markersize"])**2,
                  "marker": mmrk, "facecolors": fc,
                  "edgecolors": ec,
                  **plot_style.default_scatter_kw}
            
            for xdata, ydata, j in zip(
                    [oht_n_plot, oht_s_plot],
                    [iel_n_plot, iel_s_plot], (0, 1)):
                
                # There is apparently a glitch where matplotlib
                # (or more likely the SVG backend) outputs some
                # symbols with rounded corners despite
                # joinstyle='miter'. This is alleviated
                # (somehow) by drawing symbols one by one, hence
                # the otherwise unnecessary loop over ensemble
                # members e below:
                for e in range(len(xdata[m])):
                    axs[j].scatter(xdata[m][e], ydata[m][e],
                        label=models[m] if (j==1 and e==0)
                              else None, **kw)
    
    # Add EBM line, shading, and labels, and "CMIP6" labels, to
    # each panel.
    # 
    # The EBM line is plotted between the minimum and maximum
    # x-data values (tas for panel a, OHT for panel b), and 
    # through the ensemble-mean data point.
    # 
    # Text label offsets are set on the command line, e.g.,
    # --a-dr-ebm-text 0 0 interpreted as dx = dy = 0 in data
    # coordinates on panel a (with default values):
    for xdata, ydata, ebm_results, j in zip(
            [oht_n_plot, oht_s_plot],
            [iel_n_plot, iel_s_plot],
            [ebm_results_n, ebm_results_s], [0,1]):
        
        xlim_ebm = np.array([np.nanmin(np.concatenate(xdata)),
                             np.nanmax(np.concatenate(xdata))])
        
        xy_em = (np.nanmean(np.concatenate(xdata)),
                 np.nanmean(np.concatenate(ydata)))
    
        plot_tools.add_ebm_line(axs[j],
            ebm_results[f"diel/doht_ebm"], xlim_ebm, xy_em,
            dr_label_text=getattr(cmd,f"{'ab'[j]}_dr_ebm_text"))
    
    # Add y = 0 and x = 0 grid lines (must do this after
    # plotting data/setting any manual axes limits):
    plot_tools.add_zero_gridlines(axs, which="both")
    
    # Legend to the right of second panel:
    leg_kw = {"bbox_to_anchor": (1.02, 1.15),
              "loc": "upper left",
              "fontsize": mpl.rcParams["legend.fontsize"]+1,
              "markerscale": 1.0, "labelspacing": 0.18}
    axs[1].legend(**plot_tools._add_default_kw(
                  leg_kw, plot_style.default_legend_kw))
    
    # Set axes titles (panel labels a/b) and axis labels:
    plot_tools.add_subplot_panel_labels(axs)
    plot_tools.add_subplot_panel_titles(axs,
        ["Arctic Ocean", "Southern Ocean"])
    
    if oht_xlabel[0] == oht_xlabel[1]:
        xtxt0 = axs[0].get_position().x0
        xtxt1 = axs[1].get_position().x1
        ytxt0 = axs[0].get_position().y0
        fig.text(xtxt0 + 0.5*(xtxt1 - xtxt0), 0.28*ytxt0,
                 oht_xlabel[0], ha="center", va="center",
                 fontsize=mpl.rcParams["axes.labelsize"])
    else:
        for j in range(2):
            axs[j].set_xlabel(oht_xlabel[j], labelpad=2,
                fontsize=mpl.rcParams["axes.labelsize"]
                         -1.5*(oht_convergence[j]))
    
    axs[0].set_ylabel("Sea ice-edge change, "
                      + r"$\Delta\phi_\mathrm{i}$ (" + u"\u00b0"
                      + "N)", labelpad=3)
    
    axs[1].set_ylabel(r"$\Delta\phi_\mathrm{i}$ (" + u"\u00b0"
                      + "S)", labelpad=3)
    
    plot_tools.finish_fig(fig, savefig=cmd.savefig,
        file_name=cmd.savefigname,
        fig_metadata={"Title": cmd.savefigtitle})



if __name__ == "__main__":
    main()
