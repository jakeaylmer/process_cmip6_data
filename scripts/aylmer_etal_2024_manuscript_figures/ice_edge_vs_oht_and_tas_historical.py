"""Script to plot change in sea ice edge against change in
polar surface temperature and poleward ocean heat transport over
the historical period with observations and CMIP6 data
(manuscript Figs. 1 and 3, and supplementary information Figs.
S5 and S6).
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from utils import (ebm, observations, plot_style, plot_tools,
                   script_tools)

# Alias for long function name:
from utils.script_tools import (
    load_data_multi_model_one_hemi_two_ref_lats as ldmmohtrl)



def main():
    
    # ======================================================== #
    # Parse command-line arguments and print script parameters
    # -------------------------------------------------------- #
    prsr = script_tools.argument_parser(
        "Plot delta_iel against delta_oht and delta_tas with "
        + "observations/reanalyses")
    script_tools.add_cmip_selection_cmd_args(prsr)
    script_tools.add_ebm_fitting_cmd_args(prsr)
    script_tools.add_plotting_cmd_args(prsr)
    
    cmd, models = script_tools.get_args(prsr)
    
    n_models = len(models)
    n_yr_avg_1 = cmd.yravg1[1] - cmd.yravg1[0] + 1
    n_yr_avg_2 = cmd.yravg2[1] - cmd.yravg2[0] + 1
    
    # Determine which hemisphere to analyse from the input
    # reference latitudes:
    if cmd.reflats[0] < 0.0 and cmd.reflats[1] < 0.0:
        hemi = "s"
    else:
        hemi = "n"
    
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
     
    data = ldmmohtrl(models, diagnostics, hemi,
        ref_lats=cmd.reflats,
        time_period_1=(cmd.yravg1[0], cmd.yravg1[1]),
        time_period_2=(cmd.yravg2[0], cmd.yravg2[1]),
        experiment_ids=("historical", cmd.exp),
        verbosity=1)
    
    n_total_members = sum([len(data["aht"][m])
                           for m in range(n_models)])
    
    # Create a copy of the data needed for the plot in case
    # we have to remove any data for models to exclude when
    # EBM fitting, next step:
    tas_plot = data["tas"].copy()
    oht_plot = data["oht"].copy()
    iel_plot = data["iel"].copy()
    
    # Remove models from the data array if excluding any from
    # EBM fitting. Start at the end (reverse alphabetical order
    # by model name) so that deleting can be done on the fly:
    for x in sorted(cmd.exclude_from_ebm)[::-1]:
        j = models.index(x)
        for d in data.keys():  # i.e., every diagnostic
            del data[d][j]
    
    # Dictionary of EBM parameters
    # 
    # "name": (value, error) [OR] (None, None) [OR] "piControl"
    # 
    # corresponding to options:
    #     1) providing value directly
    #     2) calculate from input historical/future data
    #     3) calculate from pre-industrial control data
    # respectively. Then they get over-written as appropriate
    # during the ebm.fit() function:
    # 
    ebm_params = {"bc"  : (None, None),
                  "beta": (None, None),
                  "bup" : (None, None),
                  "h"   : (None, None),
                  "S"   : (None, None)}
    
    for p in ebm_params.keys():
        if p in cmd.piparams:
            ebm_params[p] = "piControl"
    
    ebm_params, ebm_results = ebm.fit(data,
        ebm_params=ebm_params, ref_lats=cmd.reflats,
        hemi=hemi, piparams_n_year_averages=n_yr_avg_1)
    
    # Convert OHT to OHTC for plotting (NB: this is approximate
    # as it does not account for land, but can be considered a
    # scaling for units/to put in a more intuitive form for what
    # is being represented, i.e., the difference in OHT across
    # the two reference latitudes. It doesn't affect the
    # validity of the underlying calculations, it's just the
    # plot/axis scale which is approximate):
    # 
    if abs(cmd.reflats[1]) < 85.0:
        
        # Approximate area of ocean between reference latitudes:
        A_approx = ebm.area_factor(cmd.reflats)
        
        # Adjust/scale relevant data and EBM results:
        for m in range(len(oht_plot)):
            oht_plot[m] *= 1.0E12 / A_approx  # TW --> W m-2
        ebm_results["diel/doht_ebm"] = (
            ebm_results["diel/doht_ebm"][0] *A_approx*1.0E-12,
            ebm_results["diel/doht_ebm"][1] *A_approx*1.0E-12)
        
        oht_xlabel = "OHT Convergence change, "
        oht_xlabel += r"$\Delta$OHTC (W m$^{-2}$)"
        oht_clabel = r"$\Delta$OHTC (W m$^{-2}$)"
        
        oht_convergence = True
    
    else:
        
        oht_xlabel = r"Poleward OHT change, $\Delta$OHT (TW)"
        oht_clabel = r"$\Delta$OHT (TW)"
        
        oht_convergence = False
    
    
    # ======================================================== #
    # Load observation data:
    # -------------------------------------------------------- #
    diel_obs = observations.get_passive_microwave_delta_iel(
        hemi, time_period_1=cmd.yravg1,time_period_2=cmd.yravg2)
    
    reanalyses = ["CFSR/CFSv2", "ERA5", "JRA-55", "MERRA-2"]
    dtas_obs = observations.get_reanalysis_delta_tas(
        reanalyses, ref_lats=cmd.reflats,
        time_period_1=cmd.yravg1, time_period_2=cmd.yravg2)
    
    doht_obs, doht_obs_err = observations.get_ecco_delta_oht(
        ref_lats=cmd.reflats, time_period_1=cmd.yravg1,
        time_period_2=cmd.yravg2,
        convergence=oht_convergence)
    
    # In the paper we quote the proportion of simulations that
    # fall strictly within the limits of the observational
    # estimates on the plot:
    for xd, yd, xo, yo, xtxt, ytxt in zip(
            (tas_plot, oht_plot, tas_plot),
            (iel_plot, iel_plot, oht_plot),
            (dtas_obs, doht_obs_err, dtas_obs),
            (diel_obs, diel_obs, doht_obs_err),
            ("tas", "oht", "tas"), ("iel", "iel", "oht")):
        n = 0
        for m in range(n_models):
            n += np.sum(
                (xd[m] > np.min(xo)) & (xd[m] < np.max(xo))
                & (yd[m] > np.min(yo)) & (yd[m] < np.max(yo)))
        
        print(f"{n} simulations within delta_{xtxt} and "
              + f"delta_{ytxt} observations "
              + f"({100.0*n/n_total_members:.1f}%)")
    
    # Repeat for the combination of all three:
    n = 0
    for m in range(n_models):
        n += np.sum(
            (tas_plot[m] > np.min(dtas_obs))
            & (tas_plot[m] < np.max(dtas_obs))
            & (iel_plot[m] > np.min(diel_obs))
            & (iel_plot[m] < np.max(diel_obs))
            & (oht_plot[m] > np.min(doht_obs_err))
            & (oht_plot[m] < np.max(doht_obs_err)))
    
    print(f"{n} simulations within delta_tas, delta_iel, and "
          + f"delta_oht observations "
          + f"({100.0*n/n_total_members:.1f}%)")
    
    
    # ======================================================= #
    # Create figure
    # ------------------------------------------------------- #
    
    fig, axs, cbar_axs = plot_tools.fig_layout_2_panel_2_cbar()
    
    plot_tools.set_axis_ticks_from_cmd(axs, cmd)
    
    # Set up discrete colormaps with boundaries set from
    # command-line arguments:
    cm_cmaps = [plot_style.default_cmap_delta_oht,
                plot_style.default_cmap_delta_tas]
    cm_norms = []
    for p in range(2):
        cats = getattr(cmd, f"{'ab'[p]}_cbar_cats")
        levs = np.arange(cats[0], cats[1]+0.5*cats[2], cats[2])
        cm_norms.append(mpl.colors.BoundaryNorm(levs,
                                                ncolors=256,
                                                extend="both"))
    
    # Plot data:
    for m in range(len(models)):
        
        if exclude_from_ebm[m]:
            
            marker_excl = plot_style.marker_exclude[
                cmd.exclude_from_ebm.index(models[m])
                % len(plot_style.marker_exclude)]
            
            for j, xdata in zip((0, 1), (tas_plot, oht_plot)):
                axs[j].scatter(xdata[m], iel_plot[m],
                    marker=marker_excl,
                    **plot_style.default_scatter_exclude_kw)
        else:
            
            xdata = [tas_plot, oht_plot]
            
            for j in range(2):
                axs[j].scatter(xdata[j][m], iel_plot[m],
                    c=xdata[(j+1)%2][m], cmap=cm_cmaps[j],
                    norm=cm_norms[j], edgecolors="none",
                    **plot_style.default_scatter_kw)
    
    # Add colorbars:
    cbars = np.array([fig.colorbar(plt.cm.ScalarMappable(
        cmap=cm_cmaps[j], norm=cm_norms[j]),
        cax=cbar_axs[j], extend="both",
        orientation="horizontal") for j in range(2)])
    
    for j, label in zip((0,1), (oht_clabel,r"$\Delta{}T$ (K)")):
        cbar_axs[j].set_xlabel(label, labelpad=2, loc="center")
        
        cticks = getattr(cmd, f"{'ab'[j]}_cbar_ticks")
        cbars[j].set_ticks(np.arange(cticks[0],
                                     cticks[1] + 0.5*cticks[2],
                                     cticks[2]))
        
        cbar_axs[j].tick_params(which="both", axis="both",
            left=False, right=False, top=False, bottom=False)
    
    
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
    for xdata, xvar, j in zip([tas_plot, oht_plot],
                              ["tas", "oht"], [0, 1]):
        
        xlim_ebm = np.array([np.nanmin(np.concatenate(xdata)),
                             np.nanmax(np.concatenate(xdata))])
        
        xy_em = (np.nanmean(np.concatenate(xdata)),
                   np.nanmean(np.concatenate(iel_plot)))
    
        plot_tools.add_ebm_line(axs[j],
            ebm_results[f"diel/d{xvar}_ebm"], xlim_ebm, xy_em,
            dr_label_text=getattr(cmd,f"{'ab'[j]}_dr_ebm_text"))
        
        plot_tools.add_text_label(axs[j], "CMIP6", xy_em,
            slope=ebm_results[f"diel/d{xvar}_data"][0],
            dr_label_text=getattr(cmd,f"{'ab'[j]}_dr_cmip6_text"),
            annotate_kw=plot_style.default_cmip6_annotate_kw)
    
    for x0, j in zip((np.mean(dtas_obs), doht_obs), (0, 1)):
        plot_tools.add_text_label(axs[j],
            "Observations/\nreanalyses",
            (x0, np.mean(diel_obs)),
            dr_label_text=getattr(cmd,f"{'ab'[j]}_dr_obs_text"),
            annotate_kw=plot_style.default_obs_annotate_kw)
    
    
    # Add observations/reanalyses `error' bars; reanalyses for
    # panel (a):
    plot_tools.add_points_along_error_bar(axs[0], dtas_obs,
        np.mean(diel_obs), labels=reanalyses,
        legend_kw={"bbox_to_anchor": (1.15, -0.16),
                   "loc": "upper right"})
    
    # OHT estimate for panel (b):
    plot_tools.add_range_as_error_bar(axs[1], doht_obs_err,
        np.mean(diel_obs), which="x")
    
    for j, xpos, xdata in zip((0, 1),
                              (np.mean(dtas_obs), doht_obs),
                              (dtas_obs, doht_obs_err)):
        
        # Ice-edge latitude estimate (both panels):
        plot_tools.add_range_as_error_bar(axs[j], diel_obs,
            xpos, which="y")
        
        # Pale-grey rectangle to highlight coincidence of each
        # observational estimate pair in each panel:
        axs[j].add_patch(mpl.patches.Rectangle(
            (np.min(xdata), np.min(diel_obs)),
            np.max(xdata) - np.min(xdata),
            np.max(diel_obs) - np.min(diel_obs),
            **plot_style.default_obs_rectangle_kw))
    
    # Add y = 0 and x = 0 grid lines (must do this after
    # plotting data/setting any manual axes limits):
    plot_tools.add_zero_gridlines(axs, which="both")
    
    # Set axes titles (panel labels a/b) and axis labels:
    plot_tools.add_subplot_panel_labels(axs)
    
    axs[0].set_xlabel("Surface temperature change, "
                      + r"$\Delta{}T$ (K)", labelpad=2)
    
    # x-axis label on panel b needs to be a bit smaller if
    # it is OHTC (the label is longer):
    axs[1].set_xlabel(oht_xlabel, labelpad=2,
        fontsize=mpl.rcParams["axes.labelsize"]
                 -1.5*(oht_convergence))
    
    axs[0].set_ylabel("Sea ice-edge change, "
                      + r"$\Delta\phi_\mathrm{i}$ (" + u"\u00b0"
                      + hemi.upper() + ")",
                      labelpad=6 if hemi=="n" else 2)
    
    axs[1].set_ylabel(r"$\Delta\phi_\mathrm{i}$ (" + u"\u00b0"
                      + hemi.upper() + ")", labelpad=2)
    
    plot_tools.finish_fig(fig, savefig=cmd.savefig,
        file_name=cmd.savefigname,
        fig_metadata={"Title": cmd.savefigtitle})



if __name__ == "__main__":
    main()
