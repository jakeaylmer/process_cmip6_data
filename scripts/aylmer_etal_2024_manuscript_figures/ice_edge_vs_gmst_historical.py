"""Script to plot change in sea ice edge against global mean
surface temperature change over the historical period with
observations and CMIP6 data (supplementary information Fig. S1).
"""

import matplotlib as mpl
import numpy as np

from utils import (maths, observations, plot_style, plot_tools,
                   script_tools)



def main():
    
    # ======================================================== #
    # Parse command-line arguments and print script parameters
    # -------------------------------------------------------- #
    prsr = script_tools.argument_parser(
        "Plot delta_iel against delta_gmst with "
        + "observations/reanalyses")
    script_tools.add_cmip_selection_cmd_args(prsr)
    script_tools.add_plotting_cmd_args(prsr, n_cbars=0,
        text_labels=["cmip6", "obs"])
    
    cmd, models = script_tools.get_args(prsr)
    
    n_models = len(models)
    n_yr_avg_1 = cmd.yravg1[1] - cmd.yravg1[0] + 1
    n_yr_avg_2 = cmd.yravg2[1] - cmd.yravg2[0] + 1
    
    
    # ======================================================== #
    # Load CMIP6 model data
    # -------------------------------------------------------- #
    data_n, data_s = script_tools.load_data_multi_model(
        models, ["iel", "tas"], lat_eval=0.0,
        time_period_1=cmd.yravg1, time_period_2=cmd.yravg2,
        experiment_ids=["historical", cmd.exp])

    # Global mean temperature was not saved directly; above,
    # load the hemisphere averages (ref_lat = 0.0 in both
    # hemispheres) and compute the global means from them:
    data_n["gmst"] = []
    data_s["gmst"] = []
    for m in range(n_models):
        gmst_m = np.mean(
            np.vstack((data_n["tas"][m], data_s["tas"][m])),
            axis=0)
        data_n["gmst"].append(gmst_m)
        data_s["gmst"].append(gmst_m)
    
    n_total_members = sum([len(data_n["gmst"][m])
                           for m in range(n_models)])
    
    
    # ======================================================== #
    # Load observation data:
    # -------------------------------------------------------- #
    diel_n_obs = observations.get_passive_microwave_delta_iel(
        "n", time_period_1=cmd.yravg1,time_period_2=cmd.yravg2)
    
    diel_s_obs = observations.get_passive_microwave_delta_iel(
        "s", time_period_1=cmd.yravg1,time_period_2=cmd.yravg2)
    
    reanalyses = ["CFSR/CFSv2", "ERA5", "JRA-55", "MERRA-2"]
    dgmst_obs = observations.get_reanalysis_delta_tas(
        reanalyses, global_mean=True, time_period_1=cmd.yravg1,
        time_period_2=cmd.yravg2)
    
    # In the paper we quote the proportion of simulations that
    # fall strictly within the limits of the observational
    # estimates on the plot:
    for yd, yo, ytxt in zip(
            (data_n["iel"], data_s["iel"]),
            (diel_n_obs, diel_s_obs),
            ("iel_n", "iel_s")):
        n = 0
        for m in range(n_models):
            n += np.sum(
                (data_n["gmst"][m] > np.min(dgmst_obs))
                & (data_n["gmst"][m] < np.max(dgmst_obs))
                & (yd[m] > np.min(yo)) & (yd[m] < np.max(yo)))
        
        print(f"{n} simulations within delta_gmst and "
              + f"delta_{ytxt} observations "
              + f"({100.0*n/n_total_members:.1f}%)")
    
    # Repeat for the combination agreeing in both hemispheres:
    n = 0
    for m in range(n_models):
        n += np.sum(
            (data_n["gmst"][m] > np.min(dgmst_obs))
            & (data_n["gmst"][m] < np.max(dgmst_obs))
            & (data_n["iel"][m] > np.min(diel_n_obs))
            & (data_n["iel"][m] < np.max(diel_n_obs))
            & (data_s["iel"][m] > np.min(diel_s_obs))
            & (data_s["iel"][m] < np.max(diel_s_obs)))
    
    print(f"{n} simulations within delta_gmst, delta_iel_n, "
          + f"and delta_iel_s observations "
          + f"({100.0*n/n_total_members:.1f}%)")
    
    
    # ======================================================= #
    # Create figure
    # ------------------------------------------------------- #
    fig, axs = plot_tools.fig_layout_2_panel()
    plot_tools.set_axis_ticks_from_cmd(axs, cmd)
    
    # Plot data:
    for m in range(len(models)):
        for j, ydata in zip((0, 1),
                            [data_n["iel"], data_s["iel"]]):
            axs[j].scatter(data_n["gmst"][m], ydata[m],
                facecolors=plot_style.default_cmip6_annotate_kw[
                    "color"], edgecolors="none",
                **plot_style.default_scatter_kw)
    
    # Add "CMIP6" labels to each panel:
    for ydata, j in zip([data_n["iel"], data_s["iel"]], [0, 1]):
        
        xy_em = (np.nanmean(np.concatenate(data_n["gmst"])),
                 np.nanmean(np.concatenate(ydata)))
        
        slope = maths.orthogonal_distance_regression(
            np.concatenate(data_n["gmst"]),
            np.concatenate(ydata))[2]
        
        plot_tools.add_text_label(axs[j], "CMIP6", xy_em,
            slope=slope,
            dr_label_text=getattr(cmd,f"{'ab'[j]}_dr_cmip6_text"),
            annotate_kw=plot_style.default_cmip6_annotate_kw)
    
    for y0, j in zip((np.mean(diel_n_obs), np.mean(diel_s_obs)),
                     (0, 1)):
        
        plot_tools.add_text_label(axs[j],
            "Observations/\nreanalyses",
            (np.mean(dgmst_obs), y0),
            dr_label_text=getattr(cmd,f"{'ab'[j]}_dr_obs_text"),
            annotate_kw=plot_style.default_obs_annotate_kw)
        
        # Add observations/reanalyses 'error' bars:
        plot_tools.add_points_along_error_bar(axs[j], dgmst_obs,
            y0, labels=reanalyses, add_legend=j==0,
            legend_kw={"loc": "upper left",
                       "bbox_to_anchor": (-0.03,1.03)})
    
    for j, ydata in zip((0, 1), (diel_n_obs, diel_s_obs)):
        
        plot_tools.add_range_as_error_bar(axs[j], ydata,
            np.mean(dgmst_obs), which="y")
        
        # Pale-grey rectangle to highlight coincidence of each
        # observational estimate pair in each panel:
        axs[j].add_patch(mpl.patches.Rectangle(
            (np.min(dgmst_obs), np.min(ydata)),
            np.max(dgmst_obs) - np.min(dgmst_obs),
            np.max(ydata) - np.min(ydata),
            **plot_style.default_obs_rectangle_kw))
    
    # Add y = 0 and x = 0 grid lines (must do this after
    # plotting data/setting any manual axes limits):
    plot_tools.add_zero_gridlines(axs, which="both")
    
    # Set subplot panel labels (a/b) and titles:
    plot_tools.add_subplot_panel_labels(axs)
    plot_tools.add_subplot_panel_titles(axs,
        ["Arctic Ocean", "Southern Ocean"])
    
    # Add shared x-axis label, centered half way between the
    # left edge of the left panel and right edge of the right
    # panel:
    xtxt0 = axs[0].get_position().x0
    xtxt1 = axs[1].get_position().x1
    ytxt0 = axs[0].get_position().y0
    fig.text(xtxt0 + 0.5*(xtxt1 - xtxt0), 0.28*ytxt0,
             "Global-mean surface air temperature change (K)",
             ha="center", va="center",
             fontsize=mpl.rcParams["axes.labelsize"])
    
    axs[0].set_ylabel("Sea ice-edge change, "
                      + r"$\Delta\phi_\mathrm{i}$ (" + u"\u00b0"
                      + "N)", labelpad=6)
    
    axs[1].set_ylabel(r"$\Delta\phi_\mathrm{i}$ (" + u"\u00b0"
                      + "S)", labelpad=2)
    
    plot_tools.finish_fig(fig, savefig=cmd.savefig,
        file_name=cmd.savefigname,
        fig_metadata={"Title": cmd.savefigtitle})



if __name__ == "__main__":
    main()
