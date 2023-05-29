import numpy as np

from process_cmip6_data.src import metadata as md
from process_cmip6_data.src import netcdf as nf
from process_cmip6_data.src import qc
from process_cmip6_data.src import script_tools

from my_python_utilities.data_tools import nc_tools as nct


def main():
    
    prsr = script_tools.qc_argument_parser()
    prsr.add_argument("--normalised", action="store_true",
                      help="Plot normalised (W m-2) diagnostics")
    cmd = prsr.parse_args()
    
    if cmd.normalised:
        units = "_wm-2"
        units_label = nf.field_units["heatflux"]
        dfactor = 1.0
    else:
        units = ""
        units_label = "TW"
        dfactor = 1000.0
    
    # Can't reasonably plot all members for multiple methods
    # on one plot:
    if cmd.ensmemberid < 0:
        cmd.ensmemberid = 0
        print("WARNING: ignoring --ensmemberid (plotting first "
              + "member only)")
    
    # For scatter plots, plot cmd.yravgs anomalies, and for
    # historical or ssp* simulations plot time differences:
    if cmd.experiment != "piControl":
        print("WARNING: ignoring --ensmemberid for scatter "
              + f"plots; plotting deltas ({cmd.yravg1[0]}-"
              + f"{cmd.yravg1[1]} to {cmd.yravg2[0]}-"
              + f"{cmd.yravg2[1]}) for all members" + "\n"
              + "Use --yravg1 Y1 Y2 and --yravg2 Y1 Y2 to "
              + "change time-averaging periods")
    
    script_tools._check_qc_experiment_args(cmd)
    
    try:
        ftype = md.ocn_prognostic_temperature[cmd.model]
    except KeyError:
        ftype = md.ocn_prognostic_temperature_when_does_not_exist
    
    gd_kw = {
        "model_id": cmd.model,
        "experiment_1": cmd.experiment,
        "experiment_2": cmd.x2,
        "lat_eval": cmd.latitude
    }
    
    # Load from aliases (defined in api._diagnostic_definitions)
    dhdt_resid_hfbasin_n, dhdt_resid_hfbasin_s, year, r_vals, \
        ripf, lat_eval_n, lat_eval_s, data_path_resid_hfbasin \
        = qc.get_data(f"dhdt_{ftype}_residual_hfbasin_cc{units}",
                      **gd_kw)
    
    dhdt_resid_hfx_n, dhdt_resid_hfx_s, _, _, _, _, _, \
        data_path_resid_hfx = qc.get_data(
            f"dhdt_{ftype}_residual_hfx_cc{units}", **gd_kw)
    
    dhdt_direct_n, dhdt_direct_s, _, _, _, _, _, \
        data_path_direct = qc.get_data(
            f"dhdt_{ftype}_direct_cc{units}", **gd_kw)
    
    data_paths = qc.shorten_data_path_labels(
        [data_path_resid_hfbasin, data_path_resid_hfx,
         data_path_direct])
    
    e, ens_label = qc.ens_member_label(cmd.ensmemberid,
                                       r_vals, ripf)
    
    n_ens = len(ripf)
    
    # ------------------------------------------------------- #
    
    exp_label = cmd.experiment
    exp_label += "" if cmd.x2 == "none" else f"_{cmd.x2}"
    
    descr_kw = {
        "model_id"                    : cmd.model,
        "experiment_id"               : exp_label,
        "ensemble_members_description": ens_label,
        "diagnostic_description"      :
            f"o{ftype}temptend{units}_compare_methods",
        "plot_type"                   : "time_series_"
                                        + f"{lat_eval_n:.0f}n"
    }
    
    # Figs. 1-2: time series, showing full yearly data and low-
    # pass (moving-average) filtered data for each method
    # -------------------------------------------------------- #
    fig1, ax1 = qc.start_figure(**descr_kw)
    
    descr_kw["plot_type"] =f"time_series_{abs(lat_eval_s):.0f}s"
    fig2, ax2 = qc.start_figure(**descr_kw)
    
    # Put data in TW (assume I have not changed raw data from
    # being in PW by default...)
    ydata_n = [dfactor*dhdt_resid_hfbasin_n,
               dfactor*dhdt_resid_hfx_n, dfactor*dhdt_direct_n]
    ydata_s = [dfactor*dhdt_resid_hfbasin_s,
               dfactor*dhdt_resid_hfx_s, dfactor*dhdt_direct_s]
    
    year_ma = qc.moving_average(year, ydata_n[0], cmd.movavg)[0]
    
    ydata_ma_n = [qc.moving_average(year, x, cmd.movavg)[1]
                  for x in ydata_n]
    ydata_ma_s = [qc.moving_average(year, x, cmd.movavg)[1]
                  for x in ydata_s]
    
    labels = ["residual (hfbasin)", "residual (hfx/hfy)",
              "direct"]
    colors = ["tab:blue", "tab:orange", "tab:green"]
    
    for ax, ydata, ydata_ma, lat_eval, hemi in zip(
            [ax1, ax2],
            [ydata_n, ydata_s],
            [ydata_ma_n, ydata_ma_s],
            [lat_eval_n, lat_eval_s],
            ["n", "s"]
        ):
        
        # Full time series:
        qc.plot_data(ax, xdata=year,
            ydata=ydata, ens_index=e, ens_label=ens_label,
            labels=[None]*3, colors=colors, linewidths=0.5,
            ripf_labels=ripf
        )
        
        # cmd.movavg filtered time series:
        qc.plot_data(ax, xdata=year_ma,
            ydata=ydata_ma, ens_index=e, ens_label=ens_label,
            labels=labels, colors=colors, linewidths=1.5,
            ripf_labels=ripf
        )
        
        ax.set_title(f"{cmd.model} {cmd.experiment} "
            + ("" if cmd.x2 == "none" else f"+ {cmd.x2}") + ", "
            + f"{ens_label} o{ftype}temptend (various methods) "
            + "\ntime series at "
            + qc.lat_eval_label(lat_eval,
                getattr(nf,f"nc_ref_lat_{hemi}_attrs")["units"])
            )
    
        ax.set_xlabel("Year")
        ax.set_ylabel(units_label)
        ax.legend(fontsize=10)
    
    
    # Figs. 3-4: scatter plots of anomalies (piControl) or
    # deltas between time periods (otherwise) for each pair of
    # methods
    # ------------------------------------------------------- #
    descr_kw["plot_type"] = f"scatter_{lat_eval_n:.0f}n"
    if cmd.experiment != "piControl":
        descr_kw["ensemble_members_description"] = \
            f"{n_ens}_members"
    fig3, axs3 = qc.start_figure(ncols=3, **descr_kw)
    
    descr_kw["plot_type"] = f"scatter_{abs(lat_eval_s):.0f}s"
    fig4, axs4 = qc.start_figure(ncols=3, **descr_kw)
    
    for j in range(len(labels)):
        labels[j] += f" ({units_label})"
    
    for axs, datasets, lat_eval, hemi in zip(
            [axs3, axs4], [ydata_n, ydata_s],
            [lat_eval_n, lat_eval_s], ["n", "s"]
        ):
        
        if cmd.experiment == "piControl":
            datasets = qc.contiguous_timestep_averages(
                datasets, cmd.yravgs)
            
            # Need to select member axis for qc.scatter_data_3
            # (if piControl, only one anyway):
            datasets = [x[:,0] for x in datasets]
        
        else:
            datasets = qc.time_average_period_differences(
                year, datasets, cmd.yravg1, cmd.yravg2)
        
        qc.scatter_data_3(axs, datasets=datasets,
            ax_labels=labels,
            legend_label_type="correlation_slope",
            colors=["tab:grey"], ax_label_colors=colors)
        
        for ax in axs:
            ax.axline([0,0], slope=1.0, color="tab:red",
                label=r"1:1")
            leg = ax.legend(fontsize=10)
            leg.set_zorder(1000)
    
    fig3.suptitle(f"{cmd.model} {cmd.experiment}"
        + ("" if cmd.x2 == "none" else f" + {cmd.x2}") + ", "
        + f"o{ftype}temptend{units} (various methods)" + "\nat "
        + qc.lat_eval_label(lat_eval,
            getattr(nf,f"nc_ref_lat_n_attrs")["units"])
        + (f" ({cmd.yravgs} year averages, anomalies)"
           if cmd.experiment == "piControl" else
           f" ({n_ens} members, difference between "
           + f"{cmd.yravg2[0]}" + u"\u2013" + f"{cmd.yravg2[1]}"
           + f" and {cmd.yravg1[0]}" + u"\u2013"
           + f"{cmd.yravg1[1]} averages)"),
        fontsize=10.5)
    
    fig4.suptitle(f"{cmd.model} {cmd.experiment}"
        + ("" if cmd.x2 == "none" else f" + {cmd.x2}") + ", "
        + f"o{ftype}temptend{units} (various methods)" + "\nat "
        + qc.lat_eval_label(lat_eval,
            getattr(nf,f"nc_ref_lat_n_attrs")["units"])
        + (f" ({cmd.yravgs} year averages, anomalies)"
           if cmd.experiment == "piControl" else
           f" ({n_ens} members, difference between "
           + f"{cmd.yravg2[0]}" + u"\u2013" + f"{cmd.yravg2[1]}"
           + f" and {cmd.yravg1[0]}" + u"\u2013"
           + f"{cmd.yravg1[1]} averages)"),
        fontsize=10.5)
    
    
    qc.finish_figures([fig1, fig2, fig3, fig4],
        "\n".join(data_paths), savefig=cmd.savefigs,
        subplots_adjust_kw={"top":0.86, "bottom":0.23})
    



if __name__ == "__main__":
    main()
