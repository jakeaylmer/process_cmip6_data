from api import model_diagnostics as mdiags
from src import netcdf as nf
from src import qc
from src import script_tools



def main():
    
    cmd = script_tools.parse_qc_cmd_args()
    
    gd_kw = {
        "model_id": cmd.model,
        "experiment_1": cmd.experiment,
        "experiment_2": cmd.x2,
        "lat_eval": cmd.latitude
    }
    
    data_paths = ["", ""]
    
    # Load from aliases (defined in
    # api._diagnostic_definitions):
    hfds_n, hfds_s, year, r_vals, ripf, lat_eval_n, lat_eval_s,\
        data_paths[0] = qc.get_data("hfds", **gd_kw)
    hfds_wm2_n, hfds_wm2_s, year, r_vals, ripf, lat_eval_n, \
        lat_eval_s, data_paths[1] = \
            qc.get_data("hfds_wm-2", **gd_kw)
    
    e, ens_label = qc.ens_member_label(cmd.ensmemberid,
                                       r_vals, ripf)
    
    full_diag_names = [mdiags.get_diagnostic_from_alias(x)
                       for x in ["hfds", "hfds_wm-2"]]
    
    # ------------------------------------------------------- #
    
    exp_label = cmd.experiment
    exp_label += "" if cmd.x2 == "none" else f"_{cmd.x2}"
    
    descr_kw = {
        "model_id"                    : cmd.model,
        "experiment_id"               : exp_label,
        "ensemble_members_description": ens_label,
        "diagnostic_description"      : full_diag_names[0],
        "plot_type"                   : "time_series_"
                                        + f"{lat_eval_n:.0f}n_"
                                        + f"{-lat_eval_s:.0f}s"
    }
    
    fig1, ax1 = qc.start_figure(**descr_kw)
    
    descr_kw["diagnostic_description"] = full_diag_names[1]
    fig2, ax2 = qc.start_figure(**descr_kw)
    
    year_ma, hfds_ma_n = qc.moving_average(year, hfds_n,
                                           cmd.movavg)
    _, hfds_ma_s = qc.moving_average(year, hfds_s, cmd.movavg)
    _, hfds_wm2_ma_n = qc.moving_average(year, hfds_wm2_n, cmd.movavg)
    _, hfds_wm2_ma_s = qc.moving_average(year, hfds_wm2_s, cmd.movavg)
    
    for fig, ax, ydata, ydata_ma, j, ylabel in zip(
            [fig1, fig2], [ax1, ax2],
            [[hfds_n, hfds_s], [hfds_wm2_n, hfds_wm2_s]],
            [[hfds_ma_n, hfds_ma_s],
                [hfds_wm2_ma_n, hfds_wm2_ma_s]],
            [0, 1],
            [nf.field_units["heattransport"],
                nf.field_units["heatflux"]]):
        
        # Full time series data:
        qc.plot_data(ax, xdata=year, ydata=ydata, ens_index=e,
            ens_label=ens_label, labels=[None]*2,
            linewidths=0.5, ripf_labels=ripf)
        
        # Moving-average filtered data:
        qc.plot_data(ax, xdata=year_ma, ydata=ydata_ma,
            ens_index=e, ens_label=ens_label, linewidths=1.5,
            labels=["hfds_n", "hfds_s"], ripf_labels=ripf)
            
        ax.set_title(f"{cmd.model} {cmd.experiment}"
            + ("" if cmd.x2 == "none" else f" + {cmd.x2}")
            + ", " + ens_label + ", " + full_diag_names[j] +"\n"
            + "time series at " + qc.lat_eval_label(lat_eval_n)
            + " / " + qc.lat_eval_label(lat_eval_s))
        
        ax.set_xlabel("Year")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=8 if "all" in ens_label else 12)
    
        qc.finish_figures([fig], str(data_paths[j]),
            savefig=cmd.savefigs,
            subplots_adjust_kw={"top":0.86, "bottom":0.18})



if __name__ == "__main__":
    main()
