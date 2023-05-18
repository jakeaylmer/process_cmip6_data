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
        "lat_eval": cmd.latitude,
        "select_diag": cmd.seldiag
    }
    
    data_paths = ["", "", ""]
    
    # Load from keywords (defined in
    # api._diagnostic_definitions):
    dhdt_n, dhdt_s, year, r_vals, ripf, lat_n_dhdt, lat_s_dhdt,\
        data_paths[0] = qc.get_data("dhdt", **gd_kw)
    
    hfds_n, hfds_s, year, r_vals, ripf, lat_n_hfds, lat_s_hfds,\
        data_paths[1] = qc.get_data("hfds", **gd_kw)
    
    oht_n, oht_s, year, r_vals, ripf, lat_n_oht, lat_s_oht,\
        data_paths[2] = qc.get_data("oht", **gd_kw)
    
    budget_n = -dhdt_n + hfds_n + oht_n
    budget_s = -dhdt_s + hfds_s - oht_s
    
    e, ens_label = qc.ens_member_label(cmd.ensmemberid,
                                       r_vals, ripf)
    
    data_paths = qc.shorten_data_path_labels(data_paths)
    
    full_diag_names = [
        mdiags.get_diagnostic_from_keyword(x, cmd.model,
                                           cmd.seldiag)
                       for x in ["dhdt", "hfds", "oht"]]
    
    # ------------------------------------------------------- #
    
    exp_label = cmd.experiment
    exp_label += "" if cmd.x2 == "none" else f"_{cmd.x2}"
    
    descr_kw = {
        "model_id"                    : cmd.model,
        "experiment_id"               : exp_label,
        "ensemble_members_description": ens_label,
        "diagnostic_description"      : "budget",
        "plot_type"                   : "time_series_"
                                        + f"{lat_n_oht:.0f}n"
    }
    
    fig1, ax1 = qc.start_figure(**descr_kw)
    descr_kw["plot_type"] = f"time_series_{-lat_s_oht:.0f}s"
    fig2, ax2 = qc.start_figure(**descr_kw)
    
    year_ma, dhdt_ma_n = qc.moving_average(year, dhdt_n,
                                           cmd.movavg)
    _, dhdt_ma_s = qc.moving_average(year, dhdt_s, cmd.movavg)
    _, hfds_ma_n = qc.moving_average(year, hfds_n, cmd.movavg)
    _, hfds_ma_s = qc.moving_average(year, hfds_s, cmd.movavg)
    _, oht_ma_n = qc.moving_average(year, oht_n, cmd.movavg)
    _, oht_ma_s = qc.moving_average(year, oht_s, cmd.movavg)
    _, budget_ma_n = qc.moving_average(year, budget_n, cmd.movavg)
    _, budget_ma_s = qc.moving_average(year, budget_s, cmd.movavg)
    
    ydata_n = [dhdt_n, hfds_n, oht_n]
    ydata_s = [dhdt_s, hfds_s, oht_s]
    ydata_ma_n = [dhdt_ma_n, hfds_ma_n, oht_ma_n]
    ydata_ma_s = [dhdt_ma_s, hfds_ma_s, oht_ma_s]
    
    labels = ["dhdt", "hfds", "oht"]
    colors = ["tab:orange", "tab:green", "tab:blue"]
    
    for ax, ydata, ydata_ma, lat_eval in zip(
            [ax1, ax2],
            [ydata_n, ydata_s],
            [ydata_ma_n, ydata_ma_s],
            [lat_n_oht, lat_s_oht]):
        
        # Full time series data:
        qc.plot_data(ax, xdata=year, ydata=ydata, ens_index=e,
            ens_label=ens_label, labels=[None]*3, colors=colors,
            linewidths=0.5, ripf_labels=ripf)
        
        # Moving-average filtered data:
        qc.plot_data(ax, xdata=year_ma, ydata=ydata_ma,
            ens_index=e, ens_label=ens_label, linewidths=1.5,
            colors=colors, labels=labels, ripf_labels=ripf)
            
        ax.set_title(f"{cmd.model} {cmd.experiment}"
            + ("" if cmd.x2 == "none" else f" + {cmd.x2}")
            + ", " + ens_label + ", ocean heat budget terms\n"
            + "time series at " + qc.lat_eval_label(lat_eval))
        
        ax.set_xlabel("Year")
        ax.set_ylabel(nf.field_units["heattransport"])
        ax.legend(fontsize=8 if "all" in ens_label else 12)
    
    descr_kw["plot_type"] = \
        f"closure_error_time_series_{lat_n_oht:.0f}n"
    fig3, ax3 = qc.start_figure(**descr_kw)
    
    qc.plot_data(ax3, xdata=year, ydata=[budget_n, budget_s],
        ens_index=e, ens_label=ens_label, labels=[None]*2,
            linewidths=0.5, ripf_labels=ripf)
    
    qc.plot_data(ax3, xdata=year_ma,
        ydata=[budget_ma_n, budget_ma_s], ens_index=e,
        ens_label=ens_label, labels=["n", "s"], linewidths=1.5,
        ripf_labels=ripf)
    
    ax3.set_xlabel("Year")
    ax3.set_ylabel(nf.field_units["heattransport"])
    ax3.legend(fontsize=12)
    
    ax3.set_title(f"{cmd.model} {cmd.experiment}"
        + ("" if cmd.x2 == "none" else f" + {cmd.x2}")
        + ", " + ens_label + ", ocean heat budget\n"
        + "closure error, time series at "
        + qc.lat_eval_label(lat_n_oht) + " / "
        + qc.lat_eval_label(lat_s_oht))
    
    qc.finish_figures([fig1, fig2, fig3],
        "\n".join([str(x) for x in data_paths]),
        savefig=cmd.savefigs,
        subplots_adjust_kw={"top":0.86, "bottom":0.22})



if __name__ == "__main__":
    main()
