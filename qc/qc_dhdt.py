from process_cmip6_data.api import model_diagnostics as mdiags
from process_cmip6_data.src import netcdf as nf
from process_cmip6_data.src import qc
from process_cmip6_data.src import script_tools


def main():
    
    cmd = script_tools.parse_qc_cmd_args()
    
    gd_kw = {
        "model_id": cmd.model,
        "experiment_1": cmd.experiment,
        "experiment_2": cmd.x2,
        "lat_eval": cmd.latitude,
        "select_diag": cmd.seldiag
    }
    
    data_paths = ["", ""]
    
    # Load from keyword (def. in api._diagnostic_definitions):
    dhdt_n, dhdt_s, year, r_vals, ripf, lat_eval_n, lat_eval_s,\
        data_paths[0] = qc.get_data("dhdt", **gd_kw)
    dhdt_wm2_n, dhdt_wm2_s, _, _, _, _, _, data_paths[1] = \
            qc.get_data("dhdt_wm-2", **gd_kw)
    
    data_paths = qc.shorten_data_path_labels(data_paths)
    
    e, ens_label = qc.ens_member_label(cmd.ensmemberid,
                                       r_vals, ripf)
    
    full_diag_names = [mdiags.get_diagnostic_from_keyword(x,
                           cmd.model, cmd.seldiag)
                       for x in ["dhdt", "dhdt_wm-2"]]
    
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
    
    year_ma, dhdt_ma_n = qc.moving_average(year, dhdt_n,
                                           cmd.movavg)
    _, dhdt_ma_s = qc.moving_average(year, dhdt_s, cmd.movavg)
    
    fig_list = [fig1]
    ax_list = [ax1]
    ydata_list = [[dhdt_n, dhdt_s]]
    ydata_ma_list = [[dhdt_ma_n, dhdt_ma_s]]
    j_list = [0]
    ylabel_list = [nf.field_units["heattransport"]]
    
    if len(data_paths) == 3:
        # has W m-2 data
        descr_kw["diagnostic_description"] = full_diag_names[1]
        fig2, ax2 = qc.start_figure(**descr_kw)
        
        fig_list.append(fig2)
        ax_list.append(ax2)
        j_list.append(1)
        ylabel_list.append(nf.field_units["heatflux"])
        ydata_list.append([dhdt_wm2_n, dhdt_wm2_s])
        
        _, dhdt_wm2_ma_n = qc.moving_average(year, dhdt_wm2_n, cmd.movavg)
        _, dhdt_wm2_ma_s = qc.moving_average(year, dhdt_wm2_s, cmd.movavg)
        
        ydata_ma_list.append([dhdt_wm2_ma_n, dhdt_wm2_ma_s])
    
    for fig, ax, ydata, ydata_ma, j, ylabel in zip(
            fig_list, ax_list, ydata_list, ydata_ma_list,
            j_list, ylabel_list):
        
        # Full time series data:
        qc.plot_data(ax, xdata=year, ydata=ydata, ens_index=e,
            ens_label=ens_label, labels=[None]*2,
            linewidths=0.5, ripf_labels=ripf)
        
        # Moving-average filtered data:
        qc.plot_data(ax, xdata=year_ma, ydata=ydata_ma,
            ens_index=e, ens_label=ens_label, linewidths=1.5,
            labels=["dhdt_n", "dhdt_s"], ripf_labels=ripf)
            
        ax.set_title(f"{cmd.model} {cmd.experiment}"
            + ("" if cmd.x2 == "none" else f" + {cmd.x2}")
            + ", " + ens_label + ", " + full_diag_names[j] +"\n"
            + "time series at " + qc.lat_eval_label(lat_eval_n)
            + " / " + qc.lat_eval_label(lat_eval_s))
        
        ax.set_xlabel("Year")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=8 if "all" in ens_label else 12)
    
        qc.finish_figures([fig],
            "\n".join([str(data_paths[k]) for k in [0,j+1]]),
            savefig=cmd.savefigs,
            subplots_adjust_kw={"top":0.86, "bottom":0.18})



if __name__ == "__main__":
    main()
