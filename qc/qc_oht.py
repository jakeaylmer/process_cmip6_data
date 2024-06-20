"""Quality control plot: ocean heat transport (OHT). Plots
a time series at a specified latitude (north and south) and a
heat transport profile (i.e., specified time average as a
function of latitude).
"""

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
        "lat_eval": None,
        "select_diag": cmd.seldiag}
    
    # Load from keyword (defined in api._diagnostic_definitions)
    oht_n, oht_s, year, r_vals, ripf, lat_n, lat_s, data_path \
        = qc.get_data("oht", **gd_kw)
    
    e, ens_label = qc.ens_member_label(cmd.ensmemberid,
                                       r_vals, ripf)
    
    full_diag_name = mdiags.get_diagnostic_from_keyword(
        "oht", cmd.model, cmd.seldiag)
    
    data_paths = qc.shorten_data_path_labels(data_path)
    
    y1y2, oht_pr_n = qc.profile_mean(year, oht_n, cmd.yravg)
    _, oht_pr_s = qc.profile_mean(year, oht_s, cmd.yravg)
    
    oht_ts_n, lat_eval_n = \
        qc.evaluate_at_latitude(oht_n, lat_n, abs(cmd.latitude))
    oht_ts_s, lat_eval_s = \
        qc.evaluate_at_latitude(oht_s, lat_s,-abs(cmd.latitude))
    
    # ------------------------------------------------------- #
    
    exp_label = cmd.experiment
    exp_label += "" if cmd.x2 == "none" else f"_{cmd.x2}"
    
    descr_kw = {
        "model_id"                    : cmd.model,
        "experiment_id"               : exp_label,
        "ensemble_members_description": ens_label,
        "diagnostic_description"      : full_diag_name,
        "plot_type"                   : "time_series_"
            + f"{lat_eval_n:.0f}n_{abs(lat_eval_s):.0f}s"}
    
    fig1, ax1 = qc.start_figure(**descr_kw)
    
    year_ma, oht_ts_ma_n = qc.moving_average(year, oht_ts_n,
                                             cmd.movavg)
    oht_ts_ma_s = qc.moving_average(year, oht_ts_s,
                                    cmd.movavg)[1]
    
    # Full time series data:
    qc.plot_data(ax1, xdata=year, ydata=[oht_ts_n, -oht_ts_s],
        ens_index=e, ens_label=ens_label, labels=[None]*2,
        linewidths=0.5, ripf_labels=ripf)
    
    # Moving-average filtered data:
    qc.plot_data(ax1, xdata=year_ma,
        ydata=[oht_ts_ma_n, -oht_ts_ma_s], ens_index=e,
        ens_label=ens_label, labels=["oht_n", u"\u2212oht_s"],
        linewidths=1.5, ripf_labels=ripf)
    
    ax1.set_title(f"{cmd.model} {cmd.experiment}"
        + ("" if cmd.x2 == "none" else f" + {cmd.x2}") + ", "
        + ens_label + ", " + full_diag_name +"\ntime series at "
        + qc.lat_eval_label(lat_eval_n) + " / "
        + qc.lat_eval_label(lat_eval_s), fontsize=11)
    
    ax1.set_xlabel("Year")
    ax1.set_ylabel(nf.field_units["heattransport"])
    ax1.legend(fontsize=8 if "all" in ens_label else 12)
    
    
    # Latitude profile
    descr_kw["plot_type"] = (f"latitude_profile_y{y1y2[0]}-"
                            + f"{y1y2[1]}")
    fig2, ax2 = qc.start_figure(**descr_kw)
    
    xdata = [lat_n]
    ydata = [oht_pr_n.T]
    labels = ["oht"]
    colors = ["k"]
    if "hfbasin" not in full_diag_name:
        xdata += [lat_s]
        ydata += [oht_pr_s.T]
        labels = ["oht_n", "oht_s"]
        colors = ["tab:red", "tab:blue"]
    
    qc.plot_data(ax2, xdata=xdata, ydata=ydata,
        ens_index=e, ens_label=ens_label, labels=labels,
        colors=colors, ripf_labels=ripf)
    
    ax2.set_title(f"{cmd.model} {cmd.experiment}"
        + ("" if cmd.x2 == "none" else f" + {cmd.x2}") + ", "
        + ens_label + ", " + full_diag_name +"\nprofile "
        + f"averaged over years {y1y2[0]}" + u"\u2013"
        + f"{y1y2[1]}", fontsize=11)
    
    qc.global_latitude_xaxis(ax2)
    ax2.set_ylabel(nf.field_units["heattransport"])
    ax2.legend(fontsize=8 if "all" in ens_label else 12)
    
    qc.finish_figures([fig1, fig2],
        "\n".join([str(x) for x in data_paths]),
        savefig=cmd.savefigs,
        subplots_adjust_kw={"top":0.86, "bottom":0.18})


if __name__ == "__main__":
    main()
