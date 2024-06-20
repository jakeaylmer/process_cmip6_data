"""Quality control plot: atmospheric vertical heat fluxes. Plots
time series of each term for a specified reference latitude.
"""

from process_cmip6_data.api import model_diagnostics as mdiags
from process_cmip6_data.src import netcdf as nf
from process_cmip6_data.src import qc
from process_cmip6_data.src import script_tools


def main():
    
    diagnostic_choices = ["gn_interp", "gn", "cc"]
    
    cmd = script_tools.parse_qc_cmd_args(diagnostic_choices)
    
    f_star = ["f_down", "f_olr", "f_sw_surf", "f_sw_toa","f_up"]
    
    ydata_n = []
    ydata_s = []
    
    gd_kw = {
        "model_id": cmd.model,
        "experiment_1": cmd.experiment,
        "experiment_2": cmd.x2,
        "lat_eval": cmd.latitude}
    
    for fn in f_star:
        
        # Load from aliases (defined in
        # api._diagnostic_definitions):
        f_star_n, f_star_s, year, r_vals, ripf, lat_eval_n, \
            lat_eval_s, data_path = qc.get_data(
                f"{fn}_{cmd.diagnostic}", **gd_kw)
        
        ydata_n.append(f_star_n.copy())
        ydata_s.append(f_star_s.copy())
    
    e, ens_label = qc.ens_member_label(cmd.ensmemberid,
                                       r_vals, ripf)
    
    full_diag_name = mdiags.get_diagnostic_from_alias(
        f"{f_star[-1]}_{cmd.diagnostic}").replace(f_star[-1],
                                                  "f_*")
    
    # ------------------------------------------------------- #
    
    exp_label = cmd.experiment
    exp_label += "" if cmd.x2 == "none" else f"_{cmd.x2}"
    
    descr_kw = {
        "model_id"                    : cmd.model,
        "experiment_id"               : exp_label,
        "ensemble_members_description": ens_label,
        "diagnostic_description"      : full_diag_name.replace(
                                            "*", "star"),
        "plot_type"                   : "time_series_"
                                        + f"{lat_eval_n:.0f}n"
    }
    
    fig1, ax1 = qc.start_figure(**descr_kw)
    descr_kw["plot_type"] =f"time_series_{abs(lat_eval_s):.0f}s"
    fig2, ax2 = qc.start_figure(**descr_kw)
    
    year_ma = qc.moving_average(year, ydata_n[0], cmd.movavg)[0]
    ydata_ma_n = [qc.moving_average(year, y, cmd.movavg)[1]
                  for y in ydata_n]
    ydata_ma_s = [qc.moving_average(year, y, cmd.movavg)[1]
                  for y in ydata_s]
    
    colors = ["tab:pink", "tab:blue", "tab:orange", "tab:olive",
              "tab:green"]
    
    for ax, ydata, ydata_ma, lat_eval in zip(
            [ax1, ax2], [ydata_n, ydata_s],
            [ydata_ma_n, ydata_ma_s],
            [lat_eval_n, lat_eval_s]):
        
        # Full time series data:
        qc.plot_data(ax, xdata=year, ydata=ydata,
            ens_index=e, ens_label=ens_label, colors=colors,
            labels=[None]*len(ydata), linewidths=0.5,
            ripf_labels=ripf)
        
        # Moving-average filtered data:
        qc.plot_data(ax, xdata=year_ma, ydata=ydata_ma,
            ens_index=e, ens_label=ens_label, colors=colors,
            labels=f_star, linewidths=1.5, ripf_labels=ripf)
        
        ax.set_title(f"{cmd.model} {cmd.experiment}"
            + ("" if cmd.x2 == "none" else f" + {cmd.x2}")
            + ", " + ens_label + ", " + full_diag_name + "\n"
            + "time series at " + qc.lat_eval_label(lat_eval))
        
        ax.set_xlabel("Year")
        ax.set_ylabel(nf.field_units["heatflux"])
        ax.legend(fontsize=8 if "all" in ens_label else 12)
    
    qc.finish_figures([fig1, fig2],
        str(data_path).replace(f_star[-1], "f_*"),
        savefig=cmd.savefigs,
        subplots_adjust_kw={"top":0.86, "bottom":0.18})


if __name__ == "__main__":
    main()
