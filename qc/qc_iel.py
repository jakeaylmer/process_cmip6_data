import numpy as np

from api import model_diagnostics as mdiags
from src import netcdf as nf
from src import qc
from src import script_tools



def main():
    
    cmd = script_tools.parse_qc_cmd_args()
    
    gd_kw = {"model_id": cmd.model, "select_diag": cmd.seldiag,
        "experiment_1": cmd.experiment, "experiment_2": cmd.x2}
    
    data_paths = [""]*2
    
    # Load from keywords (defined in
    # api._diagnostic_definitions)
    iel_n, iel_s, year, r_vals, ripf, _, _, data_paths[0] = \
        qc.get_data("iel", **gd_kw)
    iel_mon_n, iel_mon_s, _, r_vals, ripf, _,_,data_paths[1] = \
        qc.get_data("iel_mon", **gd_kw)
    
    e, ens_label = qc.ens_member_label(cmd.ensmemberid,
                                       r_vals, ripf)
    
    full_diag_names = [mdiags.get_diagnostic_from_keyword(x,
                           model_id=cmd.model,
                           alias_select=cmd.seldiag)
                       for x in ["iel", "iel_mon"]]
    
    # Compute climatological seasonal cycles:
    y1y2, iel_clim_n = qc.seasonal_climatology(year, iel_mon_n,
                                               cmd.yravg)
    iel_clim_s = qc.seasonal_climatology(year, iel_mon_s,
                                         cmd.yravg)[1]
    
    # ------------------------------------------------------- #
    
    exp_label = cmd.experiment
    exp_label += "" if cmd.x2 == "none" else f"_{cmd.x2}"
    
    descr_kw = {
        "model_id"                    : cmd.model,
        "experiment_id"               : exp_label,
        "ensemble_members_description": ens_label,
        "diagnostic_description"      : full_diag_names[0],
        "plot_type"                   : "time_series"}
    
    fig1, ax1 = qc.start_figure(**descr_kw)
    
    year_ma, iel_ma_n = qc.moving_average(year,iel_n,cmd.movavg)
    iel_ma_s = qc.moving_average(year, iel_s, cmd.movavg)[1]
    
    # Full time series data:
    qc.plot_data(ax1, xdata=year, ydata=[iel_n, -iel_s],
        ens_index=e, ens_label=ens_label, labels=[None]*2,
        linewidths=0.5, ripf_labels=ripf)
    
    # Moving-average filtered data:
    qc.plot_data(ax1, xdata=year_ma,
        ydata=[iel_ma_n, -iel_ma_s], ens_index=e,
        ens_label=ens_label, labels=["iel_n", u"\u2212iel_s"],
        linewidths=1.5, ripf_labels=ripf)
    
    ax1.set_title(f"{cmd.model} {cmd.experiment}"
        + ("" if cmd.x2 == "none" else f" + {cmd.x2}") + ", "
        + ens_label + " " + full_diag_names[0] + " time series")
    
    ax1.set_xlabel("Year")
    ax1.set_ylabel(nf.field_units["seaiceedge"])
    ax1.legend(fontsize=8 if "all" in ens_label else 12)
    
    descr_kw["diagnostic_description"] = full_diag_names[1]
    descr_kw["plot_type"] =f"seasonal_clim_y{y1y2[0]}-{y1y2[1]}"
    fig2, ax2 = qc.start_figure(**descr_kw)
    
    qc.plot_data(ax2, xdata=np.arange(1,13,1),
        ydata=[iel_clim_n, -iel_clim_s], markers="o",
        labels=["iel_n", u"\u2212iel_s"], ripf_labels=ripf)
    
    ax2.set_title(f"{cmd.model} {cmd.experiment}"
        + ("" if cmd.x2 == "none" else f" + {cmd.x2}") + ", "
        + ens_label + " " + full_diag_names[1]
        + "\nmonthly climatology, mean "
        + f"over years {y1y2[0]}" + u"\u2013" + f"{y1y2[1]}")
    
    qc.monthly_xaxis(ax2)
    ax2.set_ylabel(nf.field_units["seaiceedge"])
    ax2.legend(fontsize=8 if "all" in ens_label else 12)
    
    qc.finish_figures([fig1], data_paths[0],
        savefig=cmd.savefigs,
        subplots_adjust_kw={"top":0.9, "bottom":0.18})
    
    qc.finish_figures([fig2], data_paths[1],
        savefig=cmd.savefigs,
        subplots_adjust_kw={"top":0.86, "bottom":0.18})



if __name__ == "__main__":
    main()
