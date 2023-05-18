import numpy as np

from api import model_diagnostics as mdiags
from src import netcdf as nf
from src import qc
from src import script_tools



def main():
    
    cmd = script_tools.parse_qc_cmd_args()
    
    gd_kw = {"model_id": cmd.model,
        "experiment_1": cmd.experiment, "experiment_2": cmd.x2}
    
    data_paths = [""]*4
    
    # Load from aliases (defined in api._diagnostic_definitions)
    sia_n, sia_s, year, r_vals, ripf, _, _, data_paths[0] = \
        qc.get_data("sia", **gd_kw)
    sia_mon_n, sia_mon_s, _, r_vals, ripf, _,_,data_paths[1] = \
        qc.get_data("sia_mon", **gd_kw)
    sie_n, sie_s, year, r_vals, ripf, _, _, data_paths[2] = \
        qc.get_data("sie", **gd_kw)
    sie_mon_n, sie_mon_s, _, r_vals, ripf, _,_,data_paths[3] = \
        qc.get_data("sie_mon", **gd_kw)
    
    e, ens_label = qc.ens_member_label(cmd.ensmemberid,
                                       r_vals, ripf)
    
    full_diag_names = [mdiags.get_diagnostic_from_alias(x)
                       for x in ["sia", "sia_mon", "sie",
                                 "sie_mon"]]
    
    # Compute climatological seasonal cycles:
    y1y2, sia_clim_n = qc.seasonal_climatology(year, sia_mon_n,
                                               cmd.yravg)
    sia_clim_s = qc.seasonal_climatology(year, sia_mon_s,
                                         cmd.yravg)[1]
    sie_clim_n = qc.seasonal_climatology(year, sie_mon_n,
                                         cmd.yravg)[1]
    sie_clim_s = qc.seasonal_climatology(year, sie_mon_s,
                                         cmd.yravg)[1]
    
    # ------------------------------------------------------- #
    
    exp_label = cmd.experiment
    exp_label += "" if cmd.x2 == "none" else f"_{cmd.x2}"
    
    descr_kw = {
        "model_id"                    : cmd.model,
        "experiment_id"               : exp_label,
        "ensemble_members_description": ens_label,
        "diagnostic_description"      : "sia_sie_yr",
        "plot_type"                   : "time_series"}
    
    fig1, ax1 = qc.start_figure(**descr_kw)
    
    year_ma, sia_ma_n = qc.moving_average(year,sia_n,cmd.movavg)
    sia_ma_s = qc.moving_average(year, sia_s, cmd.movavg)[1]
    sie_ma_n = qc.moving_average(year, sie_n, cmd.movavg)[1]
    sie_ma_s = qc.moving_average(year, sie_s, cmd.movavg)[1]
    
    # Full time series data:
    qc.plot_data(ax1, xdata=year,
        ydata=[sia_n, sia_s, sie_n, sie_s],
        ens_index=e, ens_label=ens_label, labels=[None]*4,
        linewidths=0.5, linestyles=["-", "-", ":", ":"],
        ripf_labels=ripf)
    
    # Moving-average filtered data:
    qc.plot_data(ax1, xdata=year_ma,
        ydata=[sia_ma_n, sia_ma_s, sie_ma_n, sie_ma_s],
        ens_index=e, ens_label=ens_label,
        labels=["sia_n", "sia_s", "sie_n", "sie_s"],
        linewidths=1.5, linestyles=["-", "-", ":", ":"],
        ripf_labels=ripf)
    
    ax1.set_title(f"{cmd.model} {cmd.experiment}"
        + ("" if cmd.x2 == "none" else f" + {cmd.x2}") + ", "
        + ens_label + " " + full_diag_names[0] + " / "
        + full_diag_names[2] + " time series")
    
    ax1.set_xlabel("Year")
    ax1.set_ylabel(nf.field_units["seaicearea"])
    ax1.legend(fontsize=8 if "all" in ens_label else 12)
    
    descr_kw["diagnostic_description"] = "sia_sie_mon"
    descr_kw["plot_type"] =f"seasonal_clim_y{y1y2[0]}-{y1y2[1]}"
    fig2, ax2 = qc.start_figure(**descr_kw)
    
    qc.plot_data(ax2, xdata=np.arange(1,13,1), markers="o",
        ydata=[sia_clim_n, sia_clim_s, sie_clim_n, sie_clim_s],
        labels=["sia_n", "sia_s", "sie_n", "sie_s"],
        linestyles=["-", "-", ":", ":"], ripf_labels=ripf)
    
    ax2.set_title(f"{cmd.model} {cmd.experiment}"
        + ("" if cmd.x2 == "none" else f" + {cmd.x2}") + ", "
        + ens_label + " " + full_diag_names[1] + " / "
        + full_diag_names[3] + "\nmonthly climatology, mean "
        + f"over years {y1y2[0]}" + u"\u2013" + f"{y1y2[1]}")
    
    qc.monthly_xaxis(ax2)
    ax2.set_ylabel(nf.field_units["seaicearea"])
    ax2.legend(fontsize=8 if "all" in ens_label else 12)
    
    qc.finish_figures([fig1],
        str(data_paths[0]).replace("sia", "si?"),
        savefig=cmd.savefigs,
        subplots_adjust_kw={"top":0.9, "bottom":0.18})
    
    qc.finish_figures([fig2],
        str(data_paths[1]).replace("sia", "si?"),
        savefig=cmd.savefigs,
        subplots_adjust_kw={"top":0.86, "bottom":0.18})



if __name__ == "__main__":
    main()
