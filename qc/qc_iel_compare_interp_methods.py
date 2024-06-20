"""Quality control plot: sea ice-edge latitude. Plots zonal mean
annual mean time series and climatological seasonal cyclone for
each hemisphere, comparing multiple interpolation methods and/or
resolutions.
"""

import numpy as np

from process_cmip6_data.api import model_diagnostics as mdiags
from process_cmip6_data.src import netcdf as nf
from process_cmip6_data.src import qc
from process_cmip6_data.src import script_tools


def main():
    
    prsr = script_tools.qc_argument_parser()
    prsr.add_argument("--remapmethod", type=str, default="bil",
        choices=["bil", "con", "con2", "cony", "cony2", "dis",
                 "nn"],
        help="Re-gridding method (\"bil\" for bilinear, etc.)")
    cmd = prsr.parse_args()
    
    # ======================================================== #
    # NOTE: currently remap method is ignored (no separate
    # ----- keywords/aliases defined for anything except for
    #       bilinear ("bil")
    # ======================================================== #
    
    # Can't reasonably plot all members for multiple methods
    # on one plot:
    if cmd.ensmemberid < -1:
        cmd.ensmemberid = 0
        print("WARNING: ignoring --ensmemberid (plotting first "
              + "member only)")
    
    gd_kw = {"model_id": cmd.model,
        "experiment_1": cmd.experiment, "experiment_2": cmd.x2}
    
    methods = [f"{x}_{cmd.remapmethod}" for x in 
        ["4deg", "2deg", "1deg", "05deg", "025deg"]]
    res_labels = [u"4\u00b0", u"2\u00b0", u"1\u00b0",
                 u"\u00bd\u00b0", u"\u00bc\u00b0"]
    data_paths = []
    data_paths_mon = []
    
    iel_n = []
    iel_s = []
    iel_mon_n = []
    iel_mon_s = []
    
    for method in methods:
        # Load from aliases (defined in
        # api._diagnostic_definitions)
        xn, xs, year, r_vals, ripf, _, _, data_path = \
            qc.get_data(f"iel_{method}", **gd_kw)
        xmn, xms, _, _, _, _, _, data_path_mon = \
            qc.get_data(f"iel_{method}_mon", **gd_kw)
        
        iel_n.append(xn)
        iel_s.append(-xs)
        iel_mon_n.append(xmn)
        iel_mon_s.append(-xms)
        
        data_paths.append(data_path)
        data_paths_mon.append(data_path_mon)
        
    e, ens_label = qc.ens_member_label(cmd.ensmemberid,
                                       r_vals, ripf)
    
    full_diag_names = [
        mdiags.get_diagnostic_from_alias(x).replace(
            methods[0], f"*deg_{cmd.remapmethod}")
        for x in [f"iel_{methods[0]}", f"iel_{methods[0]}_mon"]]
    
    # Compute climatological seasonal cycles:
    y1y2 = qc.seasonal_climatology(year, iel_mon_n[0],
                                   cmd.yravg)[0]
    
    iel_clim_n = [qc.seasonal_climatology(year, y, cmd.yravg)[1]
                  for y in iel_mon_n]
    iel_clim_s = [qc.seasonal_climatology(year, y, cmd.yravg)[1]
                  for y in iel_mon_s]
    
    # ------------------------------------------------------- #
    
    exp_label = cmd.experiment
    exp_label += "" if cmd.x2 == "none" else f"_{cmd.x2}"
    
    descr_kw = {
        "model_id"                    : cmd.model,
        "experiment_id"               : exp_label,
        "ensemble_members_description": ens_label,
        "diagnostic_description"      : \
            full_diag_names[0].replace("*", "star"),
        "plot_type"                   : "time_series"}
    
    fig1, ax1 = qc.start_figure(**descr_kw)
    
    year_ma = qc.moving_average(year, iel_n[0], cmd.movavg)[0]
    iel_ma_n = [qc.moving_average(year, y, cmd.movavg)[1]
                for y in iel_n]
    iel_ma_s = [qc.moving_average(year, y, cmd.movavg)[1]
                for y in iel_s]
    
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red",
              "tab:purple"]*2
    linestyles = ["-"]*len(methods) + ["--"]*len(methods)
    
    # Full time series data:
    qc.plot_data(ax1, xdata=year, ydata=iel_n + iel_s,
        ens_index=e, ens_label=ens_label,
        labels=res_labels + ["\"\" " + u"\u2212s"]
               + [None]*(len(methods)-1), colors=colors,
        linestyles=linestyles, linewidths=1.5, ripf_labels=ripf)
    
    # Moving-average filtered data:
    qc.plot_data(ax1, xdata=year_ma, ydata=iel_ma_n + iel_ma_s,
        ens_index=e, ens_label=ens_label,
        labels=[None]*2*len(methods), colors=colors,
        linestyles=linestyles, linewidths=0.5, ripf_labels=ripf)
    
    ax1.set_title(f"{cmd.model} {cmd.experiment}"
        + ("" if cmd.x2 == "none" else f" + {cmd.x2}") + ", "
        + ens_label + " " + full_diag_names[0] + " time series")
    
    ax1.set_xlabel("Year")
    ax1.set_ylabel(nf.field_units["seaiceedge"])
    ax1.legend(fontsize=10)
    
    descr_kw["diagnostic_description"] = \
        full_diag_names[1].replace("*", "star")
    descr_kw["plot_type"] =f"seasonal_clim_y{y1y2[0]}-{y1y2[1]}"
    fig2, ax2 = qc.start_figure(**descr_kw)
    
    qc.plot_data(ax2, xdata=np.arange(1,13,1),
        ydata=iel_clim_n +iel_clim_s, markers="o",
        labels=res_labels + ["\"\" " + u"\u2212s"]
               + [None]*(len(methods)-1), colors=colors,
        linestyles=linestyles, ripf_labels=ripf)
    
    ax2.set_title(f"{cmd.model} {cmd.experiment}"
        + ("" if cmd.x2 == "none" else f" + {cmd.x2}") + ", "
        + ens_label + " " + full_diag_names[1]
        + "\nmonthly climatology, mean "
        + f"over years {y1y2[0]}" + u"\u2013" + f"{y1y2[1]}")
    
    qc.monthly_xaxis(ax2)
    ax2.set_ylabel(nf.field_units["seaiceedge"])
    ax2.legend(fontsize=10)
    
    qc.finish_figures([fig1], data_paths[0].replace(
            methods[0], f"*deg_{cmd.remapmethod}"),
        savefig=cmd.savefigs,
        subplots_adjust_kw={"top":0.9, "bottom":0.18})
    
    qc.finish_figures([fig2], data_paths[0].replace(
            methods[0], f"*deg_{cmd.remapmethod}"),
        savefig=cmd.savefigs,
        subplots_adjust_kw={"top":0.86, "bottom":0.18})


if __name__ == "__main__":
    main()
