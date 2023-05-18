import numpy as np

from api import model_diagnostics as mdiags
from src import netcdf as nf
from src import qc
from src import script_tools



def main():
    
    # Reverse order (alphanumeric order gives "cc_approx" first)
    cmd = script_tools.parse_qc_cmd_args(
        mdiags.get_valid_aliases_from_keyword("aht")[::-1])
    
    gd_kw = {
        "diagnostic": cmd.diagnostic,
        "model_id": cmd.model,
        "experiment_1": cmd.experiment,
        "experiment_2": cmd.x2,
        "lat_eval": None
    }
    
    # Load from aliases (defined in api._diagnostic_definitions)
    aht_n, aht_s, year, r_vals, ripf, lat_n, lat_s, data_path \
        = qc.get_data(**gd_kw)
    
    e, ens_label = qc.ens_member_label(cmd.ensmemberid,
                                       r_vals, ripf)
    
    n_ens = len(ripf)
    
    full_diag_name = mdiags.get_diagnostic_from_alias(
        cmd.diagnostic)
    
    y1y2, aht_pr_n = qc.profile_mean(year, aht_n, cmd.yravg)
    _, aht_pr_s = qc.profile_mean(year, aht_s, cmd.yravg)
    
    aht_ts_n, lat_eval_n = \
        qc.evaluate_at_latitude(aht_n, lat_n, abs(cmd.latitude))
    aht_ts_s, lat_eval_s = \
        qc.evaluate_at_latitude(aht_s, lat_s,-abs(cmd.latitude))
    
    # ------------------------------------------------------- #
    
    exp_label = cmd.experiment
    exp_label += "" if cmd.x2 == "none" else f"_{cmd.x2}"
    
    descr_kw = {
        "model_id"                    : cmd.model,
        "experiment_id"               : exp_label,
        "ensemble_members_description": ens_label,
        "diagnostic_description"      : full_diag_name,
        "plot_type"                   : "time_series_"
            + f"{lat_eval_n:.0f}n_{abs(lat_eval_s):.0f}s"
    }
    
    fig1, ax1 = qc.start_figure(**descr_kw)
    
    year_ma, aht_ts_ma_n = qc.moving_average(year, aht_ts_n,
                                             cmd.movavg)
    aht_ts_ma_s = qc.moving_average(year, aht_ts_s,
                                    cmd.movavg)[1]
    
    # Full time series data:
    qc.plot_data(ax1, xdata=year, ydata=[aht_ts_n, -aht_ts_s],
        ens_index=e, ens_label=ens_label, labels=[None]*2,
        linewidths=0.5, ripf_labels=ripf)
    
    # Moving-average filtered data:
    qc.plot_data(ax1, xdata=year_ma,
        ydata=[aht_ts_ma_n, -aht_ts_ma_s], ens_index=e,
        ens_label=ens_label, labels=["aht_n", u"\u2212aht_s"],
        linewidths=1.5, ripf_labels=ripf)
    
    ax1.set_title(f"{cmd.model} {cmd.experiment}"
        + ("" if cmd.x2 == "none" else f" + {cmd.x2}") + ", "
        + ens_label + ", " + full_diag_name +"\ntime series at "
        + qc.lat_eval_label(lat_eval_n) + " / "
        + qc.lat_eval_label(lat_eval_s))
    
    ax1.set_xlabel("Year")
    ax1.set_ylabel(nf.field_units["heattransport"])
    ax1.legend(fontsize=8 if "all" in ens_label else 12)
    
    
    # Latitude profile
    descr_kw["plot_type"] = (f"latitude_profile_y{y1y2[0]}-"
                            + f"{y1y2[1]}")
    fig2, ax2 = qc.start_figure(**descr_kw)
    
    qc.plot_data(ax2, xdata=[lat_n, lat_s],
        ydata=[aht_pr_n.T, aht_pr_s.T],
        ens_index=e, ens_label=ens_label,
        labels=["aht_n", "aht_s"],
        ripf_labels=ripf)
    
    ax2.set_title(f"{cmd.model} {cmd.experiment}"
        + ("" if cmd.x2 == "none" else f" + {cmd.x2}") + ", "
        + ens_label + ", " + full_diag_name +"\nprofile "
        + f"averaged over years {y1y2[0]}" + u"\u2013"
        + f"{y1y2[1]}")
    
    ax2.set_xlabel(nf.nc_ref_lat_n_name + " ("
                   + nf.nc_ref_lat_n_attrs["units"] + ")")
    ax2.set_ylabel(nf.field_units["heattransport"])
    ax2.legend(fontsize=8 if "all" in ens_label else 12)
    
    qc.finish_figures([fig1, fig2], str(data_path),
        savefig=cmd.savefigs,
        subplots_adjust_kw={"top":0.86, "bottom":0.18})



if __name__ == "__main__":
    main()
