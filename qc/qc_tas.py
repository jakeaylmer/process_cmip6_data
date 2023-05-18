import numpy as np

from api import model_diagnostics as mdiags
from src import netcdf as nf
from src import qc
from src import script_tools



def main():
    
    # Reverse order (alphanumeric order gives "cc_approx" first)
    cmd = script_tools.parse_qc_cmd_args(
        mdiags.get_valid_aliases_from_keyword("tas")[::-1])
    
    gd_kw = {
        "diagnostic": cmd.diagnostic,
        "model_id": cmd.model,
        "experiment_1": cmd.experiment,
        "experiment_2": cmd.x2,
        "lat_eval": cmd.latitude
    }
    
    # Load from aliases (defined in api._diagnostic_definitions)
    tas_n, tas_s, year, r_vals, ripf, lat_eval_n, lat_eval_s, \
        data_path = qc.get_data(**gd_kw)
    
    e, ens_label = qc.ens_member_label(cmd.ensmemberid,
                                       r_vals, ripf)
    
    n_ens = len(ripf)
    
    full_diag_name = mdiags.get_diagnostic_from_alias(
        cmd.diagnostic)
    
    # ------------------------------------------------------- #
    
    exp_label = cmd.experiment
    exp_label += "" if cmd.x2 == "none" else f"_{cmd.x2}"
    
    descr_kw = {
        "model_id"                    : cmd.model,
        "experiment_id"               : exp_label,
        "ensemble_members_description": ens_label,
        "diagnostic_description"      : cmd.diagnostic,
        "plot_type"                   : "time_series_"
            + f"{lat_eval_n:.0f}n_{abs(lat_eval_s):.0f}s"
    }
    
    fig1, ax1 = qc.start_figure(**descr_kw)
    
    year_ma, tas_ma_n = qc.moving_average(year, tas_n,
                                          cmd.movavg)
    tas_ma_s = qc.moving_average(year, tas_s, cmd.movavg)[1]
    
    # Full time series data:
    qc.plot_data(ax1, xdata=year, ydata=[tas_n, tas_s],
        ens_index=e, ens_label=ens_label, labels=[None]*2,
        linewidths=0.5, ripf_labels=ripf)
    
    # Moving-average filtered data:
    qc.plot_data(ax1, xdata=year_ma, ydata=[tas_ma_n, tas_ma_s],
        ens_index=e, ens_label=ens_label,
        labels=["tas_n", "tas_s"], linewidths=1.5,
        ripf_labels=ripf)
    
    ax1.set_title(f"{cmd.model} {cmd.experiment}"
        + ("" if cmd.x2 == "none" else f" + {cmd.x2}") + ", "
        + ens_label + ", " + full_diag_name +"\ntime series at "
        + qc.lat_eval_label(lat_eval_n) + " / "
        + qc.lat_eval_label(lat_eval_s))
    
    ax1.set_xlabel("Year")
    ax1.set_ylabel(nf.field_units["temperature"])
    ax1.legend(fontsize=8 if "all" in ens_label else 12)
    
    qc.finish_figures([fig1], str(data_path),
        savefig=cmd.savefigs,
        subplots_adjust_kw={"top":0.86, "bottom":0.18})



if __name__ == "__main__":
    main()
