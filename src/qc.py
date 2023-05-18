"""Routines for data in-place processing, plot creation and
formatting, for the purposes of quality control (qc).
"""

from datetime import datetime as dt, timezone as tz
from pathlib import Path

from netCDF4 import num2date
import numpy as np
import matplotlib.pyplot as plt

from api import load_data
from src.metadata import dir_out_nc_data, members, year_range
from src.netcdf import nc_time_units, nc_calendar

from my_python_utilities.data_tools.nc_tools import ncdump

global _ncdump_count
_ncdump_count = 0

# ============================================================ #
# Data post-processing
# ============================================================ #

def profile_mean(t, x, t_avg=(1, 21)):
    """Get mean over specified years from a profile.
    """
    jt = (t >= t_avg[0]) & (t <= t_avg[1])
    return (t[jt][0], t[jt][-1]), np.mean(x[jt], axis=0)



def evaluate_at_latitude(x, lat, lat_eval=65.0):
    """Evaluate profile at specified coordinate."""
    jlat = np.argmin(abs(lat - lat_eval))
    return x[:,:,jlat], lat[jlat]



def seasonal_climatology(t, x, t_avg=(1, 21)):
    """Get seasonal cycle climatology from monthly input data,
    averaging over resulting year indices j_t_avg.
    """
    jt = (t >= t_avg[0]) & (t <= t_avg[1])
    return (t[jt][0], t[jt][-1]), np.nanmean(
        np.reshape(x, (len(x)//12, 12, *np.shape(x)[1:]))[jt],
        axis=0)



def moving_average(t, x, n_ma=21):
    """Moving average filter on x(t), of width n_ma time steps.
    """
    
    nt_ma = len(t) - n_ma
    
    t_ma = np.zeros(nt_ma).astype(t.dtype)
    x_ma = np.zeros((nt_ma, *np.shape(x)[1:])).astype(x.dtype)
    
    for k in range(nt_ma):
        t_ma[k] = np.mean(t[k:k+n_ma])
        x_ma[k] = np.mean(x[k:k+n_ma], axis=0)
    
    return t_ma, x_ma



def contiguous_timestep_averages(x, n_avg=21, anomaly=True):
    """Get contiguous n_avg time-step averages, discarding
    remainder time steps at the end of the time series, as
    anomalies (default).
    """
    
    nolist = type(x) not in [list]
    x_use = [x.copy()] if nolist else x.copy()
    
    for j in range(len(x_use)):
        n_avgs = len(x_use[j]) // n_avg
        
        # If the number of time steps is not an exact multiple
        # of the averaging length, discard the remainder at the
        # end of the time series:
        n_discard = len(x_use[j]) % n_avg
        if n_discard > 0:
            x_use[j] = x_use[j][:-n_discard]
        
        if anomaly:
            x_use[j] -= np.nanmean(x_use[j], axis=0)
        
        x_use[j] = np.nanmean(
            np.reshape(x_use[j],
                       (n_avgs, n_avg, *np.shape(x_use[j])[1:])
                      ), axis=1)
    
    if nolist:
        x_use = x_use[0]
    
    return x_use



def time_average_period_differences(t, x, t_avg_1=(1980, 2000),
                                    t_avg_2=(2001, 2021)):
    """Compute differences between two time averaging periods.
    """
    jt1 = (t >= t_avg_1[0]) & (t <= t_avg_1[1])
    jt2 = (t >= t_avg_2[0]) & (t <= t_avg_2[1])
    
    if type(x) in [list]:
        xret = []
        for xj in x:
            xret.append(np.nanmean(xj[jt2], axis=0)
                        - np.nanmean(xj[jt1], axis=0))
    else:
        xret = np.nanmean(x[jt2], axis=0) \
                - np.nanmean(x[jt1], axis=0)
    
    return xret


# ============================================================ #
# Metadata, labels, etc.
# ============================================================ #

def evaluate_time_member_coords(time, r_vals, i_vals, p_vals,
        f_vals, experiment="piControl"):
    """Returns year integer array for time axes, datetime array
    for time axes, and ripf strings.
    """
    
    ripf = [f"r{r_vals[j]}i{i_vals[j]}p{p_vals[j]}f{f_vals[j]}"
        for j in range(len(r_vals))]
    
    # Assumes no modifications to the default time units and
    # calendar since data was generated:
    date = num2date(time, units=nc_time_units[experiment],
                    calendar=nc_calendar[experiment])
    
    year = np.arange(date[0].year, date[-1].year+1, 1)
    
    return year, date, ripf



def ens_member_label(code, r_vals, ripf):
    """Process command-line argument --ensmemberid (code) into
    an index value for member axes of data arrays (if
    applicable) and string label for plot file names.
    """
    
    if code <= -2:  # show all
        if len(r_vals) > 1:
            e = None
            ens_label = "all_members"
        else:
            e = 0
            ens_label = ripf[0]
    elif code == -1:  # show ensemble mean
        if len(r_vals) > 1:
            e = None
            ens_label = "ensemble_mean"
        else:
            e = 0
            ens_label = ripf[0]
    elif code == 0:  # show first member
        e = 0
        ens_label = ripf[0]
    else:  # code interpreted as specified realisation index
        e = np.nonzero(np.array(r_vals) == code)[0][0]
        ens_label = ripf[e]
    
    return e, ens_label



def lat_eval_label(lat_eval, units="degrees_north"):
    """Default latitude units are degrees_north even for
    southern hemisphere data -- however, it is easier for plot
    titles to have these in degrees_south.
    """
    
    if units == "degrees_north" and lat_eval < 0.0:
        x = f"{abs(lat_eval):.2f} degrees_south"
    else:
        x = f"{lat_eval:.2f} {units}"
    return x



def shorten_data_path_labels(data_paths):
    """Take a (list of) data directory paths and shorten them by
    removing the root output directory (md.dir_out_nc_data).
    Returns a one-longer list where the first is this removed
    part, and the rest are the shortened file paths. These can
    be printed on figures using '\n'.join(data_paths) in the
    call to finish_figures().
    """
    
    x = [data_paths] if type(data_paths) == str else data_paths
    
    xret = ["./ = " + str(dir_out_nc_data)]
    kstart = len(str(dir_out_nc_data)) + 1
    
    xret += ["./" + j[kstart:] for j in x if j != ""]
    
    return xret


# ============================================================ #
# Loading data wrappers
# ============================================================ #

def get_data(diagnostic="dhdt_pot_residual_hfbasin",
        model_id="IPSL-CM6A-LR", experiment_1="piControl",
        experiment_2="none", lat_eval=65.0, select_diag=0,
        ncdump_max=1
    ):
    """Get data and required metadata for specified diagnostic,
    model, experiment(s), and (if applicable) latitude for
    evaluation.
    
    Can load data for one experiment (experiment_1) or two
    concatenated (experiment_1 + experiment_2) if experiment_2
    is not "none".
    
    Paramter diagnostic can be either a keyword (in which case
    the integer diag_select can be used to choose a non-default
    for the specified model), an alias, or the full name (see
    api._diagnostic_definitions).
    
    """
    
    global _ncdump_count
    
    try:
        
        kw = {"model_id"    : model_id,
              "lat_eval"    : None if lat_eval is None
                              else (lat_eval, -lat_eval),
              "select_diags": select_diag}
        
        if experiment_2 == "none":
            
            date, ripf, lat_n, lat_s, diags_n, diags_s, \
                data_paths = load_data.one_experiment(
                    [diagnostic], experiment_id=experiment_1,
                    **kw)
            
            data_path = data_paths[0]
            
        else:
            date, ripf, lat_n, lat_s, diags_n, diags_s, \
                data_paths = load_data.two_experiments(
                    [diagnostic],
                    experiment_ids=(experiment_1, experiment_2),
                    **kw)
            
            data_path = data_paths[0][0]
            data_path = str(data_path).replace(
                experiment_1, "{" + experiment_1 + ","
                              + experiment_2 + "}")
            
        if lat_eval is None:
            # Full coordinates for profile (extract arrays):
            lat_n = lat_n[0]
            lat_s = lat_s[0]
        else:
            # Evaluated at specified latitude (get single
            # values for each hemisphere):
            if experiment_2 == "none":
                lat_n = lat_n[0]
                lat_s = lat_s[0]
            else:
                lat_n = lat_n[0,0]
                lat_s = lat_s[0,0]
        
        # These are the same regardless of one vs two exps.:
        r_vals =np.array([int(x[1:x.index("i")]) for x in ripf])
        year   = np.array([x.year for x in date[0]])
        data_n = diags_n[0]  # load_data returns list of diags
        data_s = diags_s[0]  # (but we only give it one here)
        
        if _ncdump_count < ncdump_max:
            ncdump(str(data_path).replace(
                "{" + experiment_1 + "," + experiment_2 + "}",
                experiment_2))
            _ncdump_count += 1

    
    except FileNotFoundError:
        
        # Both load_data.one_experiment() and
        # load_data.two_experiments() will attempt to load the
        # data for the specified diagnostic, model, and
        # experiments without checking they exist first. So the
        # first point of failure will be a FileNotFoundError.
        # 
        # Here, we may wish to continue anyway if we are loading
        # multiple diagnostics for a case where there are
        # multiple choices (read: oht, dhdt) and we wish to
        # examine a subset of them for a particular model.
        # 
        # Therefore, set the data to be "missing", the metadata
        # can be determined from the metadata module (funnily
        # enough!) and continue anyway, printing a warning:
        # ---------------------------------------------------- #
        
        print("WARNING: data for diagnostic \'" + diagnostic
              + "\' does not exist for " + model_id + " "
              + experiment_1 + ("" if experiment_2 == "none"
                  else (" and/or " + experiment_2)))
        
        if experiment_2 == "none":
            yr_s, yr_e = year_range[experiment_1][model_id]
            ripf = members[model_id][experiment_1]
        else:
            yr_s, _ = year_range[experiment_1][model_id]
            _, yr_e = year_range[experiment_2][model_id]
            ripf = sorted(list(set(
                members[model_id][experiment_1] +
                members[model_id][experiment_2]
            )))
        
        n_t = yr_e - yr_s + 1
        n_ens = len(ripf)
        
        data_n = np.nan*np.ones((n_t, n_ens), dtype=np.float64)
        data_s = np.nan*np.ones((n_t, n_ens), dtype=np.float64)
        year = np.arange(yr_s, yr_e+1, 1)
        r_vals = np.arange(n_ens, dtype=int)
        
        if lat_eval is None:
            lat_n = np.full(1, np.nan)
            lat_s = np.full(1, np.nan)
        else:
            lat_n = np.nan  # load_data.*() functions return NaN
            lat_s = np.nan  # if not applicable/missing anyway
        data_path = ""
    
    return data_n, data_s, year, r_vals, ripf, lat_n, \
           lat_s, str(data_path)



# ============================================================ #
# Plot elements
# ============================================================ #

def global_latitude_xaxis(ax):
    """Set single Axes instance x-axis to global latitude.
    """
    ax.set_xlim([-95.0, 95.0])
    ax.xaxis.set_ticks(np.arange(-90.0, 90.001, 15.0))


def monthly_xaxis(ax):
    """Set single Axes instance x-axis to monthly labels for
    climatology plots (assumes x-axis data is set to 1 for Jan.,
    2 for Feb., ..., 12 for Dec., etc., e.g., np.arange(1,13,1).
    """
    ax.set_xlim([0.5, 12.5])
    ax.xaxis.set_ticks(np.arange(1, 13, 1))
    ax.set_xlabel("Month")



def plot_data(ax, xdata=[], ydata=[],
        ens_index=0, ens_label="ensemble_mean",
        labels=[],
        ripf_labels=[],
        colors=["tab:red", "tab:blue"],
        linestyles=["-"],
        linewidths=[1.5],
        markers=None,
        ens_cmap=None
    ):
    """General line plotting on a single Axes instance. Allows
    multiple x/y data each with different colours, linestyles,
    labels, and selects from the member axes of datasets
    according to the standard ens_label string.
    """
    
    if type(xdata) not in [list, tuple]:
        xdata = [xdata]
    if type(ydata) not in [list, tuple]:
        ydata = [ydata]
    
    n_data = len(ydata)
    n_ens = len(ripf_labels)
    
    if type(colors) not in [list, tuple]:
        colors = [colors]
    if type(linestyles) not in [list, tuple]:
        linestyles = [linestyles]
    if type(linewidths) not in [list, tuple]:
        linewidths = [linewidths]
    if type(markers) not in [list, tuple]:
        markers = [markers]
    
    if ens_label == "all_members":
        
        if len(linestyles) == 1 and len(ydata) == 2:
            linestyles = ["-", "--"]
        
        if ens_cmap is None:
            # Choose default colormap based on
            # number of ensemble members:
            if n_ens <= 10:
                cols = plt.cm.get_cmap("tab10")(
                    np.arange(0.05, 1.0, 0.1))
            elif n_ens <= 20:
                cols = plt.cm.get_cmap("tab20")(
                    np.arange(0.025, 1.0, 0.05))
            else:
                cols = plt.cm.get_cmap("turbo")(
                    np.arange(0.5/n_ens, 1.0, 1.0/n_ens))
        else:
            cols = plt.cm.get_cmap(ens_cmap)(
                np.arange(0.5/n_ens, 1.0, 1.0/n_ens))
        
        for m in range(n_ens):
            for d in range(n_data):
                ax.plot(xdata[d%len(xdata)], ydata[d][:,m],
                    color=cols[m],
                    linewidth=linewidths[d%len(linewidths)],
                    linestyle=linestyles[d%len(linestyles)],
                    marker=markers[d%len(markers)],
                    mfc=cols[m]
                        if linestyles[d%len(linestyles)] == "-"
                        else "none",
                    label=None if labels[d] is None else
                          f"{labels[d]} ({ripf_labels[m]})")
            
    elif ens_label == "ensemble_mean":
        
        for d in range(n_data):
            ax.plot(xdata[d%len(xdata)],
                np.mean(ydata[d], axis=1),
                color=colors[d%len(colors)],
                linewidth=linewidths[d%len(linewidths)],
                linestyle=linestyles[d%len(linestyles)],
                marker=markers[d%len(markers)],
                mfc=colors[d%len(colors)]
                    if linestyles[d%len(linestyles)] == "-"
                    else "none",
                label=labels[d])
        
    else:
        
        for d in range(n_data):
            ax.plot(xdata[d%len(xdata)],
                ydata[d][:,ens_index],
                color=colors[d%len(colors)],
                linewidth=linewidths[d%len(linewidths)],
                linestyle=linestyles[d%len(linestyles)],
                marker=markers[d%len(markers)],
                mfc=colors[d%len(colors)]
                    if linestyles[d%len(linestyles)] == "-"
                    else "none",
                label=labels[d])



def scatter_data_3(axs, datasets=[],
        ax_labels=[],
        legend_label_type="correlation_slope",
        colors=["tab:grey"],
        ax_label_colors=["k"]
    ):
    """Takes three datasets and creates a scatter plot of each
    pair."""
    
    leg_label = [""]*3
    
    if "correlation" in legend_label_type:
        x = "$r={:.3f}$"
        leg_label = [
            x.format(np.corrcoef(datasets[0],datasets[1])[0,1]),
            x.format(np.corrcoef(datasets[0],datasets[2])[0,1]),
            x.format(np.corrcoef(datasets[1],datasets[2])[0,1])
        ]
    
    if "slope" in legend_label_type:
        x = "slope$={:.2f}$"
        if len(leg_label[0]) > 0:
            x = "\n" + x
        
        try:
            m0 = np.polyfit(datasets[0], datasets[1], 1)[0]
        except np.linalg.LinAlgError:
            m0 = np.nan
        try:
            m1 = np.polyfit(datasets[0], datasets[2], 1)[0]
        except np.linalg.LinAlgError:
            m1 = np.nan
        
        try:
            m2 = np.polyfit(datasets[1], datasets[2], 1)[0]
        except np.linalg.LinAlgError:
            m2 = np.nan
        
        leg_label[0] += x.format(m0)
        leg_label[1] += x.format(m1)
        leg_label[2] += x.format(m2)
    
    axs[0].scatter(datasets[0], datasets[1], label=leg_label[0],
        c=colors[0%len(colors)], zorder=10)
    axs[1].scatter(datasets[0], datasets[2], label=leg_label[1],
        c=colors[1%len(colors)], zorder=10)
    axs[2].scatter(datasets[1], datasets[2], label=leg_label[2],
        c=colors[2%len(colors)], zorder=10)
    
    axs[0].set_xlabel(ax_labels[0],
        color=ax_label_colors[0%len(ax_label_colors)])
    axs[1].set_xlabel(ax_labels[0],
        color=ax_label_colors[0%len(ax_label_colors)])
    axs[2].set_xlabel(ax_labels[1],
        color=ax_label_colors[1%len(ax_label_colors)])
    
    axs[0].set_ylabel(ax_labels[1],
        color=ax_label_colors[1%len(ax_label_colors)])
    axs[1].set_ylabel(ax_labels[2],
        color=ax_label_colors[2%len(ax_label_colors)])
    axs[2].set_ylabel(ax_labels[2],
        color=ax_label_colors[2%len(ax_label_colors)])



# ============================================================ #
# Figure creation, properties
# ============================================================ #

def start_figure(ncols=1, nrows=1,
        model_id="unknown_model",
        experiment_id="unknown_experiment",
        ensemble_members_description="unknown_member",
        diagnostic_description="unknown_diagnostic",
        plot_type="profile",
        extra_info=""
    ):
    """Generate figure, axes -- main layout, with metadata set
    in the window title (used for automatic filename generation
    when saving.
    """
    
    fig, ax = plt.subplots(figsize=(9, 4.5),
                           ncols=ncols, nrows=nrows)
    
    fig.canvas.manager.set_window_title(
        f"qc_{model_id}_{experiment_id}_"
        + f"{ensemble_members_description}_"
        + f"{diagnostic_description}_{plot_type}"
        + ("" if extra_info == "" else f"_{extra_info}")
    )
    
    if nrows > 1 or ncols > 1:
        for a in ax.flatten():
            a.grid(color=[0.9]*3)
    else:
        ax.grid(color=[0.9]*3)
    
    return fig, ax



def fig_metadata(fig, data_path,
        x_dt=0.01, y_dt=0.99,
        x_fp=0.01, y_fp=0.01,
        fig_text_kw={}
    ):
    """Add metadata to figure [file path(s) at the bottom,
    timestamp at the top].
    """
    
    # Default properties if not provided:
    if "fontsize" not in fig_text_kw.keys():
        fig_text_kw["fontsize"] = 7
    if "color" not in fig_text_kw.keys():
        fig_text_kw["color"] = [0.6]*3
    
    fig.text(x_dt, y_dt,
             dt.now(tz.utc).strftime("%H:%M UTC %d %b %Y"),
             ha="left", va="top", **fig_text_kw)
    
    fig.text(x_fp, y_fp, str(data_path), ha="left", va="bottom",
             **fig_text_kw)



def finish_figures(fig_list, data_path,
        savefig=False,
        savefig_directory=Path("/home/users/gb919150/phd/",
                               "process_cmip6_data/qc/qc_out"),
        subplots_adjust_kw={"top": 0.9, "bottom":0.15},
        fig_text_kw={}
    ):
    """Apply figure layouts and metadata, and either show them
    interactively (default) or save to file. If saving, the
    figure window title is used as the file name [set in the
    routine start_figure()].
    """
    
    for fig in fig_list:
        fig.tight_layout()
        fig_metadata(fig, data_path, fig_text_kw=fig_text_kw)
        fig.subplots_adjust(**subplots_adjust_kw)
    
    if savefig:
        
        Path(savefig_directory).mkdir(parents=True,
                                      exist_ok=True)
        
        for fig in fig_list:
            save_file = Path(savefig_directory,
                fig.canvas.manager.get_window_title() + ".png")
            fig.savefig(save_file, dpi=200)
            print(f"Saved: {str(save_file)}")
    
    else:
        for fig in fig_list:
            fig.show()
