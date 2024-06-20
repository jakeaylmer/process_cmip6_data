"""Save northward ocean heat transport data computed as a direct
model output (hfbasin) as a yearly-averaged time series, from
input monthly data. Need to run save_areacello.py
first (as the coordinates are obtained from there).
"""

from pathlib import Path
import numpy as np

from process_cmip6_data.src import (
    diagnostics as diags,
    load_processed_data as lpd,
    metadata as md,
    netcdf as nf,
    script_tools,
    utils
)

diag_name = "oht_from_hfbasin"
nc_var_name = "oht"

nc_long_name = ("Northward ocean heat transport calculated "
    + "on the native grid (hfbasin){}")

# Note in this case there are no separate "n" and "s"
# variables as hfbasin is an equivalent approximation
# for the transport in either direction at a given
# latitude (i.e., no accumulation of errors)
nc_var_attrs = dict()
nc_var_attrs["units"] = nf.field_units["heattransport"]
nc_var_attrs["standard_name"] = "northward_ocean_heat_transport"
nc_var_attrs["cell_methods"] = (
    f"{nf.nc_time_name}: mean " +
    f"{nf.nc_ref_lat_single_name}: point")


# Short description added to netCDF "title" attribute (need
# not be completely accurate/detailed here):
nc_title_str = "ocean heat transport"


def process_member(member_id, model_id, experiment_id):
    """Load data for hfbasin for one member, model, and
    experiment. The hfbasin data is somewhat unique to the other
    input fields in that it has a "basin" (domain of zonal
    integration) dimension that is not consistent among models
    with regards to which is global (all we want here) and the 
    name of the latitude variable. This information is coded in
    the metadata module, dictionary hfbasin_metadata, a
    dictionary where the keys are the model_ids and the values
    are lists [<int; global basin index>, <str; latitude
    variable name].
    
    Returns
    -------
    lat : 1D array of float, shape (nlat,)
        Latitude coordinates.
    
    hfbasin : 2D array of float, shape (nyears, nlat)
        Yearly-averaged time series of global hfbasin as a
        function of time (year) and latitude.
    
    """
    data_dir = Path(md.dir_raw_nc_data, model_id, experiment_id,
                    "hfbasin")
    
    data_files_in = sorted([str(x) for x in data_dir.glob(
        f"hfbasin*{model_id}*{experiment_id}*{member_id}*.nc")])
    
    lat, hfbasin = utils.nc_get_arrays(data_files_in,
        [md.hfbasin_metadata[model_id][1]], ["hfbasin"])
    
    lat = np.squeeze(lat)
    hfbasin = np.squeeze(hfbasin)
    
    # Select basin:
    hfbasin = hfbasin[:,md.hfbasin_metadata[model_id][0],:]
    
    # Reset missing values:
    hfbasin = np.where(
        hfbasin >= md.default_original_missing_value,
        md.default_new_missing_value, hfbasin)
    
    # Compute annual means:
    hfbasin_ym = diags.year_mean_1D(hfbasin, nsteps_year=12,
                                    keep_nan=True)
    
    return lat, hfbasin_ym



def main():
    
    cmd = script_tools.default_cmd_args()
    
    yr_s, yr_e = md.year_range[cmd.experiment][cmd.model]
    ens_members = md.members[cmd.model][cmd.experiment]
    
    nt = yr_e - yr_s + 1
    n_ens = len(ens_members)
    
    ens_members = md.members[cmd.model][cmd.experiment]
    n_ens = len(ens_members)
    
    ref_lats = md.default_ref_lats_oht_n
    
    pm_kw = {
        'model_id': cmd.model,
        'experiment_id': cmd.experiment}
    
    # Process first member to get n_lat:
    print(f"Processing member: {ens_members[0]} "
        + f"(1 of {n_ens})")
    
    lat, hfbasin_0 = process_member(ens_members[0], **pm_kw)
    
    _, nlat = np.shape(hfbasin_0)
    
    hfbasin = np.zeros((nt, n_ens, nlat))
    hfbasin[:,0,:] = hfbasin_0
    
    # Remaining members:
    if n_ens > 1:
        for m in range(1, n_ens):
            print(f"Processing member: {ens_members[m]} "
                + f"({m+1} of {n_ens})")
            hfbasin[:,m,:] = process_member(ens_members[m],
                **pm_kw)[1]
    
    hfbasin /= 1.0E15  # in PW
    
    hfbasin_interp = diags.interpolate_to_ref_latitudes(
        lat, hfbasin, ref_lats)
    
    # ------------------------------------------------------- #
    
    print("Saving to NetCDF...")
    
    save_nc_kw = {
        "model_id": cmd.model,
        "member_ids": ens_members,
        "experiment_id": cmd.experiment,
        "year_range": (yr_s, yr_e),
        "nc_title_str": nc_title_str}
    
    diag_kw = {"name": diag_name,
        "time_methods": nf.diag_nq_yearly,
        "other_methods": nf.diag_nq_native}
    
    nc_var_name_kw = {"name": nc_var_name,
        "time_methods": nf.nc_var_nq_yearly,
        "other_methods": nf.nc_var_nq_native}
    
    nf.save_yearly_ref_lat_single(hfbasin, lat,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(hemi="", **nc_var_name_kw),
        nc_field_type=hfbasin.dtype,
        unlimited_ref_lat_dim=False,
        nc_field_attrs={
            "long_name": nc_long_name.format(""),
            **nc_var_attrs},
        **save_nc_kw)
    
    diag_kw["other_methods"] = nf.diag_nq_native_interp
    nc_var_name_kw["other_methods"] = nf.nc_var_nq_native_interp
    
    nf.save_yearly_ref_lat_single(hfbasin_interp, ref_lats,
        nf.diag_name(**diag_kw),
        nf.nc_var_name(hemi="", **nc_var_name_kw),
        nc_field_type=hfbasin_interp.dtype,
        nc_field_attrs={
            "long_name": nc_long_name.format(
                ", interpolated to reference latitudes"
            ),
            **nc_var_attrs},
        **save_nc_kw)


if __name__ == "__main__":
    main()
