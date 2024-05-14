from argparse import ArgumentParser
from pathlib import Path
import numpy as np
import netCDF4 as nc

from process_cmip6_data.src import (
    diagnostics as diags,
    load_raw_data as lrd,
    metadata as md,
    netcdf as nf,
    script_tools
)

from my_python_utilities.data_tools import nc_tools as nct

# Basic diagnostic name (details such as time averaging, zonal
# mean, etc., and the corresponding netCDF variable name, are
# added automatically as specified in the netcdf module):
diag_name = "tas"

nc_long_name = ("Near-surface air temperature, averaged "
    + "monthly and everywhere {}ward of {}")

nc_var_attrs_n = {
    "units"        : nf.field_units["temperature"],
    "standard_name": "air_temperature",
    "cell_methods" : f"{nf.nc_time_name}: mean "
                     + f"{nf.nc_ref_lat_n_name}: mean "
                     + "(area-weighted)"
}

nc_var_attrs_s = nc_var_attrs_n.copy()
nc_var_attrs_s["cell_methods"] = (f"{nf.nc_time_name}: "
    + f"mean {nf.nc_ref_lat_s_name}: mean (area-weighted)")


def save_monthly_ref_lat(fint_n, fint_s, ref_lat_n, ref_lat_s,
        diagnostic_id, nc_field_name_n, nc_field_name_s,
        model_id, member_ids, experiment_id,
        year_range=(1, 500),
        nc_field_type_n=np.float64,
        nc_field_type_s=np.float64,
        ref_lat_n_bnds=None, ref_lat_s_bnds=None,
        unlimited_ref_lat_dim=True,
        nc_field_attrs_n={},
        nc_field_attrs_s={},
        nc_global_attrs={},
        save_dir=None, file_name=None
    ):
    """Save a monthly data for integrations/averages above a
    reference latitude (one variable for north, one for south).
    """
    
    save_dir, file_name = \
        nf.prepare_to_save(save_dir, file_name, diagnostic_id,
                        experiment_id, model_id)
    
    # Prepare the global attributes: ------------------------ #
    nc_file_attrs = nf.default_nc_file_attrs.copy()
    nf.set_nc_title_attr(nc_file_attrs, model_id, experiment_id)
    nf.set_nc_history_attr_timestamp(nc_file_attrs)
    nf.set_nc_source_attr(nc_file_attrs, model_id)
    nf.set_nc_experiment_attr(nc_file_attrs, experiment_id)
    
    for extra_attr in nc_global_attrs.keys():
        nc_file_attrs[extra_attr] = nc_global_attrs[extra_attr]
    # ------------------------------------------------------- #
    
    with nc.Dataset(Path(save_dir, file_name), "w") as ncout:
        
        nf._set_attributes(ncout, nc_file_attrs)
        
        nf.create_time_variables(ncout,
                                 nf.nc_time_units[experiment_id],
                                 nf.nc_calendar[experiment_id])
        nf.set_time_data_monthly(ncout, year_range)
        nf.set_ripf_data(ncout, member_ids)
        
        ncout.createDimension(nf.nc_ref_lat_n_dim_name,
            None if unlimited_ref_lat_dim else len(ref_lat_n))
        ncout.createDimension(nf.nc_ref_lat_s_dim_name,
            None if unlimited_ref_lat_dim else len(ref_lat_s))
        
        if (ref_lat_n_bnds is not None
                or ref_lat_s_bnds is not None
            ):
            ncout.createDimension(nf.nc_vert_dim_name, 2)
        
        nf.set_ref_lat_data(ncout, ref_lat_n, ref_lat_s,
                            ref_lat_n_bnds, ref_lat_s_bnds)
        
        for f_name, f_type, dim_name, data, attrs in zip(
                [nc_field_name_n, nc_field_name_s],
                [nc_field_type_n, nc_field_type_s],
                [nf.nc_ref_lat_n_dim_name, nf.nc_ref_lat_s_dim_name],
                [fint_n, fint_s],
                [nc_field_attrs_n, nc_field_attrs_s]
            ):
            
            ncout.createVariable(f_name,
                f_type, (nf.nc_time_dim_name,
                    nf.nc_member_dim_name, dim_name),
                fill_value=nf.nc_fill_value, zlib=True
            )
            ncout.variables[f_name][:,:,:] = data
            nf._set_attributes(ncout.variables[f_name],
                               attrs)
    
    print(f"Saved: {str(Path(save_dir, file_name))}")


def main():
    
    prsr = ArgumentParser()
    
    prsr.add_argument("-r", "--reanalysis", type=str,
                      default=md.default_reanalysis,
                      choices=md.defined_reanalyses,
                      help="Reanalysis name (case sensitive)")
    
    prsr.add_argument("-a", "--approx", action="store_true",
                      help="Do approximate version ("
                           + f"{nf.diag_nq_cell_center_approx})"
                           + " anyway")
    
    cmd = prsr.parse_args()
    
    yr_s, yr_e = md.reanalyses_year_range[cmd.reanalysis]
    
    ref_lats_n = md.default_ref_lats_n
    ref_lats_s = md.default_ref_lats_s
    
    ref_lats_n_bnds = np.array([[ref_lats_n[k], 90.0]
        for k in range(len(ref_lats_n))])
    ref_lats_s_bnds = np.array([[-90.0, ref_lats_s[k]]
        for k in range(len(ref_lats_s))])
    
    # Load the processed areacella data for reanalyses using
    # the raw-data loading routines (which are more generic
    # than the processed-data loading routines and can be
    # adapted here for the reanalysis data):
    lon, lat, lat_bnds, areacella = lrd._load_coordinate_arrays(
        Path(md.dir_out_nc_data_reanalyses, "areacella",
             f"areacella_{cmd.reanalysis}.nc"),
        [nf.nc_lon_1d_name, nf.nc_lat_1d_name,
         f"{nf.nc_lat_1d_name}_bnds", "areacella"])
    
    # Load standardised-raw tas data (manually, as there is no
    # generic routine in src):
    print("Loading raw data...")
    diag_name_load = nf.diag_name(name="tas",
        time_methods=nf.diag_nq_monthly)
    
    nc_file = Path(md.dir_out_nc_data_reanalyses,
                   diag_name_load,
                   nf.nc_file_name(diag_name_load, "reanalysis",
                                   cmd.reanalysis))
    
    tas, = nct.get_arrays([nc_file], [],
        [nf.nc_var_name(name="tas",
                        time_methods=nf.nc_var_nq_monthly)])
    
    # For saving:
    diag_kw = {"name": diag_name,
        "time_methods": nf.diag_nq_monthly,
        "space_methods": nf.diag_nq_area_mean}
    
    nc_var_kw = {"name": diag_name,
        "time_methods": nf.nc_var_nq_monthly,
        "space_methods": nf.nc_var_nq_area_mean}
    
    save_nc_kw = {
        "model_id"       : cmd.reanalysis,
        "member_ids"     : ["r1i1p1f1"],
        "experiment_id"  : "reanalysis",
        "year_range"     : (yr_s, yr_e),
        "nc_global_attrs": {
            "comment": "Atmospheric reanalysis diagnostics for "
                       + "analysis in the work, \'Modulation "
                       + "of future sea ice loss by ocean heat "
                       + "transport\'. This dataset contains "
                       + "one diagnostic for one reanalysis "
                       + "(source).",
            "source": md.reanalyses_long_names[cmd.reanalysis],
            "title": "Atmospheric reanalysis diagnostics: "
                + cmd.reanalysis,
            "references": md.reanalyses_references[cmd.reanalysis]
        }
    }
    
    if np.ndim(lat) == 1:
        
        true_lats_n, true_lats_n_bnds, true_tas_avg_n = \
            diags.integrate_horizontally_exact(tas, lat_bnds,
                areacella, hemi="n", normalise=True)
        
        true_lats_s, true_lats_s_bnds, true_tas_avg_s = \
            diags.integrate_horizontally_exact(tas, lat_bnds,
                areacella, hemi="s", normalise=True)
        
        interp_tas_avg_n = diags.interpolate_to_ref_latitudes(
            true_lats_n, true_tas_avg_n, ref_lats_n)
        interp_tas_avg_s = diags.interpolate_to_ref_latitudes(
            true_lats_s, true_tas_avg_s, ref_lats_s)
        
        # --------------------------------------------------- #
        
        print("Saving to NetCDF...")
        
        diag_kw["other_methods"] = nf.diag_nq_native
        nc_var_kw["other_methods"] = nf.nc_var_nq_native
        
        save_monthly_ref_lat(
            true_tas_avg_n, true_tas_avg_s,
            true_lats_n, true_lats_s,
            nf.diag_name(**diag_kw),
            nf.nc_var_name(hemi="n", **nc_var_kw),
            nf.nc_var_name(hemi="s", **nc_var_kw),
            ref_lat_n_bnds=true_lats_n_bnds,
            ref_lat_s_bnds=true_lats_s_bnds,
            unlimited_ref_lat_dim=False,
            nc_field_attrs_n={
                "long_name": nc_long_name.format("north",
                    "native-grid latitude bounds"),
                **nc_var_attrs_n},
            nc_field_attrs_s={
                "long_name": nc_long_name.format("south",
                    "native-grid latitude bounds"),
                **nc_var_attrs_s},
            save_dir=Path(md.dir_out_nc_data_reanalyses,
                          nf.diag_name(**diag_kw)),
            **save_nc_kw
        )
        
        diag_kw["other_methods"] = nf.diag_nq_native_interp
        nc_var_kw["other_methods"] = nf.nc_var_nq_native_interp
        
        save_monthly_ref_lat(
            interp_tas_avg_n, interp_tas_avg_s,
            ref_lats_n, ref_lats_s,
            nf.diag_name(**diag_kw),
            nf.nc_var_name(hemi="n", **nc_var_kw),
            nf.nc_var_name(hemi="s", **nc_var_kw),
            ref_lat_n_bnds=ref_lats_n_bnds,
            ref_lat_s_bnds=ref_lats_s_bnds,
            nc_field_attrs_n={
                "long_name": nc_long_name.format("north",
                    "native-grid latitude bounds, interpolated "
                    + "to reference latitudes"),
                **nc_var_attrs_n},
            nc_field_attrs_s={
                "long_name": nc_long_name.format("south",
                    "native-grid latitude bounds, interpolated "
                    + "to reference latitudes"),
                **nc_var_attrs_s},
            save_dir=Path(md.dir_out_nc_data_reanalyses,
                          nf.diag_name(**diag_kw)),
            **save_nc_kw
        )
        
    
    if np.ndim(lat) != 1 or cmd.approx:
        
        # Have to do it the "old fashioned" way
        # (integrate cells poleward of reference latitude,
        # based on cell center coordinates)
        # 
        # Either that, or want the approximate version anyway
        # in which case need to create mesh grid latitude if
        # it is 1D:
        if np.ndim(lat) == 1:
            lon, lat = np.meshgrid(lon, lat)
        
        # This function returns both the integral and mean,
        # so just take second return value:
        tas_avg_n = diags.integrate_horizontally_multi_members(
            tas, areacella, lat, ref_lats=ref_lats_n,
            hemi="n")[1]
        
        tas_avg_s = diags.integrate_horizontally_multi_members(
            tas, areacella, lat, ref_lats=ref_lats_s,
            hemi="s")[1]
        
        # --------------------------------------------------- #
        
        print("Saving to NetCDF...")
        
        diag_kw["other_methods"] = \
            nf.diag_nq_cell_center_approx
        nc_var_kw["other_methods"] = \
            nf.nc_var_nq_cell_center_approx
        
        save_monthly_ref_lat(tas_avg_n, tas_avg_s,
            ref_lats_n, ref_lats_s,
            nf.diag_name(**diag_kw),
            nf.nc_var_name(hemi="n", **nc_var_kw),
            nf.nc_var_name(hemi="s", **nc_var_kw),
            ref_lat_n_bnds=ref_lats_n_bnds,
            ref_lat_s_bnds=ref_lats_s_bnds,
            nc_field_attrs_n={
                "long_name": nc_long_name.format("north",
                    "reference latitudes based on cell-center "
                    + "latitudes"),
                **nc_var_attrs_n},
            nc_field_attrs_s={
                "long_name": nc_long_name.format("south",
                    "reference latitudes based on cell-center "
                    + "latitudes"),
                **nc_var_attrs_s},
            save_dir=Path(md.dir_out_nc_data_reanalyses,
                          nf.diag_name(**diag_kw)),
            **save_nc_kw
        )



if __name__ == "__main__":
    main()
