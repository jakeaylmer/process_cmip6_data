# This script was run on the originally-processed data, i.e.,
# prior to submission of associated manuscript, on 18 May 2024
# to update the global netCDF attributes prior to archiving of
# the data, upon recent acceptance of said manuscript.
# 
# The updates are:
#     1) Removed title of paper (which has changed anyway) in
#        "comment" attribute
#     2) "coauthors" attribute format
#     3) a) "source" attribute renamed to "source_id" (refers to
#           the model name)
#        b) "source" attribute now used for raw data source
#           (with CMIP6 data citation)
#     4) "references" updated to associated manuscript (
#        accepted in principle at the time of running this
#        script, so no DOI etc. yet), model description paper,
#        and CMIP6/ESGF data reference.
#     5) Rename "institution" attribute to "author_institution"
#     6) Reformat (and prepend) "history" attribute
# 
# These apply to all processed data.
# 
# The default values of and routines which set these attributes
# have also been updated in the module src/netcdf.py (so that if
# processing data again this script would not be needed).
# 
# The script modifies files in place, and is mostly safe against
# errors (but I currently copy files externally before modifying
# with this script anyway).
# ============================================================ #

from argparse import ArgumentParser
from datetime import datetime as dt, UTC
import os
import netCDF4 as nc

from process_cmip6_data.src import metadata as md, netcdf as nf


def update_title(nc_attrs, fname, model_id, experiment_id):
    """Update/set "title" attribute of nc_attrs (dict) based on
    file name fname, (sufficient to uniquely identify what title
    should be).
    
    Originally, "title" was set to:
    
        "CMIP6 diagnostics: <model name>, <experiment name>"
    
    Now, with the update, it is set to:
    
        "CMIP6 processed diagnostics for <model name>,
        <experiment name>: <nc_title_str>"
    
    where nc_title_str is a short description of the diagnostics
    in the file.
    """
    
    # The following dictionary gives filePATH patterns (keys)
    # that determine the nc_title_str as defined above (values):
    # 
    # The patterns must therefore be unique (tested this works
    # fine, there is only one case where it is mistaken for
    # another diagnostic, but including the / part of the file
    # path distinguishes it; see comment below):
    nc_title_str_cases = {
        "ahtc_from"       : "atmospheric heat transport convergence",
        "aht_from"        : "atmospheric heat transport",
        "areacella"       : "areacella",  # does not need updating
        "areacello"       : "areacello",  # does not need updating
        "f_down_area_mean_yr_": "surface downwelling longwave radiation polar-cap averages",
        "f_olr_area_mean_yr_": "outgoing longwave radiation polar-cap averages",
        "f_sw_surf_area_mean_yr_": "surface net downward shortwave radiation polar-cap averages",
        "f_sw_toa_area_mean_yr_": "top of atmosphere net shortwave radiation polar-cap averages",
        "f_up_area_mean_yr_": "surface upward heat flux polar-cap averages",
        "f_down_yr": "surface downwelling longwave radiation",
        "f_olr_yr": "outgoing longwave radiation",
        "f_sw_surf_yr": "surface net downward shortwave radiation",
        "f_sw_toa_yr": "top of atmosphere net shortwave radiation",
        "f_up_yr": "surface upward heat flux",
        "hfds_area_mean_yr_": "heat flux into ocean surface polar-cap averages",
        "hfds_hor_int_yr_": "heat flux into ocean surface polar-cap integrals",
        "hfds_yr": "heat flux into ocean surface",
        "iel_mon_": "sea ice-edge latitude",
        "iel_zm_": "sea ice-edge latitude zonal mean",
        "temptend_ver_int_area_mean_yr_": "ocean column heat content tendency polar-cap averages",
        "temptend_ver_int_from_hfds_hfbasin_hor_int_yr_": "ocean column heat content tendency polar-cap integrals",
        "temptend_ver_int_from_hfds_hfx_hfy_area_mean_yr_": "ocean column heat content tendency polar-cap averages",
        "temptend_ver_int_from_hfds_hfx_hfy_hor_int_yr_": "ocean column heat content tendency polar-cap integrals",
        "temptend_ver_int_hor_int_yr_": "ocean column heat content tendency polar-cap integrals",
        "temptend_ver_int_yr/": "ocean column heat content tendency",  # "/" very important! Distinguishes from OHT diagnotstic, one case where there is not a unique match...
        "ohtc_from_": "ocean heat transport convergence",
        "oht_from_": "ocean heat transport",
        "sia_": "sea ice area",
        "sie_": "sea ice extent",
        "tas_area_mean_": "near-surface air temperature polar-cap averages",
        "tas_yr": "near-surface air temperature"
    }
    
    # Search for the match:
    for j in nc_title_str_cases.keys():
        if j in fname:
            nc_title_str = nc_title_str_cases[j]
            break
    else:
        raise Exception("Unknown title attribute conversion, "
                        + f"file name is {fname}")
    
    # Title does not need updating for areacella/areacello:
    if nc_title_str not in ["areacella", "areacello"]:
        nc_attrs["title"] = ("CMIP6 processed diagnostics for "
            + f"{model_id}, {experiment_id}: {nc_title_str}")


def get_model_experiment_ids_from_filename(fname):
    """If needed, determine the model and experiment from the
    file name (in the main code, preferrably these are obtained
    from the original file attributes; this function just comes
    to the rescue if something has gone wrong). Note that fname
    is the whole path
    """
    
    # Determine the model. In all cases fname must contain the
    # model. Need to be careful of, e.g., "CanESM5" and
    # "CanESM5-CanOE", both of which match for CanESM5-CanOE:
    matches_m = [int(m in fname) for m in md.defined_models]
    n_matches_m = sum(matches_m)
    
    if n_matches_m == 0:
        raise Exception("Something has REALLY gone wrong: "
                        + "cannot determine model (source_id) "
                        + f"for file: {fname}")
    elif n_matches_m == 1:
        model_id = md.defined_models[matches_m.index(1)]
    else:
        # n_matches > 1; choose the longest matching name
        matches_m_names = []  # store names that match
        for mj in range(len(matches_m)):
            if matches_m[mj] == 1:
                matches_m_names += [md.defined_models[mj]]
        
        # get string lengths of each name that matches:
        len_matches_m_names = [len(mn) for mn in
                               matches_m_names]
        
        # choose the longest
        model_id = matches_m_names[len_matches_m_names.index(
                                   max(len_matches_m_names))]
    # -------------------------------------------------------- #
    
    # Determine the experiment. For all outputs, fname contains
    # the experiment except for areacella and areacello.
    # Fortunately, the experiments for those are defined in
    # metadata = md module:
    if "areacell" in fname:
        if "areacella" in fname:
            experiment_id = md.areacella_file_kw[model_id][0]
        elif "areacello" in fname:
            experiment_id = md.areacello_file_kw[model_id][0]
    else:
        for x in md.defined_experiments:
            if x in fname:
                experiment_id = x
                break
        else:
            raise Exception("Something has REALLY gone wrong: "
                            + "cannot determine experiment_id "
                            + f"for file: {fname}")
    
    return model_id, experiment_id


def main():
    
    prsr = ArgumentParser()
    prsr.add_argument("-i", "--infiles", type=str, default=[""],
                      nargs="*")
    cmd = prsr.parse_args()
    
    if any(["*" in x for x in cmd.infiles]):
        # this means pattern matching on the command line has
        # failed (e.g. due to typo)
        cmd.infiles = []
    
    for j in range(len(cmd.infiles)):
        
        # Ensure file can be edited:
        os.chmod(cmd.infiles[j], 0o644)
        
        with nc.Dataset(cmd.infiles[j], "a") as ncdat:
            
            # Get all current global attributes and values:
            nc_global_attrs = {k: ncdat.getncattr(k)
                for k in ncdat.ncattrs()}
            
            # Identify model and experiment
            # 
            # Originally used "source" for model ID; below, that
            # is being changed to "source_id" and "source" is
            # being used for a more descriptive version and
            # citing the data
            # 
            # Experiment was and remains "experiment_id"
            # 
            # If something goes wrong all attributes are removed
            # but all is not lost, as these can be recovered
            # from the file name:
            # 
            if "source_id" in nc_global_attrs.keys():
                # Must be running the code again (e.g., to
                # correct something), as "source_id" was never
                # included. 
                model_id = nc_global_attrs["source_id"]
                experiment_id = nc_global_attrs["experiment_id"]
            
            elif "source" in nc_global_attrs.keys():
                model_id = nc_global_attrs["source"]
                experiment_id = nc_global_attrs["experiment_id"]
            
            else:
                model_id, experiment_id = \
                    get_model_experiment_ids_from_filename(
                        str(cmd.infiles[j]))
            
            # Now delete the attributes from ncdat [this is so
            # they can be (re)assigned alphabetically, whereas
            # modifying each one in place changes the order.
            # This doesn't really matter but is just
            # aesthetically pleasing this way]:
            for k in nc_global_attrs.keys():
                ncdat.delncattr(k)
            
            # Reset "coauthors" attribute:
            nc_global_attrs["coauthors"] = \
                nf.default_nc_file_attrs["coauthors"]
            
            # Change "institution" attribute to
            # "author_institution" (to remove any ambiguity with
            # institution of data production; the original/
            # 'source' data was produced at a different
            # institution to the 'final' data in these files):
            if "institution" in nc_global_attrs.keys():
                del nc_global_attrs["institution"]
            
            nc_global_attrs["author_institution"] = \
                nf.default_nc_file_attrs["author_institution"]
            
            # Update the title. It can just be set/no danger of
            # overwriting.
            # 
            # Originally did not include any specific
            # information about what diagnostics were in the
            # file. It's obvious from the variable attributes,
            # but CF guideline/convention is to include a short
            # summary of what the dataset contains.
            # 
            # The saving routines in the src/netcdf module have
            # been modified to accept a string for such
            # purposes. See docs of update_title() function:
            update_title(nc_global_attrs,
                str(cmd.infiles[j]), model_id,
                experiment_id)
            
            # Update comment (see netcdf = nf module).
            # 
            # Originally it referred to an earlier draft of the
            # associated manuscript, but the title has changed
            # and anyway it was poorly phrased.
            # 
            # This is the same for all diagnostics except
            # areacella and areacello. Here, just copy from nf
            # and areacella/o scripts.
            # 
            # First part is the same, regardless:
            nc_global_attrs["comment"] = ("Climate "
                + "model diagnostics derived from "
                + "outputs in the CMIP6 archive for "
                + "the analysis presented in Aylmer et "
                + "al. 2024 [1]. This dataset contains ")
            
            if "areacell" in str(cmd.infiles[j]):
                
                if "areacella" in str(cmd.infiles[j]):
                    domain = "atmosphere"
                elif "areacello" in str(cmd.infiles[j]):
                    domain = "ocean"
                
                nc_global_attrs["comment"] += (domain
                    + " grid data for one model "
                    + f"({nf.nc_file_attrs_model_name}), "
                    + "taken from one experiment ("
                    + nf.nc_file_attrs_experiment_name
                    + ") and ensemble member ("
                    + f"{nf.nc_file_attrs_member_name}).")
            
            else:
                nc_global_attrs["comment"] += ("one "
                    + "diagnostic for one model "
                    + f"({nf.nc_file_attrs_model_name}) "
                    + "and one experiment ("
                    + nf.nc_file_attrs_experiment_name
                    + ") [2,3].")
            
            ref_count = nf._set_nc_references_attr(
                nc_global_attrs, model_id, experiment_id)
            
            # Reset source/source_id attributes:
            # (now source refers to data source with citation,
            # and source_id is just the model name):
            # 
            # Ad-hoc check:
            if experiment_id.startswith("esm-"):
                # it's either "esm-piControl" or "esm-hist",
                # in both cases the activity, needed for the
                # source attribute, is "CMIP". I did not want to
                # include those in metadata as
                # md.defined_experiments since it only occurs
                # twice (the other attributes required for
                # references/metadata for these cases are
                # included specifically for GFLD-ESM4 in the
                # metadata module, and it doesn't cause issues
                # for any other global attribute set in the
                # netcdf = nf module):
                md.activity[experiment_id] = "CMIP"
            
            nf._set_nc_source_attr(nc_global_attrs, model_id,
                experiment_id, reference_number=ref_count)
            
            nf._set_nc_source_id_attr(nc_global_attrs, model_id)
            
            # Add to history audit trail:
            if "history" in nc_global_attrs.keys():
                
                if "preparation for archiving" in nc_global_attrs["history"]:
                    
                    # In this case I must be running the code on
                    # the same file again. So recover the original
                    # history string so that the "preparation for
                    # archiving" part does not appear repeatedly:
                    hist_strs = nc_global_attrs["history"].split("\n")
                    nc_global_attrs["history"] = hist_strs[-1]
                    # (this assumes no other modifications to
                    # the history attribute have been made apart
                    # from what is in this file)
                
                else:
                    # Running from originally saved data. The
                    # created part of this was originally in the
                    # format: "created %H:%M UTC %d %b %Y"; for
                    # consistency, switch this around:
                    hist_strs = nc_global_attrs["history"].split(" ")
                    nc_global_attrs["history"] = \
                        " ".join(hist_strs[1:]) + ": created;"
                
                nc_global_attrs["history"] = (
                    dt.now(UTC).strftime("%H:%M UTC %d %b %Y")
                    + ": updated global attributes in "
                    + "preparation for archiving;\n"
                    + nc_global_attrs["history"])
            
            # Set attributes:
            for k in sorted(list(nc_global_attrs.keys())):
                ncdat.setncattr(k, nc_global_attrs[k])


if __name__ == "__main__":
    main()
