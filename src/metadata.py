"""Main settings and model metadata, to be filled in manually.
"""

from pathlib import Path
import numpy as np

# =========================================================== #

default_model_id             = "CESM2"
default_member_id            = "r1i1p1f1"
default_experiment_id        = "piControl"
default_future_experiment_id = "ssp370"

default_original_missing_value = np.float32(1.0E20)
default_new_missing_value      = np.nan

# For sea ice (siconc thresholds in units of fraction):
default_sie_threshold = 0.15
default_iel_threshold = 0.15

# For oht/aht:
default_ref_lats_oht_n = np.arange(-80.0, 90.01, 1.0)
default_ref_lats_oht_s = np.arange(-80.0, 90.01, 1.0)

default_ref_lats_aht_n = np.arange(-90.0, 90.01, 1.0)
default_ref_lats_aht_s = np.arange(-90.0, 90.01, 1.0)

# For area integrals:
default_ref_lats_n = np.concatenate(([0], np.arange(50.0, 85.01, 1.)))
default_ref_lats_s = np.concatenate(([0], np.arange(-50.0, -85.01, -1.)))

# NOTE: for opottemptend derived as a residual from hfbasin and
# ----- hfds, need to have the oht and area-integral reference
# latitudes somewhat matching up for consistency (i.e., all lats
# defined for area integrals should be in that for oht, 
# otherwise a "nearest-neighbour" approach is used when adding
# the two quantities to get the residual in this case).

defined_experiments = ["piControl", "historical", "ssp370",
                       "ssp585"]

activity = {
    "piControl"    : "CMIP",
    "historical"   : "CMIP",
    "ssp370"       : "ScenarioMIP",
    "ssp585"       : "ScenarioMIP"
}

dir_raw_nc_data = Path("/storage", "basic", "cpom", "gb919150",
                       "CMIP6")
dir_out_nc_data = Path("/storage", "silver", "cpom", "gb919150",
                       "phd", "cmip6_processed_data")

# =========================================================== #

# Start and end years. For PI control this is arbitrary and
# defined differently across models. Here, I just set them
# all to start at the year 1:
year_range = {
    "piControl" : {
        "AWI-CM-1-1-MR": (1, 500),
        "CanESM5"      : (1, 1051),
        "CanESM5-CanOE": (1, 501),
        "CESM2"        : (1, 1200),
        "CESM2-FV2"    : (1, 150),
        "CESM2-WACCM"  : (1, 499),
        "CNRM-CM6-1"   : (1, 500),
        "CNRM-CM6-1-HR": (1, 300),
        "CNRM-ESM2-1"  : (1, 500),
        "GFDL-ESM4"    : (1, 500),
        "GISS-E2-2-G"  : (1, 296),
        "IPSL-CM6A-LR" : (1, 1000),
        "MIROC6"       : (1, 500),
        "MPI-ESM1-2-HR": (1, 500),
        "MPI-ESM1-2-LR": (1, 1000),
        "MRI-ESM2-0"   : (1, 701),
        "NorESM2-LM"   : (1, 501),
        "NorESM2-MM"   : (1, 500),
        "UKESM1-0-LL"  : (1, 1880),
        "UKESM1-1-LL"  : (1, 462)
    }
}

defined_models = list(year_range["piControl"].keys())

# For historical, all models have the same year range, but set
# for each model anyway to have a consistent data structure:
year_range["historical"] = dict.fromkeys(defined_models,
                                         (1850, 2014))

# Some models extend ssp245 and ssp585 (and others?) beyond 2100
# but for simplicity I just use up to 2100 (note that this is
# not accounted for when loading raw data, so if extension files
# are downloaded automatically, they must be deleted manually to
# avoid "shape mismatch along axis 0" errors):
for exp in ["ssp370", "ssp585"]:
    year_range[exp] = dict.fromkeys(defined_models,
                                    (2015, 2100))

# =========================================================== #
# ENSEMBLE MEMBERS
# ----------------
# 
# List of ensemble members per model, per experiment
# =========================================================== #

members = {
    "AWI-CM-1-1-MR": {
        "piControl" : ["r1i1p1f1"],
        "historical": [f"r{x}i1p1f1" for x in range(1, 6)],
        "ssp370"    : [f"r{x}i1p1f1" for x in range(1, 6)],
        "ssp585"    : ["r1i1p1f1"]
    },
    "CanESM5": {
        "piControl" : ["r1i1p2f1"],
        "historical": [f"r{x}i1p2f1" for x in range(1, 26)],
        "ssp370"    : [f"r{x}i1p2f1" for x in range(1, 26)],
        "ssp585"    : [f"r{x}i1p2f1" for x in range(1, 26)]
    },
    "CanESM5-CanOE": {
        "piControl" : ["r1i1p2f1"],
        "historical": ["r1i1p2f1", "r2i1p2f1", "r3i1p2f1"],
        "ssp370"    : ["r1i1p2f1", "r2i1p2f1", "r3i1p2f1"],
        "ssp585"    : ["r1i1p2f1", "r2i1p2f1", "r3i1p2f1"]
    },
    "CESM2": {
        "piControl" : ["r1i1p1f1"],
        "historical": ["r4i1p1f1", "r10i1p1f1", "r11i1p1f1"],
        "ssp370"    : ["r4i1p1f1", "r10i1p1f1", "r11i1p1f1"],
        "ssp585"    : ["r4i1p1f1", "r10i1p1f1", "r11i1p1f1"]
    },
    "CESM2-FV2":
        dict.fromkeys(defined_experiments, ["r1i2p2f1"]),
    "CESM2-WACCM": {
        "piControl" : ["r1i1p1f1"],
        "historical": ["r1i1p1f1", "r2i1p1f1", "r3i1p1f1"],
        "ssp370"    : ["r1i1p1f1"],
        "ssp585"    : ["r1i1p1f1", "r2i1p1f1", "r3i1p1f1"]
    },
    "CNRM-CM6-1": {
        "piControl": ["r1i1p1f2"],
        "historical": [f"r{x}i1p1f2" for x in range(1, 7)],
        "ssp370": [f"r{x}i1p1f2" for x in range(1, 7)],
        "ssp585": [f"r{x}i1p1f2" for x in range(1, 7)]
    },
    "CNRM-CM6-1-HR":
        dict.fromkeys(defined_experiments, ["r1i1p1f1"]),
    "CNRM-ESM2-1":
        dict.fromkeys(defined_experiments, ["r1i1p1f2"]),
    "GFDL-ESM4":
        dict.fromkeys(defined_experiments, ["r1i1p1f1"]),
    "GISS-E2-2-G": {
        "piControl" : ["r1i1p3f1"],
        "historical": [f"r{x}i1p3f1" for x in range(1,6)],
        "ssp370"    : [f"r{x}i1p3f1" for x in range(1,6)],
        "ssp585"    : [f"r{x}i1p3f1" for x in range(1,6)]
    },
    "IPSL-CM6A-LR": {
        "piControl" : ["r1i1p1f1"],
        "historical": [f"r{x}i1p1f1" for x in range(1, 33)],
        "ssp370"    : [f"r{x}i1p1f1" for x in range(1,11)]
                    + ["r14i1p1f1"],
        "ssp585"    : [f"r{x}i1p1f1" for x in range(1,5)]
                    + ["r6i1p1f1", "r14i1p1f1"]
    },
    "MIROC6":
        dict.fromkeys(defined_experiments, ["r1i1p1f1"]),
    "MPI-ESM1-2-HR": {
        "piControl" : ["r1i1p1f1"],
        "historical": [f"r{x}i1p1f1" for x in range(1, 11)],
        "ssp370"    : [f"r{x}i1p1f1" for x in range(1, 11)],
        "ssp585"    : ["r1i1p1f1", "r2i1p1f1"]
    },
    "MPI-ESM1-2-LR": {
        "piControl" : ["r1i1p1f1"],
        "historical": [f"r{x}i1p1f1" for x in range(1, 11)],
        "ssp370"    : [f"r{x}i1p1f1" for x in range(1, 11)],
        "ssp585"    : [f"r{x}i1p1f1" for x in range(1, 11)]
    },
    "MRI-ESM2-0": {
        "piControl" : ["r1i1p1f1"],
        "historical": [f"r{x}i1p1f1" for x in range(1, 6)],
        "ssp370"    : [f"r{x}i1p1f1" for x in range(1, 6)],
        "ssp585"    : [f"r{x}i1p1f1" for x in range(1, 6)]
    },
    "NorESM2-LM": {
        "piControl" : ["r1i1p1f1"],
        "historical": [f"r{x}i1p1f1" for x in range(1, 4)],
        "ssp370"    : ["r1i1p1f1"],
        "ssp585"    : ["r1i1p1f1"]
    },
    "NorESM2-MM":
        dict.fromkeys(defined_experiments, ["r1i1p1f1"]),
    "UKESM1-0-LL" : {
        "piControl" : ["r1i1p1f2"],
        "historical": [f"r{x}i1p1f2" for x in range(1, 4+1)]
                    + [f"r{x}i1p1f3" for x in range(5, 7+1)]
                    + [f"r{x}i1p1f2" for x in range(8, 12+1)]
                    + [f"r{x}i1p1f2" for x in range(16, 19+1)],
        "ssp370"    : [f"r{x}i1p1f2" for x in range(1, 12+1)]
                    + [f"r{x}i1p1f2" for x in range(16, 19+1)],
        "ssp585"    : [f"r{x}i1p1f2" for x in range(1, 4+1)]
                    + ["r8i1p1f2"]
    },
    "UKESM1-1-LL" : {
        "piControl" : ["r1i1p1f2"],
        "historical": ["r1i1p1f2"],
        "ssp370"    : ["r1i1p1f2"]
    }
}


# =========================================================== #
# DATA REFERENCES
# ---------------
# Information for data references (could include one dictionary
# and single entry for each model/experiment, but some of the
# information is useful elsewhere (e.g., for constructing netCDF
# global metadata attributes) as well as for the data
# references.
# 
# Here, need the following which depend on model and in most
# cases also on experiment:
#     > institution_id        [model_id] [experiment_id]
#     > model_long_name       [model_id]
#     > data_version          [model_id] [experiment_id]
#     > data_doi              [model_id] [experiment_id]
#     > data_publication_year [model_id] [experiment_id]
#     > data_author_list      [model_id] [experiment_id]
# 
# Version labels are most up-to-date as of 14 May 2024 (and
# unlikely to change at this stage...?)
# 
# In each case above which depend on experiment_id as well,
# include an empty string "" for each defined experiment if
# there is no data for that experiment_id and model_id (this is
# so that the automatic reference generating doesn't need to get
# to complicated with checking, although that would probably be
# more robust).
# 
# After these are set, the dictionary "data_reference" is
# constructed (here) automatically, so that:
# 
# meta_data.data_reference[model_id][experiment_id] gives the
# formatted ESGF data reference for model_id, experiment_id.
# 
# =========================================================== #
institution_id = {
    "AWI-CM-1-1-MR": dict.fromkeys(defined_experiments, "AWI"),
    "CanESM5"      : dict.fromkeys(defined_experiments, "CCCma"),
    "CanESM5-CanOE": dict.fromkeys(defined_experiments, "CCCma"),
    "CESM2"        : dict.fromkeys(defined_experiments, "NCAR"),
    "CESM2-FV2"    : dict.fromkeys(defined_experiments, "NCAR"),
    "CESM2-WACCM"  : dict.fromkeys(defined_experiments, "NCAR"),
    "CNRM-CM6-1"   : dict.fromkeys(defined_experiments, "CNRM-CERFACS"),
    "CNRM-CM6-1-HR": dict.fromkeys(defined_experiments, "CNRM-CERFACS"),
    "CNRM-ESM2-1"  : dict.fromkeys(defined_experiments, "CNRM-CERFACS"),
    "IPSL-CM6A-LR" : dict.fromkeys(defined_experiments, "IPSL"),
    "GFDL-ESM4"    : dict.fromkeys(defined_experiments
                                   + ["esm-piControl", "esm-hist"],
                                   "NOAA-GFDL"),
    "GISS-E2-2-G"  : dict.fromkeys(defined_experiments, "NASA-GISS"),
    "MIROC6"       : dict.fromkeys(defined_experiments, "MIROC"),
    "MPI-ESM1-2-HR": {
        "piControl" : "MPI-M",
        "historical": "MPI-M",
        "ssp370"    : "DKRZ",
        "ssp585"    : "DKRZ"
    },
    "MPI-ESM1-2-LR": dict.fromkeys(defined_experiments, "MPI-M"),
    "MRI-ESM2-0"   : dict.fromkeys(defined_experiments, "MRI"),
    "NorESM2-LM"   : dict.fromkeys(defined_experiments, "NCC"),
    "NorESM2-MM"   : dict.fromkeys(defined_experiments, "NCC"),
    "UKESM1-0-LL"  : dict.fromkeys(defined_experiments, "MOHC"),
    "UKESM1-1-LL"  : {
        "piControl" : "MOHC",
        "historical": "MOHC",
        "ssp370"    : "MOHC",
        "ssp585"    : "N/A"
    }
}

# These are only used in the dataset references (and should
# match what the model is called in the ESGF citations; see
# their DOIs):
model_long_name = {
    "AWI-CM-1-1-MR": "AWI-CM-1.1-MR",
    "CanESM5"      : "CanESM5",
    "CanESM5-CanOE": "CanESM5-CanOE",
    "CESM2"        : "CESM2",
    "CESM2-FV2"    : "CESM2-FV2",
    "CESM2-WACCM"  : "CESM2-WACCM",
    "CNRM-CM6-1"   : "CNRM-CM6-1",
    "CNRM-CM6-1-HR": "CNRM-CM6-1-HR",
    "CNRM-ESM2-1"  : "CNRM-ESM2-1",
    "IPSL-CM6A-LR" : "IPSL-CM6A-LR",
    "GFDL-ESM4"    : "GFDL-ESM4",
    "GISS-E2-2-G"  : "GISS-E2-2-G",
    "MIROC6"       : "MIROC6",
    "MPI-ESM1-2-HR": "MPI-ESM1.2-HR",
    "MPI-ESM1-2-LR": "MPI-ESM1.2-LR",
    "MRI-ESM2-0"   : "MRI-ESM2.0",
    "NorESM2-LM"   : "NorESM2-LM",
    "NorESM2-MM"   : "NorESM2-MM",
    "UKESM1-0-LL"  : "UKESM1.0-LL",
    "UKESM1-1-LL"  : "UKESM1.1-LL"
}


# data_version are the YYYYMMDD-formatted versions used and
# recorded in the data citations. It is assumed that the dataset
# versions are the same for every member and raw field
# (impractical to double check whether that is true and highly
# unlikely to make any difference):
data_version = {
    "AWI-CM-1-1-MR": {
        "piControl" : "20191015",
        "historical": "20200720",
        "ssp370"    : "20190529",
        "ssp585"    : "20190529"
    },
    "CanESM5"      : dict.fromkeys(defined_experiments, "20190429"),
    "CanESM5-CanOE": dict.fromkeys(defined_experiments, "20190429"),
    "CESM2"        : {
        "piControl" : "20190320",
        "historical": "20190308",
        "ssp370"    : "20200528",
        "ssp585"    : "20200528"
    },
    "CESM2-FV2"    : dict.fromkeys(defined_experiments, "20220915"),
    "CESM2-WACCM"  : {
        "piControl" : "20190320",
        "historical": "20190227",
        "ssp370"    : "20190815",
        "ssp585"    : "20200702"
    },
    "CNRM-CM6-1"   : {
        "piControl" : "20180814",
        "historical": "20180917",
        "ssp370"    : "20190219",
        "ssp585"    : "20190219"
    },
    "CNRM-CM6-1-HR": {
        "piControl" : "20191021",
        "historical": "20191021",
        "ssp370"    : "20200127",
        "ssp585"    : "20191202"
    },
    "CNRM-ESM2-1"  : {
        "piControl" : "20181115",
        "historical": "20181206",
        "ssp370"    : "20191021",
        "ssp585"    : "20191021"
    },
    "GFDL-ESM4"    : {
        "piControl"    : "20180701",
        "esm-piControl": "20180701",
        "historical"   : "20190726",
        "esm-hist"     : "20180701",
        "ssp370"       : "20180701",
        "ssp585"       : "20180701"
    },
    "GISS-E2-2-G"  : {
        "piControl" : "20211002",
        "historical": "20211020",
        "ssp370"    : "20211015",
        "ssp585"    : "20211015"
    },
    "IPSL-CM6A-LR" : {
        "piControl" : "20200326",
        "historical": "20180803",
        "ssp370"    : "20190119",
        "ssp585"    : "20190903"
    },
    "MIROC6"       : {
        "piControl" : "20181212",
        "historical": "20181212",
        "ssp370"    : "20190627",
        "ssp585"    : "20190627"
    },
    "MPI-ESM1-2-HR": dict.fromkeys(defined_experiments, "20190710"),
    "MPI-ESM1-2-LR": dict.fromkeys(defined_experiments, "20190710"),
    "MRI-ESM2-0"   : {
        "piControl" : "20190222",
        "historical": "20190222",
        "ssp370"    : "20190222",
        "ssp585"    : "20191108"
    },
    "NorESM2-LM"   : {
        "piControl" : "20210118",
        "historical": "20190815",
        "ssp370"    : "20191108",
        "ssp585"    : "20191108"
    },
    "NorESM2-MM"   : {
        "piControl" : "20191108",
        "historical": "20191108",
        "ssp370"    : "20191108",
        "ssp585"    : "20230616"
    },
    "UKESM1-0-LL"  : {
        "piControl" : "20200828",
        "historical": "20190406",
        "ssp370"    : "20190510",
        "ssp585"    : "20190507"
    },
    "UKESM1-1-LL"  : {
        "piControl" : "20220505",
        "historical": "20220512",
        "ssp370"    : "20220512",
        "ssp585"    : ""
    }
}

# DOIs; this could be made more efficient since they all start
# with "10.22033/ESGF/CMIP6.":
data_doi = {
    "AWI-CM-1-1-MR": {
        "piControl" : "10.22033/ESGF/CMIP6.2777",
        "historical": "10.22033/ESGF/CMIP6.2686",
        "ssp370"    : "10.22033/ESGF/CMIP6.2803",
        "ssp585"    : "10.22033/ESGF/CMIP6.2817"
    },
    "CanESM5"      : {
        "piControl" : "10.22033/ESGF/CMIP6.3673",
        "historical": "10.22033/ESGF/CMIP6.3610",
        "ssp370"    : "10.22033/ESGF/CMIP6.3690",
        "ssp585"    : "10.22033/ESGF/CMIP6.3696"
    },
    "CanESM5-CanOE": {
        "piControl" : "10.22033/ESGF/CMIP6.10266",
        "historical": "10.22033/ESGF/CMIP6.10260",
        "ssp370"    : "10.22033/ESGF/CMIP6.10271",
        "ssp585"    : "10.22033/ESGF/CMIP6.10276"
    },
    "CESM2"        : {
        "piControl" : "10.22033/ESGF/CMIP6.7733",
        "historical": "10.22033/ESGF/CMIP6.7627",
        "ssp370"    : "10.22033/ESGF/CMIP6.7753",
        "ssp585"    : "10.22033/ESGF/CMIP6.7768"
    },
    "CESM2-FV2"    : {
        "piControl" : "10.22033/ESGF/CMIP6.11301",
        "historical": "10.22033/ESGF/CMIP6.11297",
        "ssp370"    : "",
        "ssp585"    : ""
    },
    "CESM2-WACCM"  : {
        "piControl" : "10.22033/ESGF/CMIP6.10094",
        "historical": "10.22033/ESGF/CMIP6.10071",
        "ssp370"    : "10.22033/ESGF/CMIP6.10102",
        "ssp585"    : "10.22033/ESGF/CMIP6.10115"
    },
    "CNRM-CM6-1"   : {
        "piControl" : "10.22033/ESGF/CMIP6.4163",
        "historical": "10.22033/ESGF/CMIP6.4066",
        "ssp370"    : "10.22033/ESGF/CMIP6.4197",
        "ssp585"    : "10.22033/ESGF/CMIP6.4224"
    },
    "CNRM-CM6-1-HR": {
        "piControl" : "10.22033/ESGF/CMIP6.4164",
        "historical": "10.22033/ESGF/CMIP6.4067",
        "ssp370"    : "10.22033/ESGF/CMIP6.4198",
        "ssp585"    : "10.22033/ESGF/CMIP6.4225"
    },
    "CNRM-ESM2-1"  : {
        "piControl" : "10.22033/ESGF/CMIP6.4165",
        "historical": "10.22033/ESGF/CMIP6.4068",
        "ssp370"    : "10.22033/ESGF/CMIP6.4199",
        "ssp585"    : "10.22033/ESGF/CMIP6.4226"
    },
    "GFDL-ESM4"    : {
        "piControl"    : "10.22033/ESGF/CMIP6.8669",
        "esm-piControl": "10.22033/ESGF/CMIP6.8536",
        "historical"   : "10.22033/ESGF/CMIP6.8597",
        "esm-hist"     : "10.22033/ESGF/CMIP6.8522",
        "ssp370"       : "10.22033/ESGF/CMIP6.8691",
        "ssp585"       : "10.22033/ESGF/CMIP6.8706"
    },
    "GISS-E2-2-G"  : {
        "piControl" : "10.22033/ESGF/CMIP6.7382",
        "historical": "10.22033/ESGF/CMIP6.7129",
        "ssp370"    : "10.22033/ESGF/CMIP6.11873",
        "ssp585"    : "10.22033/ESGF/CMIP6.11889"
    },
    "IPSL-CM6A-LR" : {
        "piControl" : "10.22033/ESGF/CMIP6.5251",
        "historical": "10.22033/ESGF/CMIP6.5195",
        "ssp370"    : "10.22033/ESGF/CMIP6.5265",
        "ssp585"    : "10.22033/ESGF/CMIP6.5271"
    },
    "MIROC6"       : {
        "piControl" : "10.22033/ESGF/CMIP6.5711",
        "historical": "10.22033/ESGF/CMIP6.5603",
        "ssp370"    : "10.22033/ESGF/CMIP6.5752",
        "ssp585"    : "10.22033/ESGF/CMIP6.5771"
    },
    "MPI-ESM1-2-HR": {
        "piControl" : "10.22033/ESGF/CMIP6.6674",
        "historical": "10.22033/ESGF/CMIP6.6594",
        "ssp370"    : "10.22033/ESGF/CMIP6.4399",
        "ssp585"    : "10.22033/ESGF/CMIP6.4403"
    },
    "MPI-ESM1-2-LR": {
        "piControl" : "10.22033/ESGF/CMIP6.6675",
        "historical": "10.22033/ESGF/CMIP6.6595",
        "ssp370"    : "10.22033/ESGF/CMIP6.6695",
        "ssp585"    : "10.22033/ESGF/CMIP6.6705"
    },
    "MRI-ESM2-0"   : {
        "piControl" : "10.22033/ESGF/CMIP6.6900",
        "historical": "10.22033/ESGF/CMIP6.6842",
        "ssp370"    : "10.22033/ESGF/CMIP6.6915",
        "ssp585"    : "10.22033/ESGF/CMIP6.6929"
    },
    "NorESM2-LM"   : {
        "piControl" : "10.22033/ESGF/CMIP6.8217",
        "historical": "10.22033/ESGF/CMIP6.8036",
        "ssp370"    : "10.22033/ESGF/CMIP6.8268",
        "ssp585"    : "10.22033/ESGF/CMIP6.8319"
    },
    "NorESM2-MM"   : {
        "piControl" : "10.22033/ESGF/CMIP6.8221",
        "historical": "10.22033/ESGF/CMIP6.8040",
        "ssp370"    : "10.22033/ESGF/CMIP6.8270",
        "ssp585"    : "10.22033/ESGF/CMIP6.8321"
    },
    "UKESM1-0-LL"  : {
        "piControl" : "10.22033/ESGF/CMIP6.6298",
        "historical": "10.22033/ESGF/CMIP6.6113",
        "ssp370"    : "10.22033/ESGF/CMIP6.6347",
        "ssp585"    : "10.22033/ESGF/CMIP6.6405"
    },
    "UKESM1-1-LL"  : {
        "piControl" : "10.22033/ESGF/CMIP6.16823",
        "historical": "10.22033/ESGF/CMIP6.16797",
        "ssp370"    : "10.22033/ESGF/CMIP6.16845",
        "ssp585"    : ""
    }
}

# Only provide URL if corresponding DOI does not exist:
data_url = {
    "CESM2-FV2": {
        "ssp370": "http://cera-www.dkrz.de/WDCC/meta/CMIP6/CMIP6.ScenarioMIP.NCAR.CESM2-FV2.ssp370",
        "ssp585": "http://cera-www.dkrz.de/WDCC/meta/CMIP6/CMIP6.ScenarioMIP.NCAR.CESM2-FV2.ssp585"
    }
}

# "Year" entry of references (does not necessarily match
# dataset versions or date accessed):
data_publication_year = {
    "AWI-CM-1-1-MR": {
        "piControl" : "2018",
        "historical": "2018",
        "ssp370"    : "2019",
        "ssp585"    : "2019"
    },
    "CanESM5"      : dict.fromkeys(defined_experiments, "2019"),
    "CanESM5-CanOE": dict.fromkeys(defined_experiments, "2019"),
    "CESM2"        : dict.fromkeys(defined_experiments, "2019"),
    "CESM2-FV2"    : {
        "piControl" : "2019",
        "historical": "2019",
        "ssp370"    : "2023",
        "ssp585"    : "2023"
    },
    "CESM2-WACCM"  : dict.fromkeys(defined_experiments, "2019"),
    "CNRM-CM6-1"   : {
        "piControl" : "2018",
        "historical": "2018",
        "ssp370"    : "2019",
        "ssp585"    : "2019"
    },
    "CNRM-CM6-1-HR": {
        "piControl" : "2019",
        "historical": "2019",
        "ssp370"    : "2020",
        "ssp585"    : "2019"
    },
    "CNRM-ESM2-1"  : {
        "piControl" : "2018",
        "historical": "2018",
        "ssp370"    : "2019",
        "ssp585"    : "2019"
    },
    "GFDL-ESM4"    : dict.fromkeys(defined_experiments
                                   + ["esm-piControl", "esm-hist"],
                                   "2018"),
    "GISS-E2-2-G"  : {
        "piControl" : "2019",
        "historical": "2019",
        "ssp370"    : "2021",
        "ssp585"    : "2021"
    },
    "IPSL-CM6A-LR" : {
        "piControl" : "2018",
        "historical": "2018",
        "ssp370"    : "2019",
        "ssp585"    : "2019"
    },
    "MIROC6"       : {
        "piControl" : "2018",
        "historical": "2018",
        "ssp370"    : "2019",
        "ssp585"    : "2019"
    },
    "MPI-ESM1-2-HR": dict.fromkeys(defined_experiments, "2019"),
    "MPI-ESM1-2-LR": dict.fromkeys(defined_experiments, "2019"),
    "MRI-ESM2-0"   : dict.fromkeys(defined_experiments, "2019"),
    "NorESM2-LM"   : dict.fromkeys(defined_experiments, "2019"),
    "NorESM2-MM"   : dict.fromkeys(defined_experiments, "2019"),
    "UKESM1-0-LL"  : dict.fromkeys(defined_experiments, "2019"),
    "UKESM1-1-LL"  : dict.fromkeys(defined_experiments, "2022")
}


# Authors for the dataset citations (note these differ per
# experiment, and are not the same as the model description
# references which come later):
data_author_list = {
    "AWI-CM-1-1-MR": dict.fromkeys(defined_experiments,
                     "Semmler, T., S. Danilov, T. Rackow, D. "
                     + "Sidorenko, D. Barbi, J. Hegewald, "
                     + "and coauthors"),
    "CanESM5"      : dict.fromkeys(defined_experiments,
                     "Swart, N. C., J. N. S. Cole, V. V. "
                     + "Kharin, M. Lazare, J. F. Scinocca, N. "
                     + "P. Gillett, and coauthors"),
    "CanESM5-CanOE": dict.fromkeys(defined_experiments,
                     "Swart, N. C., J. N. S. Cole, V. V. "
                     + "Kharin, M. Lazare, J. F. Scinocca, N. "
                     + "P. Gillett, and coauthors"),
    "CESM2"        : {
        "piControl" : "Danabasoglu, G., D. Lawrence, K. "
                      + "Lindsay, W. Lipscomb, and G. Strand",
        "historical": "Danabasoglu, G.",
        "ssp370"    : "Danabasoglu, G.",
        "ssp585"    : "Danabasoglu, G."
    },
    "CESM2-FV2"    : dict.fromkeys(defined_experiments,
                     "Danabasoglu, G."),
    "CESM2-WACCM"  : dict.fromkeys(defined_experiments,
                     "Danabasoglu, G."),
    "CNRM-CM6-1"   : dict.fromkeys(defined_experiments,
                     "Voldoire, A."),
    "CNRM-CM6-1-HR": dict.fromkeys(defined_experiments,
                     "Voldoire, A."),
    "CNRM-ESM2-1"  : {
        "piControl" : "Séférian, R.",
        "historical": "Séférian, R.",
        "ssp370"    : "Voldoire, A.",
        "ssp585"    : "Voldoire, A."
    },
    "GFDL-ESM4"    : {
        "piControl"    : "Krasting, J. P., J. G. John, C. "
                         + "Blanton, C. McHugh, S. Nikonov, A. "
                         + "Radhakrishnan, and coauthors",
        "esm-piControl": "Krasting, J. P., J. G. John, C. "
                         + "Blanton, C. McHugh, S. Nikonov, A. "
                         + "Radhakrishnan, and coauthors",
        "historical"   : "Krasting, J. P., J. G. John, C. "
                         + "Blanton, C. McHugh, S. Nikonov, A. "
                         + "Radhakrishnan, and coauthors",
        "esm-hist"     : "Krasting, J. P., J. G. John, C. "
                         + "Blanton, C. McHugh, S. Nikonov, A. "
                         + "Radhakrishnan, and coauthors",
        "ssp370"       : "John, J. G., C. Blanton, C. McHugh, "
                         + "A. Radhakrishnan, K. Rand, H. "
                         + "Vahlenkamp, and coauthors",
        "ssp585"       : "John, J. G., C. Blanton, C. McHugh, "
                         + "A. Radhakrishnan, K. Rand, H. "
                         + "Vahlenkamp, and coauthors"
    },
    "GISS-E2-2-G"  : dict.fromkeys(defined_experiments,
                     "NASA Goddard Institute for Space Studies "
                     + "(NASA/GISS)"),
    "IPSL-CM6A-LR" : dict.fromkeys(defined_experiments,
                     "Boucher, O., S. Denvil, G. Levavasseur, "
                     + "A. Cozic, A. Caubel, M.-A. Foujols, "
                     + "and coauthors"),
    "MIROC6"       : {
        "piControl" : "Tatebe, H. and M. Watanabe",
        "historical": "Tatebe, H. and M. Watanabe",
        "ssp370"    : "Shiogama, H., M. Abe, and H. Tatebe",
        "ssp585"    : "Shiogama, H., M. Abe, and H. Tatebe"
    },
    "MPI-ESM1-2-HR": {
        "piControl" : "Jungclaus, J., M. Bittner, K.-H. "
                     + "Wieners, F. Wachsmann, M. Schupfner, "
                     + "S. Legutke, and coauthors",
        "historical": "Jungclaus, J., M. Bittner, K.-H. "
                     + "Wieners, F. Wachsmann, M. Schupfner, "
                     + "S. Legutke, and coauthors",
        "ssp370"    : "Schupfner, M., K.-H. Wieners, F. "
                      + "Wachsmann, C. Steger, M. Bittner, J. "
                      + "Jungclaus, and coauthors",
        "ssp585"    : "Schupfner, M., K.-H. Wieners, F. "
                      + "Wachsmann, C. Steger, M. Bittner, J. "
                      + "Jungclaus, and coauthors"
    },
    "MPI-ESM1-2-LR": dict.fromkeys(defined_experiments,
                     "Wieners, K.-H., M. Giorgetta, J. "
                     + "Jungclaus, C. Reick, M. Esch, M. "
                     + "Bittner, and coauthors"),
    "MRI-ESM2-0"   : dict.fromkeys(defined_experiments,
                     "Yukimoto, S., T. Koshiro, H. Kawai, N. "
                     + "Oshima, K. Yoshida, S. Urakawa, "
                     + "and coauthors"),
    "NorESM2-LM"   : dict.fromkeys(defined_experiments,
                     "Seland, Ø., M. Bentsen, D. J. L. "
                     + "Oliviè, T. Toniazzo, A. Gjermundsen, "
                     + "L. S. Graff, and coauthors"),
    "NorESM2-MM"   : dict.fromkeys(defined_experiments,
                     "Bentsen, M., D. J. L. Oliviè, Ø. "
                     + "Seland, T. Toniazzo, A. Gjermundsen, "
                     + "Ada; L. S. Graff, and coauthors"),
    "UKESM1-0-LL"  : {
        "piControl" : "Tang, Y., S. Rumbold, R. Ellis, D. "
                      + "Kelley, J. Mulcahy, A. Sellar, "
                      + "and coauthors",
        "historical": "Tang, Y., S. Rumbold, R. Ellis, D. "
                      + "Kelley, J. Mulcahy, A. Sellar, "
                      + "and coauthors",
        "ssp370"    : "Good, P., A. Sellar, Y. Tang, S. "
                      + "Rumbold, R. Ellis, D. Kelley, and "
                      + "T. Kuhlbrodt",
        "ssp585"    : "Good, P., A. Sellar, Y. Tang, S. "
                      + "Rumbold, R. Ellis, D. Kelley, and "
                      + "T. Kuhlbrodt"
    },
    "UKESM1-1-LL" : {
        "piControl" : "Mulcahy, J., S. Rumbold, Y. Tang, J. "
                      + "Walton, C. Hardacre, M. Stringer, "
                      + "and coauthors",
        "historical": "Mulcahy, J., S. Rumbold, Y. Tang, J. "
                      + "Walton, C. Hardacre, M. Stringer, "
                      + "and coauthors",
        "ssp370"    : "Walton, J., J. Mulcahy, Y. Tang, S. "
                      + "Rumbold, C. Hardacre, M. Stringer, "
                      + "and coauthors",
        "ssp585"    : ""
    }
}

# Generate dataset references:
# ----------------------------
# 
# {author_list}, {publication_year}: {institution_id}
# {model_long_name} model output prepared for CMIP6 {activity}
# {experiment}, version {data_version}, Earth System Grid
# Federation, {data_link}
# 
dataset_reference_fmt_str = ("{}, {}: {} {} model output "
                             + "prepared for CMIP6 {} {}, "
                             + "version {}, Earth System Grid "
                             + "Federation, {}")

data_reference = {}
for m in defined_models:
    data_reference[m] = {}
    for x in defined_experiments:
        
        # Work out data link (either doi or URL):
        if (m in data_doi.keys()
                and x in data_doi[m].keys()
                and data_doi[m][x] != ""):
            
            fmt_link = f"doi:{data_doi[m][x]}"
        
        elif (m in data_url.keys()
                and x in data_url[m].keys()
                and data_url[m][x] != ""):
            
            fmt_link = f"URL: {data_url[m][x]}"
        
        else:
            fmt_link = ""
        
        data_reference[m][x] = dataset_reference_fmt_str.format(
            data_author_list[m][x], data_publication_year[m][x],
            institution_id[m][x], model_long_name[m],
            activity[x], x, data_version[m][x], fmt_link)

# Add any ad-hoc cases:
m = "GFDL-ESM4"
for x in ["esm-piControl", "esm-hist"]:
    data_reference[m][x] = dataset_reference_fmt_str.format(
        data_author_list[m][x], data_publication_year[m][x],
            institution_id[m][x], model_long_name[m],
            "CMIP", x, data_version[m][x],
            f"doi:{data_doi[m][x]}")


# =========================================================== #
# MODEL DESCRIPTION REFERENCES:
# -----------------------------
# One standard reference per model describing the model
# components, configuration, and/or general performance
# evaluation. These are just 'hard-coded' ('hard-formatted')
# for simplicity (they are all in different journals, the
# authors and other metadata are different to those of the above
# data references, and in any case are not in a consistent
# format like the data references -- i.e., there's no benefit to
# separating out the reference elements here).
# =========================================================== #
model_reference = {
    "AWI-CM-1-1-MR": "Semmler, T., S. Danilov, P. Gierz, H. F. "
                     + "Goessling, J. Hegewald, C. Hinrichs, "
                     + "and coauthors, 2020: Simulations for "
                     + "CMIP6 with the AWI Climate Model "
                     + "AWI-CM-1-1, J. Adv. Model. Earth Syst.,"
                     + " 12, e2019MS002009, "
                     + "doi:10.1029/2019MS002009",
    "CanESM5"      : "Swart, N. C., J. N. S. Cole, V. V. "
                     + "Kharin, M. Lazare, J. F. Scinocca, N. "
                     + "P. Gillett, and coauthors, 2019: The "
                     + "Canadian Earth System Model version 5 "
                     + "(CanESM5), Geosci. Model Dev., 12, "
                     + "4823-4873, "
                     + "doi:10.5194/gmd-12-4823-2019",
    "CESM2"        : "Danabasoglu, G., J. F. Lamarque, J. "
                     + "Bacmeister, D. A. Bailey, A. K. "
                     + "DuVivier, J. Edwards, and coauthors, "
                     + "2020: The Community Earth System Model "
                     + "version 2 (CESM2), J. Adv. Model. "
                     + "Earth Syst., 12, e2019MS001916, "
                     + "doi:10.1029/2019MS001916",
    "CNRM-CM6-1"   : "Voldoire, A., D. Saint-Martin, S. Sénési,"
                     + " B. Decharme, A. Alias, M. Chevallier, "
                     + "and coauthors, 2019: Evaluation of "
                     + "CMIP6 DECK experiments with CNRM-CM6-1,"
                     + " J. Adv. Model. Earth Syst., 11, "
                     + "2177-2213, doi:10.1029/2019MS001683",
    "CNRM-ESM2-1"  : "Séférian, R., P. Nabat, M. Michou, D. "
                     + "Saint-Martin, A. Voldoire, J. Colin, "
                     + "and coauthors, 2019: Evaluation of CNRM"
                     + " Earth System Model, CNRM-ESM2-1: role "
                     + "of Earth system processes in present-"
                     + "day and future climate, J. Adv. Model. "
                     + "Earth Syst., 11, 4182-4227, "
                     + "doi:10.1029/2019MS001791",
    "GFDL-ESM4"    : "Dunne, J. P., L. W. Horowitz, A. J. "
                     + "Adcroft, P. Ginoux, M. Held, J. G. "
                     + "John, and coauthors, 2020: The GFDL "
                     + "Earth System Model version 4.1 "
                     + "(GFDL-ESM 4.1): overall coupled model "
                     + "description and simulation "
                     + "characteristics, J. Adv. Model. Earth "
                     + "Syst., 12, e2019MS002015, "
                     + "doi:10.1029/2019MS002015",
    "GISS-E2-2-G"  : "Rind, D., C. Orbe, J. Jonas, L. "
                     + "Nazarenko, T. Zhou, M. Kelley, and "
                     + "coauthors, 2020: GISS model E2.2: a "
                     + "climate model optimized for the middle "
                     + "atmosphere - model structure, "
                     + "climatology, variability, and climate "
                     + "sensitivity, J. Geophys. Res., 125, "
                     + "e2019JD032204, "
                     + "doi:10.1029/2019JD032204",
    "IPSL-CM6A-LR" : "Boucher, O., J. Servonnat, A. L. "
                     + "Albright, O. Aumont, Y. Balkanski, V. "
                     + "Bastrikov, and coauthors, 2020: "
                     + "Presentation and evaluation of the "
                     + "IPSL-CM6A-LR climate model, J. Adv. "
                     + "Model. Earth Syst., 12, e2019MS002010, "
                     + "doi:10.1029/2019MS002010",
    "MIROC6"       : "Tatebe, H., T. Ogura, T. Nitta, Y. "
                     + "Komuro, K. Ogochi, T. Takemura, and "
                     + "coauthors, 2019: Description and basic "
                     + "evaluation of simulated mean state, "
                     + "internal variability, and climate "
                     + "sensitivity in MIROC6, Geosci. Model "
                     + "Dev., 12, 2727-2765, "
                     + "doi:10.5194/gmd-12-2727-2019",
    "MPI-ESM1-2-HR": "Müller, W. A., J. H. Jungclaus, T. "
                     + "Mauritsen, J. Baehr, M. Bittner, R. "
                     + "Budich, and coauthors, 2018: A higher-"
                     + "resolution version of the Max Planck "
                     + "Institute Earth System Model "
                     + "(MPI-ESM1.2-HR), J. Adv. Model. Earth "
                     + "Syst., 10, 1383-1413, "
                     + "doi:10.1029/2017MS001217",
    "MPI-ESM1-2-LR": "Mauritsen, T., J. Bader, T. Becker, J. "
                     + "Behrens, M. Bittner, R. Brokopf, and "
                     + "coauthors, 2019: Developments in the "
                     + "MPI-M Earth System Model version 1.2 "
                     + "(MPI-ESM1.2) and its response to "
                     + "increasing CO2, J. Adv. Model. Earth "
                     + "Syst., 11, 998-1038, "
                     + "doi:10.1029/2018MS001400",
    "MRI-ESM2-0"   : "Yukimoto, S., H. Kawai, T. Koshiro, N. "
                     + "Oshima, K. Yoshida, S. Urakawa, and "
                     + "coauthors, 2019: The Meteorological "
                     + "Research Institute Earth System Model "
                     + "version 2.0, MRI-ESM2.0: description "
                     + "and basic evaluation of the physical "
                     + "component, J. Meteorol. Soc. Japan, "
                     + "97, 931-965, doi:10.2151/jmsj.2019-051",
    "NorESM2-LM"   : "Seland, Ø., M. Bentsen, D. Olivié, T. "
                     + "Toniazzo, A. Gjermundsen, L. S. Graff, "
                     + "and coauthors, 2020: Overview of the "
                     + "Norwegian Earth System Model (NorESM2) "
                     + "and key climate response of CMIP6 "
                     + "DECK, historical, and scenario "
                     + "simulations, Geosci. Model Dev., 13, "
                     + "6165-6200, "
                     + "doi:10.5194/gmd-13-6165-2020",
    "UKESM1-0-LL"  : "Sellar, A. A., C. G. Jones, J. P. "
                     + "Mulcahy, Y. Tang, A. Yool, A. "
                     + "Wiltshire, and coauthors, 2019: "
                     + "UKESM1: description and evaluation of "
                     + "the U. K. Earth System Model, J. Adv. "
                     + "Model. Earth Syst., 11, 4513-4588, "
                     + "doi:10.1029/2019MS001739"
}

# Some model reference papers apply to different variations of
# the same core model:
model_reference["CanESM5-CanOE"] = model_reference["CanESM5"]
model_reference["CESM2-FV2"] = model_reference["CESM2"]
model_reference["CESM2-WACCM"] = model_reference["CESM2"]
model_reference["CNRM-CM6-1-HR"] = model_reference["CNRM-CM6-1"]
model_reference["NorESM2-MM"] = model_reference["NorESM2-LM"]
model_reference["UKESM1-1-LL"] = model_reference["UKESM1-0-LL"]
# =========================================================== #


# =========================================================== #
# GRID PROPERTIES
# ---------------
# Dimensions, coordinate netCDF names, location of cell area
# and vertices, etc.
# 
# =========================================================== #

# Atmosphere grid dimensions (ny_atm, nx_atm):
grid_dims_atm = {
    "AWI-CM-1-1-MR": (192, 384),
    "CanESM5"      : (64, 128),
    "CanESM5-CanOE": (64, 128),
    "CESM2"        : (192, 288),
    "CESM2-FV2"    : (96, 144),
    "CESM2-WACCM"  : (192, 288),
    "CNRM-CM6-1"   : (128, 256),
    "CNRM-CM6-1-HR": (360, 720),
    "CNRM-ESM2-1"  : (128, 256),
    "IPSL-CM6A-LR" : (143, 144),
    "GFDL-ESM4"    : (180, 288),
    "GISS-E2-2-G"  : (90, 144),
    "MIROC6"       : (128, 256),
    "MPI-ESM1-2-HR": (192, 384),
    "MPI-ESM1-2-LR": (96, 192),
    "MRI-ESM2-0"   : (160, 320),
    "NorESM2-LM"   : (96, 144),
    "NorESM2-MM"   : (192, 288),
    "UKESM1-0-LL"  : (144, 192),
    "UKESM1-1-LL"  : (144, 192)
}


# Atmosphere grid cell areas, absolute paths:
areacella_file_kw = {
    "AWI-CM-1-1-MR": ("piControl", "r1i1p1f1", "gn"),
    "CanESM5"      : ("piControl", "r1i1p2f1", "gn"),
    "CanESM5-CanOE": ("piControl", "r1i1p2f1", "gn"),
    "CESM2"        : ("piControl", "r1i1p1f1", "gn"),
    "CESM2-FV2"    : ("piControl", "r1i1p1f1", "gn"),
    "CESM2-WACCM"  : ("piControl", "r1i1p1f1", "gn"),
    "CNRM-CM6-1"   : ("piControl", "r1i1p1f2", "gr"),
    "CNRM-CM6-1-HR": ("piControl", "r1i1p1f2", "gr"),
    "CNRM-ESM2-1"  : ("piControl", "r1i1p1f2", "gr"),
    "GFDL-ESM4"    : ("piControl", "r1i1p1f1", "gr1"),
    "GISS-E2-2-G"  : ("piControl", "r1i1p1f1", "gn"),
    "IPSL-CM6A-LR" : ("piControl", "r1i1p1f1", "gr"),
    "MIROC6"       : ("piControl", "r1i1p1f1", "gn"),
    "MPI-ESM1-2-HR": ("piControl", "r1i1p1f1", "gn"),
    "MPI-ESM1-2-LR": ("piControl", "r1i1p1f1", "gn"),
    "MRI-ESM2-0"   : ("piControl", "r1i1p1f1", "gn"),
    "NorESM2-LM"   : ("piControl", "r1i1p1f1", "gn"),
    "NorESM2-MM"   : ("piControl", "r1i1p1f1", "gn"),
    "UKESM1-0-LL"  : ("piControl", "r1i1p1f2", "gn"),
    "UKESM1-1-LL"  : ("piControl", "r1i1p1f2", "gn")
}

def areacella_file(model_id):
    return Path(dir_raw_nc_data, model_id,
        areacella_file_kw[model_id][0],
        f"areacella_fx_{model_id}_"
        + f"{areacella_file_kw[model_id][0]}_"
        + f"{areacella_file_kw[model_id][1]}_"
        + f"{areacella_file_kw[model_id][2]}.nc"
    )


# Grid coordinate names in raw netCDF data (lon, lat):
lonlat_atm_nc_names = {
    "AWI-CM-1-1-MR": ("lon", "lat"),
    "CanESM5"      : ("lon", "lat"),
    "CanESM5-CanOE": ("lon", "lat"),
    "CESM2"        : ("lon", "lat"),
    "CESM2-FV2"    : ("lon", "lat"),
    "CESM2-WACCM"  : ("lon", "lat"),
    "CNRM-CM6-1"   : ("lon", "lat"),
    "CNRM-CM6-1-HR": ("lon", "lat"),
    "CNRM-ESM2-1"  : ("lon", "lat"),
    "GFDL-ESM4"    : ("lon", "lat"),
    "GISS-E2-2-G"  : ("lon", "lat"),
    "IPSL-CM6A-LR" : ("lon", "lat"),
    "MIROC6"       : ("lon", "lat"),
    "MPI-ESM1-2-HR": ("lon", "lat"),
    "MPI-ESM1-2-LR": ("lon", "lat"),
    "MRI-ESM2-0"   : ("lon", "lat"),
    "NorESM2-LM"   : ("lon", "lat"),
    "NorESM2-MM"   : ("lon", "lat"),
    "UKESM1-0-LL"  : ("lon", "lat"),
    "UKESM1-1-LL"  : ("lon", "lat")
}

lonlat_bnds_atm_nc_names = {
    "AWI-CM-1-1-MR": ("lon_bnds", "lat_bnds"),
    "CanESM5"      : ("lon_bnds", "lat_bnds"),
    "CanESM5-CanOE": ("lon_bnds", "lat_bnds"),
    "CESM2"        : ("lon_bnds", "lat_bnds"),
    "CESM2-FV2"    : ("lon_bnds", "lat_bnds"),
    "CESM2-WACCM"  : ("lon_bnds", "lat_bnds"),
    "CNRM-CM6-1"   : None,
    "CNRM-CM6-1-HR": None,
    "CNRM-ESM2-1"  : None,
    "GFDL-ESM4"    : ("lon_bnds", "lat_bnds"),
    "GISS-E2-2-G"  : ("lon_bnds", "lat_bnds"),
    "IPSL-CM6A-LR" : None,
    "MIROC6"       : ("lon_bnds", "lat_bnds"),
    "MPI-ESM1-2-HR": ("lon_bnds", "lat_bnds"),
    "MPI-ESM1-2-LR": ("lon_bnds", "lat_bnds"),
    "MRI-ESM2-0"   : ("lon_bnds", "lat_bnds"),
    "NorESM2-LM"   : ("lon_bnds", "lat_bnds"),
    "NorESM2-MM"   : ("lon_bnds", "lat_bnds"),
    "UKESM1-0-LL"  : ("lon_bnds", "lat_bnds"),
    "UKESM1-1-LL"  : ("lon_bnds", "lat_bnds")
}

def lonlat_bnds_atm_file(model_id):
    if lonlat_bnds_atm_nc_names[model_id] is None:
        raise FileNotFoundError
    else:
        return areacella_file(model_id)



# Ocean grid dimensions (ny_ocn, nx_ocn):
grid_dims_ocn = {
    "AWI-CM-1-1-MR": (),
    "CanESM5"      : (291, 360),
    "CanESM5-CanOE": (291, 360),
    "CESM2"        : (384, 320),
    "CESM2-FV2"    : (384, 320),
    "CESM2-WACCM"  : (384, 320),
    "CNRM-CM6-1"   : (294, 362),
    "CNRM-CM6-1-HR": (1050, 1442),
    "CNRM-ESM2-1"  : (294, 362),
    "IPSL-CM6A-LR" : (332, 362),
    "GFDL-ESM4"    : (576, 720),
    "GISS-E2-2-G"  : (180, 288),  # has independent lon/lat
    "MIROC6"       : (256, 360),
    "MPI-ESM1-2-HR": (404, 802),
    "MPI-ESM1-2-LR": (220, 256),
    "MRI-ESM2-0"   : (363, 360),
    "NorESM2-LM"   : (385, 360),
    "NorESM2-MM"   : (385, 360),
    "UKESM1-0-LL"  : (330, 360),
    "UKESM1-1-LL"  : (330, 360)
}
    

# Ocean grid cell areas, absolute paths:
areacello_file_kw = {
    "AWI-CM-1-1-MR": ("piControl", "r1i1p1f1", "gn"),
    "CanESM5"      : ("piControl", "r1i1p2f1", "gn"),
    "CanESM5-CanOE": ("piControl", "r1i1p2f1", "gn"),
    "CESM2"        : ("piControl", "r1i1p1f1", "gn"),
    "CESM2-FV2"    : ("piControl", "r1i2p2f1", "gn"),
    "CESM2-WACCM"  : ("piControl", "r1i1p1f1", "gn"),
    "CNRM-CM6-1"   : ("piControl", "r1i1p1f2", "gn"),
    "CNRM-CM6-1-HR": ("piControl", "r1i1p1f2", "gn"),
    "CNRM-ESM2-1"  : ("piControl", "r1i1p1f2", "gn"),
    "GFDL-ESM4"    : ("piControl", "r1i1p1f1", "gn"),
    "GISS-E2-2-G"  : ("piControl", "r1i1p1f1", "gn"),
    "IPSL-CM6A-LR" : ("piControl", "r1i1p1f1", "gn"),
    "MIROC6"       : ("piControl", "r1i1p1f1", "gn"),
    "MPI-ESM1-2-HR": ("piControl", "r1i1p1f1", "gn"),
    "MPI-ESM1-2-LR": ("piControl", "r1i1p1f1", "gn"),
    "MRI-ESM2-0"   : ("piControl", "r1i1p1f1", "gn"),
    "NorESM2-LM"   : ("piControl", "r1i1p1f1", "gn"),
    "NorESM2-MM"   : ("piControl", "r1i1p1f1", "gn"),
    "UKESM1-0-LL"  : ("piControl", "r1i1p1f2", "gn"),
    "UKESM1-1-LL"  : ("piControl", "r1i1p1f2", "gn")
}

def areacello_file(model_id):
    return Path(dir_raw_nc_data, model_id,
        areacello_file_kw[model_id][0],
        f"areacello_Ofx_{model_id}_"
        + f"{areacello_file_kw[model_id][0]}_"
        + f"{areacello_file_kw[model_id][1]}_"
        + f"{areacello_file_kw[model_id][2]}.nc"
    )


# Grid coordinate names in raw netCDF data (lon, lat):
lonlat_ocn_nc_names = {
    "AWI-CM-1-1-MR": ("lon", "lat"),
    "CanESM5"      : ("longitude", "latitude"),
    "CanESM5-CanOE": ("longitude", "latitude"),
    "CESM2"        : ("lon", "lat"),
    "CESM2-FV2"    : ("lon", "lat"),
    "CESM2-WACCM"  : ("lon", "lat"),
    "CNRM-CM6-1"   : ("lon", "lat"),
    "CNRM-CM6-1-HR": ("lon", "lat"),
    "CNRM-ESM2-1"  : ("lon", "lat"),
    "GFDL-ESM4"    : ("lon", "lat"),
    "GISS-E2-2-G"  : ("lon", "lat"),
    "IPSL-CM6A-LR" : ("nav_lon", "nav_lat"),
    "MIROC6"       : ("longitude", "latitude"),
    "MPI-ESM1-2-HR": ("longitude", "latitude"),
    "MPI-ESM1-2-LR": ("longitude", "latitude"),
    "MRI-ESM2-0"   : ("longitude", "latitude"),
    "NorESM2-LM"   : ("longitude", "latitude"),
    "NorESM2-MM"   : ("longitude", "latitude"),
    "UKESM1-0-LL"  : ("longitude", "latitude"),
    "UKESM1-1-LL"  : ("longitude", "latitude")
}

lonlat_bnds_ocn_nc_names = {
    "AWI-CM-1-1-MR": ("lon_bnds", "lat_bnds"),
    "CanESM5"      : ("vertices_longitude","vertices_latitude"),
    "CanESM5-CanOE": ("vertices_longitude","vertices_latitude"),
    "CESM2"        : ("lon_bnds", "lat_bnds"),
    "CESM2-FV2"    : ("lon_bnds", "lat_bnds"),
    "CESM2-WACCM"  : ("lon_bnds", "lat_bnds"),
    "CNRM-CM6-1"   : ("bounds_lon", "bounds_lat"),
    "CNRM-CM6-1-HR": ("bounds_lon", "bounds_lat"),
    "CNRM-ESM2-1"  : ("bounds_lon", "bounds_lat"),
    "GFDL-ESM4"    : ("lon_bnds", "lat_bnds"),
    "GISS-E2-2-G"  : ("lon_bnds", "lat_bnds"),
    "IPSL-CM6A-LR" : ("bounds_nav_lon", "bounds_nav_lat"),
    "MPI-ESM1-2-HR": ("vertices_longitude", "vertices_latitude"),
    "MIROC6"       : ("vertices_longitude", "vertices_latitude"),
    "MPI-ESM1-2-LR": ("vertices_longitude", "vertices_latitude"),
    "MRI-ESM2-0"   : ("vertices_longitude", "vertices_latitude"),
    "NorESM2-LM"   : ("vertices_longitude", "vertices_latitude"),
    "NorESM2-MM"   : ("vertices_longitude", "vertices_latitude"),
    "UKESM1-0-LL"  : ("vertices_longitude", "vertices_latitude"),
    "UKESM1-1-LL"  : ("vertices_longitude", "vertices_latitude")
}

lonlat_bnds_ocn_file = areacello_file



# ========================================================== #
# SPECIFIC FIELD METADATA
# ========================================================== #

variable_domain = {
    "hfbasin"     : "ocn",
    "hfds"        : "ocn",
    "hfls"        : "atm",
    "hfss"        : "atm",
    "hfx"         : "ocn",
    "hfy"         : "ocn",
    "rlds"        : "atm",
    "rlus"        : "atm",
    "rlut"        : "atm",
    "rsds"        : "atm",
    "rsdt"        : "atm",
    "rsus"        : "atm",
    "rsut"        : "atm",
    "ocontemptend": "ocn",
    "opottemptend": "ocn",
    "siconc"      : "ocn",
    "siconca"     : "atm",
    "ta"          : "atm",
    "tas"         : "atm"
}


def degree_C_to_K(x):
    return x + 273.15

def K_to_degree_C(x):
    return x - 273.15


# Store the index of hfbasin that corresponds to global basin,
# and the name of the latitude coordinate axis:
# 
hfbasin_metadata = {
    "CanESM5"      : [2, "lat"],
    "CanESM5-CanOE": [2, "lat"],
    "GFDL-ESM4"    : [2, "y"],
    "GISS-E2-2-G"  : [2, "lat"],
    "IPSL-CM6A-LR" : [0, "nav_lat"],
    "MIROC6"       : [2, "lat"],
    "MPI-ESM1-2-LR": [0, "lat"],
    "MPI-ESM1-2-HR": [0, "lat"],
    "MRI-ESM2-0"   : [2, "lat"],
    "NorESM2-LM"   : [3, "lat"],
    "NorESM2-MM"   : [3, "lat"],
    "UKESM1-0-LL"  : [1, "lat"],
    "UKESM1-1-LL"  : [1, "lat"]
}


# "pot" for potential temperature, "con" for conservative:
# This is used to decide whether to load "opottemptend" or
# "ocontemptend" (3D tendency of sea water temperature
# expressed as heat content):
# 
ocn_prognostic_temperature = {
    "AWI-CM-1-1-MR"  : "pot",
    "CanESM5"        : "pot",
    "CanESM5-CanOE"  : "pot",
    "CESM2"          : "pot",
    "CESM2-FV2"      : "pot",
    "CESM2-WACCM"    : "pot",
    "IPSL-CM6A-LR"   : "con",
    "MIROC6"         : "pot",
    "MPI-ESM1-2-HR"  : "pot",
    "MPI-ESM1-2-LR"  : "pot",
    "MRI-ESM2-0"     : "pot"
}

# For the residual diagnostic of ocean temperature tendency
# (it doesn't matter because the quantity is heat content,
# not temperature, anyway):
ocn_prognostic_temperature_when_does_not_exist = "pot"


# ========================================================== #
# REANALYSES METADATA
# ========================================================== #

dir_raw_nc_data_reanalyses = Path("/storage", "basic", "cpom",
    "gb919150", "monthly_near_surface_air_temperature")

dir_out_nc_data_reanalyses = Path("/storage", "silver", "cpom",
                                  "gb919150", "phd",
                                  "reanalyses")

reanalysis_long_name = {
    "CFSR"   : "Climate Forecast System Reanalysis (CFSR)",
    "CFSv2"  : "Climate Forecast System Reanalysis Version "
               + "2 (CFSv2)",
    "ERA5"   : "European Centre for Medium-Range Weather "
               + "Forecasts (ECMWF) Reanalysis v5 (ERA5)",
    "JRA-55" : "Japanese 55-year Reanalysis (JRA-55)",
    "MERRA-2": "Modern-Era Retrospective analysis for "
               + "Research and Applications, Version 2 "
               + "(MERRA-2)"
}

defined_reanalyses = list(reanalysis_long_name.keys())
default_reanalysis = defined_reanalyses[0]


reanalysis_data_reference = {
    "CFSR"   : "Saha, S., and coauthors, 2010: NCEP Climate "
               + "Forecast System Reanalysis (CFSR) Monthly "
               + "Products, January 1979 to December 2010, "
               + "Research Data Archive at the National Center "
               + "for Atmospheric Research, Computational and "
               + "Information Systems Laboratory, "
               + "doi:10.5065/D6DN438J "
               + "[Accessed 23 February 2024]",
    "CFSv2"  : "Saha, S., and coauthors, 2012: NCEP Climate "
               + "Forecast System Version 2 (CFSv2) Monthly "
               + "Products, Research Data Archive at the "
               + "National Center for Atmospheric Research, "
               + "Computational and Information Systems "
               + "Laboratory. doi:10.5065/D69021ZF "
               + "[Accessed 23 February 2024]",
    "ERA5"   : "Hersbach, H., B. Bell, P. Berrisford, G. "
               + "Biavati, A. Horányi, J. Muñoz Sabater, and "
               + "coauthors, 2023: ERA5 monthly averaged data "
               + "on single levels from 1940 to present. "
               + "Copernicus Climate Change Service (C3S) "
               + "Climate Data Store (CDS), "
               + "doi:10.24381/cds.f17050d7 "
               + "[Accessed 23 February 2024]",
    "JRA-55" : "Japan Meteorological Agency/Japan, 2013: "
               + "JRA-55: Japanese 55-year Reanalysis, Monthly "
               + "Means and Variances, Research Data Archive "
               + "at the National Center for Atmospheric "
               + "Research, Computational and Information "
               + "Systems Laboratory, doi:10.5065/D60G3H5B "
               + "[Accessed 23 February 2024]",
    "MERRA-2": "Global Modeling and Assimilation Office "
               + "(GMAO), 2015: MERRA-2 instM_2d_asm_Nx: 2d,Mon"
               + "thly mean,Single-Level,Assimilation,Single-Le"
               + "vel Diagnostics V5.12.4, Greenbelt, MD, USA, "
               + "Goddard Earth Sciences Data and Information "
               + "Services Center (GES DISC), "
               + "doi:10.5067/5ESKGQTZG7FO "
               + "[Accessed 23 February 2024]"
}


reanalysis_reference = {
    "CFSR"   : "Saha, S. et al., 2010: The NCEP Climate "
               + "Forecast System Reanalysis, Bull. Am. "
               + "Meteorol. Soc., 91(8), 1015-1058, "
               + "doi:10.1175/2010BAMS3001.1",
    "CFSv2"  : "Saha, S. et al., 2014: The NCEP Climate "
               + "Forecast System Version 2, J. Clim., 27(6), "
               + "2185-2208, doi:10.1175/JCLI-D-12-00823.1",
    "ERA5"   : "Hersbach, H. et al., 2020: The ERA5 global "
               + "reanalysis, Q. J. R. Meteorol. Soc., "
               + "146(730), 1999-2049, doi:10.1002/qj.3803",
    "JRA-55" : "Kobayashi, S. et al., 2015: The JRA-55 "
               + "Reanalysis: General Specifications and Basic "
               + "Characteristics, J. Meteorol. Soc. Japan, "
               + "93(1), 5-48, doi:10.2151/jmsj.2015-001",
    "MERRA-2": "Gelaro, R. et al., 2017: The Modern-Era "
               + "Retrospective Analysis for Research and "
               + "Applications, Version 2 (MERRA-2), J. Clim., "
               + "30(14), 5419-5454, "
               + "doi:10.1175/JCLI-D-16-0758.1"
}

# Raw files must be put in the same format for the processing
# scripts under bespoke/atmospheric_reanalyses to work:
reanalysis_nc_file_fmt = {"tas": "t2.global.monthly.{}.nc"}

reanalysis_year_range = {
    "CFSR"      : (1979, 2010),
    "CFSv2"     : (2011, 2023),
    "ERA5"      : (1979, 2023),
    "JRA-55"    : (1958, 2023),
    "MERRA-2"   : (1980, 2023)
}

reanalysis_grid_dims_atm = {
    "CFSR"      : (576, 1152),
    "CFSv2"     : (880, 1760),
    "ERA5"      : (720, 1440),
    "JRA-55"    : (320, 640),
    "MERRA-2"   : (361, 576)
}

# For ERA5, these are the names after remapping with CDO (see
# shell script for why that is done); in the raw data they are
# "longitude" and "latitude":
reanalysis_nc_coord_names = {
    "CFSR"   : ["lon", "lat"],
    "CFSv2"  : ["lon", "lat"],
    "ERA5"   : ["lon", "lat"],
    "JRA-55" : ["g4_lon_2", "g4_lat_1"],
    "MERRA-2": ["lon", "lat"]
}

# Variable names are not standardised across reanalyses:
reanalysis_nc_names = {
    "tas": {
        "CFSR"   : "TMP_L103_Avg",
        "CFSv2"  : "TMP_L103_Avg",
        "ERA5"   : "t2m",
        "JRA-55" : "TMP_GDS4_HTGL_S113",
        "MERRA-2": "T2M"
    }
}


# ========================================================== #
# PASSIVE MICROWAVE SEA ICE CONCENTRATION METADATA
# ========================================================== #

# Two datasets are used, obtained from the National Snow and
# Ice Data Center (NSIDC), for passive microwave observations
# of sea ice concentration.

# Output data directory:
dir_out_nc_data_nsidc = Path("/storage", "silver", "cpom",
                             "gb919150", "phd",
                             "passive_microwave")
# No input directory needed as this data is processed from a 
# bash script where the paths are passed in explicitly.


# Dimensions of raw data in each hemisphere ("n"orthern and
# "s"outhern) on the NSIDC 25 km polar stereographic grid:
nsidc_grid_dims = {"n": (448, 304), "s": (332, 316)}

nsidc_lonlat_file = {
    "n": Path("/storage", "basic", "cpom", "gb919150", "NSIDC",
              "NSIDC-0771_polar_stereographic_ancilliary_v1",
              "NSIDC0771_LatLon_PS_N25km_v1.0.nc"),
    "s": Path("/storage", "basic", "cpom", "gb919150", "NSIDC",
              "NSIDC-0771_polar_stereographic_ancilliary_v1",
              "NSIDC0771_LatLon_PS_S25km_v1.0.nc")
}

nsidc_areacell_file = {
    "n": Path("/storage", "basic", "cpom", "gb919150", "NSIDC",
              "NSIDC-0771_polar_stereographic_ancilliary_v1",
              "NSIDC0771_CellArea_PS_N25km_v1.0.nc"),
    "s": Path("/storage", "basic", "cpom", "gb919150", "NSIDC",
              "NSIDC-0771_polar_stereographic_ancilliary_v1",
              "NSIDC0771_CellArea_PS_S25km_v1.0.nc")
}

# Valid ice mask for the northern hemisphere; one file per
# month. There is no equivalent valid ice mask for the southern
# hemisphere. This must be coded as a string so that it can
# be .format()'ed with the month 01, 02, ..., etc.:
nsidc_valid_ice_mask_n_nc_file_fmt = \
    str(Path("/storage", "basic", "cpom", "gb919150", "NSIDC",
        "NSIDC-0622_north_polar_stereographic_valid_ice_mask_v1",
        "NIC_valid_ice_mask.N25km.{:02}.1972-2007.nc"))

nsidc_nc_time_units = "days since 1978-01-01"

nsidc_source = {
    "NSIDC-0051": "Sea ice concentrations from Nimbus-7 SMMR "
                  + "and DMSP SSM/I-SSMIS passive microwave "
                  + "data, version 2",
    "NSIDC-0079": "Bootstrap sea ice concentrations from "
                  + "Nimbus-7 SMMR and DMSP SSM/I-SSMIS, "
                  + "version 4"
}

# Dataset references:
nsidc_date_accessed = "15 April 2024"
nsidc_data_reference = {
    "NSIDC-0051": "DiGirolamo, N., C. L. Parkinson, D. J. "
                  + "Cavalieri, P. Gloersen, and H. J. "
                  + "Zwally, 2022: Sea ice concentrations "
                  + "from Nimbus-7 SMMR and DMSP SSM/I-SSMIS "
                  + "passive microwave data, version 2 [data "
                  + "set], Boulder, Colorado USA, NASA "
                  + "National Snow and Ice Data Center "
                  + "Distributed Active Archive Center, "
                  + "doi:10.5067/MPYG15WAA4WX [Accessed "
                  + nsidc_date_accessed + "]",
    "NSIDC-0079": "Comiso, J. C., 2023: Bootstrap sea ice "
                  + "concentrations from Nimbus-7 SMMR and "
                  + "DMSP SSM/I-SSMIS, version 4 [data set], "
                  + "Boulder, Colorado USA, NASA National "
                  + "Snow and Ice Data Center Distributed "
                  + "Active Archive Center, "
                  + "doi:10.5067/X5LG68MH013O [Accessed "
                  + nsidc_date_accessed + "]",
    "NSIDC-0622": "Meier, W. N., J. Stroeve, F. Fetterer, M. "
                  + "Savoie, and H. Wilcox, 2015: Polar "
                  + "stereographic valid ice masks derived "
                  + "from National Ice Center monthly sea ice "
                  + "climatologies, version 1 [data set], "
                  + "Boulder, Colorado USA, NASA National Snow "
                  + "and Ice Data Center Distributed Active "
                  + "Archive Center, doi:10.5067/M4PUJAQRI2DS "
                  + f"[Accessed {nsidc_date_accessed}]",
    "NSIDC-0771": "Stewart, J. S., W. N. Meier, and D. J. "
                  + "Scott, 2022: Polar stereographic "
                  + "ancillary grid information, version 1 "
                  + "[data set], Boulder, Colorado USA, "
                  + "National Snow and Ice Data Center, "
                  + "doi:10.5067/N6INPBT8Y104 "
                  + f"[Accessed {nsidc_date_accessed}]"
}

defined_nsidc_datasets = list(nsidc_data_reference.keys())
