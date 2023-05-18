"""Main settings and model metadata, to be filled in manually.
"""

from pathlib import Path
import numpy as np

# =========================================================== #

default_model_id             = "CESM2-FV2"
default_member_id            = "r1i2p2f1"
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

defined_experiments = ["piControl", "historical", "ssp126",
                       "ssp245", "ssp370", "ssp585"]

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
        "CanESM5"      : (1, 1051),
        "CanESM5-CanOE": (1, 501),
        "CESM2"        : (1, 1200),
        "CESM2-FV2"    : (1, 150),
        "CESM2-WACCM"  : (1, 499),
        "CNRM-CM6-1"   : (1, 500),
        "CNRM-CM6-1-HR": (1, 300),
        "CNRM-ESM2-1"  : (1, 500),
        "IPSL-CM6A-LR" : (1, 1000),
        "MPI-ESM1-2-HR": (1, 500),
        "MPI-ESM1-2-LR": (1, 1000),
        "MRI-ESM2-0"   : (1, 701),
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
for exp in ["ssp126", "ssp245", "ssp370", "ssp585"]:
    year_range[exp] = dict.fromkeys(defined_models,
                                    (2015, 2100))

# =========================================================== #
# ENSEMBLE MEMBERS
# ----------------
# 
# List of ensemble members per model, per experiment
# (this is the most painful part!):
# =========================================================== #

members = {
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
    "CESM2-FV2": {
        "piControl" : ["r1i2p2f1"],
        "historical": ["r1i2p2f1"],
        "ssp370"    : ["r1i2p2f1"],
        "ssp585"    : ["r1i2p2f1"]
    },
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
    "CNRM-CM6-1-HR": {
        "piControl": ["r1i1p1f2"],
        "historical": ["r1i1p1f2"],
        "ssp370": ["r1i1p1f2"],
        "ssp585": ["r1i1p1f2"]
    },
    "CNRM-ESM2-1" : {
        "piControl": ["r1i1p1f2"],
        "historical": ["r1i1p1f2"],
        "ssp370": ["r1i1p1f2"],
        "ssp585": ["r1i1p1f2"]
    },
    "IPSL-CM6A-LR": {
        "piControl" : ["r1i1p1f1"],
        "historical": [f"r{x}i1p1f1" for x in range(1, 33)],
        "ssp370"    : [f"r{x}i1p1f1" for x in range(1,11)]
                    + ["r14i1p1f1"],
        "ssp585"    : [f"r{x}i1p1f1" for x in range(1,5)]
                    + ["r6i1p1f1", "r14i1p1f1"]
    },
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
        "ssp585"    : ["r1i1p1f1"]
    },
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
# GRID PROPERTIES
# ---------------
# Dimensions, coordinate netCDF names, location of cell area
# and vertices, etc.
# 
# =========================================================== #

# Atmosphere grid dimensions (ny_atm, nx_atm):
grid_dims_atm = {
    "CanESM5"      : (),
    "CanESM5-CanOE": (64, 128),
    "CESM2"        : (),
    "CESM2-FV2"    : (96, 144),
    "CESM2-WACCM"  : (),
    "CNRM-CM6-1"   : (),
    "CNRM-CM6-1-HR": (),
    "CNRM-ESM2-1"  : (),
    "IPSL-CM6A-LR" : (143, 144),
    "MPI-ESM1-2-HR": (),
    "MPI-ESM1-2-LR": (),
    "MRI-ESM2-0"   : (),
    "UKESM1-0-LL"  : ()
}


# Atmosphere grid cell areas, absolute paths:
areacella_file_kw = {
    "CanESM5-CanOE": ("piControl", "r1i1p2f1", "gn"),
    "CESM2-FV2"    : ("piControl", "r1i1p1f1", "gn"),
    "IPSL-CM6A-LR" : ("piControl", "r1i1p1f1", "gr")
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
    "CanESM5"      : (),
    "CanESM5-CanOE": ("lon", "lat"),
    "CESM2"        : (),
    "CESM2-FV2"    : ("lon", "lat"),
    "CESM2-WACCM"  : (),
    "CNRM-CM6-1"   : (),
    "CNRM-CM6-1-HR": (),
    "CNRM-ESM2-1"  : (),
    "IPSL-CM6A-LR" : ("lon", "lat"),
    "MPI-ESM1-2-HR": (),
    "MPI-ESM1-2-LR": (),
    "MRI-ESM2-0"   : (),
    "UKESM1-0-LL"  : ()
}

lonlat_bnds_atm_nc_names = {
    "CanESM5-CanOE": ("lon_bnds", "lat_bnds"),
    "CESM2-FV2"    : ("lon_bnds", "lat_bnds"),
    "IPSL-CM6A-LR" : None
}

def lonlat_bnds_atm_file(model_id):
    if lonlat_bnds_atm_nc_names[model_id] is None:
        raise FileNotFoundError
    else:
        return areacella_file(model_id)



# Ocean grid dimensions (ny_ocn, nx_ocn):
grid_dims_ocn = {
    "CanESM5"      : (291, 360),
    "CanESM5-CanOE": (291, 360),
    "CESM2"        : (384, 320),
    "CESM2-FV2"    : (384, 320),
    "CESM2-WACCM"  : (384, 320),
    "CNRM-CM6-1"   : (294, 362),
    "CNRM-CM6-1-HR": (1050, 1442),
    "CNRM-ESM2-1"  : (294, 362),
    "IPSL-CM6A-LR" : (332, 362),
    "MPI-ESM1-2-HR": (404, 802),
    "MPI-ESM1-2-LR": (220, 256),
    "MRI-ESM2-0"   : (363, 360),
    "UKESM1-0-LL"  : (330, 360)
}
    

# Ocean grid cell areas, absolute paths:
areacello_file_kw = {
    "CanESM5"      : ("piControl", "r1i1p2f1", "gn"),
    "CanESM5-CanOE": ("piControl", "r1i1p2f1", "gn"),
    "CESM2"        : ("piControl", "r1i1p1f1", "gn"),
    "CESM2-FV2"    : ("piControl", "r1i2p2f1", "gn"),
    "CESM2-WACCM"  : ("piControl", "r1i1p1f1", "gn"),
    "CNRM-CM6-1"   : ("piControl", "r1i1p1f2", "gn"),
    "CNRM-CM6-1-HR": ("piControl", "r1i1p1f2", "gn"),
    "CNRM-ESM2-1"  : ("piControl", "r1i1p1f2", "gn"),
    "IPSL-CM6A-LR" : ("piControl", "r1i1p1f1", "gn"),
    "MPI-ESM1-2-HR": ("piControl", "r1i1p1f1", "gn"),
    "MPI-ESM1-2-LR": ("piControl", "r1i1p1f1", "gn"),
    "MRI-ESM2-0"   : ("piControl", "r1i1p1f1", "gn"),
    "UKESM1-0-LL"  : ("piControl", "r1i1p1f2", "gn")
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
    "CanESM5"      : ("longitude", "latitude"),
    "CanESM5-CanOE": ("longitude", "latitude"),
    "CESM2"        : ("lon", "lat"),
    "CESM2-FV2"    : ("lon", "lat"),
    "CESM2-WACCM"  : ("lon", "lat"),
    "CNRM-CM6-1"   : ("lon", "lat"),
    "CNRM-CM6-1-HR": ("lon", "lat"),
    "CNRM-ESM2-1"  : ("lon", "lat"),
    "IPSL-CM6A-LR" : ("nav_lon", "nav_lat"),
    "MPI-ESM1-2-HR": ("longitude", "latitude"),
    "MPI-ESM1-2-LR": ("longitude", "latitude"),
    "MRI-ESM2-0"   : ("longitude", "latitude"),
    "UKESM1-0-LL"  : ("longitude", "latitude")
}

lonlat_bnds_ocn_nc_names = {
    "CanESM5-CanOE": ("vertices_longitude","vertices_latitude"),
    "CESM2-FV2"    : ("lon_bnds", "lat_bnds"),
    "CNRM-ESM2-1"  : ("bounds_lon", "bounds_lat"),
    "IPSL-CM6A-LR" : ("bounds_nav_lon", "bounds_nav_lat")
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
    "IPSL-CM6A-LR" : [0, "nav_lat"],
    "MPI-ESM1-2-LR": [0, "lat"],
    "MPI-ESM1-2-HR": [0, "lat"],
    "MRI-ESM2-0"   : [2, "lat"],
    "UKESM1-1-LL"  : [1, "lat"]
}


# "pot" for potential temperature, "con" for conservative:
# This is used to decide whether to load "opottemptend" or
# "ocontemptend" (3D tendency of sea water temperature
# expressed as heat content):
# 
ocn_prognostic_temperature = {
    "CanESM5"        : "pot",
    "CanESM5-CanOE"  : "pot",
    "CESM2"          : "pot",
    "CESM2-FV2"      : "pot",
    "CESM2-WACCM"    : "pot",
    "CESM2-WACCM-FV2": "pot",
    "CNRM-CM6-1"     : "pot",
    "CNRM-CM6-1-HR"  : "pot",
    "CNRM-ESM2-1"    : "pot",
    "IPSL-CM6A-LR"   : "con",
    "MPI-ESM1-2-HR"  : "pot",
    "MPI-ESM1-2-LR"  : "pot",
    "MRI-ESM2-0"     : "pot",
    "UKESM1-0-LL"    : "pot"
}

# For the residual diagnostic of ocean temperature tendency
# (it doesn't matter because the quantity is heat content,
# not temperature, anyway):
ocn_prognostic_temperature_when_does_not_exist = "pot"