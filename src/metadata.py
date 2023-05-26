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
# (this is the most painful part!):
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
