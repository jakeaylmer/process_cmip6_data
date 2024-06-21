"""Lists of models available for each experiment or combination
of historical and future experiment."""

available_experiments = ["piControl", "historical",
                         "ssp370", "ssp585"]

# All available models (keys), and which experiments they
# provide (index matching available_experiments):
experiments_available_per_model = {
    "AWI-CM-1-1-MR": [True,  True,  True,  True ],
    "CESM2"        : [True,  True,  True,  True ],
    "CESM2-FV2"    : [True,  True,  True,  True ],
    "CESM2-WACCM"  : [True,  True,  True,  True ],
    "CNRM-CM6-1"   : [False, True,  True,  True ],
    "CNRM-CM6-1-HR": [True,  True,  True,  True ],
    "CNRM-ESM2-1"  : [False, True,  True,  True ],
    "CanESM5"      : [True,  True,  True,  True ],
    "CanESM5-CanOE": [True,  True,  True,  True ],
    "GFDL-ESM4"    : [False,  True,  True,  True ],
    "GISS-E2-2-G"  : [True,  True,  True,  True ],
    "IPSL-CM6A-LR" : [True,  True,  True,  True ],
    "MIROC6"       : [True,  True,  True,  True ],
    "MPI-ESM1-2-HR": [True,  True,  True,  True ],
    "MPI-ESM1-2-LR": [True,  True,  True,  True ],
    "MRI-ESM2-0"   : [True,  True,  True,  True ],
    "NorESM2-LM"   : [True,  True,  True,  True ],
    "NorESM2-MM"   : [True,  True,  True,  True ],
    "UKESM1-0-LL"  : [True,  True,  True,  True ],
    "UKESM1-1-LL"  : [True,  True,  True,  False]
}

available_models = sorted(list(
    experiments_available_per_model.keys()))

by_experiment = {}
for x in range(len(available_experiments)):
    by_experiment[available_experiments[x]] = \
        [k for k in available_models
         if experiments_available_per_model[k][x]]

# Add joint experiments:
by_experiment["historical+ssp370"] = \
    [k for k in available_models if (
        experiments_available_per_model[k][1] and
        experiments_available_per_model[k][2])]

by_experiment["historical+ssp585"] = \
    [k for k in available_models if (
        experiments_available_per_model[k][1] and
        experiments_available_per_model[k][3])]
