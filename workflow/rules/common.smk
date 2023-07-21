workdir: config['workdir']

## read in sample list
import pandas as pd

## read in sample and corresponding fq files talbe
SAMPLES = (
    pd.read_csv(config["samples"], sep="\t")
    .set_index("sample_id", drop=False)
    .sort_index()
)


## read in samples for aggregation
if config["aggreate"]:
    SAMPLES_AGGR = (
        pd.read_csv(config["samples_aggr"], sep="\t")
        .set_index("sample_id", drop=False)
        .sort_index()
    )
else:
    SAMPLES_AGGR = SAMPLES     ## must be defined

## paths for pipeline and/or reference data
wd = config["workdir"]
pipe_dir = config["pipeline_dir"]
env_dir = config["pipeline_env"]  #${CONDA_PREFIX} dosen't work
#umi_list = config["umi_list"]


## read in refrence files' info
REF = pd.read_csv(config["ref_files"], sep="\t", header = None, index_col = 0)
blacklist = REF.loc["blacklist"][1]   ## ENCODE blacklist


#############################################
## get taget outputs based on the config file
## either for individual samples or aggregate
## all samples listed in sample_aggr.tsv !!!
#############################################

def get_rule_all_input():
    ## ensure extra env installed
    extra_env = "extra_env/all_extra_env_installed",

    ## fixed outputs
    #meth_qc = "aggregated/meth_qc.txt",
    return  extra_env

############################
## other functions for input
#############################

###############################
##  get corresponding bwa_index
def get_cellranger_arc_ref():
    if config["PDX"]:
        #return reference for PDX data
        return config["pdx_idx"]
    else:
        return REF.loc["cellranger_arc_ref"][1]
