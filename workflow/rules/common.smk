workdir: config['work_dir']

## read in sample list
import pandas as pd

## read in sample and corresponding fq files talbe


## paths for pipeline and/or reference data
work_dir = config["work_dir"]
pipe_dir = config["pipe_dir"]


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
