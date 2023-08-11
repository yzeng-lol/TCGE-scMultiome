workdir: config['work_dir']

## read in sample list
import os
import pandas as pd

## paths for pipeline and/or reference data
work_dir = config["work_dir"]
pipe_dir = config["pipe_dir"]
env_dir = os.getenv("CONDA_PREFIX")

######################################################
## read in sample and corresponding library file table
SAMPLES = (
    pd.read_csv(config["samples"], sep="\t")
    .set_index("sample_id", drop=False)
    .sort_index()
)


#############################################
## get taget outputs based on the config file
## either for individual samples or aggregate
## all samples listed in sample_aggr.tsv !!!
#############################################

def get_rule_all_input():
    ## ensure extra env installed
    extra_env = "extra_env/all_extra_env_installed",
    arc_out = expand("arc_count/{samples}/outs/atac_fragments.tsv.gz", samples = SAMPLES["sample_id"]),
    main_out = expand("main_seurat/{samples}.RDS", samples = SAMPLES["sample_id"]),
    qc_report = expand("main_seurat/{samples}_scMultiome_QC_and_Primary_Results_Report.html", samples = SAMPLES["sample_id"]),
    ## fixed outputs
    #meth_qc = "aggregated/meth_qc.txt",
    return  extra_env + arc_out + main_out + qc_report

############################
## other functions for input
#############################

###############################
##  get corresponding bwa_index
def get_cellranger_arc_ref():
    if config["pdx"]:
        #return reference for PDX data
        return config["pdx_ref"]
    else:
        return config["arc_ref"]

def get_library(wildcards):
    return SAMPLES.loc[wildcards.sample]["sample_library"]
