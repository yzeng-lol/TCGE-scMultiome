
args = commandArgs(trailingOnly=TRUE)

## assigning commandArgs
qc_rmd      = args[1]                       ## renamed qc_report.rmd for each sample
project_id  = args[2]
pipe_dir    = args[3]                       ## for DoMultiBarHeatmap
integr_dir  = args[4]

suppressMessages(library(rmarkdown))        ## for HTML QC report


##########################
## Generate HTML QC report
##########################
{
render(qc_rmd, output_dir = "./integration",
       params = list(project_id = project_id, pipe_dir = pipe_dir, integr_dir = integr_dir))
}


## under testing
if(FALSE){
# cp /cluster/home/yzeng/snakemake/iSHARC/workflow/scripts/qc_report.Rmd /cluster/projects/tcge/scMultiome/iSHARC_test/main_seurat/Lung_scMultiome_QC_Report.Rmd
# Rscript --vanilla /cluster/home/yzeng/snakemake/iSHARC/workflow/scripts/qc_report.R Lung /cluster/projects/tcge/scMultiome/iSHARC_test/main_seurat/Lung.RDS /cluster/projects/tcge/scMultiome/iSHARC_test/main_seurat/Lung_scMultiome_QC_Report.Rmd


### testing on H4H

# cd  /cluster/projects/tcge/scMultiome/processed_data/TCGE-scMOS-PTLD
# conda activate  /cluster/home/yzeng/miniconda3/envs/iSHARC_extra_env/ebcd8f410d62ed43e71d1cebab0b3296_    ## new
# R

library(rmarkdown)

#### testing PDAC_PDA_87784
render("/cluster/home/yzeng/snakemake/iSHARC/workflow/scripts/integration_report.Rmd",
       output_dir = "./integration",  params = list(project_id = "RLN",
                                         pipe_dir = "/cluster/home/yzeng/snakemake/iSHARC",
                                         integr_dir = "/cluster/projects/tcge/scMultiome/processed_data/TCGE-scMOS-PTLD/integration"))



## testing locally
rm (list = ls())
setwd("/Users/yong/OneDrive - UHN/Projects/snakemake/iSHARC_test/PRCA")

library(rmarkdown)

#### testing PDAC_PDA_87784
render("/Users/yong/OneDrive - UHN/Projects/snakemake/iSHARC_test/PRCA/integration_report.rmd",
       output_dir = "./",  params = list(readin = "/Users/yong/OneDrive - UHN/Projects/snakemake/iSHARC_test/PRCA/integration/RNA_ATAC_integrated_by_WNN.RDS",
                                         project_id = "PRCA",
                                         pipe_dir = "/Users/yong/OneDrive - UHN/Projects/snakemake/iSHARC",
                                         integr_dir = "/Users/yong/OneDrive - UHN/Projects/snakemake/iSHARC_test/PRCA/integration"))


## testing sample Lung
#render("/Users/yong/OneDrive_UHN/Projects/snakemake/iSHARC/workflow/scripts/qc_and_results_report.rmd",
#       output_dir = "main_seurat",  params = list(readin = "/Users/yong/OneDrive_UHN/Projects/snakemake/iSHARC_test/main_seurat/Lung.RDS", sample_id = "Lung"))

## testing sample pbmc
#render("/Users/yong/OneDrive_UHN/Projects/snakemake/iSHARC/workflow/scripts/qc_and_results_report.rmd",
#       output_dir = "main_seurat",  params = list(readin = "/Users/yong/OneDrive - UHN/Projects/TCGE/scMultimoics/pbmc_granulocyte_sorted_3k/scMultiome/scMultiome.RDS", sample_id = "PBMC"))


}
