
args = commandArgs(trailingOnly=TRUE)

## assigning commandArgs
sample_id   = args[1]
sample_rds  = args[2]
qc_rmd      = args[3]                       ## renamed qc_report.rmd for each sample
pipe_dir    = args[4]                       ## for DoMultiBarHeatmap

suppressMessages(library(rmarkdown))        ## for HTML QC report


##########################
## Generate HTML QC report
##########################
{
render(qc_rmd, output_dir = "main_seurat",
       params = list(readin = sample_rds, sample_id = sample_id, pipe_dir = pipe_dir))
}


## under testing
if(FALSE){
# cp /cluster/home/yzeng/snakemake/iSHARC/workflow/scripts/qc_report.Rmd /cluster/projects/tcge/scMultiome/iSHARC_test/main_seurat/Lung_scMultiome_QC_Report.Rmd
# Rscript --vanilla /cluster/home/yzeng/snakemake/iSHARC/workflow/scripts/qc_report.R Lung /cluster/projects/tcge/scMultiome/iSHARC_test/main_seurat/Lung.RDS /cluster/projects/tcge/scMultiome/iSHARC_test/main_seurat/Lung_scMultiome_QC_Report.Rmd

## testing locally
rm (list = ls())
setwd("/Users/yong/OneDrive_UHN/Projects/snakemake/iSHARC_test/main_seurat")

library(rmarkdown)

#### testing PDAC_PDA_87784
render("/Users/yong/OneDrive_UHN/Projects/snakemake/iSHARC/workflow/scripts/qc_and_results_report.rmd",
       output_dir = "./",  params = list(readin = "/Users/yong/OneDrive_UHN/Projects/snakemake/iSHARC_test/main_seurat/PDAC_PDA_87784.RDS",
                                         sample_id = "PDAC_PDA_87784",
                                         pipe_dir = "/Users/yong/OneDrive - UHN/Projects/snakemake/iSHARC"))

## testing sample Lung
#render("/Users/yong/OneDrive_UHN/Projects/snakemake/iSHARC/workflow/scripts/qc_and_results_report.rmd",
#       output_dir = "main_seurat",  params = list(readin = "/Users/yong/OneDrive_UHN/Projects/snakemake/iSHARC_test/main_seurat/Lung.RDS", sample_id = "Lung"))

## testing sample pbmc
#render("/Users/yong/OneDrive_UHN/Projects/snakemake/iSHARC/workflow/scripts/qc_and_results_report.rmd",
#       output_dir = "main_seurat",  params = list(readin = "/Users/yong/OneDrive - UHN/Projects/TCGE/scMultimoics/pbmc_granulocyte_sorted_3k/scMultiome/scMultiome.RDS", sample_id = "PBMC"))


}
