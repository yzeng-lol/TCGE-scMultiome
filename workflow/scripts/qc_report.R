
args = commandArgs(trailingOnly=TRUE)

## assigning commandArgs
sample_id   = args[1]
sample_rds  = args[2]
qc_rmd      = args[3]                       ## renamed qc_report.rmd for each sample

suppressMessages(library(rmarkdown))        ## for HTML QC report


##########################
## Generate HTML QC report
##########################
{
render(qc_rmd, output_dir = "main_seurat",
       params = list(readin = sample_rds, sample_id = sample_id))
}


## under testing
if(FALSE){
# cp /cluster/home/yzeng/snakemake/iSHARC/workflow/scripts/qc_report.Rmd /cluster/projects/tcge/scMultiome/iSHARC_test/main_seurat/Lung_scMultiome_QC_Report.Rmd
# Rscript --vanilla /cluster/home/yzeng/snakemake/iSHARC/workflow/scripts/qc_report.R Lung /cluster/projects/tcge/scMultiome/iSHARC_test/main_seurat/Lung.RDS /cluster/projects/tcge/scMultiome/iSHARC_test/main_seurat/Lung_scMultiome_QC_Report.Rmd
}
