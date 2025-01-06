################################################################################
## Generate the HTML report for the QC and primary results per each individually
## sample
##
## Contact : Yong Zeng <yong.zeng@uhn.ca>
##
## Copyright: ...
################################################################################


##############################
### parse and assign arguments
##############################
{
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

## adding parameters
## by default ArgumentParser will add an help option
## run "Rscript main_for_individual_sample.R -h" for help info
parser$add_argument("-s", "--sample_id", required=TRUE,
                    help = "Unique sample ID")
parser$add_argument("-exo", "--extended_analyses_seurat_object", required=TRUE,
                    help = "the seurat object with extended analyses")
parser$add_argument("-pipe", "--pipe_dir", required=TRUE,
                    help = "The PATH to iSHARC pipeline, which local dependences included")
parser$add_argument("-rmd", "--report_rmd_file", required=TRUE,
                     help = "The RMD file for QC and primary results HTML report")

## assigning passing arguments
args <- parser$parse_args()
print(args)

sample_id <- args$sample_id
exo_file <- args$extended_analyses_seurat_object
pipe_dir <- args$pipe_dir  ## for DoMultiBarHeatmap
rmd_file <- args$report_rmd_file

## output dir
out_dir <- paste0("individual_samples/", sample_id)
}

##########################
## Generate HTML QC report
##########################
{
suppressMessages(library(rmarkdown))        ## for HTML QC report

render(rmd_file, output_dir = out_dir,
       params = list(readin = exo_file, sample_id = sample_id, pipe_dir = pipe_dir))
}


################################
## testing
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


}
