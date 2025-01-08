################################################################################
## Generate the HTML report for QC and primary results of each individual sample
##
## Contact : Yong Zeng <yong.zeng@uhn.ca>
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
