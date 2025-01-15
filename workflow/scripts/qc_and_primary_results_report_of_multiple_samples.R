################################################################################
## Generate the HTML report for QC and primary results of integrated
## multiple samples.
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
# parser$add_argument("-p", "--project_id", required=TRUE,
#                    help = "Unique project ID for multiple samples")
parser$add_argument("-dir", "--integration_dir", required=TRUE,
                    help = "the dirctory for the integrated outputs")
parser$add_argument("-pipe", "--pipe_dir", required=TRUE,
                    help = "The PATH to iSHARC pipeline, which local dependences included")
parser$add_argument("-rmd", "--report_rmd_file", required=TRUE,
                     help = "The RMD file for QC and primary results for multiple HTML report")

## assigning passing arguments
args <- parser$parse_args()
print(args)

#project_id <- args$project_id
integr_dir <- args$integration_dir          ## ensure to be the full path, without forward slash at the end
pipe_dir <- args$pipe_dir                   ## for DoMultiBarHeatmap
rmd_file <- args$report_rmd_file

suppressMessages(library(rmarkdown))        ## for HTML QC report

}

##########################
## Generate HTML QC report
##########################
{
render(rmd_file, output_dir = "./integrated_samples",
       params = list(integr_dir = integr_dir, pipe_dir = pipe_dir))
}
