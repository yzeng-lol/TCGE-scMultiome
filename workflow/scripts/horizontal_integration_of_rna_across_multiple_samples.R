################################################################################
##  horizontal_integration_of_rna_across_multiple_sample.R is a function for
##  horizontal integration of RNA across multiple sample by:
##      * simply merging
##      * integrating with Harmony
##      * integrating with Seurat anchor strategy
##
## Contact : Yong Zeng <yong.zeng@uhn.ca>
################################################################################


#######################################
### parse and assign arguments
### reading sample list for aggregation
#######################################
{
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

## adding parameters
## by default ArgumentParser will add an help option
parser$add_argument("-si", "--samples_integration", required=TRUE,
                    help = "list of sample IDs to be aggregated in TSV format")

parser$add_argument("-knn_k", "--knn_k_param", type = "integer", default = 20,
                    help = "k for the k-nearest neighbor algorithm")
parser$add_argument("-dims_n", "--dimentions_n", type = "integer", default = 50,
                    help = "number of reduced dimentions (e.g., PCs) for functions: RunUMAP, FindNeighbors, FindMultiModalNeighbors")
parser$add_argument("-comm_res", "--community_resolution", type = "double", default = 0.8,
                    help = "Value above (below) 1.0 if you want to obtain a larger (smaller) number of communities")

parser$add_argument("-t", "--threads", type = "integer", default = 12,
                     help = "Number of cores for the parallelization")
parser$add_argument("-fgm", "--future_globals_maxSize", type = "integer", default = 12,
                    help = "Maximum memory in GB for the future parallelization global variables")

## assigning passing arguments
args <- parser$parse_args()
print(args)

knn_k <- args$knn_k_param
dims_n <- args$dimentions_n
comm_res <- args$community_resolution


## output dir
out_dir <- paste0(getwd(), "/integrated_samples/rna/") ## with forward slash at the end

## read in aggregation list
sample_list <- read.table(args$samples_integration , header = T)
sample_ids <- sample_list$sample_id

}

##################################################
## loading required packages
##################################################
{
  suppressMessages(library(Seurat))
  suppressMessages(library(Signac))
  suppressMessages(library(dplyr))
  suppressMessages(library(ggplot2))
  suppressMessages(library(harmony))

  ## enable the Parallelization with the future packages
  suppressMessages(library(future))
  plan("multicore", workers = args$threads)
  options(future.globals.maxSize = args$future_globals_maxSize * 1024^3)

}


##########################################
## read in seurat objects per sample list
###########################################
{
seurat_object_list <- list()

for(i in 1: length(sample_ids))
{
  seurat_object_list[[i]] <- readRDS(paste0(getwd(), "/individual_samples/", sample_ids[i], "/", sample_ids[i], "_extended_seurat_object.RDS"))

  DefaultAssay(seurat_object_list[[i]]) <- 'RNA'

  ## adding sample id to meta table
  seurat_object_list[[i]][["sample_id"]] <- sample_ids[i]

  ## to remove assays not used to save space
  ## seurat_object_list[[i]][["atac_arc"]] <- NULL
  seurat_object_list[[i]][["ATAC"]] <- NULL
}

}

#################################################
## integrate snRNA-seq data from multiple samples
#################################################
{

  #############################
  ## merge without integration
  {
    ## "RNA" has been set as default assay
    ## merge raw data per seurat project
    rna_merged <- merge(x = seurat_object_list[[1]],
                        y = seurat_object_list[2:length(seurat_object_list)],
                        #add.cell.ids = sample_ids,  ## will problematic for atac integration
                        merge.data = TRUE           ## Merge the data slots as well
                        )

    ## RNA_merged differ from existing SCT
    rna_merged <- SCTransform(rna_merged, assay = 'RNA', new.assay.name = 'RNA_integrated')
    rna_merged <- RunPCA(rna_merged, npcs = dims_n, verbose = FALSE)
    rna_merged <- RunUMAP(rna_merged, dims = 1:dims_n, verbose = FALSE)
    rna_merged <- FindNeighbors(rna_merged, dims = 1:dims_n, reduction = "pca", k.param = knn_k) %>% FindClusters(resolution = comm_res)

    ## generalized the resolution
    rna_merged$RNA_integrated_clusters <- rna_merged[[paste0("RNA_integrated_snn_res.", comm_res)]]

    saveRDS(rna_merged, paste0(out_dir, "RNA_integrated_by_merging.RDS"))
  }

  ######################################
  ## merge and integrated with Harmoney
  ## sample id as the group.by.vars  by default

  ##   workflow recommended by Harmony developers :
  ##      1) Merge the raw Seurat objects for all samples to integrate;
  ##      2) then perform normalization, variable feature selection and PC calculation on this merged object
  {
    rna_harmonized <- RunHarmony(rna_merged,
                                 group.by.vars = c("sample_id"),
                                 reduction = "pca",
                                 assay.use = "RNA_integrated",
                                 reduction.save = "harmony")

    rna_harmonized <- RunUMAP(rna_harmonized, reduction = "harmony", dims = 1:dims_n, verbose = FALSE)
    rna_harmonized <- FindNeighbors(rna_harmonized, reduction = "harmony", dims = 1:dims_n, k.param = knn_k) %>%
                      FindClusters(resolution = comm_res)        ## clusters list in  "RNA_integrated_snn_res.0.8", which is differ from rna_merged$RNA_integrated_snn_res.0.8

    ## generalized the resolution
    rna_harmonized$RNA_integrated_clusters <- rna_harmonized[[paste0("RNA_integrated_snn_res.", comm_res)]]

    saveRDS(rna_harmonized, paste0(out_dir, "RNA_integrated_by_harmony.RDS"))
  }


  #########################
  ## integrated with Seurat
  ## The Reciprocal PCA (RPCA) and  Reference-based integration strategies were
  ## applied to improve efficiency and run times
  ## https://satijalab.org/seurat/archive/v4.3/integration_large_datasets
  {
    for(i in 1: length(sample_ids))
    {
    ## using the normalized assay per each Seurat object
    ## the features of SCT@data is less than RNA@data
    ## https://github.com/satijalab/seurat/issues/4082

    DefaultAssay(seurat_object_list[[i]]) <- 'SCT'
    }

    ## features selected for integration
    ## SCTransform-normalized datasets (by default 2000 features will be selected)
    rna_features <- SelectIntegrationFeatures(object.list = seurat_object_list,
                                              assay = rep("SCT", length(seurat_object_list)))

    ## updating pca reduction using the rna_features across the object
    seurat_object_list <- lapply(X = seurat_object_list,
                                FUN = function(x) {
                                  x <- ScaleData(x, features = rna_features, verbose = FALSE)
                                  x <- RunPCA(x, features = rna_features, verbose = FALSE)}
                                )

    ## Find anchors, and using the first sample as reference, and "rpca" for reduction
    ## reference: A vector specifying the object/s to be used as a reference during integration.
    ## reference: If NULL (default), all pairwise anchors are found (no reference/s).
    ## Could using the first object as reference to save computational time, but might affect the outcomes as well
    rna_anchors <- FindIntegrationAnchors(object.list = seurat_object_list,
                                          #reference = c(1),
                                          reduction = "rpca",
                                          dims = 1:dims_n)

    ## Adding "SCT_Integrated" Assay
    rna_anchored <- IntegrateData(anchorset = rna_anchors, new.assay.name = "RNA_integrated", dims = 1:dims_n)
    rna_anchored <- ScaleData(rna_anchored, verbose = FALSE)
    rna_anchored <- RunPCA(rna_anchored, npcs = dims_n, verbose = FALSE)
    rna_anchored <- RunUMAP(rna_anchored, dims = 1:dims_n, verbose = FALSE)
    rna_anchored <- FindNeighbors(rna_anchored, dims = 1:dims_n, reduction = "pca", k.param = knn_k) %>% FindClusters(resolution = comm_res)

    ## generalized the resolution
    rna_anchored$RNA_integrated_clusters <- rna_anchored[[paste0("RNA_integrated_snn_res.", comm_res)]]

    saveRDS(rna_anchored, paste0(out_dir, "RNA_integrated_by_anchors.RDS"))
  }


  ####################################################
  ##  UMAPs for merged, harmonized and integrated RNA
  ####################################################
  {
  ## UMAP plots group by sample_id
  p1 <- DimPlot(rna_merged,  reduction = "umap", group.by = "sample_id",
                label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("RNA_merged") + NoLegend()

  p2 <- DimPlot(rna_harmonized,  reduction = "umap", group.by = "sample_id",
                label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("RNA_harmonized") + NoLegend()

  p3 <- DimPlot(rna_anchored,  reduction = "umap", group.by = "sample_id",
                label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("RNA_anchored")

  g <- p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(out_dir, "Merged_Harmonized_Anchored_RNA_UMAPs_labeled_by_sample.pdf"), width = 12, height = 4)


  ## UMAP plots group by clusters (different from assay to assay)
  p1 <- DimPlot(rna_merged,  reduction = "umap", group.by = "RNA_integrated_clusters",
                label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("RNA_merged")

  p2 <- DimPlot(rna_harmonized,  reduction = "umap", group.by = "RNA_integrated_clusters",
                label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("RNA_harmonized")

  p3 <- DimPlot(rna_anchored,  reduction = "umap", group.by = "RNA_integrated_clusters",
                label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("RNA_anchored")

  g <- p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(out_dir, "Merged_Harmonized_Anchored_RNA_UMAPs_labeled_by_cluster.pdf"), width = 15, height = 4)
  }

  print("The horizontal integration of RNA across multiple samples has been successfully completed!!")

}
