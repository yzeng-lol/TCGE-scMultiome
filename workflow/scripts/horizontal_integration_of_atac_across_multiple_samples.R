################################################################################
##  horizontal_integration_of_atac_across_multiple_sample.R is a function for
##  horizontal integration of ATAC across multiple sample by:
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
parser$add_argument("-pipe", "--pipe_dir", required=TRUE,
                    help = "The PATH to iSHARC pipeline, which local dependences included")

## assigning passing arguments
args <- parser$parse_args()
print(args)

## output dir
out_dir <- paste0(getwd(), "/integrated_samples/atac/") ## with forward slash at the end

## read in aggregation list
sample_list <- read.table(args$samples_integration , header = T)
sample_ids <- sample_list$sample_id

pipe_dir <- args$pipe_dir

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
  suppressMessages(library(IRanges))
}


##########################################
## read in seurat objects per sample list
###########################################
{
seurat_object_list <- list()
for(i in 1: length(sample_ids))    # for testing
{
  seurat_object_list[[i]] <- readRDS(paste0(getwd(), "/individual_samples/", sample_ids[i], "/", sample_ids[i], "_extended_seurat_object.RDS"))

  DefaultAssay(seurat_object_list[[i]]) <- 'ATAC'

  ## adding sample id to meta talbe
  seurat_object_list[[i]][["sample_id"]] <- sample_ids[i]

  ## to remove assays not used to save space
  ## seurat_object_list[[i]][["atac_arc"]] <- NULL
  seurat_object_list[[i]][["RNA"]] <- NULL
  seurat_object_list[[i]][["SCT"]] <- NULL
}

## reference data for scATAC-seq data
anno_rds <- paste0(pipe_dir, "/workflow/dependencies/EnsDb.Hsapiens.v86_2UCSC_hg38.RDS")
anno_gene <- readRDS(anno_rds)
genome_info <- seqinfo(anno_gene)

}

###########################################################################
## integrate matched scATAC-seq data from multiple samples
###########################################################################
{
  ##########################################
  ## creating common peak set across samples
  {
    peaks_merged <- vector()

    for(i in 1: length(sample_ids))
    {
      ## merged as list
      peaks_merged  <- c(peaks_merged , seurat_object_list[[i]]@assays$ATAC@ranges)
    }

    ## appending granges
    peaks_merged  <-  unlist(as(peaks_merged, "GRangesList"))

    # Reduce and Filter out bad peaks based on length
    # peaks_merged <- disjoin(x = peaks_merged)
    peaks_merged <- reduce(x = peaks_merged)     ## less number of peaks

    peaks_width <- width(peaks_merged)
    peaks_merged <- peaks_merged[peaks_width < 10000 & peaks_width > 20]
    length(peaks_merged)
  }

  ################
  ## merge samples
  {
    ## Function Re-quantify integrated peaks in each sample
    ## on H4H
    if(TRUE){
      for(i in 1: length(sample_ids))
      {
        ##
        peaks_count <- FeatureMatrix(
                        fragments = Fragments(seurat_object_list[[i]]),    ## for cluster
                        features  = peaks_merged,
                        cells = colnames(seurat_object_list[[i]]))

        seurat_object_list[[i]][["ATAC_integrated"]] <- CreateChromatinAssay(counts = peaks_count,
                                                                             genome = genome_info,
                                                                             fragments = Fragments(seurat_object_list[[i]]),
                                                                             min.cells = 0,
                                                                             min.features = -1,        ##  set as to negative to ensure same number of cells!!
                                                                             annotation = anno_gene)
      }
    }

    ## Function Re-quantify integrated peaks in each sample
    ## for local testing
    if(FALSE){
      for(i in 1: length(sample_ids))
      {

        fpath_local <- paste0(getwd(), "/arc_count/", sample_ids[i], "/outs/atac_fragments.tsv.gz" )
        fragments_local <- CreateFragmentObject(path = fpath_local, cells = colnames(seurat_object_list[[i]]))

        peaks_count <- FeatureMatrix(
                        #fragments = Fragments(x),    ## for cluster
                        fragments = fragments_local,        ##  for local testing
                        features  = peaks_merged,
                        cells = colnames(seurat_object_list[[i]]))

        seurat_object_list[[i]][["ATAC_integrated"]] <- CreateChromatinAssay(counts = peaks_count,
                                                                            genome = genome_info,
                                                                            fragments = fragments_local,
                                                                            min.cells = 0,
                                                                            min.features = -1,        ##  set as to negative to ensure same number of cells!!
                                                                            annotation = anno_gene)
      }
    }

    ###########################
    ## merge samples
    ## switch to ATAC_integrated
    for(i in 1: length(sample_ids))
    {
      DefaultAssay(seurat_object_list[[i]]) <- 'ATAC_integrated'
    }

    atac_merged <- merge(x = seurat_object_list[[1]],
                         y = seurat_object_list[2:length(seurat_object_list)])

    ## re analysis merged atac-seq
    atac_merged <- FindTopFeatures(atac_merged, min.cutoff = 'q0')    ## q0 -> 100% ; 95 -> 95% cells
    atac_merged <- RunTFIDF(atac_merged)
    atac_merged <- RunSVD(atac_merged)
    atac_merged <- RunUMAP(atac_merged, reduction = "lsi", dims = 2:50)

    atac_merged <- FindNeighbors(atac_merged , reduction = "lsi", dims = 2:50) %>%
                   FindClusters(algorithm = 3)

    saveRDS(atac_merged, paste0(out_dir, "/ATAC_integrated_by_merging.RDS"))

  }

  #####################
  ## harmonized samples
  {
    ## sample id as the group.by.vars  by default
    atac_harmonized <-  RunHarmony(atac_merged,
                                   group.by.vars = c("sample_id"),
                                   reduction = "lsi",
                                   assay.use = "ATAC_integrated",
                                   project.dim = F,            # https://github.com/immunogenomics/harmony/issues/86
                                   reduction.save = "harmony")

    atac_harmonized  <- RunUMAP(atac_harmonized, reduction = "harmony", dims = 2:50, verbose = FALSE)
    atac_harmonized  <- FindNeighbors(atac_harmonized, reduction = "harmony", dims = 2:50) %>%
                        FindClusters(algorithm = 3)

    saveRDS(atac_harmonized, paste0(out_dir, "/ATAC_integrated_by_harmony.RDS"))
  }

  ############################################
  ## integrated samples on top of atac_merged
  {
    ## switch to ATAC_integrated
    for(i in 1: length(sample_ids))
    {
      DefaultAssay(seurat_object_list[[i]]) <- 'ATAC_integrated'
    }

    ## redo lsi for ATAC_merged per object, gonna replace previous lsi using ATAC assay for FindIntegrationAnchors
    seurat_object_list <- lapply(X = seurat_object_list,
                                 FUN = function(x) {
                                   x <- FindTopFeatures(x, min.cutoff = 'q0')
                                   x <- RunTFIDF(x)
                                   x <- RunSVD(x)
                                   }
    )

    ## identify ATAC integration features and anchors
    # atac_features <- SelectIntegrationFeatures(object.list = seurat_object_list,
    #                                          assay = rep("ATAC_merged", length(seurat_object_list)))
    ## might cause "number of items to replace is not a multiple of replacement length"

    ## reference: A vector specifying the object/s to be used as a reference during integration.
    ## reference: If NULL (default), all pairwise anchors are found (no reference/s).
    ## Could using the first object as reference to save computational time, but might affect the outcomes as well

    atac_anchors <- FindIntegrationAnchors(
                              object.list = seurat_object_list,
                              #reference = c(1),
                              #anchor.features = atac_features,
                              anchor.features = rownames(rownames(atac_merged)),
                              reduction = "rlsi",
                              k.filter = 200,    # default 200; How many neighbors (k) to use when filtering anchors
                              dims = 2:50)

    # integrate LSI embeddings
    # https://github.com/satijalab/seurat/issues/6341
    # Error in idx[i, ] <- res[[i]][[1]] : number of items to replace is not a multiple of replacement length
    atac_anchored <- IntegrateEmbeddings(anchorset = atac_anchors,
                                         reductions = atac_merged[["lsi"]],
                                         new.reduction.name = "lsi_integrated",
                                         k.weight = 100      # default 100; Number of neighbors to consider when weighting anchors
                                         )

    # create a new UMAP using the integrated embedings
    atac_anchored  <- RunUMAP(atac_anchored, reduction = "lsi_integrated", dims = 2:50)

    atac_anchored  <- FindNeighbors(atac_anchored, reduction = "lsi_integrated", dims = 2:50) %>%
                      FindClusters(algorithm = 3)

    saveRDS(atac_anchored, paste0(out_dir, "/ATAC_integrated_by_anchors.RDS"))

  }


  ####################################################
  ##  UMAPs for merged, harmonized and integrated RNA
  {
    ## UMAP plots group by sample_id
    p1 <- DimPlot(atac_merged,  reduction = "umap", group.by = "sample_id",
                  label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("ATAC_merged") + NoLegend()

    p2 <- DimPlot(atac_harmonized,  reduction = "umap", group.by = "sample_id",
                label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("ATAC_harmonized") # + NoLegend()

    p3 <- DimPlot(atac_anchored,  reduction = "umap", group.by = "sample_id",
                  label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("ATAC_anchored")

    g <- p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste0(out_dir, "/Merged_Harmonized_Anchored_ATAC_UMAPs_labeled_by_sample.pdf"), width = 12, height = 4)


    ## UMAP plots group by clusters (different from assay to assay)

    p1 <- DimPlot(atac_merged,  reduction = "umap", group.by = "ATAC_integrated_snn_res.0.8",
                  label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("ATAC_merged")

    p2 <- DimPlot(atac_harmonized ,  reduction = "umap", group.by = "ATAC_integrated_snn_res.0.8",
                  label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("ATAC_harmonized")

    p3 <- DimPlot(atac_anchored,  reduction = "umap", group.by = "ATAC_integrated_snn_res.0.8",
                  label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("ATAC_anchored")

    g <- p1 + p2 + p3  & theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste0(out_dir, "/Merged_Harmonized_Anchored_ATAC_UMAPs_labeled_by_cluster.pdf"), width = 15, height = 4)
  }
}
