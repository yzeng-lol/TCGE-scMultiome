args = commandArgs(trailingOnly=TRUE)

## assigning commandArgs
integr_meth  = args[1]     ## anchor: seurat anchor by default; harmony: by harmony

integrated_rna = paste0(getwd(), "/", args[2])
integrated_atac = paste0(getwd(), "/", args[3])

#setwd(paste0(getwd(), "/integration"))

out_dir = paste0(getwd(), "/integration/wnn/", integr_meth)

###############################
### for testing and fine-tuning
###############################
if(FALSE){
#########
## on H4H
# cd  /cluster/projects/tcge/scMultiome/iSHARC_test/main_seurat
# conda activate  /cluster/home/yzeng/miniconda3/envs/iSHARC_extra_env/ebcd8f410d62ed43e71d1cebab0b3296_    ## new
#R

rm(list = ls())
setwd("/cluster/projects/tcge/scMultiome/processed_data/TCGE-scMOS-PTLD/integration")

integr_meth  = harmony
out_dir = paste0(getwd(), "/harmony")

####################################
## testing locally with multiple RDS
rm(list = ls())
setwd("/Users/yong/OneDrive - UHN/Projects/snakemake/iSHARC_test/PRCA/integration")

out_dir = getwd()
integr_meth  = "harmony"


}

##################################################
## The main workflow for scMultiome data analysis
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

#####################################################
## creat a seurat object with integrated RNA and ATAC
## integrate RNA and ATAC using WNN
#####################################################
if(FALSE){
if(integr_meth == "anchor"){

  ## using harmonized for testing
  #scMultiome <- readRDS("RNA_integrated_by_anchors.RDS")
  #atac_integrated <- readRDS("ATAC_integrated_by_anchors.RDS")

  scMultiome <- readRDS(integrated_rna)
  atac_integrated <- readRDS(integrated_atac)

  scMultiome[["ATAC_integrated"]] <- atac_integrated[["ATAC_integrated"]]   ## add assay
  scMultiome[["lsi_integrated"]] <- atac_integrated@reductions$lsi        ## add reduction
  scMultiome[["umap_atac"]]  <- atac_integrated@reductions$umap             ## add reduction

  scMultiome[["ATAC_integrated_snn_res.0.8"]]  <- atac_integrated$ATAC_integrated_snn_res.0.8    ## add meta data
  scMultiome[["RNA_ATAC_perSample_wsnn_res.0.8"]]  <- scMultiome$wsnn_res.0.8

  scMultiome <- FindMultiModalNeighbors(scMultiome, reduction.list = list("pca", "lsi_integrated"), dims.list = list(1:50, 2:50))

} else if (integr_meth == "harmony") {

  ## using harmonized for testing
  scMultiome <- readRDS("RNA_integrated_by_harmony.RDS")
  atac_integrated <- readRDS("ATAC_integrated_by_harmony.RDS")

  scMultiome[["ATAC_integrated"]] <- atac_integrated[["ATAC_integrated"]]   ## add assay
  scMultiome[["harmony_atac"]] <- atac_integrated@reductions$harmony        ## add reduction
  scMultiome[["umap_atac"]]  <- atac_integrated@reductions$umap             ## add reduction

  scMultiome[["ATAC_integrated_snn_res.0.8"]]  <- atac_integrated$ATAC_integrated_snn_res.0.8    ## add meta data
  scMultiome[["RNA_ATAC_perSample_wsnn_res.0.8"]]  <- scMultiome$wsnn_res.0.8

  scMultiome <- FindMultiModalNeighbors(scMultiome, reduction.list = list("harmony", "harmony_atac"), dims.list = list(1:50, 2:50))
}
}

{
## integerated data
scMultiome <- readRDS(integrated_rna)
atac_integrated <- readRDS(integrated_atac)

scMultiome[["ATAC_integrated"]] <- atac_integrated[["ATAC_integrated"]]   ## add assay
scMultiome[["lsi_integrated"]] <- atac_integrated@reductions$lsi        ## add reduction
scMultiome[["umap_atac"]]  <- atac_integrated@reductions$umap             ## add reduction

scMultiome[["ATAC_integrated_snn_res.0.8"]]  <- atac_integrated$ATAC_integrated_snn_res.0.8    ## add meta data
scMultiome[["RNA_ATAC_perSample_wsnn_res.0.8"]]  <- scMultiome$wsnn_res.0.8

scMultiome <- FindMultiModalNeighbors(scMultiome, reduction.list = list("pca", "lsi_integrated"), dims.list = list(1:50, 2:50))


## Apply WNN
scMultiome <- RunUMAP(scMultiome, nn.name = "weighted.nn", reduction.name = "umap_wnn", )
scMultiome <- FindClusters(scMultiome, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

scMultiome[["RNA_ATAC_integrated_wsnn_res.0.8"]]  <- scMultiome$wsnn_res.0.8   ## wsnn_res.0.8 has been overwritten
scMultiome[["wsnn_res.0.8"]] <- NULL

saveRDS(scMultiome, paste0(out_dir, "/RNA_ATAC_integrated_by_WNN.RDS"))

}


###################################
## UMAP and cell labeling plots
##################################
# scMultiome <- readRDS(paste0(out_dir, "/RNA_ATAC_integrated_by_WNN.RDS"))

{
  ## QC result split by samples
  {

  plot_features <- c("nCount_RNA","percent.mt", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal")

  ## QC metrics per samples
  VlnPlot(scMultiome, features = plot_features,
          group.by = "sample_id", split.by = "sample_id", ncol = 5) & theme(axis.title.x = element_blank())

  ggsave(paste0(out_dir, "/QC_metrics_grouped_and_split_by_sample.pdf"), width = 16, height = 4)

  ## QC metrics per clusters
  for(pf in plot_features)
  {
  VlnPlot(scMultiome, features = pf,group.by = "WNN_SingleR_anno", split.by = "sample_id") & theme(axis.title.x = element_blank())

  ggsave(paste0(out_dir, "/QC_metrics_", pf, "_grouped_by_SingleR_anno_and_split_by_sample.pdf"), width = 10, height = 4)
  }

  }


  ## UMAPs group by clusters
  {
  p1 <- DimPlot(scMultiome, reduction = "umap",  group.by = "RNA_integrated_snn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA_integrated")
  p2 <- DimPlot(scMultiome, reduction = "umap_atac", group.by = "ATAC_integrated_snn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC_integrated")
  p3 <- DimPlot(scMultiome, reduction = "umap_wnn",  group.by = "RNA_ATAC_integrated_wsnn_res.0.8",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN_RNA_ATAC")
  g <- p1 + p2 + p3  & theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(out_dir, "/RNA_ATAC_integrated_WNN_UMAPs_grouped_by_cluster.pdf"), width = 12, height = 4)

  ## UMAPs group by sample_ids
  p1 <- DimPlot(scMultiome, reduction = "umap",  group.by = "sample_id", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA_integrated") + NoLegend()
  p2 <- DimPlot(scMultiome, reduction = "umap_atac", group.by = "sample_id", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC_integrated") + NoLegend()
  p3 <- DimPlot(scMultiome, reduction = "umap_wnn",  group.by = "sample_id",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN_RNA_ATAC")
  g <- p1 + p2 + p3  & theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(out_dir, "/RNA_ATAC_integrated_WNN_UMAPs_grouped_by_sample.pdf"), width = 12, height = 4)

  ## morality weights
  p1 <- VlnPlot(scMultiome, features = "RNA_integrated.weight", group.by = 'RNA_ATAC_integrated_wsnn_res.0.8', sort = F, pt.size = 0.1)
  p1 <- p1 + ggtitle("RNA weights") & theme(axis.text.x = element_blank(), axis.title.x = element_blank())  + NoLegend()

  p2 <- VlnPlot(scMultiome, features = "ATAC_integrated.weight", group.by = 'RNA_ATAC_integrated_wsnn_res.0.8', sort = F, pt.size = 0.1)
  p2 <- p2 + ggtitle("ATAC weights") + labs(x = "WNN_clusters") + NoLegend()

  g <- p1 + p2
  ggsave(paste0(out_dir, "/RNA_ATAC_integrated_WNN_weights_for_clusters.pdf"), width = 12, height = 5)
  }


  #########################################
  ## cell annotation and clusters per sample
  {

    ## cancer vs no cancer cells
    dat  <- scMultiome[[c("sample_id","copykat_anno")]]
    dat_sum <- dat %>%
      group_by(sample_id, copykat_anno) %>%
      summarise(cnt = n())

    p1 <- ggplot(dat_sum, aes( y = cnt, x = sample_id, fill = copykat_anno))
    p1 <- p1 + geom_bar(position="fill", stat="identity")
    p1 <- p1 + labs(y = "Fraction of cells",  x = "")
    p1 <- p1 + guides(fill=guide_legend(title="CopyKat prediction")) + theme_classic()

    ## cell annotation based single sample WNN clustering based prediction
    dat  <- scMultiome[[c("sample_id","WNN_SingleR_anno")]]
    dat_sum <- dat %>%
      group_by(sample_id, WNN_SingleR_anno) %>%
      summarise(cnt = n())

    p2 <- ggplot(dat_sum, aes( y = cnt, x = sample_id, fill = WNN_SingleR_anno))
    p2 <- p2 + geom_bar(position="fill", stat="identity")
    p2 <- p2 + labs(y = "Fraction of cells",  x = "")
    p2 <- p2 + guides(fill=guide_legend(title="SingleR annotation")) + theme_classic()

    ## cell annotation based single sample WNN clustering based prediction
    dat  <- scMultiome[[c("sample_id","RNA_ATAC_integrated_wsnn_res.0.8")]]
    dat_sum <- dat %>%
      group_by(sample_id, RNA_ATAC_integrated_wsnn_res.0.8) %>%
      summarise(cnt = n())

    p3 <- ggplot(dat_sum, aes( y = cnt, x = sample_id, fill = RNA_ATAC_integrated_wsnn_res.0.8))
    p3 <- p3 + geom_bar(position="fill", stat="identity")
    p3 <- p3 + labs(y = "Fraction of cells",  x = "")
    p3 <- p3 + guides(fill=guide_legend(title="WNN_clusters")) + theme_classic()

    g <- p1 + p2 + p3
    ggsave(paste0(out_dir, "/Fraction_of_cells_types_based_indiviudal_sample_labeling.pdf"), width = 15, height = 6)

  }

  ###################################################
  ## annotated cell types per integrated WNN clusters
  {
    dat  <- scMultiome[[c("RNA_ATAC_integrated_wsnn_res.0.8","WNN_SingleR_anno")]]
    dat_sum <- dat %>%
      group_by(RNA_ATAC_integrated_wsnn_res.0.8, WNN_SingleR_anno) %>%
      summarise(cnt = n())

    p1 <- ggplot(dat_sum, aes( y = cnt, x = RNA_ATAC_integrated_wsnn_res.0.8, fill = WNN_SingleR_anno))
    p1 <- p1 + geom_bar(position="fill", stat="identity")
    p1 <- p1 + labs(y = "Fraction of cells",  x = "WNN clusters based integrated RNA and ATAC")
    p1 <- p1 + guides(fill=guide_legend(title="SingleR annotated per sample")) + theme_classic()
    ggsave(paste0(out_dir, "/Fraction_of_SingleR_annoated_cells_types_based_indiviudal_sample_across_WNN_clusters.pdf"), width = 10, height = 5)

  }


}
