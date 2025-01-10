################################################################################
##  vertical_integration_of_individual_sample.R is a function for vertical
##  integration of matched RNA and ATAC for each individual sample, and the
##  analyses include :
##      * Second found QC filtering (optional)
##      * cell cycle correction for RNA (optional)
##      * normalization and dimensionalities reduction
##      * integration of RNA and ATAC using WNN
##      * clustering using the integrated or separated modalities.
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
parser$add_argument("-iso", "--initial_seurat_object", required=TRUE,
                    help = "the initialized seurat object with RNA and ATAC assays")
parser$add_argument("-srf", "--second_round_filter", type = "logical", default = TRUE,
                    help = "whether perform second round cell filtering based on seleted QC metircs' MADs and suggested cutoffs")
parser$add_argument("-min_RNA", "--min_nCount_RNA", type = "integer", default = 500,
                    help = "Minimal nCount_RNA for second round QC filtering")
parser$add_argument("-max_RNA", "--max_nCount_RNA", type = "integer", default = 25000,
                    help = "Maximum nCount_RNA for second round QC filtering")
parser$add_argument("-min_ATAC", "--min_nCount_ATAC", type = "integer", default = 1000,
                    help = "Minimal nCount_ATAC for second round QC filtering")
parser$add_argument("-max_ATAC", "--max_nCount_ATAC", type = "integer", default = 70000,
                    help = "Maximum nCount_ATAC for second round QC filtering")
parser$add_argument("-max_MT", "--max_pct_MT", type = "integer", default = 20,
                    help = "Maximum percentage of MT for second round QC filtering")
parser$add_argument("-mim_TSS", "--min_TSS_Enrichment", type = "integer", default = 1,
                    help = "Minimal TSS enrichment score for second round QC filtering")
parser$add_argument("-mim_NS", "--max_Nucleosome_Signal", type = "integer", default = 1,
                    help = "Maximun nucleosome signal score for second round QC filtering")
parser$add_argument("-rcc", "--regress_cell_cycle", type = "logical", default = TRUE,
                    help = "whether regress out the effects of cell cycyle for RNA assays")
parser$add_argument("-knn_k", "--knn_k_param", type = "integer", default = 20,
                    help = "k for the k-nearest neighbor algorithm")
parser$add_argument("-dims_n", "--dimentions_n", type = "integer", default = 50,
                    help = "number of reduced dimentions (e.g., PCs) for functions: RunUMAP, FindNeighbors, FindMultiModalNeighbors")
parser$add_argument("-comm_res", "--community_resolution", type = "float", default = 0.8,
                    help = "Value above (below) 1.0 if you want to obtain a larger (smaller) number of communities")
parser$add_argument("-t", "--threads", type = "integer", default = 12,
                     help = "Number of cores for the parallelization")
parser$add_argument("-fgm", "--future_globals_maxSize", type = "integer", default = 12,
                    help = "Maximum memory in GB for the future parallelization global variables")

## assigning passing arguments
args <- parser$parse_args()
print(args)

sample_id <- args$sample_id
iso_file <- args$initial_seurat_object
second_round_filter <- args$second_round_filter
regress_cell_cycle <- args$regress_cell_cycle
knn_k <- args$knn_k_param
dims_n <- args$dimentions_n
comm_res <- args$community_resolution


## output dir
out_dir <- paste0(getwd(), "/individual_samples/", sample_id, "/") ## with forward slash at the end

}

############################
### loading required packages
#############################
{
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))          ## for sankey plot
suppressMessages(library(tidyr))

## enable the Parallelization with the future packages
suppressMessages(library(future))
plan("multicore", workers = args$threads)
options(future.globals.maxSize = args$future_globals_maxSize * 1024^3)

}

####################################
### second round cell filtering
###################################
{
## read in seurat object
scMultiome <- readRDS(iso_file)

if(second_round_filter){
## selected QC metircs for second round cell filtering
qc_df <- data.matrix(scMultiome [[c("nCount_RNA", "pct_MT", "nCount_ATAC",
                                    "TSS_Enrichment", "Nucleosome_Signal")]])

## MADs for excluding outlier
upper_3mad <- apply(qc_df, 2, function(x) median(x, na.rm = TRUE) + 3 * mad(x, na.rm = TRUE))
lower_3mad <- apply(qc_df, 2, function(x) median(x, na.rm = TRUE) - 3 * mad(x, na.rm = TRUE))

## taking both MADs and preset cutoffs into account
#rna_upper <- upper_3mad["nCount_RNA"]
#rna_lower <- lower_3mad["nCount_RNA"]
#atac_uppper <- upper_3mad["nCount_ATAC"]
#atac_lower <- lower_3mad["nCount_ATAC"]

rna_upper <-  min(args$max_nCount_RNA, upper_3mad["nCount_RNA"])
rna_lower <-  max(args$min_nCount_RNA, lower_3mad["nCount_RNA"])
atac_uppper <- min(args$max_nCount_ATAC, upper_3mad["nCount_ATAC"])
atac_lower  <- max(args$min_nCount_ATAC, lower_3mad["nCount_ATAC"])
mt_upper  <- min(args$max_MT, upper_3mad["pct_MT"])
tss_lower <- max(args$min_TSS_Enrichment, lower_3mad["TSS_Enrichment"])
ns_upper  <- min(args$max_Nucleosome_Signal, upper_3mad["Nucleosome_Signal"])

cells_before <- nrow(scMultiome@meta.data)
## subsetting the scMultiome project
scMultiome <- subset(
                x = scMultiome,
                subset = nCount_RNA < rna_upper &
                         nCount_RNA > rna_lower &
                         nCount_ATAC < atac_uppper &
                         nCount_ATAC > atac_lower &
                         pct_MT < mt_upper &
                         TSS_Enrichment > tss_lower &
                         Nucleosome_Signal < ns_upper
                 )

## cells after the second round filtering
cells_after <- nrow(scMultiome@meta.data)

print(paste0(cells_after, " of ", cells_before, " cell remained after second round filtering !!"))
}

}



########################################################################
## RNA alone Normalization, Dimension Reduction, Clustering and Embedding
########################################################################
{
DefaultAssay(scMultiome) <- "RNA"

################
## normalization
## SCTransform :: an alternative to the NormalizeData, FindVariableFeatures, ScaleData workflow
## BiocManager::install("glmGamPoi")
# install sctransform from Github
## devtools::install_github("satijalab/sctransform", ref = "develop")

##########################################
## also need to regress out the cell cycle
## https://github.com/satijalab/seurat/issues/1679
## however, in some cases, we’ve found that this can negatively impact downstream analysis,
## particularly in differentiating processes (like murine hematopoiesis)
## we suggest regressing out the difference between the G2M and S phase scores
## https://satijalab.org/seurat/articles/cell_cycle_vignette.html

if(regress_cell_cycle)
{
  ## by three steps
  # 1) normalize data with SCTransform()
  scMultiome <- SCTransform(
                  scMultiome,
                  assay = 'RNA',
                  new.assay.name = 'SCT',
                  vars.to.regress = c('pct_MT'), #, 'nFeature_RNA', 'nCount_RNA')
                  #vst.flavor = "v2"                 ## will fail the CellCycleScoreing
                )

  # 2) perform cell cycle analysis based on normalized RNA
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes

  scMultiome <- CellCycleScoring(
                  scMultiome,
                  s.features = s.genes,
                  g2m.features = g2m.genes,
                  assay = 'SCT',
                  set.ident = TRUE
                )

  # 3)normalize again but this time including also the cell cycle scores
  scMultiome <- SCTransform(
                  scMultiome,
                  assay = 'RNA',
                  new.assay.name = 'SCT',
                  vars.to.regress = c('pct_MT','S.Score', 'G2M.Score'),
                  vst.flavor = "v2", verbose = FALSE
                )
} else {
  scMultiome <- SCTransform(
                  scMultiome,
                  assay = 'RNA',
                  new.assay.name = 'SCT',
                  vars.to.regress = c('pct_MT'),
                  vst.flavor = "v2", verbose = FALSE
                )
}

## !! the default assay has been shifted to SCT
## dimensional reduction and clustering
scMultiome <- RunPCA(scMultiome) %>%
              RunUMAP(dims = 1:dims_n, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

## RNA alone based clustering
scMultiome  <- FindNeighbors(scMultiome, reduction = "pca", dims = 1:dims_n, k.param = knn_k) %>%
               FindClusters(resolution = comm_res)

}


###########################################################################
## ATAC alone Normalization, Dimension Reduction, Clustering and Embedding
###########################################################################
{

## Focusing the ATAC assay with MACS2 called peaks
if(FALSE){
# exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(scMultiome) <- "ATAC_ARC"
scMultiome <- RunTFIDF(scMultiome)            # TF-IDF normalization

scMultiome <- FindTopFeatures(scMultiome, min.cutoff = 'q0')    ## q0 -> 100% ; 95 -> 95% cells
scMultiome <- RunSVD(scMultiome)
scMultiome <- RunUMAP(scMultiome, reduction = 'lsi', dims = 2:dims_n, reduction.name = "umap.atac.arc", reduction.key = "atacArcUMAP_")

## ATAC alone based clustering
scMultiome  <- FindNeighbors(scMultiome, reduction = "lsi", dims = 2:dims_n,  k.param = knn_k) %>%
                FindClusters(algorithm = 3, resolution = comm_res)
}

# exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(scMultiome) <- "ATAC"
scMultiome <- RunTFIDF(scMultiome)            # TF-IDF normalization

scMultiome <- FindTopFeatures(scMultiome, min.cutoff = 'q0')    ## q0 -> 100% ; 95 -> 95% cells
scMultiome <- RunSVD(scMultiome)
scMultiome <- RunUMAP(scMultiome, reduction = 'lsi', dims = 2:dims_n, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

## ATAC alone based clustering
scMultiome  <- FindNeighbors(scMultiome, reduction = "lsi", dims = 2:dims_n,  k.param = knn_k) %>%
               FindClusters(algorithm = 3, resolution = comm_res)


## note the reduction lsi will be overwritten using assay ATAC !!!
}


##############################################
## vertically integrate RNA and ATAC using WNN
##############################################
{
## run WNN
scMultiome <- FindMultiModalNeighbors(scMultiome, reduction.list = list("pca", "lsi"), dims.list = list(1:dims_n, 2:dims_n))
scMultiome <- RunUMAP(scMultiome, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_")
scMultiome <- FindClusters(scMultiome, graph.name = "wsnn", algorithm = 3, resolution = comm_res)

## the default resolution (=0.8)
scMultiome$RNA_clusters <- scMultiome[[paste0("SCT_snn_res.", comm_res)]]
scMultiome$ATAC_clusters <- scMultiome[[paste0("ATAC_snn_res.", comm_res)]]
scMultiome$WNN_clusters <- scMultiome[[paste0("wsnn_snn_res.", comm_res)]]

## reorder cluster levels
N_levels <- length(unique(scMultiome$RNA_clusters))
scMultiome$RNA_clusters<- factor(scMultiome$RNA_clusters,
                                  levels = 0:(N_levels-1))

N_levels <- length(unique(scMultiome$ATAC_clusters))
scMultiome$ATAC_clusters<- factor(scMultiome$ATAC_clusters,
                                 levels = 0:(N_levels-1))

N_levels <- length(unique(scMultiome$WNN_clusters))
scMultiome$WNN_clusters <- factor(scMultiome$WNN_clusters,
                                  levels = 0:(N_levels-1))

saveRDS(scMultiome, file = paste0(out_dir, sample_id, "_vertically_integrated_seurat_object.RDS"))

write.csv(scMultiome@meta.data,  file = paste0(out_dir, sample_id, "_vertically_integrated_meta_data.csv"))


###############################################
## UMAPs for ATAC, RAN separatly and integrated
p1 <- DimPlot(scMultiome, reduction = "umap.atac", group.by = "ATAC_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p2 <- DimPlot(scMultiome, reduction = "umap.rna",  group.by = "RNA_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p3 <- DimPlot(scMultiome, reduction = "umap.wnn", group.by = "WNN_clusters",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

g <- p1 + p2 + p3  & theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(out_dir, sample_id, "_ATAC_RNA_WNN_clustering_UMAPs.pdf"), width = 12, height = 4)

########################################################
### modality weight per integrated clusters
p1 <- VlnPlot(scMultiome, features = "ATAC.weight", group.by = 'WNN_clusters', sort = FALSE, pt.size = 0.1) +  ggtitle("ATAC weights")
p1 <- p1 + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + theme(legend.position = 'none')

p2 <- VlnPlot(scMultiome, features = "SCT.weight", group.by = 'WNN_clusters', sort = FALSE, pt.size = 0.1)
p2 <- p2 + ggtitle("RNA weights") + theme(legend.position = 'none')
g <- p1 + p2   & theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(out_dir, sample_id, "_ATAC_RNA_WNN_weights.pdf"), width = 12, height = 6)

####################################################################################
### Sankey plot for clustering changes across individually and integrated modalities
Cell_ID <- rownames(scMultiome@meta.data)
RNA_Clusters <- paste0("RNA_", scMultiome@meta.data$RNA_clusters)
ATAC_Clusters <- paste0("ATAC_", scMultiome@meta.data$ATAC_clusters)
WNN_Clusters <- paste0("WNN_", scMultiome@meta.data$WNN_clusters)

df <- data.frame(Cell_ID, RNA_Clusters, ATAC_Clusters, WNN_Clusters)
df <- as_tibble(df)

## unique nodes and indexing
nodes <- df %>%
         pivot_longer(-Cell_ID, values_to = "name_node") %>%
         distinct(name_node) %>%
         arrange(name_node) %>%
         mutate(idx = (1:n()) - 1)

## links with source, target and values (count)
links <- bind_rows(df %>% dplyr::select(source = ATAC_Clusters, target = RNA_Clusters),
                   df %>% dplyr::select(source = RNA_Clusters, target = WNN_Clusters)) %>%
         group_by(source, target) %>% dplyr::count(target) %>%
         mutate(value = n) %>% dplyr::select(!n) %>%  ungroup()

links$source <- nodes$idx[match(links$source, nodes$name_node)]
links$target <- nodes$idx[match(links$target, nodes$name_node)]

g <- plot_ly(
     type = "sankey",
     orientation = "h",
     node = list(label = nodes$name_node, pad = 15, thickness = 15),
     link = as.list(links))
## export the interactive HTML files
html_name <- paste0(out_dir, sample_id, "_ATAC_RNA_WNN_clusters_Sankey_plot.html")
htmlwidgets::saveWidget(as_widget(g), html_name)

print("The vertical integration of individual sample has been successfully completed!!")

}
