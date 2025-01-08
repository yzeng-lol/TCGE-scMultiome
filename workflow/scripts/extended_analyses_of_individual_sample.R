################################################################################
## extended analyses based on integrated RNA and ATAC for each individual sample
## The analyses includes:
##      * cell type annotation
##      * Identify cluster specific DEGs and enriched functions
##      * Identify cluster specific DARs and enriched motif/TFs
##      * Gene regulatory network analysis
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
parser$add_argument("-vio", "--vertically_integrated_seurat_object", required=TRUE,
                    help = "the seurat object with vertically integrated ATAC and RNA using WNN")
parser$add_argument("-pipe", "--pipe_dir", required=TRUE,
                    help = "The PATH to iSHARC pipeline, which local dependences included")

## assigning passing arguments
args <- parser$parse_args()
print(args)

sample_id <- args$sample_id
vio_file <- args$vertically_integrated_seurat_object
pipe_dir <- args$pipe_dir


## output dir
out_dir <- paste0(getwd(), "/individual_samples/", sample_id, "/") ## with forward slash at the end

}


############################
### loading required packages
#############################
{
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(SingleR))
suppressMessages(library(clusterProfiler))
suppressMessages(library(JASPAR2020))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(dplyr))
suppressMessages(library(tidygraph))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggraph))
suppressMessages(library(devtools))
suppressMessages(library(TFBSTools))

## enable the Parallelization with the future packages
suppressMessages(library(future))
plan("multicore", workers = 12)


#########################################
## loading packages from the dependencies

##
## loading the KEGG.db from the local, since it will be problematic for clusterProfiler's defaultst KEGG analysis require internet access
## https://github.com/YuLab-SMU/createKEGGdb/tree/master
if(!require("KEGG.db"))  install.packages(paste0(pipe_dir, "/workflow/dependencies/KEGG.db_1.0.tar.gz"))
library(KEGG.db)

## install "dlm", which is required for copykat
if(!require("dlm"))  install.packages(paste0(pipe_dir, "/workflow/dependencies/dlm_1.1-6.tar.gz"))
library(dlm)

## load copykat (v1.1.0)
if(!require("copykat"))  devtools::load_all(paste0(pipe_dir, "/workflow/dependencies/copykat"))
library(copykat)

## load Pando
## require seuratObject < 5.0.0;
if(!require("Pando"))  install.packages(paste0(pipe_dir, "/workflow/dependencies/Pando_1.0.4.tar.gz"))  ## latest compatible verion
library(Pando)

suppressMessages(library(doParallel))       ## parallelization for pando


## intall TFBSTools from local instead
#if(!require("TFBSTools"))  install.packages(paste0(pipe_dir, "/workflow/dependencies/TFBSTools_1.44.0.tar.gz"), INSTALL_opts = '--no-lock')
#library(TFBSTools)

}


##############################################################
## annotating the cell clusters with integrated ATAC and RNA
##############################################################
{
## readin vertically integrated seurat object
scMultiome <- readRDS(vio_file)

## ensure the active.ident is WNN_Clusters
scMultiome  <- SetIdent(scMultiome , value = scMultiome @meta.data$WNN_clusters)

##############################################################
## auto annotation using publicly available reference datasets
{
## anno_ref <-  BlueprintEncodeData()       ## form package celldex, internet required
ref_rds <- paste0(pipe_dir, "/workflow/dependencies/BlueprintEncodeData.RDS")
anno_ref <- readRDS(ref_rds)

## fetch SCT normalized GEX matrix
expr <- GetAssayData(object = scMultiome, assay = "SCT", slot = "data")

### using ENCODE
expr_anno <- SingleR(test = expr, ref = anno_ref, labels = anno_ref$label.main, clusters =  Idents(scMultiome))

## match cluster labels and annotated labels
idx_m <- match(Idents(scMultiome), rownames(expr_anno))

## add labels scMultiome object
scMultiome[["WNN_clusters_singler_annot"]] <- expr_anno$labels[idx_m]
}

############################################################
## Distinguish the tumor cells from the normal cells
if(TRUE){
## Using copyKAT  to predicts tumor and normal cells
## RNA-seq data based

setwd(out_dir)                                     ## ensure output copykat related results to desired folder

expr_raw <- as.matrix(scMultiome@assays$RNA@counts)

copykat_res <- copykat(rawmat = expr_raw, sam.name = sample_id , id.type = "S", ngene.chr = 5, win.size = 25,
                        KS.cut = 0.1,  distance = "euclidean", norm.cell.names = "", output.seg = "FLASE",
                        plot.genes = "TRUE", genome = "hg20", n.cores = 1)

## adding copykat predict labels to the scMultiome metatable
idx_s <- match(rownames(scMultiome@meta.data),copykat_res$prediction$cell.names)

## no
cell_type <- rep("not.predicted", length(idx_s))    ## there are cells will excluded for copykat prediction
cell_type[!is.na(idx_s)] <- copykat_res$prediction$copykat.pred[idx_s[!is.na(idx_s)]]

scMultiome[["WNN_clusters_copykat_annot"]] <- cell_type

setwd("../../")  ## cd back to workdir
}

print("The annotation has been successfully completed!!")
}


##########################################
# identify cluster-specific genes (DEGs)
## plus one vs other :: cluster markers
## results were add to seuratObject@misc
########################################
{
clusters <- levels(scMultiome)
L <- length(clusters)
cluster_cnt <- table(scMultiome$seurat_clusters)

DefaultAssay(scMultiome) <- "SCT"
sct_deg <- list()
sct_deg_names <- vector()

## top 5 degs per clusters
deg_list <- vector()

################
## Identify DEGs
for (i in 1:L)
{
  if(cluster_cnt[i] <=3 ){

  sct_deg[[i]] <- list()

  } else {
  ## one vs all others: prefiltering
  sct_deg[[i]] <- FindMarkers(scMultiome, ident.1 = clusters[i], ident.2 = NULL,
                              min.pct = 0.5,                ## detected at least 50% frequency in either ident.
                              logfc.threshold = log(2),     ## at least two-fold change between the average expression of comparisons
                              min.diff.pct = 0.25,          ## Pre-filter features whose detection percentages across the two groups are similar
                              )

  #sct_deg_names[i] <- paste0(clusters[i], "_specific")
  ## select top 5 upregualted degs per clusters

  idx_fc <- sct_deg[[i]][,  2] > 0
  pval <-   sct_deg[[i]][idx_fc,  1]

    if (length(pval) == 0) {
      next
      } else {
      names(pval) <- rownames(sct_deg[[i]])[idx_fc]
      pval_s <- sort(pval)
      idx_g <- min(length(pval_s), 5)
      deg_list <- c(deg_list, names(pval_s)[1:idx_g])
      }
  }
}

#########################################################
names(sct_deg) <-  clusters         ##sct_deg_names
deg_list <- unique(deg_list)        ## remove duplicates for the top 5 DEGs

## add to assay's Misc of seuratObject: can call by : scMultiome@misc$
Misc(scMultiome@assays$SCT, slot = "DEGs") <- sct_deg
Misc(scMultiome@assays$SCT, slot = "DEGs_top5") <- deg_list       ## top5 gene names

#################################################################
### Functional enrichment analysis for WNN cluster-specific genes
if(TRUE){
    ## ID convert function
    id_convert <- function(x){

      library(clusterProfiler)

      g_symbol <- rownames(x)
      if(length(g_symbol) <= 10){
        ## requiring at least 10 genes
        g_ezid <- NULL
      } else {
        g_con <- bitr(g_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
        g_ezid <- g_con$ENTREZID
      }

      return(g_ezid)
    }

    deg_list <- lapply(scMultiome@assays$SCT@misc$DEGs, id_convert)

    ## remove empty elements in the list
    deg_list <- deg_list[lapply(deg_list, length) > 0]

    ## KEGG enrichment
    compKEGG <- compareCluster(geneCluster   = deg_list,
                               fun           = "enrichKEGG",
                               pvalueCutoff  = 0.05,
                               pAdjustMethod = "BH",
                               use_internal_data =T)    ## using local build library "KEGG.db"

    if(length(compKEGG@compareClusterResult$ID) > 0){
    p1 <- dotplot(compKEGG, title = "KEGG pathway enrichment") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45))
    ggsave(paste0(out_dir, sample_id, "_WNN_clusters_specific_DEGs_KEGG.pdf"), width = 12, height = 8)
    }

     ## GO enrichment
    compGO <- compareCluster(geneCluster   = deg_list ,
                             fun           = "enrichGO",
                             OrgDb='org.Hs.eg.db',
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH")

    if(length(compGO@compareClusterResult$ID) > 0){
    p2 <- dotplot(compGO, title = "GO enrichment ") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45))
    ggsave(paste0(out_dir, sample_id, "_WNN_clusters_specific_DEGs_GO.pdf"), width = 12, height = 8)
    }
}



################################################################
## draw heat map for top 5 clusters specifically expressed genes

source(paste0(pipe_dir, "/workflow/scripts/DoMutiBarHeatmap.R"))

if(length(deg_list) > 0) {
# Show we can sort sub-bars
#DoHeatmap(scMultiome, features = deg_list, size = 4, angle = 0)

DoMultiBarHeatmap(scMultiome, features = scMultiome@misc$SCT_DEGs_top5, assay = 'SCT',
                  group.by='WNN_clusters_singler_annot', label = FALSE,
                  additional.group.by = c('WNN_clusters_copykat_annot', "Phase", 'ATAC_clusters',  'RNA_clusters', 'WNN_clusters'),
                  additional.group.sort.by = c('WNN_clusters'))

ggsave(paste0(out_dir, sample_id, "_top5_WNN_clusters_specific_DEGs_heatMap.png"), width = 16, height = 8.5)


##############################
## Linking peaks to top 5 DEGs
DefaultAssay(scMultiome) <- "ATAC"
bsgenome <- BSgenome.Hsapiens.UCSC.hg38     ## for GC correction

# first compute the GC content for each peak
scMultiome <- RegionStats(scMultiome, genome = bsgenome)

# link peaks to specified genes
## by computing the correlation between gene expression and accessibility at nearby peaks,
## and correcting for bias due to GC content, overall accessibility, and peak size,
## eg for top 5 DEGs per clusters

## will be saved to scMultiome@assays$ATAC@links
scMultiome <- LinkPeaks(
  object = scMultiome,
  peak.assay = "ATAC",
  expression.assay = "SCT",     ## all genes in SCT is time consuming !!!
  genes.use = scMultiome@misc$SCT_DEGs_top5
)
write.csv(Links(scMultiome), paste0(out_dir, sample_id, "_top5_WNN_clusters_specific_DEGs_linked_peaks.csv"))
}

print("The analysis of WNN clusters specific DEGs has been successfully completed!!")
}


############################################################
# identify cluster-specific DNA accessible regions (DARs)
## plus one vs other :: cluster markers
## results were add to seuratObject@misc
############################################################
{

# motif matrices from the JASPAR database
pfm <- getMatrixSet(x = JASPAR2020,
                      opts = list(collection = "CORE", species = "Homo sapiens"))

# add motif information to scMultiome
scMultiome <- AddMotifs(object = scMultiome,
                          genome = BSgenome.Hsapiens.UCSC.hg38,
                          pfm = pfm)

### DARs and motif enrichment
clusters <- levels(scMultiome)
L <- length(clusters)

DefaultAssay(scMultiome) <- "ATAC"
atac_dar <- list()
atac_dar_motif <- list()
#top_dar <- list()         ## top DARs for motif enrichment analysis
enriched_motifs <- vector()     ## top 5 per cluster
atac_dar_names <- c()

## identify DARs
for (i in 1:L)
{

  if(cluster_cnt[i] <=3 ){

    atac_dar[[i]] <- list()
    atac_dar_motif[[i]] <- list()

  } else {

  ## one vs all others: prefiltering
  atac_dar[[i]] <- FindMarkers(scMultiome, ident.1 = clusters[i], ident.2 = NULL,
                              min.pct = 0.05,                ## detected at least 5% frequency in either ident.
                              logfc.threshold = log(2),     ## at least two-fold change between the average expression of comparisons
                              min.diff.pct = 0.25,          ## Pre-filter features whose detection percentages across the two groups are similar
                              test.use = 'LR',              ##  using logistic regression
                              #latent.vars = 'peak_region_fragments'     #meta data missing  ## mitigate the effect of differential sequencing depth
                              )

  ## top DARs for motif enrichment analysis
  idx_top <- atac_dar[[i]]$p_val_adj < 0.005   ## requiring at least 10 RegionStats

  if(nrow(atac_dar[[i]]) == 0 | sum(idx_top) < 10){
    atac_dar_motif[[i]] <- list()
  } else{

  ## only examine top DARs for motif enrichment analysis
  #idx_top <- atac_dar[[i]]$p_val_adj < 0.005
  motif_res <- FindMotifs(object = scMultiome, features = rownames(atac_dar[[i]])[idx_top])

  ## add p.adjust using BH correction, which might miss for early version of FindMotifs
  p.adjust <- p.adjust(motif_res$pvalue, method = "BH")
  motif_res <- data.frame(motif_res, p.adjust)

  motif_en <- motif_res[motif_res$p.adjust < 0.05, ]
  atac_dar_motif[[i]]  <- motif_en

  ## top 5 enriched motifs
  if(nrow(motif_en) > 0){
    idx_top5 <- min(nrow(motif_en), 5)
    enriched_motifs <- c(enriched_motifs, motif_en$motif.name[1:idx_top5])
  }

  }

  }

}

names(atac_dar) <- names(atac_dar_motif) <-  clusters
enriched_motifs <- unique(enriched_motifs)        ## enriched_motif name

#################################################
### pull out top 5 motifs per cluster for heatmap
### -log10(p-adj)
{
  motif_hm <- matrix(0,  length(clusters), length(enriched_motifs))
  colnames(motif_hm) <- enriched_motifs
  rownames(motif_hm) <- clusters

  for (i in 1:L)
  {
    idx_motif <- match(colnames(motif_hm), atac_dar_motif[[i]]$motif.name)
    if(sum(is.na(idx_motif)) == length(idx_motif)){
      next;} else {
    motif_hm[i, !is.na(idx_motif)] <- -log10(atac_dar_motif[[i]]$p.adjust[idx_motif[!is.na(idx_motif)]])
    }
  }

  ## rm rows are all 0s
  motif_hm[motif_hm < -log10(0.05)] <- 0     ## ensuring motifs are significant
  idx_0 <- rowSums(motif_hm) == 0

  motif_hm <- motif_hm[!idx_0, ]

  ##

}

## unlist(lapply(atac_dar, nrow))
Misc(scMultiome@assays$ATAC, slot = "DARs") <- atac_dar
Misc(scMultiome@assays$ATAC, slot = "DARs_motif") <- atac_dar_motif
Misc(scMultiome@assays$ATAC, slot = "DARs_motif_hm") <- motif_hm

print("The analysis of WNN clusters specific DARs has been successfully completed!!")
}

##############################################
##  GRN: TF-gene gene regulatory network (GRN)
##  Pando: https://quadbio.github.io/Pando/
##############################################
if(TRUE){

#########################################
## packages installed from "dependencies"
#  data('phastConsElements20Mammals.UCSC.hg38')  ## in /Pando-1.0.0/data
## pando curated motifs in /Pando-1.0.0/data

data('motifs')      ## PFM for 1590 TFs
data('motif2tf')
data('SCREEN.ccRE.UCSC.hg38')    ## regions are contrained to SCREEN candidate cREs

##################################################
## infer GRN based on combined cluster-specific DEGs
## regions are contained to SCREEN candidate cREs
{
  ## combined and unique DEGs
  clusters <- levels(scMultiome)
  L <- length(clusters)
  all_degs <- vector()
  combined_grn_list <- list()  ## graph list for culster specific DEGs

  for(i in 1:L)
  {

    ## limited to cluster specific DEGs
    all_degs <- c(all_degs, rownames(scMultiome@assays$SCT@misc$DEGs[[i]]))
  }
  all_degs <- unique(all_degs)

  ## for testing
  ## all_degs <- rownames(scMultiome@assays$RNA@counts)[1:1000]

  #################
  ## initialization
  scMultiome_GRN <- initiate_grn(
    scMultiome,
    rna_assay = 'SCT',
    peak_assay = 'ATAC',
    regions = SCREEN.ccRE.UCSC.hg38)

  ################
  ## Scaning motif
  # Will select candidate regulatory regions near genes

  #  constrain to given gene list
  if(FALSE){
    patterning_genes <- read_tsv('patterning_genes.tsv')
    pattern_tfs <- patterning_genes %>%
      filter(type=='Transcription factor') %>%
      pull(symbol)
    motif2tf_use <- motif2tf %>%
      filter(tf %in% pattern_tfs)
    motifs_use <- motifs[unique(motif2tf_use$motif)]
    motif2tf_use
  }


  scMultiome_GRN <- find_motifs(
    scMultiome_GRN,
    pfm = motifs,
    motif_tfs = motif2tf,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )

  ################
  ## infering GRN
  registerDoParallel(12)

  scMultiome_GRN  <- infer_grn(
    scMultiome_GRN ,
    peak_to_gene_method = 'Signac',
    ## limited to cluster specific DARs as target genes
    genes = all_degs,  # Combined DEGs
    parallel = T
  )

  ## summary of inferred GRN
  n_fit <- nrow(GetNetwork(scMultiome_GRN)@fit)

if(n_fit >= 1) {
    ## modules discovery
    scMultiome_GRN <- find_modules(
      scMultiome_GRN,
      p_thresh = 0.05,
      nvar_thresh = 2,
      min_genes_per_module = 5,
      rsq_thresh = 0.05
    )

  modules <- NetworkModules(scMultiome_GRN)@meta
  write.csv(modules, paste0(out_dir, sample_id, "_combined_WNN_clusters_specific_DEGs_GRN_Modules.csv"))

  featurs <- NetworkFeatures(scMultiome_GRN)

    ##
  modules_use <- NetworkModules(scMultiome_GRN)@meta %>%
      filter(target%in%featurs, tf%in%featurs)

  n_modules <- nrow(modules_use)

    if(n_modules > 2) {
    scMultiome_GRN <- get_network_graph(scMultiome_GRN,
                                            rna_assay = "SCT",       ## using SCT normalized RNA
                                            umap_method = "none",    ## without UMAP embedding
                                            graph_name = "Combined_DEGs_GRN")

    g <- plot_network_graph(scMultiome_GRN , layout='fr') + ggtitle ("Combined_Clusters_Specific_DEGs")
    ggsave(paste0(out_dir, sample_id, "_combined_WNN_clusters_specific_DEGs_GRN.pdf"))

    ## customize the GRN
    Combined_DEGs_GRN <-  NetworkGraph(scMultiome_GRN , graph='Combined_DEGs_GRN')  ## graph object

    ## get edges or nodes by:
    ## edges <- Combined_DEGs_GRN %>% activate(edges) %>% data.frame()
    ## nodes <- Combined_DEGs_GRN %>% activate(nodes) %>% data.frame()
    ## could modify tbl_graph separately and combine them afterward using:
    ## GRN_new <- tbl_graph(nodes = new_nodes, edges = new_edges)

    edge_color = c('-1'='darkgrey', '1'='orange')  ## darkgrey for repression and orange for activate
    node_fill_color = c("0" = "lightgrey", "1" = "Brown")
    node_color = c("0" = "white", "1" = "Black")

    tf_list <- unique(modules$tf)
    target_list <- unique(modules$target)
    target_tf <- intersect(target_list, tf_list)

    Combined_DEGs_GRN_mut <- Combined_DEGs_GRN %>%
                             activate(edges) %>%
                             mutate(dir = sign(estimate)) %>%
                             activate(nodes) %>%
                             mutate(name_label = ifelse(name %in% target_tf, name, ""), ## label tf_target
                                    isTF = ifelse(name %in% tf_list, "1", "0"),
                                    isTargetTF = ifelse(name %in% target_tf, "1", "0"))

    saveRDS(Combined_DEGs_GRN_mut, paste0(out_dir, sample_id, "_combined_WNN_clusters_specific_DEGs_GRN.RDS"))
    combined_grn_list[[1]] <- Combined_DEGs_GRN_mut

    g <- ggraph(Combined_DEGs_GRN_mut, layout='fr') +
         geom_edge_diagonal(aes(color=factor(dir)), width = 0.2) + scale_edge_color_manual(values=edge_color) +
         geom_node_point(aes(size = centrality, fill = factor(isTF), color = factor(isTargetTF)), shape=21) +
         scale_fill_manual(values = node_fill_color) + scale_color_manual(values = node_color) +
         geom_node_text(aes(label = name_label, size = centrality*0.75), repel=T) +
         theme_void() + theme(legend.position = "none") + ggtitle("Combined_DEGs_GRN")

    ggsave(paste0(out_dir, sample_id, "_combined_WNN_clusters_specific_DEGs_GRN.pdf"))

    }
 }
}

#######################################################
## infer GRN based on cluster-specific DEGs per cluster
## regions are contrained to SCREEN candidate cREs
{
  #################
  ## initialization
  ####################################
  ## subsetting by WNN based clusters

  clusters <- levels(scMultiome)
  L <- length(clusters)
  cluster_cnt <- table(scMultiome$seurat_clusters)
  grn_list <- list()  ## graph list for culster specific DEGs
  k = 1

  for(i in 1:L)
  {
    print(i)

    ## limited to cluster specific DEGs
    n_deg <- nrow(scMultiome@assays$SCT@misc$DEGs[[i]])         ## number of DEGs

    ## requiring cluster with at least 30 samples
    ## based on https://leanscape.io/the-importance-of-identifying-the-right-sample-size-for-business-improvement/#:~:text=Why%20is%2030%20the%20minimum,as%20the%20sample%20size%20increases.
    #if(cluster_cnt[i] < 30 | length(dar_regions) < 2 ){
    if(length(n_deg) == 0){          ## for n_deg == Null

      next

    } else if (n_deg < 2){           ##  require at least more than 2 DEGs

      next

      } else {

      ## cluster subseting
      scMultiome_sub <- subset(scMultiome, idents = clusters[i])

      scMultiome_GRN <- initiate_grn(
        scMultiome_sub,
        rna_assay = 'SCT',
        peak_assay = 'ATAC',
        regions = SCREEN.ccRE.UCSC.hg38)

      # Will select candidate regulatory regions near genes
      scMultiome_GRN <- find_motifs(
        scMultiome_GRN,
        pfm = motifs,
        motif_tfs = motif2tf,
        genome = BSgenome.Hsapiens.UCSC.hg38
      )

      ################
      ## infering GRN
      registerDoParallel(12)

      scMultiome_GRN  <- infer_grn(
        scMultiome_GRN ,
        peak_to_gene_method = 'Signac',
        ## limited to cluster specific DARs as target genes
        genes = rownames(scMultiome@assays$SCT@misc$DEGs[[i]]),  # target genes to consider
        parallel = T
      )

      ## summary of inferred GRN
      n_fit <- nrow(GetNetwork(scMultiome_GRN)@fit)

      if(n_fit < 1) {
        next
      } else {
        ## modules discovery
        scMultiome_GRN <- find_modules(
          scMultiome_GRN,
          p_thresh = 0.1,
          nvar_thresh = 2,
          min_genes_per_module = 1,
          rsq_thresh = 0.05
        )

        modules <- NetworkModules(scMultiome_GRN)@meta
        write.csv(modules, paste0(out_dir, sample_id, "_WNN_cluster_", clusters[i], "_specific_DEGs_GRN_Modules.csv"))

        featurs <- NetworkFeatures(scMultiome_GRN)

        ##
        modules_use <- NetworkModules(scMultiome_GRN)@meta %>%
          filter(target%in%featurs, tf%in%featurs)

        n_modules <- nrow(modules_use)
        if(n_modules <= 2) {
          next
        } else {

          scMultiome_GRN <- get_network_graph(scMultiome_GRN,
                                              rna_assay = "SCT",       ## using SCT normalized RNA
                                              umap_method = "none",    ## without UMAP embedding
                                              graph_name = "Specific_DEGs_GRN")

          g <-  plot_network_graph(scMultiome_GRN, layout='fr')  + ggtitle(paste0("Cluster_",clusters[i], "_specific_DEGs"))
          ggsave(paste0(out_dir, sample_id, "_WNN_cluster_",clusters[i], "_specific_DEGs_GRN.pdf"))

          ## customize the GRN
          Specific_DEGs_GRN <-  NetworkGraph(scMultiome_GRN , graph='Specific_DEGs_GRN')  ## graph object

          edge_color = c('-1'='darkgrey', '1'='orange')  ## darkgrey for repression and orange for activate
          node_fill_color = c("0" = "lightgrey", "1" = "Brown")
          node_color = c("0" = "white", "1" = "Black")

          tf_list <- unique(modules$tf)
          target_list <- unique(modules$target)
          target_tf <- intersect(target_list, tf_list)

          Specific_DEGs_GRN_mut <- Specific_DEGs_GRN %>%
            activate(edges) %>%
            mutate(dir = sign(estimate)) %>%
            activate(nodes) %>%
            mutate(name_label = ifelse(name %in% target_tf, name, ""), ## label tf_target
                   isTF = ifelse(name %in% tf_list, "1", "0"),
                   isTargetTF = ifelse(name %in% target_tf, "1", "0"))

          g <- ggraph(Specific_DEGs_GRN_mut, layout='fr') +
            geom_edge_diagonal(aes(color=factor(dir)), width = 0.2) + scale_edge_color_manual(values=edge_color) +
            geom_node_point(aes(size = centrality, fill = factor(isTF), color = factor(isTargetTF)), shape=21) +
            scale_fill_manual(values = node_fill_color) + scale_color_manual(values = node_color) +
            geom_node_text(aes(label = name, size = centrality*0.75), repel=T) +
            theme_void() + theme(legend.position = "none") + ggtitle(paste0("Cluster_", clusters[i], "_specific_DEGs_GRN"))

          ggsave(paste0(out_dir, sample_id, "_WNN_cluster_", clusters[i], "_specific_DEGs_GRN.pdf"))

          grn_list[[k]] <- Specific_DEGs_GRN_mut
          names(grn_list)[k] <- paste0("WNN_cluster_", clusters[i], "_specific_DEGs_GRN")
          k = k + 1
        }

      }

    }

  }

  # saveRDS(grn_list, paste0(sample_id, "_cluster_specific_DEGs_GRN.RDS"))

}


#######################################################
## adding to scMutiome object
Misc(scMultiome, slot = "Combined_DEGs_GRN") <- combined_grn_list  ## list to save as tbl_graph object
Misc(scMultiome, slot = "Cluster_DEGs_GRN") <- grn_list

print("The gene regulatory network has been successfully completed!!")

}

################
## output Rdata
################
saveRDS(scMultiome, file = paste0(out_dir, sample_id, "_extended_seurat_object.RDS"))
write.csv(scMultiome@meta.data,  file = paste0(out_dir, sample_id, "_extended_meta_data.csv"))

print("The main seurat has been successfully executed !!")
