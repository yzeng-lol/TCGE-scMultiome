args = commandArgs(trailingOnly=TRUE)

## assigning commandArgs
sample_id   = args[1]
sample_rds  = args[2]

out_dir = paste0(getwd(), "/GRN/")         ## with forward slash at the end!!
                   
###############################
### for testing and fine-tuning
###############################
if(FALSE){
#########
## on H4H
# cd  /cluster/projects/tcge/scMultiome/iSHARC_test/main_seurat
# conda activate  /cluster/home/yzeng/miniconda3/envs/iSHARC_extra_env/ebcd8f410d62ed43e71d1cebab0b3296_    ## new

#R
sample_id   =  "PDAC_PDA_87784"
sample_rds  =  "PDAC_PDA_87784.RDS"

out_dir <- "/cluster/projects/tcge/scMultiome/iSHARC_test/main_seurat/"
scMultiome <- readRDS(sample_rds)


## update Matrix version 
#install.packages("/cluster/home/yzeng/snakemake/iSHARC/workflow/dependencies/Matrix_1.5-0.tar.gz")

if(!require("Pando"))  devtools::load_all("/cluster/home/yzeng/snakemake/iSHARC/workflow/dependencies/Pando-1.0.0")
library("Pando")



###########################
## testing locally with RDS
rm(list = ls())
setwd("/Users/yong/OneDrive - UHN/Projects/snakemake/iSHARC_test/main_seurat")
# scMultiome <- readRDS("Lung.RDS")
# sample_id   = "Lung"

sample_id   =  "PDAC_PDA_87784"
sample_rds  =  "PDAC_PDA_87784.RDS"
scMultiome <- readRDS(sample_rds)

out_dir <- "/Users/yong/OneDrive - UHN/Projects/snakemake/iSHARC_test/main_seurat/"

## load Pando locally 
suppressMessages(library(devtools))

if(!require("Pando"))  devtools::load_all("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/snakemake/iSHARC_test/PRCA/Pando-1.0.0")
library(Pando)

}

##################################################
## The main workflow for scMultiome data analysis
## loading required packages
##################################################
{
suppressMessages(library(tidyverse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Signac)) 
suppressMessages(library(ggraph))
suppressMessages(library(dplyr))
suppressMessages(library(tidygraph))


#########################################
## packages installed from "dependencies"
suppressMessages(library(devtools))

## load Pando
if(!require("Pando"))  devtools::load_all(paste0(pipe_dir, "/workflow/dependencies/Pando-1.0.0"))
library(Pando)

#data('phastConsElements20Mammals.UCSC.hg38')  ## in /Pando-1.0.0/data
## pando curated motifs in /Pando-1.0.0/data
data('motifs')      ## PFM for 1590 TFs
data('motif2tf')
data('SCREEN.ccRE.UCSC.hg38')    ## regions are contrained to SCREEN candidate cREs

}

##################################################
## infer GRN based on combined cluster-specific DEGs 
## regions are contained to SCREEN candidate cREs
###################################################
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
  
  #################
  ## initialization   
  scMultiome_GRN <- initiate_grn(
    scMultiome,
    rna_assay = 'SCT',
    peak_assay = 'ATAC', 
    regions = SCREEN.ccRE.UCSC.hg38                                  ## 
    #regions = phastConsElements20Mammals.UCSC.hg38    ## pando 2244053 ranges
  )
  
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
  scMultiome_GRN  <- infer_grn(
    scMultiome_GRN ,
    peak_to_gene_method = 'Signac',
    ## limited to cluster specific DARs as target genes
    genes = all_degs,  # Combined DEGs
    parallel = F
  )
  
  ## modules discovery, apply more stringent cutoff
  scMultiome_GRN <- find_modules(
    scMultiome_GRN, 
    p_thresh = 0.05,
    nvar_thresh = 2, 
    min_genes_per_module = 5, 
    rsq_thresh = 0.05
  )
  
  modules <- NetworkModules(scMultiome_GRN)@meta
  write.csv(modules, paste0(sample_id, "_combined_DEGs_GRN_Modules.csv"))
  
  ## plot network costomized plots 
  scMultiome_GRN <- get_network_graph(scMultiome_GRN, 
                                      rna_assay = "SCT",       ## using SCT normalized RNA
                                      umap_method = "none",    ## without UMAP embedding
                                      graph_name = "Combined_DEGs_GRN")
  
  #g <- plot_network_graph(scMultiome_GRN , layout='fr') + ggtitle ("Combined_Clusters_Specific_DEGs") 
  #ggsave(paste0("Combined_Clusters_Specific_DEGs_GRN.pdf"))
  
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
  
  # saveRDS(Combined_DEGs_GRN_mut, paste0(sample_id, "_combined_DEGs_GRN.RDS")) 
  combined_grn_list[[1]] <- Combined_DEGs_GRN_mut
  
  g <- ggraph(Combined_DEGs_GRN_mut, layout='fr') +
       geom_edge_diagonal(aes(color=factor(dir)), width = 0.2) + scale_edge_color_manual(values=edge_color) +
       geom_node_point(aes(size = centrality, fill = factor(isTF), color = factor(isTargetTF)), shape=21) + 
       scale_fill_manual(values = node_fill_color) + scale_color_manual(values = node_color) +
       geom_node_text(aes(label = name_label, size = centrality*0.75), repel=T) +
       theme_void() + theme(legend.position = "none") + ggtitle("Combined_DEGs_GRN")
  
  ggsave(paste0(sample_id, "_combined_DEGs_GRN.pdf"))
 
  
}

#######################################################
## infer GRN based on cluster-specific DEGs per cluster
## regions are contrained to SCREEN candidate cREs
###################################################
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
    if(n_deg < 2){
      
      next
      
    } else {
      
      ## cluster subseting    
      scMultiome_sub <- subset(scMultiome, idents = clusters[i])
      
      scMultiome_GRN <- initiate_grn(
        scMultiome_sub,
        rna_assay = 'SCT',
        peak_assay = 'ATAC', 
        regions = SCREEN.ccRE.UCSC.hg38                                  ## 
        #regions = phastConsElements20Mammals.UCSC.hg38    ## pando 2244053 ranges
      )
      
      # Will select candidate regulatory regions near genes
      scMultiome_GRN <- find_motifs(
        scMultiome_GRN, 
        pfm = motifs, 
        motif_tfs = motif2tf,
        genome = BSgenome.Hsapiens.UCSC.hg38
      )
      
      ################
      ## infering GRN
      #library(doParallel)
      #registerDoParallel(4)
      
      scMultiome_GRN  <- infer_grn(
        scMultiome_GRN ,
        peak_to_gene_method = 'Signac',
        ## limited to cluster specific DARs as target genes
        genes = rownames(scMultiome@assays$SCT@misc$DEGs[[i]]),  # target genes to consider 
        parallel = F
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
        write.csv(modules, paste0(sample_id, "_cluster_", clusters[i], "_specific_DEGs_GRN_Modules.csv"))
        
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
          
          #g <-  plot_network_graph(scMultiome_GRN, layout='fr')  + ggtitle(paste0("Cluster_",clusters[i], "_specific_DEGs")) 
          #ggsave(paste0("Cluster_",clusters[i], "_specific_DEGs_GRN.pdf"))
          
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
          
          ggsave(paste0(sample_id, "_cluster_", clusters[i], "_specific_DEGs_GRN.pdf"))
          
          grn_list[[k]] <- Specific_DEGs_GRN_mut
          names(grn_list)[k] <- paste0("Cluster_", clusters[i], "_specific_DEGs_GRN")
          k = k + 1                         
        }
        
      }
      
    }
    
  }
  
  # saveRDS(grn_list, paste0(sample_id, "_cluster_specific_DEGs_GRN.RDS")) 
  
}

#######################################################
## adding to scMutiome object 
#######################################################
{
Misc(scMultiome, slot = "Combined_DEGs_GRN") <- combined_grn_list  ## list to save as tbl_graph object
Misc(scMultiome, slot = "Cluster_DEGs_GRN") <- grn_list
saveRDS(scMultiome, file = paste0(out_dir, sample_id, ".RDS"))
}

########################################################
## infer GRN based on cluster-specific DEGs and DRAs
## Direct linkage between DRAs and DEGs
####################################################
if(FALSE){
  #################
  ## initialization 
  
  #data('phastConsElements20Mammals.UCSC.hg38')  ## in /Pando-1.0.0/data
  ## pando curated motifs in /Pando-1.0.0/data
  data('motifs')
  data('motif2tf')
  
  ####################################
  ## subsetting by WNN based clusters
  
  clusters <- levels(scMultiome)
  L <- length(clusters)
  cluster_cnt <- table(scMultiome$seurat_clusters)
  
  for(i in 1:L)
  {
    print(i)
    
    ## limited to cluster specific DARs, at least 2 DARs
    n_dar <- nrow(scMultiome@assays$ATAC@misc$DARs[[i]])    ## number of DARs
    n_deg <- nrow(scMultiome@assays$SCT@misc$DEGs[[i]])         ## number of DEGs
    
    ## requiring cluster with at least 30 samples
    ## based on https://leanscape.io/the-importance-of-identifying-the-right-sample-size-for-business-improvement/#:~:text=Why%20is%2030%20the%20minimum,as%20the%20sample%20size%20increases.
    #if(cluster_cnt[i] < 30 | length(dar_regions) < 2 ){
    if(n_dar < 2 | n_deg < 2){
      
      next
      
    } else {
      
      dar_regions <- StringToGRanges(rownames(scMultiome@assays$ATAC@misc$DARs[[i]]))
      
      ## cluster subseting    
      scMultiome_sub <- subset(scMultiome, idents = clusters[i])
      
      scMultiome_GRN <- initiate_grn(
        scMultiome_sub,
        rna_assay = 'SCT',
        peak_assay = 'ATAC',
        regions = dar_regions
        #regions = phastConsElements20Mammals.UCSC.hg38    ## pando 2244053 ranges
      )
      
      
      ################
      ## Scaning motif
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
      
      
      # Will select candidate regulatory regions near genes
      scMultiome_GRN <- find_motifs(
        scMultiome_GRN, 
        pfm = motifs, 
        motif_tfs = motif2tf,
        genome = BSgenome.Hsapiens.UCSC.hg38
      )
      
      ################
      ## infering GRN
      #library(doParallel)
      #registerDoParallel(4)
      scMultiome_GRN  <- infer_grn(
        scMultiome_GRN ,
        peak_to_gene_method = 'Signac',
        ## limited to cluster specific DARs as target genes
        genes = rownames(scMultiome@assays$SCT@misc$DEGs[[i]]),
        parallel = F
      )
      ## summary of inferred GRN
      n_fit <- nrow(GetNetwork(scMultiome_GRN)@fit)
      
      if(n_fit < 1) {
        next
      } else {
        
        print(i)
        
        ## modules discovery 
        scMultiome_GRN <- find_modules(
          scMultiome_GRN, 
          p_thresh = 0.1,
          nvar_thresh = 2, 
          min_genes_per_module = 1, 
          rsq_thresh = 0.05
        )
        
        #modules <- NetworkModules(scMultiome_GRN)@meta
        featurs <- NetworkFeatures(scMultiome_GRN)
        
        modules_use <- NetworkModules(scMultiome_GRN)@meta %>%
          filter(target%in%featurs, tf%in%featurs)
        
        n_modules <- nrow(modules_use)
        if(n_modules == 0) {
          next
        } else {
          
          scMultiome_GRN <- get_network_graph(scMultiome_GRN, 
                                              rna_assay = "SCT")
          
          g <- plot_network_graph(scMultiome_GRN )
          ggsave(paste0("Cluster_", i, "_specific_GRN.pdf"))
        }
        
      }
      
      
      
    }
    
  }
  
}
