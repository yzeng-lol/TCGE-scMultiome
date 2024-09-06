args = commandArgs(trailingOnly=TRUE)

## assigning commandArgs
sample_id   = args[1]
filtered_h5 = args[2]
frag_file   = args[3]
macs2_dir   = args[4]
pipe_dir    = args[5]

regCellCycle = TRUE      ## set as true by defaulting
## cellcycle analysis might be problematic due to low number of cells
## or ERROR in `cut_number()`:
##     ! Insufficient data values to produce 24 bins.


out_dir = paste0(getwd(), "/main_seurat/")         ## with forward slash at the end!!

#macs2_dir   = "/cluster/home/yzeng/miniconda3/envs/iSHARC/bin/macs2"      ## NO forward slash at the end
#pipe_dir ="/cluster/home/yzeng/snakemake/iSHARC"                          ## NO forward slash at the end

###############################
### for testing and fine-tuning
###############################
if(FALSE){
#########
## on H4H
# cd  /cluster/projects/tcge/scMultiome/iSHARC_test/main_seurat
# conda activate  /cluster/home/yzeng/miniconda3/envs/iSHARC_extra_env/ebcd8f410d62ed43e71d1cebab0b3296_    ## new

## conda activate /cluster/home/yzeng/miniconda3/envs/iSHARC_extra_env/672743604e74f3edd5a4eb4b70c92872_  ## old
#R
sample_id   = "Lung"
filtered_h5 = "/cluster/projects/tcge/scMultiome/iSHARC_test/arc_count/Lung/outs/filtered_feature_bc_matrix.h5"
frag_file   = "/cluster/projects/tcge/scMultiome/iSHARC_test/arc_count/Lung/outs/atac_fragments.tsv.gz"
macs2_dir   = "/cluster/home/yzeng/miniconda3/envs/iSHARC/bin/macs2"
pipe_dir = "/cluster/home/yzeng/snakemake/iSHARC"
anno_rds    = "/cluster/home/yzeng/snakemake/iSHARC/workflow/dependencies/EnsDb.Hsapiens.v86_2UCSC_hg38.RDS"

###########################
## testing locally with RDS
rm(list = ls())
setwd("/Users/yong/OneDrive - UHN/Projects/snakemake/iSHARC_test/main_seurat")
# scMultiome <- readRDS("Lung.RDS")
# sample_id   = "Lung"
scMultiome <- readRDS("PDAC_PDA_87784.RDS")
sample_id   = "PDAC_PDA_87784"

out_dir <- "/Users/yong/OneDrive - UHN/Projects/snakemake/iSHARC_test/main_seurat/"
pipe_dir  = "/Users/yong/OneDrive - UHN/Projects/snakemake/iSHARC"

}

##################################################
## The main workflow for scMultiome data analysis
## loading required packages
##################################################
{
suppressMessages(library(hdf5r))           ## to read HDF5 files
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(SingleR))         ## auto annotation
suppressMessages(library(celldex))         ## annotation reference
suppressMessages(library(qlcMatrix))       ## for LinkPeaks
suppressMessages(library(future))          ## for paralleling
suppressMessages(library(biovizBase))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(JASPAR2020))       ## motif enrichment analysis
suppressMessages(library(TFBSTools))

suppressMessages(library(rmarkdown))        ## for HTML QC report

#########################################
## packages installed from "dependencies"
suppressMessages(library(devtools))

## install "dlm", which is required for copykat
if(!require("dlm"))  install.packages(paste0(pipe_dir, "/workflow/dependencies/dlm_1.1-6.tar.gz"))
library(dlm)

## install "ComplexHeatmap", which is required for ArchR
if(!require("ComplexHeatmap"))  install.packages(paste0(pipe_dir, "/workflow/dependencies/ComplexHeatmap_2.18.0.tar.gz"))
library(ComplexHeatmap)

## load copykat
if(!require("copykat"))  devtools::load_all(paste0(pipe_dir, "/workflow/dependencies/copykat"))
library(copykat)

## load ArchR
if(!require("ArchR"))  devtools::load_all(paste0(pipe_dir, "/workflow/dependencies/ArchR"))
library(ArchR)

}


###########################################################################
## load scMultiome data processed by cellranger-arc count
## scMultiome data with snRNA-seq and scATAC-seq data
## You can download the test dataset from the 10x Genomics website here.
## Please make sure following files were successfully generated:
##   1) Filtered feature barcode matrix (HDF5)
##   2) ATAC Per fragment information file (TSV.GZ)
##   3) ATAC Per fragment information index (TSV.GZ index)
## the 10x hdf5 file contains both data types.
###########################################################################

#######################################
## inputs and parameters initialization
#######################################
{
## hdf 5 file after joint cell calling
inputdata <- Read10X_h5(filtered_h5)

## ATAC-seq fragments file
frag.file <- frag_file

### gene anno_gene to UCSC style failed, which might be due to version of Signac and GenomeInfoDb
## tried to load locally generated anno_gene with above codes, failed as well ..
if(FALSE){
anno_gene <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)     ## signac

###### change NCBI chromosome format "1, 2, X, Y, MT" to UCSC format "chr1, chr2, chrX,Y,M"
## seqlevelsStyle(anno_gene) <- 'UCSC'       ## failed due to "cannot open URL ..."
anno_gene_v <- "hg38"
genome(anno_gene) <- anno_gene_v
}

anno_rds <- paste0(pipe_dir, "/workflow/dependencies/EnsDb.Hsapiens.v86_2UCSC_hg38.RDS")
anno_gene <- readRDS(anno_rds)
genome_info <- seqinfo(anno_gene)

## Parallelization using plan for both seurat and signac
## plan("multiprocess", workers = 4)
}

######################################
## add RNA and ATAC to a Seurat object
######################################
{
## sparse Matrix of class "dgCMatrix"
rna_counts <- inputdata$`Gene Expression`
atac_counts <- inputdata$Peaks                  ## peaks called by cellranger-arc
rm(inputdata)

# Create Seurat object and add chrM percentage to meta data
scMultiome <- CreateSeuratObject(counts = rna_counts)
scMultiome[["percent.mt"]] <- PercentageFeatureSet(scMultiome, pattern = "^MT-")
rm(rna_counts)

################################################
# Now add in the ATAC-seq data by cellranger-arc
# Only use peaks in standard chromosomes: chr1-22 + chrX + chrY
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

## add ATAC assay by cellranger-arc
## if genome is assigned as "hg38" or "GRCh38", the internet is required for download!!!
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = genome_info,
  fragments = frag.file,
  min.cells = 0,
  min.features = 0,
  annotation = anno_gene     ## gene annotation doesn't work for h4h so far
)

scMultiome[["atac_arc"]] <- chrom_assay
#rm(atac_counts, chrom_assay)

########################
## call peaks using MAC2
DefaultAssay(scMultiome) <- "atac_arc"
## need to specify the outdir

peaks <- CallPeaks(scMultiome, macs2.path = macs2_dir,
                  outdir = out_dir,  fragment.tempdir = out_dir)   ## add more parem

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

## quantify counts in each peak
## counting fragments (as ArchR)rather than cut sites by cellranger-arc
macs2_counts <- FeatureMatrix(
  fragments = Fragments(scMultiome),
  features = peaks,
  sep = c("-", "-"),                 ##  c(":", "-") will lead to invalid character indexing
  cells = colnames(scMultiome)
)

## using the same granges format
tmp <- unlist(strsplit(rownames(macs2_counts), "-"))
LL <- length(tmp)
n_name <- paste0(tmp[seq(1, LL, 3)], ":", tmp[seq(2, LL, 3)], "-", tmp[seq(3, LL, 3)])
rownames(macs2_counts) <- n_name


# create the ATAC assay using the MACS2 peak set and add it to the Seurat object
scMultiome[["ATAC"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  sep = c(":", "-"),
  genome = genome_info,
  fragments = frag.file,
  min.cells = 0,
  min.features = -1,                ##  set as to negative to ensure same number of cells!!
  annotation = anno_gene
)

rm(macs2_counts)

####################################
#####  add ATAC specific QC metircs
DefaultAssay(scMultiome) <- "ATAC"

# compute nucleosome signal score per cell
scMultiome <- NucleosomeSignal(object = scMultiome)

# compute TSS enrichment score per cell
scMultiome <- TSSEnrichment(object = scMultiome, fast = FALSE)

## add blacklist ratio and fraction of reads in peaks
## not availabe for the testing data
## scMultiome$pct_reads_in_peaks <- scMultiome$peak_region_fragments / scMultiome$passed_filters * 100
## scMultiome$blacklist_ratio <- scMultiome$blacklist_region_fragments / scMultiome$peak_region_fragments
}


#############################
## joint cell filtering again
#############################
{
## need to be customized based on distribution
g <- VlnPlot(scMultiome, features = c("nCount_RNA","percent.mt", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 5, log = TRUE, pt.size = 0, group.by = "orig.ident") + NoLegend()
g <- g & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(paste0(out_dir, sample_id, "_VlnPlot_4_QC_metrics.pdf"), width = 12, height = 4)

## Low count for a cell indicates that it may be dead/dying or an empty droplet.
## however, High count indicates that the "cell" may in fact be a doublet (or multiplet).
if(FALSE){

  ## checking the qc quantiles
  # summary(scMultiome[[c("nCount_ATAC", "nCount_RNA","percent.mt")]])
  qc_df <- data.matrix(scMultiome [[c("nCount_RNA","percent.mt", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal")]])
  probs_s <- c(0, 0.025, 0.25, 0.50, 0.75, 0.975, 1)
  qc_percentiles <- apply(qc_df, 2, function(x) quantile(x, probs = probs_s, na.rm = TRUE))
  #qc_quantiles
  qc_percentiles  <- data.frame(qc_percentiles)

  idx_RNA <- qc_df$nCount_RNA > max(1000, qc_percentiles$nCount_RNA[2]) &
             qc_df$nCount_RNA < min(25000, qc_percentiles$nCount_RNA[6])
  idx_mt <- qc_df$percent.mt < min(20, qc_percentiles$percent.mt[6])
  idx_ATAC <- qc_df$nCount_ATAC > max(5000, qc_percentiles$nCount_ATAC[2]) &
              qc_df$nCount_ATAC < min(70000, qc_percentiles$nCount_ATAC[6])
  idx_tss <- qc_df$TSS.enrichment > max(1, qc_percentiles$TSS.enrichment[2])
  idx_ns <- qc_df$nucleosome_signal < min(2, qc_percentiles$nucleosome_signal[6])


  scMultiome <- subset(
                  x = scMultiome,
                  subset = idx_RNA & idx_mt & idx_ATAC & idx_tss & idx_ns)
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
## there is version 2 : https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
## need to shift to version 2

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

if(regCellCycle)
{
  ## by three steps
  # 1) normalize data with SCTransform()
  scMultiome <- SCTransform(
    scMultiome,
    assay = 'RNA',
    new.assay.name = 'SCT',
    vars.to.regress = c('percent.mt'), #, 'nFeature_RNA', 'nCount_RNA')
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
    vars.to.regress = c('percent.mt','S.Score', 'G2M.Score'),
    #vst.flavor = "v2", verbose = FALSE
  )
} else {
  scMultiome <- SCTransform(
    scMultiome,
    assay = 'RNA',
    new.assay.name = 'SCT',
    vars.to.regress = c('percent.mt'),
    #vst.flavor = "v2", verbose = FALSE
  )
}


## !! default assay has been shifted to SCT
## dimensional reduction and clustering
scMultiome <- RunPCA(scMultiome) %>%
              RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

## RNA alone based clustering
scMultiome  <- FindNeighbors(scMultiome, reduction = "pca", dims = 1:50) %>%
               FindClusters()

}


###########################################################################
## ATAC alone Normalization, Dimension Reduction, Clustering and Embedding
###########################################################################
{
# exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(scMultiome) <- "atac_arc"
scMultiome <- RunTFIDF(scMultiome)            # TF-IDF normalization

scMultiome <- FindTopFeatures(scMultiome, min.cutoff = 'q0')    ## q0 -> 100% ; 95 -> 95% cells
scMultiome <- RunSVD(scMultiome)
scMultiome <- RunUMAP(scMultiome, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac.arc", reduction.key = "atacArcUMAP_")

## ATAC alone based clustering
scMultiome  <- FindNeighbors(scMultiome, reduction = "lsi", dims = 2:30) %>%
                FindClusters(algorithm = 3)

# exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(scMultiome) <- "ATAC"
scMultiome <- RunTFIDF(scMultiome)            # TF-IDF normalization

scMultiome <- FindTopFeatures(scMultiome, min.cutoff = 'q0')    ## q0 -> 100% ; 95 -> 95% cells
scMultiome <- RunSVD(scMultiome)
scMultiome <- RunUMAP(scMultiome, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

## ATAC alone based clustering
scMultiome  <- FindNeighbors(scMultiome, reduction = "lsi", dims = 2:30) %>%
               FindClusters(algorithm = 3)


## note the reduction lsi will be overwritten using assay ATAC !!!
}


###################################
## integrate RNA and ATAC using WNN
##################################
{
scMultiome <- FindMultiModalNeighbors(scMultiome, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
scMultiome <- RunUMAP(scMultiome, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
scMultiome <- FindClusters(scMultiome, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

## cell identifies will be assigned as WNN clusters res: scMultiome@active.ident; scMultiome$seurat_clusters

## DimPlo
## using clusters by each individual data
p1 <- DimPlot(scMultiome, reduction = "umap.rna",  group.by = "SCT_snn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(scMultiome, reduction = "umap.atac", group.by = "ATAC_snn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
#p3 <- DimPlot(scMultiome, reduction = "umap.atac.arc", group.by = "atac_arc_snn_res.0.8",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC_ARC")
p4 <- DimPlot(scMultiome, reduction = "wnn.umap", group.by = "wsnn_res.0.8",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

g <- p1 + p2 + p4  & theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(out_dir, sample_id, "_clustering_UMAPs.pdf"), width = 12, height = 4)

}

##############################################################
## auto annotation using publicly available reference datasets
##############################################################
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
scMultiome[["WNN_SingleR_anno"]] <- expr_anno$labels[idx_m]

}

#####################################################
## Distinguish the tumor cells from the normal cells
#####################################################


if(TRUE){
## Using copyKAT to predicts tumor and normal cells
## RNA-seq data based

setwd(out_dir)                                     ## ensure output copykat results to desired folder
expr_raw <- as.matrix(scMultiome@assays$RNA@counts)

copykat_res <- copykat(rawmat = expr_raw, sam.name = sample_id , id.type = "S", ngene.chr = 5, win.size = 25,
                        KS.cut = 0.1,  distance = "euclidean", norm.cell.names = "", output.seg = "FLASE",
                        plot.genes = "TRUE", genome = "hg20", n.cores = 1)

## adding copykat predict labels to the scMultiome metatable
idx_s <- match(rownames(scMultiome@meta.data),copykat_res$prediction$cell.names)

cell_type <- rep("NA", length(idx_s))
cell_type[!is.na(idx_s)] <- copykat_res$prediction$copykat.pred[idx_s[!is.na(idx_s)]]

scMultiome[["copykat_anno"]] <- cell_type

setwd("../")  ## cd back to workdir
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
deg_list <- vector() ## top 5 degs per clusters

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



names(sct_deg) <-  clusters         ##sct_deg_names
deg_list <- unique(deg_list)        ## remove duplicates

## add to assay's Misc of seuratObject: can call by : scMultiome@misc$
Misc(scMultiome@assays$SCT, slot = "DEGs") <- sct_deg
Misc(scMultiome@assays$SCT, slot = "DEGs_top5") <- deg_list       ## top5 gene names

################################################################
## draw heat map for top 5 clusters specifically expressed genes
source(paste0(pipe_dir, "/workflow/scripts/DoMutiBarHeatmap.R"))

if(length(deg_list) > 0) {
# Show we can sort sub-bars
#DoHeatmap(scMultiome, features = deg_list, size = 4, angle = 0)
DoMultiBarHeatmap(scMultiome, features = scMultiome@misc$SCT_DEGs_top5, assay = 'SCT',
                  group.by='WNN_SingleR_anno', label = FALSE,
                  additional.group.by = c('copykat_anno', "Phase", 'ATAC_snn_res.0.8',  'SCT_snn_res.0.8', 'wsnn_res.0.8'),
                  additional.group.sort.by = c('wsnn_res.0.8'))
ggsave(paste0(out_dir, sample_id, "_top5_DEGs_heatMap.png"), width = 16, height = 8.5)

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
write.csv(Links(scMultiome), paste0(out_dir, sample_id, "_top5_DEGs_linked_peaks.csv"))
}


}


############################################################
# identify cluster-specific genes regulatory regions (DARs)
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

}

##############################################
##  GRN: TF-gene gene regulatory network (GRN)
##  Pando: https://quadbio.github.io/Pando/
##############################################
{

##################################################
## loading required packages
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
{
Misc(scMultiome, slot = "Combined_DEGs_GRN") <- combined_grn_list  ## list to save as tbl_graph object
Misc(scMultiome, slot = "Cluster_DEGs_GRN") <- grn_list
#saveRDS(scMultiome, file = paste0(out_dir, sample_id, ".RDS"))
}

}

################
## output Rdata
################
saveRDS(scMultiome, file = paste0(out_dir, sample_id, ".RDS"))
write.csv(scMultiome@meta.data,  file = paste0(out_dir, sample_id, "_metaData.csv"))

print("The main seurat has been successfully executed !!")
