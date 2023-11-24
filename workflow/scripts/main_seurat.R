args = commandArgs(trailingOnly=TRUE)

## assigning commandArgs
sample_id   = args[1]
filtered_h5 = args[2]
frag_file   = args[3]
macs2_dir   = args[4]
pipe_dir    = args[5]

regCellCycle = TRUE      ## set as true by defaulting

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
setwd("/Users/yong/OneDrive - UHN/Projects/snakemake/iSHARC_test")
scMultiome <- readRDS("Lung.RDS")
out_dir <- "/Users/yong/OneDrive - UHN/Projects/snakemake/iSHARC_test/main_seurat/"
sample_id   = "Lung"
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

DefaultAssay(scMultiome) <- "SCT"
sct_deg <- list()
sct_deg_names <- vector()
deg_list <- vector() ## top 5 degs per clusters

################
## Identify DEGs
for (i in 1:L)
{
  ## one vs all others: prefiltering
  sct_deg[[i]] <- FindMarkers(scMultiome, ident.1 = clusters[i], ident.2 = NULL,
                              min.pct = 0.5,                ## detected at least 50% frequency in either ident.
                              logfc.threshold = log(2),     ## at least two-fold change between the average expression of comparisons
                              min.diff.pct = 0.25,          ## Pre-filter features whose detection percentages across the two groups are similar
                              )

  sct_deg_names[i] <- paste0(clusters[i], "_specific")

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
names(sct_deg) <-  sct_deg_names
deg_list <- unique(deg_list)        ## remove duplicates

## add to assay's Misc of seuratObject: can call by : scMultiome@misc$
Misc(scMultiome, slot = "SCT_DEGs") <- sct_deg
Misc(scMultiome, slot = "SCT_DEGs_top5") <- deg_list       ## top5 gene names

################################################################
## draw heat map for top 5 clusters specifically expressed genes
source(paste0(pipe_dir, "/workflow/scripts/DoMutiBarHeatmap.R"))

# Show we can sort sub-bars
#DoHeatmap(scMultiome, features = deg_list, size = 4, angle = 0)
DoMultiBarHeatmap(scMultiome, features = scMultiome@misc$SCT_DEGs_top5, assay = 'SCT',
                  group.by='WNN_SingleR_anno', label = FALSE,
                  additional.group.by = c('copykat_anno', "Phase", 'ATAC_snn_res.0.8',  'SCT_snn_res.0.8', 'wsnn_res.0.8'),
                  additional.group.sort.by = c('wsnn_res.0.8'))
ggsave(paste0(out_dir, sample_id, "_top5_DEGs_heatMap.png"), width = 16, height = 8.5)

##############################
## Linking peaks to top 5 DEGs
{
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


##########################################
# identify cluster-specific genes (DEGs)
## plus one vs other :: cluster markers
## results were add to seuratObject@misc
########################################
{
clusters <- levels(scMultiome)
L <- length(clusters)

DefaultAssay(scMultiome) <- "ATAC"
atac_dar <- list()
top_dar <- list()         ## top DARs for motif enrichment analysis
atac_dar_names <- c()


## identify DARs
for (i in 1:L)
{
  ## one vs all others: prefiltering
  atac_dar[[i]] <- FindMarkers(scMultiome, ident.1 = clusters[i], ident.2 = NULL,
                              min.pct = 0.05,                ## detected at least 5% frequency in either ident.
                              logfc.threshold = log(2),     ## at least two-fold change between the average expression of comparisons
                              min.diff.pct = 0.25,          ## Pre-filter features whose detection percentages across the two groups are similar
                              test.use = 'LR',              ##  using logistic regression
                              #latent.vars = 'peak_region_fragments'     #meta data missing  ## mitigate the effect of differential sequencing depth
                              )

  top_dar[[i]]  <- rownames(atac_dar[[i]][atac_dar[[i]]$p_val < 0.005, ])
  atac_dar_names[i] <- paste0(clusters[i], "_specific")
}
names(atac_dar) <- names(top_dar) <-  atac_dar_names
## unlist(lapply(atac_dar, nrow))
Misc(scMultiome, slot = "ATAC_DARs") <- atac_dar
Misc(scMultiome, slot = "ATAC_DARs_top") <- top_dar

## motif enrichment
## under tuning
if(FLASE){
# motif matrices from the JASPAR database
pfm <- getMatrixSet(x = JASPAR2020,
         opts = list(collection = "CORE", species = "Homo sapiens"))

# add motif information to scMultiome
scMultiome <- AddMotifs(object = scMultiome,
                genome = BSgenome.Hsapiens.UCSC.hg38,
                pfm = pfm)

## motif enrichment
enriched.motifs <- FindMotifs(object = scMultiome, features = top.da.peak)

}

}


################
## output Rdata
################
saveRDS(scMultiome, file = paste0(out_dir, sample_id, ".RDS"))

write.csv(scMultiome@meta.data,  file = paste0(out_dir, sample_id, "_metaData.csv"))
print("The main seurat has been successfully executed !!")


##############################################
##  GRN: TF-gene gene regulatory network (GRN)
## tools:
##  FigR: https://buenrostrolab.github.io/FigR/
##  Pando: https://quadbio.github.io/Pando/
##############################################
{


}
