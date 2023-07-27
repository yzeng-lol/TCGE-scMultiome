args = commandArgs(trailingOnly=TRUE)

## assigning commandArgs
sample_id   = args[1]
filtered_h5 = args[2]
frag_file   = args[3]
macs2_dir   = args[4]

regCellCycle = TRUE

out_dir = paste0(getwd(), "/main_seurat/")         ## with forward slash at the end!!
macs2_dir   = "/cluster/home/yzeng/miniconda3/envs/iSHARC/bin/macs2"
anno_rds    = "/cluster/home/yzeng/snakemake/iSHARC/workflow/dependencies/EnsDb.Hsapiens.v86_2UCSC_hg38.RDS"


### for testing
if(FALSE){
#conda activate /cluster/home/yzeng/miniconda3/envs/iSHARC_extra_env/9e63b7702a5c67416a77f3ad2e11f273_
#R
sample_id   = "Lung"
filtered_h5 = "/cluster/projects/tcge/scMultiome/iSHARC_test/arc_count/Lung/outs/filtered_feature_bc_matrix.h5"
frag_file   = "/cluster/projects/tcge/scMultiome/iSHARC_test/arc_count/Lung/outs/atac_fragments.tsv.gz"
macs2_dir   = "/cluster/home/yzeng/miniconda3/envs/iSHARC/bin/macs2"
anno_rds    = "/cluster/home/yzeng/snakemake/iSHARC/workflow/dependencies/EnsDb.Hsapiens.v86_2UCSC_hg38.RDS"
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
suppressMessages(library(qlcMatrix))       ## for LinkPeaks
suppressMessages(library(future))           ## for paralleling
suppressMessages(library(biovizBase))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(rmarkdown))        ## for HTML QC report

###################################################################
## packages installed from "dependencies"
suppressMessages(library(devtools))
# devtools::load_all(pkgs_path)
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

### gene anno_gene to UCSC style failed, which might due to version of Signac and GenomeInfoDb
if(FALSE){
anno_gene <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)     ## signac

###### change NCBI chromosome format "1, 2, X, Y, MT" to UCSC format "chr1, chr2, chrX,Y,M"
## seqlevelsStyle(anno_gene) <- 'UCSC'       ## failed due to "cannot open URL ..."
anno_gene_v <- "hg38"
genome(anno_gene) <- anno_gene_v
}

## tried to load locally generated anno_gene with above codes, failed as well ..
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
VlnPlot(scMultiome, features = c("nCount_RNA","percent.mt", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 5, log = TRUE, pt.size = 0) + NoLegend()
ggsave(paste0(out_dir, sample_id, "_VlnPlot_4_QC_metrics.pdf"), width = 12, height = 4)

## checking the qc quantiles
# summary(scMultiome[[c("nCount_ATAC", "nCount_RNA","percent.mt")]])
qc_f <- data.matrix(scMultiome[[c("nCount_RNA","percent.mt", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal")]])
qc_quantiles <- apply(qc_f, 2, function(x) quantile(x, probs = seq(0, 1, 0.05), na.rm = TRUE))
qc_quantiles

## Low count for a cell indicates that it may be dead/dying or an empty droplet.
## however, High count indicates that the "cell" may in fact be a doublet (or multiplet).
if(FALSE){
scMultiome <- subset(
  x = scMultiome,
  subset = nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20 &
    nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1)
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
p3 <- DimPlot(scMultiome, reduction = "umap.atac.arc", group.by = "atac_arc_snn_res.0.8",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC_ARC")
p4 <- DimPlot(scMultiome, reduction = "wnn.umap", group.by = "wsnn_res.0.8",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

g <- p1 + p2 + p3 + p4  & theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(out_dir, sample_id, "_UMAP_plots_clustering_by_self.pdf"), width = 12, height = 12)

## different UMAPs using same WNN labels
p1 <- DimPlot(scMultiome, reduction = "umap.rna", group.by = "wsnn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(scMultiome, reduction = "umap.atac", group.by = "wsnn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(scMultiome, reduction = "umap.atac.arc", group.by = "wsnn_res.0.8",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC_ARC")
p4 <- DimPlot(scMultiome, reduction = "wnn.umap", group.by = "wsnn_res.0.8",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
g <- p1 + p2 + p3 + p4  & theme(plot.title = element_text(hjust = 0.5))  ## & NoLegend()
ggsave(paste0(out_dir, sample_id, "_UMAP_plots_clustering_by_WNN.pdf"), width = 12, height = 12)


## WNN UMAPs using diff cultering labels
p1 <- DimPlot(scMultiome, reduction = "wnn.umap", group.by = "SCT_snn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN_RNA_label")
p2 <- DimPlot(scMultiome, reduction = "wnn.umap", group.by = "ATAC_snn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN_ATAC_label")
p3 <- DimPlot(scMultiome, reduction = "wnn.umap", group.by = "atac_arc_snn_res.0.8",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN_ATAC_ARC_label")
p4 <- DimPlot(scMultiome, reduction = "wnn.umap", group.by = "wsnn_res.0.8",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
g <- p1 + p2 + p3 + p4  & theme(plot.title = element_text(hjust = 0.5))  ## & NoLegend()
ggsave(paste0(out_dir, sample_id, "_UMAP_plot_WNN_clustering_by_self.pdf"), width = 12, height = 12)

}


##########################################
# identify DEGs and DARs
## all pair-wise comaprisons
## plus one vs other :: cluster markers
## apply prefiltering to save time
## results were add to seuratObject@misc
########################################
{

clusters <- levels(scMultiome)
L <- length(clusters)

##################################
### Differentially Expressed Genes
DefaultAssay(scMultiome) <- "SCT"
sct_deg <- list()
sct_deg_names <- c()
k = 1
for (i in 1:(L-1))
{
  ## all pair-wise comparison
  if(FALSE){
  for(j in (i + 1): L)
  {
    sct_deg[[k]] <- FindMarkers(scMultiome, ident.1 = clusters[i], ident.2 = clusters[j],
                                min.pct = 0.5,                ## detected at least 50% frequency in either ident.
                                logfc.threshold = log(2),     ## at least two-fold change between the average expression of comparisons
                                min.diff.pct = 0.25,          ## Pre-filter features whose detection percentages across the two groups are similar
                                )
    sct_deg_names[k] <- paste0(clusters[i], "_vs_", clusters[j])
    k = k + 1
  }
  }

  ## one vs all others: prefiltering
  sct_deg[[k]] <- FindMarkers(scMultiome, ident.1 = clusters[i], ident.2 = NULL,
                              min.pct = 0.5,                ## detected at least 50% frequency in either ident.
                              logfc.threshold = log(2),     ## at least two-fold change between the average expression of comparisons
                              min.diff.pct = 0.25,          ## Pre-filter features whose detection percentages across the two groups are similar
                              )

  sct_deg_names[k] <- paste0(clusters[i], "_vs_all_others")
  k = k + 1
}
names(sct_deg) <-  sct_deg_names
## unlist(lapply(sct_deg, nrow))

## add to assay's Misc of seuratObject: can call by : scMultiome@misc$
Misc(scMultiome, slot = "SCT_DEGs") <- sct_deg

######################################
### Differentially  Accessible Regions
DefaultAssay(scMultiome) <- "ATAC"
atac_dar <- list()
atac_dar_names <- c()
k = 1
for (i in 1:(L-1))
{
  ## all pair-wise comparison
  if(FALSE){
  for(j in (i + 1): L)
  {
    atac_dar[[k]] <- FindMarkers(scMultiome, ident.1 = clusters[i], ident.2 = clusters[j],
                                min.pct = 0.05,                ## detected at least 5% frequency in either ident.
                                logfc.threshold = log(2),     ## at least two-fold change between the average expression of comparisons
                                min.diff.pct = 0.25,          ## Pre-filter features whose detection percentages across the two groups are similar
                                test.use = 'LR',              ##  using logistic regression
                                # latent.vars = 'peak_region_fragments'   #meta data missing   ## mitigate the effect of differential sequencing depth
    )
    atac_dar_names[k] <- paste0(clusters[i], "_vs_", clusters[j])
    k = k + 1
  }
  }

  ## one vs all others: prefiltering
  atac_dar[[k]] <- FindMarkers(scMultiome, ident.1 = clusters[i], ident.2 = NULL,
                              min.pct = 0.05,                ## detected at least 5% frequency in either ident.
                              logfc.threshold = log(2),     ## at least two-fold change between the average expression of comparisons
                              min.diff.pct = 0.25,          ## Pre-filter features whose detection percentages across the two groups are similar
                              test.use = 'LR',              ##  using logistic regression
                              #latent.vars = 'peak_region_fragments'     #meta data missing  ## mitigate the effect of differential sequencing depth
  )

  atac_dar_names[k] <- paste0(clusters[i], "_vs_all_others")
  k = k + 1
}
names(atac_dar) <-  atac_dar_names
## unlist(lapply(atac_dar, nrow))
Misc(scMultiome, slot = "ATAC_DARs") <- atac_dar

}

#######################################
## Linking peaks to genes and vice versa
## limit to DEGs or DARs ???
#######################################
if(FALSE){
DefaultAssay(scMultiome) <- "ATAC"

library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38     ## for GC correction

# first compute the GC content for each peak
scMultiome <- RegionStats(scMultiome, genome = bsgenome)

# link peaks to specified genes
## by computing the correlation between gene expression and accessibility at nearby peaks,
## and correcting for bias due to GC content, overall accessibility, and peak size,
## eg for testing
gene_names <- head(rownames(scMultiome@misc$SCT_DEGs[[1]]))

## will be saved to scMultiome@assays$ATAC@links
scMultiome <- LinkPeaks(
          object = scMultiome,
          peak.assay = "ATAC",
          expression.assay = "SCT",     ## all genes in SCT is time consuming !!!
          genes.use =   gene_names
        )

## link genes to specified regions , n
## eg testing
reg_names <- rownames(scMultiome@misc$ATAC_DARs[[1]])
closest_genes_2reg <- ClosestFeature(scMultiome, regions = reg_names)
}

##############
## output data
##############
saveRDS(scMultiome, file = paste0(out_dir, sample_id, ".RDS"))


##########################
## Generate HTML QC report
##########################
if(FALSE){
render(paste0(scr_dir, "/workflow/scripts/qc_report.Rmd"), output_dir = "main_seurat",
       params = list(readin = paste0(out_dir, sample_id, ".RDS"), sample_id = sample_id))
}



print("The main seurat has been successfully executed !!")
