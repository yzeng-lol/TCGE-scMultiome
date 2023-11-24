
# Dependencies
## Installation of packages from local
The packages are unavailable/outdated through conda/bioconda/r were cloned. Remember to remove the .git folder from cloned packages to avoid committing conflict.  

1. [copyKAT](https://github.com/navinlabcode/copykat) is employed to Distinguish cancer cells from non-malignant cell types, as well as different tumor subclones.
2. [dlm](https://cran.r-project.org/web/packages/dlm/index.html) is required for copyKAT
3. [ArchR](https://github.com/GreenleafLab/ArchR) is primarily employed for ...
4. [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html) is required for ArchR

## Reference datasets
1. BlueprintEncodeData.RDS is used by SingleR for auto cell annotation
2. EnsDb.Hsapiens.v86_2UCSC_hg38.RDS: Gene annotation in UCSC format for scATAC-seq analysis  
