
# Dependencies
## Installation of local packages
The following packages, which are unavailable or outdated in Conda, were cloned from their respective GitHub repositories. If the cloned repos are not tarred source code, be sure to remove the .git folder from cloned repo to avoid potential conflict during commits.  

1. [KEGG.db](https://github.com/YuLab-SMU/createKEGGdb/tree/master): Required for running clusterProfiler in offline mode without internet access.
2. [copyKAT](https://github.com/navinlabcode/copykat): Used to distinguish cancer cells from non-cancer cell types.
3. [dlm](https://cran.r-project.org/web/packages/dlm/index.html): A dependency required for copyKAT.
4. [Pando](https://github.com/quadbio/Pando): Used for gene regulatory network (GRN) analysis.


## Reference datasets
1. BlueprintEncodeData.RDS is used by SingleR for auto cell annotation
2. EnsDb.Hsapiens.v86_2UCSC_hg38.RDS: Gene annotation in UCSC format for scATAC-seq analysis
