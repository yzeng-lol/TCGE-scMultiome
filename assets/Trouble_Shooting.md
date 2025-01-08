## libmamba Aborting While Installing Extra Environments
If you encounter the error message "libmamba Non-conda folder exists at prefix" while installing the extra environments, this issue is likely related to the Mamba version (v2) within the ISHARC environment. For more details, refer to the issue [here](https://github.com/jiarong/VirSorter2/issues/222).

You can resolve this by downgrading Mamba to version 1 with the following command:

```bash
conda activate iSHARC
mamba --version
mamba install -c conda-forge mamba=1.5.11
```

## Error During Integration of Multiple Samples Using Seurat
If you encounter the error:

```r
Error in idx[i, ] <- res[[i]][[1]] :
number of items to replace is not a multiple of replacement length
```
This issue occurs during the integration of multiple samples with Seurat's FindIntegrationAnchors and IntegrateEmbeddings functions. The error typically arises when the number of anchor cells is less than the k.weight parameter, which specifies the number of neighbours to use when weighting anchors. You might need to remove the samples with low cell number. For more solutions, refer to the issue discussion [here](https://github.com/satijalab/seurat/issues/6341).
