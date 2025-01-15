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


## Error: Exceeding the Maximum Memory for Parallelization
If you encounter the error:

```r
Error in getGlobalsAndPackages(expr, envir = envir, globals = TRUE) :
  The total size of the X globals that need to be exported for the future expression ('FUN()') is XXX GiB. This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize'). The X largest globals are ...
```
This issue occurs because the size of global variables are larger than the default memory limit for parallelization. To get around this, increase the memory limit by adjusting the `future_globals_maxSize` in your configuration file to the required memory size (XXX GiB) as suggested in the [Seurat documentation](https://satijalab.org/seurat/archive/v3.0/future_vignette.html).
