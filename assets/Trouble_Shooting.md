## libmamba Aborting While Installing Extra Environments
If you encounter the error message "libmamba Non-conda folder exists at prefix" while installing the extra environments, this issue is likely related to the Mamba version (v2) within the ISHARC environment. For more details, refer to the issue [here](https://github.com/jiarong/VirSorter2/issues/222).

You can resolve this by downgrading Mamba to version 1 with the following command:

```bash
conda activate iSHARC
mamba --version
mamba install -c conda-forge mamba=1.5.11
```
