# iSHARC: **I**ntegrating **s**cMultiome data for **h**eterogeneity **a**nd **r**egulatory analysis in **c**ancer (v1.0.0-beta)


## Intoduction
This pipeline is designed for automated end-to-end quality control (QC) and analysis of scMultiome data. The pipeline was developed by [Yong Zeng](mailto:yzeng@uhnresearch.ca) based on some prior work of Mathieu Lupien Lab.


### Features
- **Portability**: The pipeline was developed with [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html), which will automatically deploy the execution environments. Users can also specify the version of softwares to be install in YAML files. It can also be performed across different cluster engines (e.g. SLURM) or stand-alone machines.
- **Flexibility**: The pipeline can be applied to individual samples, as well as to aggregate multiple samples from large-scale profiling.

### Citation
Work-in-progress

### How it works
This schematic diagram shows you how pipeline will be working:
<img src="figures/scMultiome_workflow.png" alt="Schematic_diagram" style="width:100.0%" />


## Installation
1) Make sure that you have a Conda-based Python3 distribution(e.g.,the [Miniconda](https://docs.conda.io/en/latest/miniconda.html)). The Miniconda3-py38_23.3.1-0-Linux-x86_64.sh for Linux is prefered to avoid potential conflicts. The installation of [Mamba](https://github.com/mamba-org/mamba) is also recommended:

	```bash
	$ conda install -n base -c conda-forge mamba
	```

2) Git clone this pipeline.
	```bash
	$ cd
	$ git clone https://github.com/yzeng-lol/iSHARC
	```

3) Install pipeline\'s core enviroment
	```bash
	$ cd iSHARC
	$ conda activate base
	$ mamba env create --file conda_env.yaml
	```

4) Test run
	> **IMPORTANT**: EXTRA ENVIRONMENTS WILL BE INSTALLED, MAKE SURE YOU STILL HAVE INTERNET ACCESS.


5) Run on HPCs
  (Under construction ...)
	You can also submit this pipeline to clusters with the template ./workflow/sbatch_Snakefile_template.sh. This template is for SLURM, however, it could be modified to different resource management systems. More details about cluster configuration can be found at [here](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

	```bash
	## Test run by SLURM submission
	$ sed -i 's,/path/to,'"$PWD"',g' ./workflow/sbatch_Snakefile_template.sh    # replace PATHs for testing
	$ sbatch ./workflow/sbatch_Snakefile_template.sh
	```

## Assets and Troubleshooting
There are several scripts
