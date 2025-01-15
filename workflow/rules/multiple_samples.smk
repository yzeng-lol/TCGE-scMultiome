################################################################################
## Perform QC and exploration analyses by integrating multiple samples per study
################################################################################

## horizontally merge and integrate snRNA-seq data across samples
rule horizontal_integration_of_rna_across_multiple_samples:
    input:
        #samples_integr
        expand("individual_samples/{samples}/{samples}_extended_seurat_object.RDS", samples = SAMPLES_INTEGR["sample_id"])
    output:
        "integrated_samples/rna/{sample}_integrated_by_harmony.RDS",
        "integrated_samples/rna/{sample}_integrated_by_anchors.RDS"
    #resources:
    #    mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"],
        integr_list = config["samples_integr"],
        fgm = config["future_globals_maxSize"],
        knn_k = config["clustering_params"]["knn_k"],
        dims_n = config["clustering_params"]["dims_n"],
        comm_res = config["clustering_params"]["comm_res"]
    threads:
        config["threads"]
    log:
        "logs/{sample}_horizontally_integrated_by_harmony_and_anchors.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        "(Rscript --vanilla {params.pipe_dir}/workflow/scripts/horizontal_integration_of_rna_across_multiple_samples.R "
        "   --threads {threads} "
        "   --future_globals_maxSize {params.fgm} "
        "   --knn_k_param {params.knn_k} "
        "   --dimentions_n {params.dims_n} "
        "   --community_resolution {params.comm_res} "
        "   --samples_integration {params.integr_list}) 2> {log}"


## horizontally merge and integrate snATAC-seq data across samples
rule horizontal_integration_of_atac_across_multiple_samples:
    input:
        #samples_integr
        expand("individual_samples/{samples}/{samples}_extended_seurat_object.RDS", samples = SAMPLES_INTEGR["sample_id"])
    output:
        "integrated_samples/atac/{sample}_integrated_by_harmony.RDS",
        "integrated_samples/atac/{sample}_integrated_by_anchors.RDS"
    #resources:
    #    mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"],
        integr_list = config["samples_integr"],
        fgm = config["future_globals_maxSize"],
        knn_k = config["clustering_params"]["knn_k"],
        dims_n = config["clustering_params"]["dims_n"],
        comm_res = config["clustering_params"]["comm_res"]
    threads:
        config["threads"]
    log:
        "logs/{sample}_horizontally_integrated_by_harmony_and_anchors.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        "(Rscript --vanilla {params.pipe_dir}/workflow/scripts/horizontal_integration_of_atac_across_multiple_samples.R "
        "   --threads {threads} "
        "   --future_globals_maxSize {params.fgm} "
        "   --knn_k_param {params.knn_k} "
        "   --dimentions_n {params.dims_n} "
        "   --community_resolution {params.comm_res} "
        "   --samples_integration {params.integr_list} "
        "   --pipe_dir {params.pipe_dir}) 2> {log}"


## Verically integrate horizontally integrated (Harmonized) snRNA-seq and scATAC-seq for multiple samples using WNN
rule vertical_integration_of_multiple_harmonized_samples:
    input:
        integrated_rna = "integrated_samples/rna/RNA_integrated_by_harmony.RDS",
        integrated_atac = "integrated_samples/atac/ATAC_integrated_by_harmony.RDS",
    output:
        "integrated_samples/wnn/harmony/{sample}_integrated_by_WNN.RDS"
    #resources:
    #    mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"],
        fgm = config["future_globals_maxSize"],
        knn_k = config["clustering_params"]["knn_k"],
        dims_n = config["clustering_params"]["dims_n"],
        comm_res = config["clustering_params"]["comm_res"]
    threads:
        config["threads"]
    log:
        "logs/{sample}_harmonized_vertically_integrated_by_WNN.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        "(Rscript --vanilla {params.pipe_dir}/workflow/scripts/vertical_integration_of_multiple_samples.R "
        "   --integration_method harmony "
        "   --threads {threads} "
        "   --future_globals_maxSize {params.fgm} "
        "   --knn_k_param {params.knn_k} "
        "   --dimentions_n {params.dims_n} "
        "   --community_resolution {params.comm_res} "
        "   --integrated_rna  {input.integrated_rna} "
        "   --integrated_atac {input.integrated_atac}) 2> {log}"


## Integrate snRNA-seq and scATAC-seq for multiple samples (SEURAT anchors) with WNN
rule vertical_integration_of_multiple_anchored_samples:
    input:
        integrated_rna = "integrated_samples/rna/RNA_integrated_by_anchors.RDS",
        integrated_atac = "integrated_samples/atac/ATAC_integrated_by_anchors.RDS",
    output:
        "integrated_samples/wnn/anchor/{sample}_integrated_by_WNN.RDS"
    #resources:
    #    mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"],
        fgm = config["future_globals_maxSize"],
        knn_k = config["clustering_params"]["knn_k"],
        dims_n = config["clustering_params"]["dims_n"],
        comm_res = config["clustering_params"]["comm_res"]
    threads:
        config["threads"]
    log:
        "logs/{sample}_anchored_vertically_integrated_by_WNN.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        "(Rscript --vanilla {params.pipe_dir}/workflow/scripts/vertical_integration_of_multiple_samples.R "
        "   --integration_method anchor "
        "   --threads {threads} "
        "   --future_globals_maxSize {params.fgm} "
        "   --knn_k_param {params.knn_k} "
        "   --dimentions_n {params.dims_n} "
        "   --community_resolution {params.comm_res} "
        "   --integrated_rna  {input.integrated_rna} "
        "   --integrated_atac {input.integrated_atac}) 2> {log}"


############################################
## generate integrated QC report per project
############################################
rule html_report_of_multiple_samples:
    input:
        "integrated_samples/wnn/harmony/RNA_ATAC_integrated_by_WNN.RDS",
        "integrated_samples/wnn/anchor/RNA_ATAC_integrated_by_WNN.RDS"
    output:
        "integrated_samples/{sample}_QC_and_Primary_Results.html"      ## project_id
    #resources:
    #    mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"],
        work_dir = config["work_dir"]
    log:
        "logs/{sample}_integrated_samples_report.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        ## generating qc report named by sample id
        "(cp {params.pipe_dir}/workflow/scripts/qc_and_primary_results_report_of_multiple_samples.Rmd "
        "    {params.work_dir}/integrated_samples/{wildcards.sample}_QC_and_Primary_Results.Rmd && "
        "Rscript --vanilla {params.pipe_dir}/workflow/scripts/qc_and_primary_results_report_of_multiple_samples.R "
        "  --report_rmd_file {params.work_dir}/integrated_samples/{wildcards.sample}_QC_and_Primary_Results.Rmd "
        "  --integration_dir {params.work_dir}/integrated_samples "
        "  --pipe_dir {params.pipe_dir}) 2> {log}"
