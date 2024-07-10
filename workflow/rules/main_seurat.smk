#################################################
## Perform QC and exploration analysis per sample
#################################################
rule main_seurat:
    input:
        gex = "arc_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        atac = "arc_count/{sample}/outs/atac_fragments.tsv.gz"
    output:
        "main_seurat/{sample}.RDS"
    resources:
        mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"],
        macs2_dir = env_dir
    log:
        "logs/{sample}_main_seurat.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        "(Rscript --vanilla {params.pipe_dir}/workflow/scripts/main_seurat.R "
        "{wildcards.sample} {input.gex} {input.atac} {params.macs2_dir}/bin/macs2 "
        "{params.pipe_dir}) 2> {log}"

################################
## generate QC report per sample
################################
rule qc_and_results_report:
    input:
        "main_seurat/{sample}.RDS"
    output:
        "main_seurat/{sample}_scMultiome_QC_and_Primary_Results_Report.html"
    resources:
        mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"],
        work_dir = config["work_dir"]
    log:
        "logs/{sample}_qc_report.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        ## generating qc report named by sample id
        "(cp {params.pipe_dir}/workflow/scripts/qc_and_results_report.Rmd "
        "    {params.work_dir}/main_seurat/{wildcards.sample}_scMultiome_QC_and_Primary_Results_Report.Rmd && "
        "Rscript --vanilla {params.pipe_dir}/workflow/scripts/qc_and_results_report.R "
        "  {wildcards.sample} {params.work_dir}/{input} "
        "  {params.work_dir}/main_seurat/{wildcards.sample}_scMultiome_QC_and_Primary_Results_Report.Rmd "
        "  {params.pipe_dir}) 2> {log}"
