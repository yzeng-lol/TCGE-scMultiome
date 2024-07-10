################################################################################
## Perform QC and exploration analyses by integrating multiple samples per study
################################################################################

## merge and integrate snRNA-seq data
rule integrate_rna:
    input:
        samples_integr
    output:
        "integration/rna/{sample}_integrated_by_harmony.RDS",
        "integration/rna/{sample}_integrated_by_anchors.RDS"
    #resources:
    #    mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"]
    log:
        "logs/{sample}_integrated_by_harmony_and_anchors.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        "(Rscript --vanilla {params.pipe_dir}/workflow/scripts/integrate_rna.R {input}) 2> {log}"

## merge and integrate scATAC-seq data
rule integrate_atac:
    input:
        samples_integr
    output:
        "integration/atac/{sample}_integrated_by_harmony.RDS",
        "integration/atac/{sample}_integrated_by_anchors.RDS"
    #resources:
    #    mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"]
    log:
        "logs/{sample}_integrated_by_harmony_and_anchors.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        "(Rscript --vanilla {params.pipe_dir}/workflow/scripts/integrate_atac.R {input} "
        "{params.pipe_dir}) 2> {log}"

## Integrate snRNA-seq and scATAC-seq for multiple samples (Harmonized) with WNN
rule integrate_rna_atac_harmony:
    input:
        integrated_rna = "integration/rna/RNA_integrated_by_harmony.RDS",
        integrated_atac = "integration/atac/ATAC_integrated_by_harmony.RDS",
    output:
        "integration/wnn/harmony/{sample}_integrated_by_WNN.RDS"
    #resources:
    #    mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"]
    log:
        "logs/{sample}_harmonized_integrated_by_WNN.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        "(Rscript --vanilla {params.pipe_dir}/workflow/scripts/integrate_rna_atac.R harmony "
        "{input.integrated_rna}  {input.integrated_atac}) 2> {log}"


## Integrate snRNA-seq and scATAC-seq for multiple samples (SEURAT anchors) with WNN
rule integrate_rna_atac_anchor:
    input:
        integrated_rna = "integration/rna/RNA_integrated_by_anchors.RDS",
        integrated_atac = "integration/atac/ATAC_integrated_by_anchors.RDS",
    output:
        "integration/wnn/anchor/{sample}_integrated_by_WNN.RDS"
    #resources:
    #    mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"]
    log:
        "logs/{sample}_anchored_integrated_by_WNN.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        "(Rscript --vanilla {params.pipe_dir}/workflow/scripts/integrate_rna_atac.R anchor "
        "{input.integrated_rna}  {input.integrated_atac}) 2> {log}"


############################################
## generate integrated QC report per project
############################################
rule integration_report:
    input:
        "integration/wnn/harmony/RNA_ATAC_integrated_by_WNN.RDS",
        "integration/wnn/anchor/RNA_ATAC_integrated_by_WNN.RDS"
    output:
        "integration/{sample}_scMultiome_QC_and_Primary_Results_Report.html"      ## project_id
    #resources:
    #    mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"],
        work_dir = config["work_dir"]
    log:
        "logs/{sample}_report.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        ## generating qc report named by sample id
        "(cp {params.pipe_dir}/workflow/scripts/integration_report.Rmd "
        "    {params.work_dir}/integration/{wildcards.sample}_scMultiome_QC_and_Primary_Results_Report.Rmd && "
        "Rscript --vanilla {params.pipe_dir}/workflow/scripts/integration_report.R "
        "  {params.work_dir}/integration/{wildcards.sample}_scMultiome_QC_and_Primary_Results_Report.Rmd "
        "  {params.pipe_dir} {wildcards.sample} {params.work_dir}/integration) 2> {log}"
