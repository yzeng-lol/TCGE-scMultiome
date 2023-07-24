#######################
##
#######################
rule main_seurat:
    input:
        gex = "arc_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        atac = "arc_count/{sample}/outs/atac_fragments.tsv.gz"
    output:
        "main_seurat/{sample}.RDS"
    resources:
        mem_mb=60000
    params:
        scr_dir = config["pipe_dir"]
    log:
        "main_seurat/{sample}_main_seurat.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        "(Rscript --vanilla {params.scr_dir}/workflow/scripts/main_seurat.R "
        "{wildcards.sample} {input.gex} {input.atac}) 2> {log}"
