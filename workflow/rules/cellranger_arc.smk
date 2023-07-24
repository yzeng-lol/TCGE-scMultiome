#######################
## cellranger_arc_count
#######################
rule cellranger_arc_count:
    input:
        get_cellranger_arc_ref(),
        get_library
    output:
        "arc_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        "arc_count/{sample}/outs/atac_fragments.tsv.gz",
    params:
        arc_dir = config["arc_dir"],
        out_dir = config["work_dir"]
    threads: 12
    log:
        "logs/{sample}_cellranger_arc_count.log"
    shell:
        "(export PATH={params.arc_dir}:$PATH && "
        "cellranger-arc count --id {wildcards.sample} --reference={input[0]} "
        "--libraries={input[1]}  --localcores={threads} && "
        ## cellranger will created ID folder in work_dir
        "mv {wildcards.sample} {params.out_dir}/arc_count/) 2> {log}"
