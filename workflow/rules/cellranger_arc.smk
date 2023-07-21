#######################
## cellranger_arc_count
#######################
rule cellranger_arc_count:
    input:
        #config["bwa_index"],
        get_cellranger_arc_ref(),
        get_library
    output:
        "cellranger_arc_count/{sample}/outs/filtered_feature_bc_matrix.h5"
    threads: 12
    log:
        "logs/{sample}_cellranger_arc_count.log"
    shell:
        "(cellranger-arc count --id {} --reference= {}| "
        "--library={}) 2> {log}"
