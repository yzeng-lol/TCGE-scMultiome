#################################################
## prepare library files for cellrancer_arc_count
#################################################
rule prep_library_file:
    input:
        samples_list
    output:
        "arc_count/{sample}_library.csv"
    params:
        id = get_sample_seq_id,
        gex_path = get_gex_seq_path,
        atac_path = get_atac_seq_path
    shell:
        ## generate library.csv per sample
        "echo 'fastqs,sample,library_type' > arc_count/{wildcards.sample}_library.csv && "
        "echo '{params.gex_path},{params.id},Gene Expression' >>  arc_count/{wildcards.sample}_library.csv && "
        "echo '{params.atac_path},{params.id},Chromatin Accessibility' >>  arc_count/{wildcards.sample}_library.csv "


#######################
## cellranger_arc_count
#######################
rule cellranger_arc_count:
    input:
        get_cellranger_arc_ref(),
        "arc_count/{sample}_library.csv"
        #get_library
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
        ## output will create the {sample} folder prior the cellranger, leading to error of
        ## {sample} is not a pipestance directory for cellranger
        "(export PATH={params.arc_dir}:$PATH && "
        "cd {params.out_dir}/arc_count && "
        "cellranger-arc count --id {wildcards.sample}_res --reference={input[0]} "
        "--libraries={params.out_dir}/{input[1]}  --localcores={threads} && "
        "rm -r {wildcards.sample} && mv {wildcards.sample}_res {wildcards.sample}) 2> {log}"
