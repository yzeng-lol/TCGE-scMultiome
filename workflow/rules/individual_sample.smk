#########################################################################
## vertical integration of matched RNA and ATAC for each individual sample
## The analyses includes:
##      * re-peak calling using MACS for ATAC
##      * create seurat object with RNA, ATAC_ARC and ATAC assays
#########################################################################
rule initializaion_of_individual_sample:
    input:
        fbm = "arc_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        atac = "arc_count/{sample}/outs/atac_fragments.tsv.gz",
        pbm = "arc_count/{sample}/outs/per_barcode_metrics.csv"
    output:
        "individual_samples/{sample}/{sample}_initial_seurat_object.RDS"
    resources:
        mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"],
        macs_dir = env_dir
    log:
        "logs/{sample}_initialization.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        "(mkdir -p individual_samples/{wildcards.sample} && "
        "Rscript --vanilla {params.pipe_dir}/workflow/scripts/initialization_of_individual_sample.R "
        "  --sample_id {wildcards.sample} --feature_barcode_matrix {input.fbm} "
        "  --per_barcode_metrics {input.pbm}  --atac_file {input.atac} "
        "  --macs_dir {params.macs_dir}/bin/macs2 "
        "  --pipe_dir {params.pipe_dir}) 2> {log}"


#########################################################################
## vertical integration of matched RNA and ATAC for each individual sample
## The analyses includes:
##      * QC assessment and second round cell filltering as coustomized
##      * cell cycle correction for RNA
##      * normalizedion, dimentionality reduction
##      * integration of RNA and ATAC using WNN
##      * clustering using the integrated or sparated modalites
#########################################################################
rule vertical_integration_of_individual_sample:
    input:
        "individual_samples/{sample}/{sample}_initial_seurat_object.RDS"
    output:
        "individual_samples/{sample}/{sample}_vertically_integrated_seurat_object.RDS"
    resources:
        mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"],
        srf = config["second_round_filter"],
        rcc = config["regress_cell_cycle"],
        fgm = config["future_globals_maxSize"],
        min_RNA = config["second_round_cutoffs"]["nCount_RNA_min"],
        max_RNA = config["second_round_cutoffs"]["nCount_RNA_max"],
        min_ATAC = config["second_round_cutoffs"]["nCount_ATAC_min"],
        max_ATAC = config["second_round_cutoffs"]["nCount_ATAC_max"],
        max_MT = config["second_round_cutoffs"]["pct_MT_max"],
        min_TSS = config["second_round_cutoffs"]["TSS_Enrichment_min"],
        max_NS = config["second_round_cutoffs"]["Nucleosome_Signal_max"],
        knn_k = config["clustering_params"]["knn_k"],
        dims_n = config["clustering_params"]["dims_n"],
        comm_res = config["clustering_params"]["comm_res"]    
    threads:
        config["threads"]
    log:
        "logs/{sample}_vertical_integration.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        "(mkdir -p individual_samples/{wildcards.sample} && "
        "Rscript --vanilla {params.pipe_dir}/workflow/scripts/vertical_integration_of_individual_sample.R "
        "  --threads {threads} "
        "  --future_globals_maxSize {params.fgm} "
        "  --sample_id {wildcards.sample}  "
        "  --initial_seurat_object  {input} "
        "  --second_round_filter {params.srf} "
        "  --min_nCount_RNA {params.min_RNA} "
        "  --max_nCount_RNA {params.max_RNA} "
        "  --min_nCount_ATAC {params.min_ATAC} "
        "  --max_nCount_ATAC {params.max_ATAC} "
        "  --max_pct_MT {params.max_MT} "
        "  --min_TSS_Enrichment {params.min_TSS} "
        "  --max_Nucleosome_Signal {params.max_NS} "
        "  --knn_k_param {params.knn_k} "
        "  --dimentions_n {params.dims_n} "
        "  --community_resolution {params.comm_res} "
        "  --regress_cell_cycle {params.rcc}) 2> {log}"

################################################################################
## extended analyses based on integrated RNA and ATAC for each individual sample
## The analyses includes:
##      * cell type annotation
##      * Identify cluster specific DEGs and enriched functions
##      * Identify cluster specific DARs and enriched motif/TFs
##      * Gene regulatory network analysis
################################################################################
rule extended_analyses_of_individual_sample:
    input:
        "individual_samples/{sample}/{sample}_vertically_integrated_seurat_object.RDS"
    output:
        "individual_samples/{sample}/{sample}_extended_seurat_object.RDS"
    resources:
        mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"],
        fgm = config["future_globals_maxSize"]
    threads:
        config["threads"]
    log:
        "logs/{sample}_extended_analyses.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        "(mkdir -p individual_samples/{wildcards.sample} && "
        "Rscript --vanilla {params.pipe_dir}/workflow/scripts/extended_analyses_of_individual_sample.R "
        "  --threads {threads} "
        "  --future_globals_maxSize {params.fgm} "
        "  --sample_id {wildcards.sample}  "
        "  --vertically_integrated_seurat_object  {input} "
        "  --pipe_dir {params.pipe_dir}) 2> {log}"


##################################################################################
## generate the HTML report of QC and primairy results for each indiviudal  sample
##################################################################################
rule html_report_of_individual_sample:
    input:
        "individual_samples/{sample}/{sample}_extended_seurat_object.RDS"
    output:
        "individual_samples/{sample}/{sample}_QC_and_Primary_Results.html"
    resources:
        mem_mb=60000
    params:
        pipe_dir = config["pipe_dir"],
        work_dir = config["work_dir"]
    log:
        "logs/{sample}_html_report.log"
    conda:
        "extra_env/R_pkgs.yaml"
    shell:
        ## generating qc report named by sample id
        "(cp {params.pipe_dir}/workflow/scripts/qc_and_primary_results_of_individual_sample.Rmd "
        "    {params.work_dir}/individual_samples/{wildcards.sample}/{wildcards.sample}_QC_and_Primary_Results.Rmd && "
        "Rscript --vanilla {params.pipe_dir}/workflow/scripts/qc_and_primary_results_of_individual_sample.R "
        "  --sample_id {wildcards.sample} "
        "  --extended_analyses_seurat_object {params.work_dir}/{input} "
        "  --report_rmd_file {params.work_dir}/individual_samples/{wildcards.sample}/{wildcards.sample}_scMultiome_QC_and_Primary_Results_Report.Rmd "
        "  --pipe_dir {params.pipe_dir}) 2> {log}"
