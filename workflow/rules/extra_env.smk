################################################################################
## Check and install extra conda envs automatically,
## Need internet during the inital installation and updating
## wildcard {sample} : installed

### all extra env
rule install_all_extra_env:
    input:
        'extra_env/cellranger_arc_{sample}',
        'extra_env/R_pkgs_{sample}',
    output:
        'extra_env/all_extra_env_{sample}'    ## using wildcard works for --cluster as well
    shell:
        'touch {output}'

## individual extra env
rule install_extra_env_4_cellranger_arc:
    output:
        'extra_env/cellranger_arc_{sample}'
    shell:
        "touch {output}  && "
        "export PATH=/cluster/tools/software/centos7/cellranger-arc/2.0.2:$PATH && "
        "echo $PATH >>  {output}"

rule install_extra_env_4_seurat:
    output:
        'extra_env/R_pkgs_{sample}'
    conda:
        'extra_env/R_pkgs.yaml'
    shell:
        'touch {output}'
