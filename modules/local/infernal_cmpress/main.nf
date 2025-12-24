process INFERNAL_CMPRESS {
    tag "$cm.baseName"
    label 'process_medium'

    conda "bioconda::infernal=1.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/infernal:1.1.5--pl5321h031d066_3' :
        'quay.io/biocontainers/infernal:1.1.5--pl5321h031d066_3' }"

    input:
    path cm

    output:
    path "Rfam.cm*", emit: cm_index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def is_compressed = cm.name.endsWith('.gz')
    """
    # Decompress if needed
    if [ "${is_compressed}" = "true" ]; then
        gunzip -c ${cm} > Rfam.cm
    else
        cp ${cm} Rfam.cm
    fi

    # Index the CM database
    cmpress Rfam.cm

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infernal: \$(cmscan -h | grep -E '^# INFERNAL' | sed 's/# INFERNAL //' | sed 's/ .*//')
    END_VERSIONS
    """
}
