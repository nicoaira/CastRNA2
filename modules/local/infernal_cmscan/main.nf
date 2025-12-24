process INFERNAL_CMSCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::infernal=1.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/infernal:1.1.5--pl5321h031d066_3' :
        'quay.io/biocontainers/infernal:1.1.5--pl5321h031d066_3' }"

    input:
    tuple val(meta), path(fasta)
    path cm_index

    output:
    tuple val(meta), path("*.tblout"), emit: tblout
    tuple val(meta), path("*.out")   , emit: output
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cmscan \\
        --tblout ${prefix}.tblout \\
        --cpu ${task.cpus} \\
        ${args} \\
        Rfam.cm \\
        ${fasta} \\
        > ${prefix}.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infernal: \$(cmscan -h | grep -E '^# INFERNAL' | sed 's/# INFERNAL //' | sed 's/ .*//')
    END_VERSIONS
    """
}
