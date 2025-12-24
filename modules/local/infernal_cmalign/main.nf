process INFERNAL_CMALIGN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::infernal=1.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/infernal:1.1.5--pl5321h031d066_3' :
        'quay.io/biocontainers/infernal:1.1.5--pl5321h031d066_3' }"

    input:
    tuple val(meta), path(fasta)
    path cm_index

    output:
    tuple val(meta), path("*.sto"), emit: stockholm
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rfam_id = meta.rfam_id
    """
    # Extract the specific CM for this Rfam family
    cmfetch Rfam.cm ${rfam_id} > ${rfam_id}.cm

    # Align the sequence to the CM
    cmalign \\
        --outformat Stockholm \\
        -o ${prefix}.sto \\
        --cpu ${task.cpus} \\
        ${args} \\
        ${rfam_id}.cm \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infernal: \$(cmalign -h | grep -E '^# INFERNAL' | sed 's/# INFERNAL //' | sed 's/ .*//')
    END_VERSIONS
    """
}
