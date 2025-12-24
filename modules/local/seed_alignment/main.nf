process SEED_ALIGNMENT {
    tag "$meta.id"
    label 'process_single'

    // No container needed - uses system Python (only stdlib required)

    input:
    tuple val(meta), path(fasta)
    path seed_file

    output:
    tuple val(meta), path("*.sto"), emit: stockholm
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rfam_id = meta.rfam_id
    """
    extract_seed_alignment.py \\
        --seed ${seed_file} \\
        --fasta ${fasta} \\
        --rfam-id ${rfam_id} \\
        --seq-id "${meta.id}" \\
        --output ${prefix}.sto

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
