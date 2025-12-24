process CONSENSUS_PROJECTION {
    tag "$meta.id"
    label 'process_single'

    // No container needed - uses system Python (only stdlib required)

    input:
    tuple val(meta), path(stockholm)

    output:
    tuple val(meta), path("*.dbn"), emit: dotbracket
    tuple val(meta), path("*.json"), emit: report
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def noncanonical = params.noncanonical_pairs ?: 'keep'
    def keep_brackets = params.keep_all_bracket_types ? '--preserve-pseudoknots' : ''
    """
    consensus_projection.py \\
        --input ${stockholm} \\
        --output ${prefix}.dbn \\
        --report ${prefix}.json \\
        --noncanonical ${noncanonical} \\
        ${keep_brackets} \\
        --seq-id "${meta.id}" \\
        --rfam-id "${meta.rfam_id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
