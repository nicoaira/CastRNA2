process PARSE_CMSCAN {
    tag "$meta.id"
    label 'process_single'

    // No container needed - uses system Python (only stdlib required)

    input:
    tuple val(meta), path(tblout)

    output:
    tuple val(meta_updated), path("*.json"), emit: result, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_score = task.ext.min_score ?: params.cmscan_score
    def min_evalue = task.ext.min_evalue ?: params.cmscan_evalue
    // meta_updated will be set by the script output
    meta_updated = meta
    """
    parse_cmscan.py \\
        --tblout ${tblout} \\
        --seq-id "${meta.id}" \\
        --min-score ${min_score} \\
        --min-evalue ${min_evalue} \\
        --output ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    meta_updated = meta
    """
    echo '{"seq_id": "${meta.id}", "rfam_id": null, "status": "NO_HIT"}' > ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
