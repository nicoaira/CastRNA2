process GENERATE_SUMMARY {
    label 'process_single'

    // No container needed - uses system Python (only stdlib required)

    input:
    path reports
    path excluded_reports, stageAs: 'excluded/*'

    output:
    path "summary.tsv"   , emit: summary
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate_summary.py \\
        --reports ${reports} \\
        --excluded excluded/ \\
        --output summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
