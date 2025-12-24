/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT_CHECK Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validates input files and creates channels with metadata.
    Supports two input modes:
    1. FASTA file (all sequences go through cmscan inference)
    2. CSV samplesheet with optional Rfam IDs per sequence
----------------------------------------------------------------------------------------
*/

workflow INPUT_CHECK {

    take:
    input_file  // Path to input file (FASTA or CSV)

    main:

    // Detect input type based on file extension
    def input_path = file(input_file)
    def is_csv = input_path.name.endsWith('.csv')

    if (is_csv) {
        // Samplesheet mode: parse CSV and create channels
        ch_input = Channel
            .fromPath(input_file)
            .splitCsv(header: true)
            .map { row ->
                // Validate required fields
                if (!row.seq_id) {
                    error "ERROR: Missing 'seq_id' in samplesheet row: ${row}"
                }

                // Create meta map
                def meta = [
                    id: row.seq_id,
                    rfam_id: row.rfam_id ?: null,  // null if empty
                    single_end: true
                ]

                // Handle fasta path - if not provided, we'll need to extract from multi-fasta
                def fasta = row.fasta ? file(row.fasta, checkIfExists: true) : null

                return [meta, fasta]
            }

        // Separate sequences with and without Rfam IDs
        ch_input
            .branch {
                meta, fasta ->
                    known: meta.rfam_id != null && meta.rfam_id != ''
                    unknown: true
            }
            .set { ch_branched }

        ch_known = ch_branched.known
        ch_unknown = ch_branched.unknown

    } else {
        // FASTA mode: split multi-fasta into individual sequences
        ch_input = Channel
            .fromPath(input_file)
            .splitFasta(record: [id: true, seqString: true, desc: true])
            .map { record ->
                def meta = [
                    id: record.id,
                    rfam_id: null,  // Will be inferred
                    single_end: true
                ]

                // Create a temporary fasta file for this sequence
                // The fasta content will be passed through
                def fasta_content = ">${record.id}"
                if (record.desc) {
                    fasta_content += " ${record.desc}"
                }
                fasta_content += "\n${record.seqString}\n"

                return [meta, fasta_content]
            }

        // All sequences go to unknown (need inference)
        ch_known = Channel.empty()
        ch_unknown = ch_input
    }

    emit:
    known   = ch_known    // tuple(meta, fasta) with rfam_id set
    unknown = ch_unknown  // tuple(meta, fasta_content) without rfam_id
    is_csv  = is_csv      // boolean flag for downstream handling
}


/*
 * Helper workflow to handle multi-fasta input with a mapping CSV
 * This allows using a multi-fasta file with a separate mapping file
 */
workflow INPUT_CHECK_WITH_MAPPING {

    take:
    fasta_file   // Path to multi-fasta file
    mapping_file // Path to mapping CSV (seq_id,rfam_id)

    main:

    // Parse mapping file into a map
    ch_mapping = Channel
        .fromPath(mapping_file)
        .splitCsv(header: true)
        .map { row ->
            [row.seq_id, row.rfam_id ?: null]
        }
        .collectFile(name: 'mapping.txt', newLine: true) { id, rfam_id ->
            "${id}\t${rfam_id ?: ''}"
        }

    // Split fasta and join with mapping
    ch_sequences = Channel
        .fromPath(fasta_file)
        .splitFasta(record: [id: true, seqString: true, desc: true])
        .map { record ->
            [record.id, record]
        }

    // Create the mapping channel
    ch_map = Channel
        .fromPath(mapping_file)
        .splitCsv(header: true)
        .map { row ->
            [row.seq_id, row.rfam_id ?: null]
        }

    // Join sequences with their Rfam IDs
    ch_joined = ch_sequences
        .join(ch_map, by: 0, remainder: true)
        .map { seq_id, record, rfam_id ->
            def meta = [
                id: record.id,
                rfam_id: rfam_id,
                single_end: true
            ]

            def fasta_content = ">${record.id}"
            if (record.desc) {
                fasta_content += " ${record.desc}"
            }
            fasta_content += "\n${record.seqString}\n"

            [meta, fasta_content]
        }

    // Branch based on whether Rfam ID is known
    ch_joined
        .branch {
            meta, fasta_content ->
                known: meta.rfam_id != null && meta.rfam_id != ''
                unknown: true
        }
        .set { ch_branched }

    emit:
    known   = ch_branched.known
    unknown = ch_branched.unknown
}
