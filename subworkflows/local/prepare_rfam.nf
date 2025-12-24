/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE_RFAM Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Downloads Rfam database files if not provided, then indexes with cmpress
----------------------------------------------------------------------------------------
*/

include { RFAM_DOWNLOAD   } from '../../modules/local/rfam_download/main'
include { INFERNAL_CMPRESS } from '../../modules/local/infernal_cmpress/main'

workflow PREPARE_RFAM {

    main:

    ch_versions = Channel.empty()

    // Determine alignment method - seed alignment requires the seed file
    def use_seed_alignment = params.alignment_method == 'seed'

    // Determine if we need to download or use provided files
    if (params.rfam_cm) {
        // User provided CM file - use directly
        ch_rfam_cm = Channel.fromPath(params.rfam_cm, checkIfExists: true)
        // Clanin is optional - use if provided
        if (params.rfam_clanin) {
            ch_rfam_clanin = Channel.fromPath(params.rfam_clanin, checkIfExists: true)
        } else {
            ch_rfam_clanin = Channel.empty()
        }
        // Seed file - use if provided, otherwise empty (will be downloaded if needed and not provided)
        if (params.rfam_seed) {
            ch_rfam_seed = Channel.fromPath(params.rfam_seed, checkIfExists: true)
        } else if (use_seed_alignment) {
            error """
            ============================================================
            ERROR: Seed alignment method requires Rfam seed file.
            
            You specified --alignment_method seed but did not provide --rfam_seed.
            When using user-provided --rfam_cm, you must also provide --rfam_seed
            for seed-based alignment.
            
            Either:
              1. Provide --rfam_seed /path/to/Rfam.seed.gz
              2. Remove --rfam_cm to download both files automatically
              3. Use --alignment_method cmalign instead
            ============================================================
            """
        } else {
            ch_rfam_seed = Channel.empty()
        }
        ch_rfam_version = Channel.of("user-provided")

    } else {
        // Download from Rfam FTP
        def release = params.rfam_release ?: "CURRENT"

        RFAM_DOWNLOAD(release, use_seed_alignment)

        ch_rfam_cm = RFAM_DOWNLOAD.out.cm
        // clanin is optional - older Rfam versions don't have it
        ch_rfam_clanin = RFAM_DOWNLOAD.out.clanin.ifEmpty(Channel.empty())
        // seed is optional - only downloaded if using seed alignment method
        ch_rfam_seed = RFAM_DOWNLOAD.out.seed.ifEmpty(Channel.empty())
        ch_rfam_version = RFAM_DOWNLOAD.out.version
        ch_versions = ch_versions.mix(RFAM_DOWNLOAD.out.versions)
    }

    // Index the CM database
    INFERNAL_CMPRESS(ch_rfam_cm)
    ch_versions = ch_versions.mix(INFERNAL_CMPRESS.out.versions)

    emit:
    cm_index = INFERNAL_CMPRESS.out.cm_index  // Indexed CM database files
    seed     = ch_rfam_seed                    // Seed alignments (optional, for seed method)
    clanin   = ch_rfam_clanin                  // Clan definitions (optional)
    version  = ch_rfam_version                 // Version info
    versions = ch_versions                     // Tool versions
}
