process RFAM_DOWNLOAD {
    tag "$release"
    label 'process_single'

    // No conda/container required - uses system curl/wget
    // Most systems have curl or wget pre-installed

    input:
    val release
    val download_seed

    output:
    path "Rfam.cm.gz"   , emit: cm
    path "Rfam.seed.gz" , emit: seed, optional: true
    path "Rfam.clanin"  , emit: clanin, optional: true
    path "version.txt"  , emit: version
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def base_url = "https://ftp.ebi.ac.uk/pub/databases/Rfam/${release}"
    """
    # Download using curl (more commonly available) with wget fallback
    download_file() {
        local url="\$1"
        local output="\$2"
        if command -v curl &> /dev/null; then
            curl -fsSL "\$url" -o "\$output"
        elif command -v wget &> /dev/null; then
            wget -q "\$url" -O "\$output"
        else
            echo "ERROR: Neither curl nor wget found. Please install one of them."
            exit 1
        fi
    }

    # Try to download a file, return 0 on success, 1 on failure
    try_download_file() {
        local url="\$1"
        local output="\$2"
        if command -v curl &> /dev/null; then
            curl -fsSL "\$url" -o "\$output" 2>/dev/null && return 0 || return 1
        elif command -v wget &> /dev/null; then
            wget -q "\$url" -O "\$output" 2>/dev/null && return 0 || return 1
        else
            return 1
        fi
    }

    # Download CM file (~330 MB)
    echo "Downloading Rfam.cm.gz from ${base_url}..."
    download_file "${base_url}/Rfam.cm.gz" Rfam.cm.gz

    # Download seed alignments if requested (~1.5 GB compressed)
    if [ "${download_seed}" = "true" ]; then
        echo "Downloading Rfam.seed.gz from ${base_url}..."
        download_file "${base_url}/Rfam.seed.gz" Rfam.seed.gz
    fi

    # Try to download clanin file (optional - not available in older Rfam versions)
    echo "Downloading Rfam.clanin..."
    if try_download_file "${base_url}/Rfam.clanin" Rfam.clanin; then
        echo "Successfully downloaded Rfam.clanin"
    else
        echo "WARNING: Rfam.clanin not available for release ${release} (older versions do not include this file)"
        echo "Continuing without clan information..."
    fi

    # Record version info
    echo "release: ${release}" > version.txt
    echo "url: ${base_url}" >> version.txt
    echo "downloaded: \$(date -Iseconds)" >> version.txt
    if [ -f Rfam.clanin ]; then
        echo "clanin: available" >> version.txt
    else
        echo "clanin: not available" >> version.txt
    fi
    if [ -f Rfam.seed.gz ]; then
        echo "seed: available" >> version.txt
    else
        echo "seed: not downloaded" >> version.txt
    fi

    # Get download tool version
    if command -v curl &> /dev/null; then
        DOWNLOAD_TOOL="curl"
        DOWNLOAD_VERSION=\$(curl --version | head -n1 | awk '{print \$2}')
    else
        DOWNLOAD_TOOL="wget"
        DOWNLOAD_VERSION=\$(wget --version | head -n1 | sed 's/GNU Wget //' | sed 's/ .*//')
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \${DOWNLOAD_TOOL}: \${DOWNLOAD_VERSION}
        rfam_release: ${release}
    END_VERSIONS
    """
}
