process MERQURY {

    tag "merqury"

    container params.polishing_container

    publishDir "${params.outdir}/merqury", mode: 'copy'

    input:
        val ready
        path asm
        path r1
        path r2

    output:
        path "merqury_regions.bed", emit: bed
        path "**"

    script:
    """
    set -euo pipefail

    meryl count \
        k=21 \
        output reads.meryl \
        ${r1} \
        ${r2}
    
    export MERQURY="/opt/conda/envs/polish/share/merqury"

    merqury.sh \
        reads.meryl \
        "${asm}" \
        merqury_out

    candidate=\$(find . -maxdepth 2 -type f -name '*.bed' \
        -printf '%s\t%p\n' | sort -nr | head -n 1 | cut -f2-)
    if [ -z "\${candidate}" ]; then
        echo "ERROR: Merqury did not produce a BED track" >&2
        exit 1
    fi
    cp "\${candidate}" merqury_regions.bed
    """
}
