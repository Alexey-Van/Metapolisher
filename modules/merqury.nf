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

    cp *.bed merqury_regions.bed
    """
}