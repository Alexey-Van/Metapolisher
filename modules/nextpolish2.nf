process NP2 {

    tag "nextpolish2"

    container params.polishing_container

    publishDir "${params.outdir}/nextpolish2", mode: 'copy'

    input:
        val ready
        path asm
        path hifi_bam
        path r1
        path r2

    output:
        path "np2.vcf", emit: vcf
        path "np2.polished.fa", emit: fa

    script:
    """
    set -euo pipefail

    PREFIX=np2

    yak count \
        -b37 \
        -t ${task.cpus} \
        -k21 \
        -o short.21.yak \
        ${r1} \
        ${r2}

    yak count \
        -b37 \
        -t ${task.cpus} \
        -k31 \
        -o short.31.yak \
        ${r1} \
        ${r2}

    nextPolish2 \
        -t ${task.cpus} \
        -o np2.polished.fa \
        ${hifi_bam} \
        ${asm} \
        short.21.yak \
        short.31.yak
    

    minimap2 \
        --cs \
        -cx asm5 \
        ${asm} \
        np2.polished.fa \
        > asm.paf

    sort -k6,6 -k8,8n asm.paf > asm.sorted.paf

    paftools.js call \
        -f ${asm} \
        asm.sorted.paf \
        > np2.vcf
    """
}