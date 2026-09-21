process AUTO_POLISH {

    tag "automated_polish"

    container params.polishing_container

    publishDir "${params.outdir}/automated_polishing", mode: 'copy'

    input:
        val ready
        path draft
        path reads
        path readmers_meryl
        val reads_type

    output:
        path "*.consensus.fasta", emit: polished
        path "**"

    script:
    """
    set -euo pipefail

    THREADS=${task.cpus}
    ITERATIONS=1

    cp ${draft} iter_0.consensus.fasta

    CURRENT=iter_0.consensus.fasta

    for i in \$(seq 1 \${ITERATIONS}); do

        echo "[Iteration \$i]"

        if [ "${reads_type}" = "hifi" ]; then

            echo "[1] HiFi polishing"

            K=15

            meryl count \
                k=\${K} \
                output meryl_db \
                \${CURRENT}

            meryl print greater-than 100 meryl_db > repetitive_k15.txt

            winnowmap \
                --MD \
                -W repetitive_k15.txt \
                -ax map-pb \
                -t ${task.cpus} \
                \${CURRENT} \
                ${reads} \
                > aln.sam

        elif [ "${reads_type}" = "ont" ]; then

            echo "[1] ONT polishing"

            K=15

            meryl count \
                k=\${K} \
                output meryl_db \
                \${CURRENT}

            meryl print greater-than 100 meryl_db > repetitive_k15.txt

            winnowmap \
                --MD \
                -W repetitive_k15.txt \
                -ax map-ont \
                -t ${task.cpus} \
                \${CURRENT} \
                ${reads} \
                > aln.sam

        else

            echo "ERROR: unknown reads type: ${reads_type}"
            exit 1

        fi

        samtools view \
            -h \
            -F 0x104 \
            aln.sam \
            > filtered.sam

        racon \
            -t ${task.cpus} \
            ${reads} \
            filtered.sam \
            \${CURRENT} \
            > racon.fasta

        meryl count \
            k=21 \
            output racon.meryl \
            racon.fasta

        minimap2 \
            --cs \
            -cx asm5 \
            \${CURRENT} \
            racon.fasta \
            > asm.paf

        sort -k6,6 -k8,8n asm.paf > asm.sorted.paf

        paftools.js call \
            -f \${CURRENT} \
            asm.sorted.paf \
            > variants.vcf

        bgzip variants.vcf
        tabix -p vcf variants.vcf.gz

        merfin \
            -sequence \${CURRENT} \
            -seqmers racon.meryl \
            -readmers ${readmers_meryl} \
            -peak 106 \
            -output merfin \
            -threads ${task.cpus} \
            -filter \
            -vcf variants.vcf.gz

        awk '
        BEGIN { OFS="\t" }
        /^#/ { print; next }
        {
            key = $1 FS $2
            if (!(key in seen)) {
                seen[key] = 1
                print
            }
        }
        ' merfin.filter.vcf > merfin.fixed.vcf

        bgzip -f merfin.fixed.vcf
        tabix -f -p vcf merfin.fixed.vcf.gz

        # Normalize
        bcftools norm \
            -f \${CURRENT} \
            -m -both \
            merfin.fixed.vcf.gz \
            -Oz \
            -o merfin.norm.vcf.gz

        bcftools index -f merfin.norm.vcf.gz

        # Consensus
        bcftools consensus \
            merfin.norm.vcf.gz \
            -f \${CURRENT} \
            -H 1 \
            > iter_\${i}.consensus.fasta

        CURRENT=iter_\${i}.consensus.fasta

    done
    """
}