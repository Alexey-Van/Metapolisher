/mnt/data/vanichkinao/MetaPolisher/TGS-GapCloser/tgsgapcloserbin/tgsgapcloser  \
        --scaff  /mnt/data/vanichkinao/MetaPolisher/ragtag_without_cor/ragtag.scaffold.fasta \
        --reads  /mnt/data/vanichkinao/MetaPolisher/data/CHM13_T2T_ONT_fastq_guppy_6.3.7_hac_aa.fastq.gz.1 \
        --output test_gap_closer \
        --pilon  /mnt/data/vanichkinao/Y/envs/ragtag/share/pilon-1.24-0/pilon.jar  \
        --ngs    /mnt/data/vanichkinao/MetaPolisher/data/ngs_reads.fastq.gz  \
        --samtools /mnt/data/vanichkinao/Y/envs/ragtag/bin/samtools  \
        --java    /mnt/data/vanichkinao/Y/envs/ragtag/bin/java \
        >pipe.log 2>pipe.err