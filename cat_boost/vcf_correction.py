#!/usr/bin/env python3
"""
Convert custom R/V format to standard VCFv4.2
"""

from pathlib import Path

input_file = 'data/T2T/original_fasta/true_raw.vcf'
output_vcf = "data/T2T/original_fasta/truth.vcf"

chrom_length = {"chr1": 248255384}  # при необходимости расширить на все хромосомы

with open(input_file) as fin, open(output_vcf, "w") as fout:
    # Header
    fout.write("##fileformat=VCFv4.2\n")
    for chrom, length in chrom_length.items():
        fout.write(f"##contig=<ID={chrom},length={length}>\n")
    fout.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    fout.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of variant\">\n")
    fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for line in fin:
        if not line.strip():
            continue
        parts = line.strip().split("\t")
        if parts[0] == "V":
            chrom = parts[1]
            pos = int(parts[2])
            end = int(parts[3])
            qual = parts[5]
            ref = parts[6]
            alt = parts[7]

            info_fields = []
            if alt.startswith("<") or alt == "-":
                # Symbolic SV
                alt_val = "<DEL>" if alt in ["-", "DEL"] else alt
                info_fields.append(f"SVTYPE=DEL")
                info_fields.append(f"END={end}")
                ref_val = "N"  # референс любой, т.к. SV
            else:
                # SNV или маленький индел
                ref_val = ref
                alt_val = alt

            info = ";".join(info_fields) if info_fields else "."
            fout.write(f"{chrom}\t{pos}\t.\t{ref_val}\t{alt_val}\t{qual}\tPASS\t{info}\n")
