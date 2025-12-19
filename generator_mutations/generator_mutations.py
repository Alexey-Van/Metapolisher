from Bio import SeqIO
import random
import sys
import argparse
from Bio.Seq import Seq 

random.seed(42)

# ------------------------
# Command line arguments
# ------------------------
parser = argparse.ArgumentParser(description="Внесение мутаций в FASTA")
parser.add_argument("--fasta_in", default="data/T2T/chm13_chr1.fa", help="Входной FASTA")
parser.add_argument("--fasta_out", default="chm13_chr1_mutated.fa", help="Выходной FASTA")
parser.add_argument("--vcf_out", default="chm13_chr1_mutated.vcf", help="Выходной VCF")
parser.add_argument("--tsv_out", default="chm13_chr1_mutated.tsv", help="Выходной TSV")

parser.add_argument("--num_snvs", type=int, default=200, help="Количество SNV")
parser.add_argument("--num_indels", type=int, default=100, help="Количество маленьких INDEL")
parser.add_argument("--num_big_del", type=int, default=50, help="Количество больших делеций")

args = parser.parse_args()

# ------------------------
# Чтение FASTA
# ------------------------
record = SeqIO.read(args.fasta_in, "fasta")
seq = list(str(record.seq))
chrom = record.id

mutations = []


def write_mut(chrom, pos, ref, alt, mut_type):
    """Сохраняем информацию о мутации."""
    mutations.append([chrom, pos, ref, alt, mut_type])


# ------------------------
# 1) SNV
# ------------------------
for _ in range(args.num_snvs):
    pos = random.randint(0, len(seq) - 1)
    ref = seq[pos]
    alt = random.choice([b for b in "ACGT" if b != ref])
    seq[pos] = alt
    write_mut(chrom, pos + 1, ref, alt, "SNV")


# ------------------------
# 2) small INDEL
# ------------------------
for _ in range(args.num_indels):
    pos = random.randint(0, len(seq) - 5)
    
    if random.random() < 0.5:
        # deletion of 1–3 bp
        length = random.randint(1, 3)
        ref = "".join(seq[pos:pos+length])
        alt = seq[pos]
        del seq[pos:pos+length]
        write_mut(chrom, pos + 1, ref, alt, f"DEL_{length}")
    else:
        # insertion of 1–3 bp
        ins_len = random.randint(1, 3)
        ins = "".join(random.choice("ACGT") for _ in range(ins_len))
        ref = seq[pos]
        alt = ref + ins
        seq[pos] = alt
        write_mut(chrom, pos + 1, ref, alt, f"INS_{ins_len}")


# ------------------------
# 3) large deletions 100–5000 bp
# ------------------------
for _ in range(args.num_big_del):
    start = random.randint(0, len(seq) - 6000)
    length = random.randint(100, 5000)
    end = start + length
    ref = "".join(seq[start:end])
    alt = seq[start]
    del seq[start:end]
    write_mut(chrom, start + 1, ref[:10] + "...", alt, f"LARGE_DEL_{length}")


# ------------------------
# save of mutated FASTA
# ------------------------
from Bio.Seq import Seq
record.seq = Seq("".join(seq))
SeqIO.write(record, args.fasta_out, "fasta")
