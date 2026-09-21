from Bio import SeqIO

min_len = 6
with open("homopolymers.bed", "w") as out:
    for record in SeqIO.parse("../data/T2T/modified_fasta/chm13_chr1_mutated.fa", "fasta"):
        seq = str(record.seq)
        i = 0
        while i < len(seq):
            j = i
            while j < len(seq) and seq[j] == seq[i]:
                j += 1
            if j - i >= min_len:
                out.write(f"{record.id}\t{i}\t{j}\thomopolymer\n")
            i = j
