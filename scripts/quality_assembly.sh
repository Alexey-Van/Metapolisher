#!/usr/bin/env bash

set -e

############################################
# PARAMETERS
############################################

ASM=$1            # assembly.fasta
READS=$2          # reads.fastq.gz (или список файлов)
BUSCO_DB=$3       # busco lineage (например mammalia_odb10)
THREADS=$4        # число потоков
KMER=$5           # размер k-mer
OUTDIR=$6         # директория результатов
REF=$7            # reference genome (optional, для QUAST)

############################################
# PREPARE
############################################
 
mkdir -p ${OUTDIR}
cd ${OUTDIR}


# ############################################
# # INSTALL MERQURY + MERYL
# ############################################

# echo "Installing Merqury..."

# git clone https://github.com/marbl/merqury.git
# /mnt/data/vanichkinao/Y/micromamba create -n quast_env -c bioconda -c conda-forge meryl quast
# /mnt/data/vanichkinao/Y/micromamba create -y -n busco_env -c bioconda busco

############################################
# RUN MERQURY
############################################

echo "Running Merqury..."

/mnt/data/vanichkinao/Y/micromamba run -n quast_env meryl count \
    k=${KMER} \
    threads=${THREADS} \
    output reads.meryl \
    ${READS}

bash /mnt/data/vanichkinao/MetaPolisher/merqury/merqury.sh \
    reads.meryl \
    ${ASM} \
    merqury_results

############################################
# RUN QUAST
############################################

echo "Running QUAST..."

if [ -z "$REF" ]; then

    /mnt/data/vanichkinao/Y/micromamba run -n quast_env quast \
        ${ASM} \
        -o quast_results \
        -t ${THREADS}

else

    /mnt/data/vanichkinao/Y/micromamba run -n quast_env quast \
        ${ASM} \
        -r ${REF} \
        -o quast_results \
        -t ${THREADS}

fi

############################################
# INSTALL BUSCO
############################################

echo "Installing BUSCO..."



############################################
# RUN BUSCO
############################################

echo "Running BUSCO..."

/mnt/data/vanichkinao/Y/micromamba run -n busco_env busco \
-i ${ASM} \
-l ${BUSCO_DB} \
-m genome \
-c ${THREADS} \
-o busco_results

############################################
# DONE
############################################

echo "All analyses completed"

echo "Results directories:"
echo "Merqury -> ${OUTDIR}/merqury_results"
echo "BUSCO -> ${OUTDIR}/busco_results"
echo "QUAST -> ${OUTDIR}/quast_results"