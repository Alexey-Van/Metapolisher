#!/bin/bash
set -e

ASSEMBLY_INITAL="$1"          
ASSEMBLY_TRUE="$2"  
out_prefix="$3"     
    
minimap2  --cs -cx asm5 ${ASSEMBLY_INITAL} ${ASSEMBLY_TRUE} > ${out_prefix}.minimap2.paf
sort -k6,6 -k8,8n ${out_prefix}.minimap2.paf > ${out_prefix}.minimap2.sorted.paf
paftools.js call -f ${ASSEMBLY_INITAL} ${out_prefix}.minimap2.sorted.paf > ${out_prefix}.paftools.vcf