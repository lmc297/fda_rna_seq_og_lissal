#!/bin/bash

# loop through trimmed forward reads
for f in trimmed_reads_stranded/*_R1_val_1.fq.gz
do
# map trimmed paired-end reads to reference genome with BWA MEM
# sort with samtools sort and convert to BAM format
bwa mem -t 28 reference_genome/listeria_GCA_016306305.1.fna $f ${f%_R1_val_1.fq.gz}_R2_val_2.fq.gz | samtools sort -@ 2 -O bam - > ${f%_R1_val_1.fq.gz}.bam
# index sorted BAM with samtools
samtools index -@ 2 ${f%_R1_val_1.fq.gz}.bam
# ignore unmapped reads
# count number of mapped reads
samtools view -F 0x4 ${f%_R1_val_1.fq.gz}.bam | cut -f 1 | sort | uniq | wc -l > ${f%_R1_val_1.fq.gz}.counts
done
