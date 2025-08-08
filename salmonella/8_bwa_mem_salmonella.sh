#!/bin/bash

for f in trimmed_reads_stranded/*_R1_val_1.fq.gz
do
bwa mem -t 28 reference_genome/salmonella_GCF_000006945.2.fna $f ${f%_R1_val_1.fq.gz}_R2_val_2.fq.gz | samtools sort -@ 2 -O bam - > ${f%_R1_val_1.fq.gz}.bam
samtools index -@ 2 ${f%_R1_val_1.fq.gz}.bam
samtools view -F 0x4 ${f%_R1_val_1.fq.gz}.bam | cut -f 1 | sort | uniq | wc -l > ${f%_R1_val_1.fq.gz}.counts
done
