#!/usr/bin/python

import sys, os, glob
from Bio import SeqIO

# add annotation information to DESeq results

# create dictionary linking genes to eggnog annotations
eggnog = {}
with open("eggnog/out.emapper.annotations", "r") as infile:
	for line in infile:
		if not line.startswith("##"):
			if "#query" in line:
				eheader = line.strip()
			else:
				seqid = line.split("\t")[0].strip()
				eggnog[seqid] = line.strip() 

# create dictionary linking genes to NCBI annotations
d = {}
dlocus = {}
with open("eggnog/listeria_GCA_016306305.1.ffn", "r") as infile:
	for record in SeqIO.parse(infile, "fasta"):	
		seqdes = str(record.description).strip()
		seqid = str(record.id).strip()
		d[seqid] = seqdes
		locus = seqdes.split("locus_tag=")[1].strip()
		locus = locus.split("]")[0].strip()	
		dlocus[locus] = seqid

# loop through DESeq results TSV files
# add annotaation information and print TSV files
for f in glob.glob("./deseq_**.tsv"):
	fname = f.split("/")[-1].strip()
	with open(f, "r") as infile:
		for line in infile:
			if "baseMean" not in line:
				locus = line.split("\t")[0].strip()
				if locus in dlocus.keys():
					seqid = dlocus[locus]
				else:
					seqid = "no_ffn"
				if seqid in d.keys():
					annot = d[seqid]
				else:
					annot = "no_ffn"
				if seqid in eggnog.keys():
					egg = eggnog[seqid]
				else:
					egg = "\t".join(["no_eggnog"] * len(eheader.split("\t"))).strip()
				final = line.strip() + "\t" + annot + "\t" + egg
				with open("annot_" + fname, "a") as outfile:
					print(final, file = outfile)
			else:
				final = line.strip() + "\t" + "ncbi_annotation" + "\t" + eheader
				with open("annot_" + fname, "a") as outfile:
					print(final, file = outfile)
