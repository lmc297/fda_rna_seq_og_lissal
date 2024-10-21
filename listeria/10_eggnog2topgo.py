#!/usr/bin/python

import sys, os, glob
from collections import Counter
from Bio import SeqIO

# create ref2go, which maps eggnog genes to GO terms
ref2go = {}
with open("out.emapper.annotations", "r") as infile:
	for line in infile:
		if not line.startswith("#"):
			gene = line.split("\t")[0].strip()
			allgo = line.split("\t")[9].strip()
			allgo = allgo.split(",")
			allgo = ", ".join(allgo).strip()
			if "GO" in allgo:
				ref2go[gene] = allgo

# create eggnog.tsv and print header line
with open("eggnog.tsv", "a") as outfile:
	print("locus" + "\t" + "terms", file = outfile)

# parse eggnog queries fasta
# for each gene, get its eggnog GO terms
# if available, print gene-to-GO mappings to eggnog.tsv
with open("queries.raw", "r") as infile:
	for record in SeqIO.parse(infile, "fasta"):
		seqdes = str(record.description).strip()
		locus = seqdes.split("locus_tag=")[1].strip()
		locus = locus.split("]")[0].strip()
		seqid = str(record.id).strip()
		if seqid in ref2go.keys():
			final = locus + "\t" + ref2go[seqid]
			with open("eggnog.tsv", "a") as outfile:
				print(final, file = outfile)
