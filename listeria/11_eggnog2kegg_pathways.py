#!/usr/bin/python

import sys, os, glob
from collections import Counter
from Bio import SeqIO

# create ref2go, which maps eggnog genes to kegg pathways
ref2go = {}
with open("out.emapper.annotations", "r") as infile:
	for line in infile:
		if not line.startswith("#"):
			gene = line.split("\t")[0].strip()
			allgo = line.split("\t")[12].strip()
			allgo = allgo.split(",")
			if gene not in ref2go.keys():
				ref2go[gene] = []
			for ag in allgo:
				if ag.strip() != "-":
					ref2go[gene].append(ag.strip())

# parse queries fasta supplied to eggnog
# if a gene has a kegg pathway, print to kegg_pathways.tsv
with open("queries.raw", "r") as infile:
	for record in SeqIO.parse(infile, "fasta"):
		seqdes = str(record.description).strip()
		locus = seqdes.split("locus_tag=")[1].strip()
		locus = locus.split("]")[0].strip()
		seqid = str(record.id).strip()
		if seqid in ref2go.keys():
			vals = ref2go[seqid]
			for v in vals:
				final = v.strip() + "\t" + locus.strip()
				with open("kegg_pathways.tsv", "a") as outfile:
					print(final, file = outfile)
