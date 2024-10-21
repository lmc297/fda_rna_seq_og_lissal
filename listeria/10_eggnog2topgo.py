#!/usr/bin/python

import sys, os, glob
from collections import Counter
from Bio import SeqIO

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

with open("eggnog.tsv", "a") as outfile:
	print("locus" + "\t" + "terms", file = outfile)

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
