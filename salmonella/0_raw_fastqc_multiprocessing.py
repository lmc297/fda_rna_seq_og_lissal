#!/usr/bin/env python

import sys, os, glob
import multiprocessing

from multiprocessing.pool import ThreadPool

# run FastQC v0.12.1 on each raw read file (*.fastq.gz)
def run_fastqc(f):
  cmd = "fastqc -t 1 {0}".format(f)
  os.system(cmd)

if __name__ == "__main__":
  pool = multiprocessing.Pool(int(sys.argv[1]))
  tasks = glob.glob("raw_reads/*.fastq.gz")
  pool.map(run_fastqc, tasks)
