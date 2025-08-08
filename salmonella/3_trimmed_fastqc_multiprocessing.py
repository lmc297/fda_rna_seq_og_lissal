#!/usr/bin/env python

import sys, os, glob
import multiprocessing

from multiprocessing.pool import ThreadPool

def run_fastqc(f):
  cmd = "fastqc -t 1 {0}".format(f)
  os.system(cmd)

if __name__ == "__main__":
  pool = multiprocessing.Pool(int(sys.argv[1]))
  tasks = glob.glob("trimmed_fastp/*.fastq.gz")
  pool.map(run_fastqc, tasks)
