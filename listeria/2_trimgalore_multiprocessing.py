#!/usr/bin/env python

import sys, os, glob
import multiprocessing

from multiprocessing.pool import ThreadPool

def run_trimgalore(f1):
  f2 = f1.replace("_R1.fastq.gz", "_R2.fastq.gz")
  cmd = "trim_galore --cores 1 --fastqc --paired --gzip -o trimmed_reads {0} {1}".format(f1, f2)
  os.system(cmd)

if __name__ == "__main__":
  pool = multiprocessing.Pool(int(sys.argv[1]))
  tasks = glob.glob("raw_reads/*_R1.fastq.gz")
  pool.map(run_trimgalore, tasks)
