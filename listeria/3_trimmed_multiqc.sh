#!/bin/bash

# run MultiQC v1.16 on trimmed fastqc results
multiqc --outdir trimmed_multiqc trimmed_fastqc
