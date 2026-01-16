#!/bin/bash
SGE_Batch -c "bowtie2-build HSVd_genome.fasta HSVd" -r HSVd_bowtie2_build
