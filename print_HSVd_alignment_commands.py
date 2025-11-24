import sys
import os

usage = 'Usage: ' + sys.argv[0] + ' <bowtie 2 index name>' + '<fastq files>'

if len(sys.argv) < 3:
    print(usage)
    exit()

fastq_files = sys.argv[2:]
bowtie2_idx = sys.argv[1]


for fastq_file in fastq_files:
    #print(fastq_file)
    new_sam_suffix = os.path.basename(bowtie2_idx) + ".sam"
    #print(new_sam_suffix)
    new_sam_name = fastq_file.replace("unaligned_reads.fastq", new_sam_suffix)
    #print(new_sam_name)
    new_log_suffix = os.path.basename(bowtie2_idx) + ".log"
    new_log_name = fastq_file.replace("unaligned_reads.fastq", new_log_suffix)
    unaligned_to_HSVd = fastq_file.replace("unaligned_reads.fastq", "unaligned_to_HSVd_reads.fastq")

    command = f'hqsub "bowtie2 -x {bowtie2_idx} -U {fastq_file} -S {new_sam_name} --un {unaligned_to_HSVd}" -q "*" -r {new_log_name}'

    print(command)

