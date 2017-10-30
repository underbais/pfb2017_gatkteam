




#!/usr/bin/env python3

# Script to create a function that will sort fastq seqs and convert .sam files to .bam files.
# We will use the function as one step in our 'pipeline'.

import sys
import os
import subprocess



input_sam = sys.argv[1]
sample_name = os.path.splitext(input_sam)[0]


#def run_samtools():
def call_samtools_sort(input_sam, threads = 1, output_format = "BAM", samtools_output = "samtools_output"):
    execution = 0
    try:

        # get source file path:
        # create samtools output folder
        if not os.path.exists(samtools_output):
            os.makedirs(samtools_output)

        # run samtools_sort
#        cmd = "samtools sort -@ 4 ERR1977352.sam -O BAM -o ERR1977352.sorted.bam 2> samtools_52sort.err"
        cmd = "samtools sort -@ {} {} -O {} -o {}/{}.sorted.bam".format(threads, input_sam, output_format, samtools_output, sample_name)
        proc = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(cmd)
        # get errors
        outs, errs = proc.communicate()
        out_file = open("{}.sorted.bam.stdout".format(sample_name), "w")
        out_log = open("{}.sorted.bam.stderr.log".format(sample_name), "w")


    except:
        execution = 1
        return execution
        print("function failed")


#os.system("samtools)

#run(["samtools", "-S", "-b", fh

k=call_samtools_sort(input_sam)
print(k)
sys.exit(0)










#bowtie2.bam  bowtie2.sam.gz
