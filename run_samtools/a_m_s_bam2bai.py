#!/usr/bin/env python3

# Script to create a function that will take a aarg.mrkdup.sorted.bam file and index it.
# Function will return a 'aarg.mrkdup.sorted.bam.bai' file
# We will use the function as one step in our 'pipeline'.

import sys
import os
import subprocess



input_bam = sys.argv[1]
sample_name = (os.path.basename(input_bam)).split('.')[0]


def call_samtools_index(input_bam, picard_output = "picard_output"):
    execution = 0
    try:

#       cmd = samtools index ERR1977350.arrg.mrkdup.sorted.bam
        cmd = "samtools index {}".format(input_bam)
        proc = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        print(cmd)
        
        # get errors
        outs,errs = proc.communicate()
        out_log = open("{}/{}.mrkdup.sorted.bam.stderr.log".format(picard_output, sample_name), "w")
        out_log.close()

    except:
        execution = 1
        return execution
        print("function failed")

k=call_samtools_index(input_bam)
print(k)
sys.exit(0)

