#!/usr/bin/env python3

# Script to create a function that will take a sorted .bam file and mark/remove PCR duplicates.
# Function will return a '.mrkdup.sorted.bam' file
# We will use the function as one step in our 'pipeline'.

import sys
import os
import subprocess



input_bam = sys.argv[1]
sample_name = (os.path.basename(input_bam)).split('.')[0]


def call_picard_mrdkup(input_bam, picard_output = "picard_output"):
    execution = 0
    try:

        # create picard output folder
        if not os.path.exists(picard_output):
            os.makedirs(picard_output)

        # run samtools_sort
#        cmd = java -Djava.io.tmpdir=./ -jar /usr/local/anaconda/share/picard-2.14-0/picard.jar MarkDuplicates I=ERR1977350.sorted.bam O=ERR1977350.mrkdup.sorted.bam M=ERR1977350.mrkdup_metrics.txt
        picard_path = 'java -Djava.io.tmpdir=./ -jar /usr/local/anaconda/share/picard-2.14-0/picard.jar'
        cmd = "{} MarkDuplicates I={} O={}/{}.mrkdup.sorted.bam M={}/{}.mrkdup_metrics.txt".format(picard_path, input_bam, picard_output, sample_name, picard_output, sample_name)
        proc = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(cmd)
        
        # get errors
        outs, errs = proc.communicate()
        out_file = open("{}/{}.mrkdup.sorted.bam.stdout".format(picard_output, sample_name), "w")
        out_log = open("{}/{}.mrkdup.sorted.bam.stderr.log".format(picard_output, sample_name), "w")

    except:
        execution = 1
        return execution
        print("function failed")

k=call_picard_mrdkup(input_bam)
print(k)
sys.exit(0)










#bowtie2.bam  bowtie2.sam.gz
