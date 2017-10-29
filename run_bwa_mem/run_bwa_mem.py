#! /usr/bin/env python3
# python3 run_bwa_mem.py ./../ref/hg38.fa ./../../cunde/FASTQ/ERR1977350_1.fastq ./../../cunde/FASTQ/ERR1977350_2.fastq > test

import os, sys, re
import subprocess

usage = "\n\n\tusage:{} refgen fastq_r1 fastq_r2\n\n".format(sys.argv[0]) 
      
if len(sys.argv) <= 3:
    sys.stderr.write(usage)
    sys.exit(1)
            
#def callbwa():
    #get input files
refgen = sys.argv[1]
fastq_r1 = sys.argv[2]
fastq_r2 = sys.argv[3]
out_file_name = fastq_r1.split("_")[0]


    #run bwa mem
cmd = "bwa mem -M -t 8 {} {} {}".format(refgen,fastq_r1,fastq_r2)
proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    #get output and errors
outs,errs = proc.communicate()
out_file =  open("{}.sam".format(out_file_name), "w")
out_log =  open("{}.log".format(out_file_name), "w")

    # write output (from byte to ascii) in .sam 
outs_dec = outs.decode("ascii")
outs_errs = errs.decode("ascii")
out_file.write(outs_dec)
out_log.write(outs_errs)
out_file.close()
out_log.close()

sys.exit(0)

