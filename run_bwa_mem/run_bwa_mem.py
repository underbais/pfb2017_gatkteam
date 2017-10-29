#! /usr/bin/env python3
# python3 run_bwa_mem.py ./../ref/hg38.fa ./../../cunde/FASTQ/ERR1977350_1.fastq ./../../cunde/FASTQ/ERR1977350_2.fastq > test

import os, sys, re
import subprocess

usage = "\n\n\tusage:{} refgen fastq_r1 fastq_r2 my_path\n\n".format(sys.argv[0]) 
      
if len(sys.argv) <= 3:
    sys.stderr.write(usage)
    sys.exit(1)

#get input files
refgen = sys.argv[1]
fastq_r1 = sys.argv[2]
fastq_r2 = sys.argv[3]

def call_bwa(refgen, fastq_r1, fastq_r2, threads = 1, mark_split = "-M", bwa_output = "bwa_mapping_output", my_path = './'):
    #execution = 0
    try:            
        # create bwa output folder
        if not os.path.exists(bwa_output):
            os.makedirs(bwa_output)
        out_file_name = fastq_r1.split("_")[0]
        file_out = "{}/{}.sam".format(my_path,out_file_name)
        log_out =  "{}/{}.log".format(my_path,out_file_name)

        #run bwa mem
        cmd = "bwa mem {} -t {} {} {} {}".format(mark_split, threads, refgen,fastq_r1,fastq_r2)
        print(cmd)
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        #get output and errors
        outs,errs = proc.communicate()
        #out_file = open(os.path.join(my_path,bwa_output), "w")
        out_file =  open(file_out, "w")
        out_log =  open(log_out, "w")

        # write output (from byte to ascii) in .sam
        outs_dec = outs.decode("ascii")
        outs_errs = errs.decode("ascii")
        out_file.write(outs_dec)
        out_log.write(outs_errs)
        out_file.close()
        out_log.close()
        execution = file_out

        # check exception content
    except Exception as e:
        print(str(e))
        execution = None
    return execution  # return the proc object inclunding binary output and error log

# get function output and  error:
filepath = call_bwa(refgen, fastq_r1, fastq_r2)
if filepath is None:
    print('"something went wrong"')
    sys.exit(1)
print(filepath)
sys.exit(0)

