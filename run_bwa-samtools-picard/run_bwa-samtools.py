#! /usr/bin/env python3
# python3 run_bwa_mem.py ./../ref/hg38.fa ./../../cunde/FASTQ/ERR1977350_1.fastq ./../../cunde/FASTQ/ERR1977350_2.fastq > test

import os, sys, re
import subprocess

usage = "\n\n\tusage:{} refgen fastq_r1 fastq_r2 my_path\n\n".format(sys.argv[0]) 
      
if len(sys.argv) <= 3:
    sys.stderr.write(usage)
    sys.exit(1)

#### BWA function
#################

#get input files
refgen = sys.argv[1]
fastq_r1 = sys.argv[2]
fastq_r2 = sys.argv[3]

def call_bwa(refgen, fastq_r1, fastq_r2, threads = 1, mark_split = "-M", bwa_output = "bwa_output", my_path = './'):
    #execution = 0
    try:            
        # create bwa output folder
        if not os.path.exists(bwa_output):
            os.makedirs(bwa_output)
        out_file_name = fastq_r1.split("_")[0]
        file_out = "{}/{}/{}.sam".format(my_path,bwa_output,out_file_name)
        print(file_out)
        log_out =  "{}/{}.log".format(my_path,bwa_output)

        #run bwa mem
        cmd = "bwa mem {} -t {} {} {} {} > {}".format(mark_split, threads, refgen,fastq_r1,fastq_r2,file_out)
        #print(cmd)
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        #get output and errors
        outs,errs = proc.communicate()
        #out_file = open(os.path.join(my_path,bwa_output), "w")
        out_log =  open(log_out, "w")

        # write output (from byte to ascii) in .sam
        outs_errs = errs.decode("ascii")
        out_log.write(outs_errs)
        out_log.close()
        execution = file_out

        # check exception content
    except Exception as e:
        print(str(e))
        execution = None
    return execution  # return the proc object inclunding binary output and error log

# get function output and  error:
filepath_bwa = call_bwa(refgen, fastq_r1, fastq_r2)
if filepath_bwa is None:
    print('"something went wrong in bwa"')
    sys.exit(1)
#print('out_bwa',filepath_bwa)
#sys.exit(0)

#### SAM function
#################

input_sam = filepath_bwa
sample_name = os.path.basename(input_sam)
sample_name = os.path.splitext(sample_name)[0]
#print(sample_name)

#def run_samtools():
def call_samtools_sort(input_sam, threads = 1, output_format = "BAM", samtools_output = "samtools_output", my_path = './'):
    #execution = 0
    try:
        # create samtools output folder
        if not os.path.exists(samtools_output):
            os.makedirs(samtools_output)       
        file_out = "{}/{}.sam".format(my_path,samtools_output,sample_name)
        log_out =  "{}/{}.log".format(my_path,samtools_output)

        # run samtools_sort
        cmd = "samtools sort -@ {} {} -O {} -o {}/{}.sorted.bam".format(threads, input_sam, output_format, samtools_output, sample_name)
        proc = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #print(cmd)
        outs,errs = proc.communicate()
        out_errs = errs.decode("ascii")
        out_log =  open(log_out, "w")
        #print(type(outs))
        out_log.write(out_errs)
        out_log.close()
        execution = file_out
            
    except Exception as e2:
        print(str(e2))
        execution = None
    return execution 

# get function output and  error: 
filepath_samtool = call_samtools_sort(input_sam)
if filepath_samtool is None:
    print('"something went wrong in samtools"')
    sys.exit(1)
#print(filepath_samtool)
sys.exit(0)

#### PICARD function
####################
# Function will return a '.mrkdup.sorted.bam' file

input_bam = filepath_samtool
#print(input_bam)
sample_name = (os.path.basename(input_bam)).split('.')[0]
#print(sample_name)

def call_picard_mrkup(input_bam, picard_output = "picard_output"):
    execution = 0
    try:

        # create picard output folder
        if not os.path.exists(picard_output):
            os.makedirs(picard_output)
            
        file_out = "{}/{}/{}.mrkdup.sorted.bam.stdout".format(my_path,picard_output,sample_name)
        print(file_out)
        log_out =  "{}/{}.mrkdup.sorted.bam.stderr.log".format(my_path,picard_output)

        # run samtools_sort
        picard_path = 'java -Djava.io.tmpdir=./ -jar /usr/local/anaconda/share/picard-2.14-0/picard.jar'
        cmd = "{} MarkDuplicates I={} O={}/{}.mrkdup.sorted.bam M={}/{}.mrkdup_metrics.txt".format(picard_path, input_bam, picard_output, sample_name, picard_output, sample_name)
        proc = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(cmd)
        
        # get errors
        outs, errs = proc.communicate()
        out_errs = errs.decode("ascii")
        out_log = open(log_out, "w")
        out_log.write(out_errs)
        out_log.close()
        
    except Exception as e3:
        print(str(e3))
        #execution = None
        return execution

# get function output and  error:
filepath_picard_mrkup = call_picard_mrkup(input_bam)
if filepath_picard_mrkup is None:
    print('"something went wrong in Picard_mark_duplicates"')
    sys.exit(1)
    #print(filepath_samtool)
sys.exit(0)
            
