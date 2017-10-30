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
        #print(file_out)
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
        #print(str(e))
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
#print(input_sam)
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
        file_out = "{}/{}/{}.sorted.bam".format(my_path,samtools_output,sample_name)
        #print(file_out)
        log_out =  "{}/{}.log".format(my_path,samtools_output)

        # run samtools_sort
        cmd = "samtools sort -@ {} {} -O {} -o {}".format(threads, input_sam, output_format, file_out)
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
        #print(str(e2))
        execution = None
    return execution 

# get function output and  error: 
filepath_samtool = call_samtools_sort(input_sam)
if filepath_samtool is None:
    print('"something went wrong in samtools"')
    sys.exit(1)
#print(filepath_samtool)
#sys.exit(0)

#### PICARD MD function
####################
# Function will return a '.mrkdup.sorted.bam' file

input_bam = filepath_samtool
#print(input_bam)
sample_name = (os.path.basename(input_bam)).split('.')[0]
#print(sample_name)

#my_path = os.getcwd()
#print('current_dir', dir_path)

def call_picard_md(input_bam, picard_MD_output = "picard_output", my_path = "./"):
    execution = 0
    try:

        # create picard output folder
        if not os.path.exists(picard_MD_output):
            os.makedirs(picard_MD_output)
           
        file_out = "{}/{}/{}.mrkdup.sorted.bam".format(my_path,picard_MD_output,sample_name)
        #print(file_out)
        log_out =  "{}/picard.mrkdup.sorted.bam.log".format(my_path)
        #print(log_out)

        # run samtools_sort
        picard_path = 'java -Djava.io.tmpdir=./ -jar /usr/local/anaconda/share/picard-2.14-0/picard.jar'
        cmd = "{} MarkDuplicates I={} O={} M={}.metrics.txt".format(picard_path, input_bam, file_out, file_out)
        proc = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #print(cmd)
        
        # get errors
        outs, errs = proc.communicate()
        out_errs = errs.decode("ascii")
        out_log = open(log_out, "w")
        out_log.write(out_errs)
        out_log.close()
        execution = file_out
        
    except Exception as e3:
        print(str(e3))
        execution = None
    return execution

# get function output and  error:
filepath_picard_md = call_picard_md(input_bam)
if filepath_picard_md is None:
    print('"something went wrong in Picard_mark_duplicates"')
    sys.exit(1)
    print(filepath_picard_md)
#sys.exit(0)
            
##### PICARD function ARRG
##########################

input_picard_md = filepath_picard_md
#print(input_picard_md)
sample_name = (os.path.basename(input_picard_md)).split('.')[0]
#sample_name = os.path.basename(input_bam).split(".")
#print(sample_name)

outpath = sys.argv[2]

def call_picard_groups(input_picard_md, picard_arrg_out = "picard_arrg_out", my_path = "./"):
    execution = 0
    try:
        if not os.path.exists(picard_arrg_out):
            os.makedirs(picard_arrg_out)
        file_out = "{}/{}/{}.arrg.mrkdup.sorted.bam".format(my_path,picard_arrg_out,sample_name)
        #print(file_out)
        log_out =  "{}/picard.arrg.mrkdup.sorted.bam.log".format(my_path)
        #print(log_out)
            
	#run gatk HaplotypeCaller
        cmd = "java -jar /usr/local/anaconda/share/picard-2.14-0/picard.jar AddOrReplaceReadGroups I={} O={} RGLB=lib1 RGPL=illumina RGPU=dummy RGSM={}".format(input_picard_md, file_out, sample_name)
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        #get output and errors
        outs,errs = proc.communicate()
        
	#write output (from byte to ascii)
        out_errs = errs.decode("ascii")
        out_log = open(log_out, "w")
        out_log.write(out_errs)
        out_log.close()
        execution = file_out
        
    except Exception as e:
        print(str(e))
        execution = 1
    return execution  #return the proc object inclunding binary output and error log

filepath_picard_arrg = call_picard_groups(input_picard_md)
print(filepath_picard_arrg)

if filepath_picard_arrg is None:
    print('"something went wrong in Picard_ARRG"')
    sys.exit(1)
print(filepath_picard_arrg)
sys.exit(0)


#### Picard CSD_FAIDX
#####################

