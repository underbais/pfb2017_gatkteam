#!/usr/bin/env python3

# creating sequence dictionary with Picard

import os, sys, re
import subprocess

usage = "\n\n\tusage:{} input_ref \n\n".format(sys.argv[0])

if len(sys.argv) <= 1:
    sys.stderr.write(usage)
    sys.exit(1)

#get input  params

input_ref = sys.argv[1] #full path to ref fasta

def run_picard_CSD(input_ref):
	execution = 0
	try:
		if not os.path.exists(outpath):
			os.makedirs(outpath)

		#run Picard CSD and samtools faidx on reference fasta file
		out_file_name  = os.path.basename(input_ref).split(".") # list of split input bam file name
		out_file = "{}.dict".format(out_file_name[0])
		cmd1 = "nohup java -jar /usr/local/anaconda/share/picard-2.14-0/picard.jar CreateSequenceDictionary R={} O={}".format(input_ref, out_file)
		cmd2 = "nohup samtools faidx {}".format(input_ref)
		proc1 = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		proc2 = subprocess.Popen(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		#get output and errors
		outs,errs = proc1.communicate()
		out_log =  open("{}.log".format(out_file_name[0]), "w")

		#write output
		out_errs = errs.decode("ascii")
		out_log.write(out_errs)
		out_log.close()
	except Exception as e:
		print(str(e))
		execution = 1
	return execution  #return the proc object inclunding binary output and error log

test = run_picard_groups(input_bam, outpath)
print(test)
sys.exit(0)
