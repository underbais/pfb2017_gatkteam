#!/usr/bin/env python3

import os, sys, re
import subprocess

usage = "\n\n\tusage:{} refpath bampath conf outpath\n\n".format(sys.argv[0])

if len(sys.argv) <= 3:
    sys.stderr.write(usage)
    sys.exit(1)

#get input files
refpath = sys.argv[1]
bampath = sys.argv[2]
conf = int(sys.argv[3])
outpath = sys.argv[4]

def run_gatk(refpath, bampath, conf, outpath = "vcfout"):
	execution = 0
	try:
        # create bwa output folder
		if not os.path.exists(outpath):
			os.makedirs(vcfout)

		#run gatk HaplotypeCaller
		for bamfile in os.listdir(bampath):
			if bamfile.endswith('.bam'):
				out_file_name = bamfile.split(".ba")[0]
				out_file =  open("{}.vcf".format(out_file_name), "w")
				cmd = "gatk -T HaplotypeCaller -R {} -I {} --genotyping_mode DISCOVERY --stand_call_conf {} -o {}".format(refpath, bamfile, conf, out_file)
				proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        		    	#get output and errors
				outs,errs = proc.communicate()
				out_log =  open("{}.log".format(out_file_name), "w")

	        	 	# write output (from byte to ascii) in .sam
				outs_dec = outs.decode("ascii")
				outs_errs = errs.decode("ascii")
				out_file.write(outs_dec)
				out_log.write(outs_errs)
				out_file.close()
				out_log.close()
	except:
        	execution = 1
	return execution  # return the proc object inclunding binary output and error log

run_gatk(refpath, bampath, conf, outpath)

sys.exit(0)
