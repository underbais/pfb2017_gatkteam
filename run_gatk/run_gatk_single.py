
#!/usr/bin/env python3

import os, sys, re
import subprocess

usage = "\n\n\tusage:{} refpath bampath conf outpath\n\n".format(sys.argv[0])

if len(sys.argv) <= 4:
    sys.stderr.write(usage)
    sys.exit(1)

#get input  params
refpath = sys.argv[1]
bamfile = sys.argv[2]
conf = int(sys.argv[3])
outpath = sys.argv[4]

def run_gatk(refpath, bamfile, conf, outpath = "vcfout"):
	execution = 0
	try:
		if not os.path.exists(outpath):
			os.makedirs(outpath)

		#run gatk HaplotypeCaller
		out_file_name = os.path.basename(bamfile).split(".ba")[0]
		out_file = "{}/{}.vcf".format(outpath,out_file_name)
		cmd = "gatk -T HaplotypeCaller -R {} -I {} --genotyping_mode DISCOVERY -stand_call_conf {} -o {}".format(refpath, bamfile, conf, out_file)
		proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		#get output and errors
		outs,errs = proc.communicate()
		out_log =  open("{}.log".format(out_file), "w")

		#write output (from byte to ascii) in .sam
		out_errs = errs.decode("ascii")
		out_log.write(out_errs)
		out_log.close()
	except Exception as e:
		print(str(e))
		execution = 1
	return execution  #return the proc object inclunding binary output and error log

test = run_gatk(refpath, bamfile, conf, outpath)
print(test)
sys.exit(0)
