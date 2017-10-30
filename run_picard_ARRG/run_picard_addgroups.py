
#!/usr/bin/env python3

import os, sys, re
import subprocess

usage = "\n\n\tusage:{} input_bam \n\n".format(sys.argv[0])

if len(sys.argv) <= 1:
    sys.stderr.write(usage)
    sys.exit(1)

#get input  params

input_bam = sys.argv[1] #taking absolute path to input BAM file
outpath = sys.argv[2]

def run_picard_groups(input_bam, outpath = "picard_groups_bam_out"):
	execution = 0
	try:
		if not os.path.exists(outpath):
			os.makedirs(outpath)

		#run gatk HaplotypeCaller
		out_file_name  = os.path.basename(input_bam).split(".") # list of split input bam file name
		out_file = "{}/{}.arrg.{}.{}.{}".format(outpath, out_file_name[0], out_file_name[1], out_file_name[2], out_file_name[3])
		cmd = "nohup java -jar /usr/local/anaconda/share/picard-2.14-0/picard.jar AddOrReplaceReadGroups I={} O={} RGLB=lib1 RGPL=illumina RGPU=dummy RGSM={}".format(input_bam, out_file, out_file_name[0])
		proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		#get output and errors
		outs,errs = proc.communicate()
		out_log =  open("{}.{}.{}.log".format(out_file_name[0], out_file_name[1], out_file_name[2]), "w")

		#write output (from byte to ascii) in .sam
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
