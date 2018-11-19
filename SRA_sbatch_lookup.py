#!/usr/bin/env python3

import sys
import argparse
import os

#Input arguments
parser = argparse.ArgumentParser(description='Download reads from Short Read Archive')
parser.add_argument('-i', dest="input", metavar='Input file', type=str,
                    help='Text file with list of SRA IDs of any kind (Project, sample, run)')
parser.add_argument('-o', dest="output", metavar='Output folder', type=str,
                    help='Path to download folder')
parser.add_argument('-p', dest="partition", metavar='partition',type=str, help='Queue partition to use (daytime or project).',default = 'project')
args = parser.parse_args()


#Create output directory
if not os.path.isdir(args.output):
	os.makedirs(args.output)		#Create directory for output

sample_IDs = []

with open(args.input) as f:
	for line in f:
		line = line.rstrip('\n')
		sample_IDs.append(line)	

f.close()
run_IDs = []


#sbatch -c 6 --mem=23G --time=03:00:00 -J 'ENA_download' -p surveillance --wrap="wget -qO- \'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term='+sample+'\' | grep '+sample+' | cut -f1 -d\",\" > '+tmp_path"

run_to_sample = {}
line_count = 0
for sample in sample_IDs:
	tmp_path = os.path.join(args.output,sample+'_efetch_temp.txt')
	cmd = 'wget -qO- \'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term='+sample+'\' | grep '+sample+' | cut -f1 -d\",\" > '+tmp_path
	os.system(cmd)
	#cmd = 'wget -qO- \'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term='+sample+'\' | grep '+sample+' | cut -f1 -d\",\"'
	#os.system('cat '+tmp_path)
	with open(tmp_path) as f:
		for line in f:
			line_count += 1
			run_ID = line.rstrip('\n')
			run_IDs.append(run_ID)
			if line_count > 1:
				uniq_sample_name = sample+'_'+str(line_count)
			else:
				uniq_sample_name = sample
			run_to_sample[run_ID] = uniq_sample_name
	#os.remove(tmp_path)
	f.close()

print('Downloading read files from '+str(len(run_IDs))+' entries')


#Download reads with fastq-dump
def dl(run_ID,partition):
	cmd1 = '/tools/linuxbrew/bin/fastq-dump --gzip --split-files -O ' + args.output + ' ' + run_ID
	cmd = 'sbatch -c 4 --mem=4G --time=03:00:00 -J \"ENA_download\" -p ' + partition + ' --wrap=\"' + cmd1 + '\"'
	print(cmd)
	os.system(cmd)
	return cmd

if __name__ == '__main__':
	for run_ID in run_IDs:
		dl(run_ID,args.partition)

#for run_ID in run_IDs:
#	sample_ID = run_to_sample[run_ID]
#	r1 = os.path.join(args.output,run_ID+'_1.fastq.gz')
#	r2 = os.path.join(args.output,run_ID+'_2.fastq.gz')
#	r1_new = os.path.join(args.output,sample_ID+'_1.fastq.gz')
#	r2_new = os.path.join(args.output,sample_ID+'_2.fastq.gz')
#	if os.path.exists(r1):
#		os.rename(r1,r1_new)
#		os.rename(r2,r2_new)
	#else:
	#	print(sample_ID)

samples_with_reads = []
print('sample_ID\trun_ID')
for run_ID in run_IDs:
	sample_ID = run_to_sample[run_ID]
	r1 = os.path.join(args.output,run_ID+'_1.fastq.gz')
	r2 = os.path.join(args.output,run_ID+'_2.fastq.gz')
	if os.path.exists(r1) and os.path.exists(r2):
		samples_with_reads.append(sample_ID)
		print(sample_ID+'\t'+run_ID)

#for sample_ID in sample_IDs:
#	if sample_ID not in samples_with_reads:
#		print(sample_ID+'\tMISSING READS')

#for run in run_IDs:
#	cache_path = os.listdir('~/ncbi/public/sra',run_ID+'.sra.cache')
#	if os.path.exists(cache_path):
#		os.remove(cache_path)
#	cache_path = os.listdir('~/ncbi/public/sra',run_ID+'.sra')
#	if os.path.exists(cache_path):
#		os.remove(cache_path)

