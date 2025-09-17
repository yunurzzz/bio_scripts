#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import os
import argparse
import gzip
import time
import logging


'''
trim_galore -> STAR -> featureCounts
paired data
'''

def time_print(info):
	#print("\033[32m%s\033[0m %s"%(time.strftime('[%H:%M:%S]',time.localtime(time.time())), info))
	print("%s %s"%(time.strftime('[%H:%M:%S]',time.localtime(time.time())), info))
	#logging.info()


def get_opts():
	group = argparse.ArgumentParser()
	group.add_argument("-s", "--sample_list", help="Input sample file, .gz supported", required=True)
	#group.add_argument("--gzip", help="gzip clean data, if run in hardware, dont set this", action="store_true")
	#group.add_argument("--ncbi_data", help="If data download from ncbi, set this", action="store_true")
	# /home/dodo/mRNAseq/reference
	group.add_argument("-i", "--star_index", help="Index of STAR. Required.", required=True)
	group.add_argument("-f", "--fastq_path", help="Path of fastq data", required=True)
	group.add_argument("-a", "--annotation", help="Path of GTF/GFF annotation file", required=True)
	group.add_argument("-o", "--outdir", help="outdir. Default: outdir", default="outdir")
	#group.add_argument("--step", help="[all|aln|count], Specify which steps you want to run this script. all: run the entire pipeline (default); aln: start from alignment to the end; count: start from filtered TEs to finalizing the run.")
	#group.add_argument("--gzip", help="[0|1], gzip clean data. If run in hardware, donot set this. Default: 0", default=0)
	group.add_argument("-t", "--threads", help="Number of threads to run this script. Default: 12", default=12, type=int)
	return group.parse_args()


def rnaseq_pipeline(sample, index, fqpath, annotation, outdir, threads):
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	sample = os.path.abspath(sample)
	index = os.path.abspath(index)
	fqpath = os.path.abspath(fqpath)
	annotation = os.path.abspath(annotation)

	script_dir = sys.path[0]
	time_print("Entering: %s"%outdir)
	os.chdir(outdir)

	###### clean ######
	#trimgalore_path = "/home/dodo/TrimGalore-0.6.6"
	with open(sample, 'r') as sample_list:
		time_print("Step1: running trim_galore")
		clean_dir = "01.clean"
		if not os.path.exists(clean_dir):
			os.mkdir(clean_dir)
		else:
			time_print("trim_galore result found, rerun")
			os.system("rm -rf ./01.clean/*")
		os.chdir("01.clean")
		for line in sample_list:
			sample_data = line.strip().split()
			sample1 = sample_data[0]
			sample2 = sample_data[1]
			#out_pre = sample_list.strip()
			#cmd = "%s/trim_galore --illumina --fastqc --paired -j 4 -q 30 --length 25 --dont_gzip %s/%s %s/%s >> clean.log 2>&1"%(trimgalore_path, fqpath, sample1, fqpath, sample2)
			cmd = "trim_galore --illumina --fastqc --paired -j 4 -q 30 --length 25 --dont_gzip %s/%s %s/%s >> clean.log 2>&1"%(fqpath, sample1, fqpath, sample2)
			time_print("Running: %s"%cmd)
			os.system(cmd)
		time_print("clean done")
		os.chdir("..")
		#print(os.system('pwd'))

	###### aln ######
	time_print("Step2: running STAR")
	aln_dir = "02.aln"
	if not os.path.exists(aln_dir):
		os.mkdir(aln_dir)
	else:
		time_print("STAR result found, rerun")
		os.system("rm -rf ./02.aln/*")
	os.chdir("02.aln")
	#print(os.system('pwd'))
	os.system("ls ../01.clean/*R1_val_1.fq* > R1.txt")
	os.system("ls ../01.clean/*R2_val_2.fq* > R2.txt")
	os.system("paste R1.txt R2.txt > clean.txt")
	os.system("rm R1.txt R2.txt")
	with open('clean.txt', 'r') as clean_list:
		for line in clean_list:
			clean_data = line.strip().split()
			clean1 = clean_data[0]
			#print(clean1)
			clean2 = clean_data[1]
			out_prx_tmp = (clean1.split('/'))[-1]
			#print(out_prx_tmp)
			out_prx = out_prx_tmp[:-12]
			#print(out_prx)
			#out_pre = sample_list.strip()
			cmd = "STAR --twopassMode Basic \
					--quantMode TranscriptomeSAM GeneCounts \
					--runThreadN %s \
					--genomeDir %s \
					--alignIntronMin 20 \
					--alignIntronMax 50000 \
					--outSAMtype BAM SortedByCoordinate \
					--sjdbOverhang 150 \
					--outFilterMismatchNmax 2 \
					--outSJfilterReads Unique \
					--outSAMmultNmax 1 \
					--outFileNamePrefix %s \
					--outSAMmapqUnique 60 \
					--readFilesIn %s %s >> aln.log 2>&1"%(threads, index, out_prx, clean1, clean2)
			time_print("Running: %s"%cmd)
			os.system(cmd)
		os.system("rm clean.txt")
		time_print("align done")
		os.chdir("..")

	###### count ######
	time_print("Step3: running featureCounts")
	aln_dir = "03.count"
	if not os.path.exists(aln_dir):
		os.mkdir(aln_dir)
	else:
		time_print("featureCounts result found, rerun")
		os.system("rm -rf ./03.count/*")
	os.chdir("03.count")
	os.system("ln -s ../02.aln/*Aligned.sortedByCoord.out.bam ./")
	cmd = "featureCounts -T %s -p -t exon -g gene_id -a %s -o result.count *Aligned.sortedByCoord.out.bam >> count.log 2>&1"%(threads, annotation)
	time_print("Running: %s"%cmd)
	os.system(cmd)
	os.system("cut -f1,7- result.count > count.tsv")
	os.system("rm *.bam")
	time_print("count done")
	os.chdir("..")

	###### gzip ######
	#time_print("Step4: gzip clean fastq")
	os.chdir("01.clean")
	cmd = "pigz -p %s *.fq"%(threads)
	os.system(cmd)
	#time_print("gzip done")

	time_print("Finished")


if __name__ == "__main__":
	opts = get_opts()
	sample = opts.sample_list
	'''
	if opts.ncbi_data:
		ncbi_data = True
	else:
		ncbi_data = False
	'''
	index = opts.star_index
	fqpath = opts.fastq_path
	annotation = opts.annotation
	outdir = opts.outdir
	threads = opts.threads
	rnaseq_pipeline(sample, index, fqpath, annotation, outdir, threads)
