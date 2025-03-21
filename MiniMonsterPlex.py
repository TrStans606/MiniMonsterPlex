# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 13:49:59 2023

@author: treyd
"""
#importing packages
import os
import glob
import subprocess 
import argparse
import multiprocessing

#setting up arg parsing for the output folder
parser = argparse.ArgumentParser(
	prog='MiniMonsterPlex',
	description=(
		'A pipeline for variant call anaylsis of pathogens'
	),
	epilog=(
		'The default output file name is Out'
	),
)

#command line option for output folder
parser.add_argument(
	'-o',
	action='store',
	help=(
		'User defined name for the output file.'
	),
	default='Out'
)

#command line option for metadata
parser.add_argument(
	'-m',
	action='store',
	help=(
		'Name of metadatafile for the strain being used'
	),
	default=None
)

#command line option for alternative input folder

parser.add_argument(
	'-f',
	action='store',
    nargs='?',
	help=(
		'file path to folder containing fastq.gz files for input. Default is fastq/'
	),
	default='fastq',
	required=False
)

#command line option for filtering by specific isolates
parser.add_argument(
	'-i',
	action='store',
    nargs='+',
	help=(
		'spaced list of isolates you want included in Tree Building'
	),
	required=False
)
#command line option for filtering by specific isolates via a txt
parser.add_argument(
	'-il',
	action='store',
    nargs='?',
	help=(
		'new line seperated txt file of isolates you want in tree building. Argument should be file path'
	),
	required=False
)
#command line option for filtering by specific hosts
parser.add_argument(
	'-hf',
	action='store',
    nargs='*',
	help=(
		'spaced list of host(s) you want included in Tree Building'
	),
	required=False
)

#command line option for filtering by specific hosts
parser.add_argument(
	'-hfl',
	action='store',
    nargs='?',
	help=(
		'new line seperated txt file of host(s) you want in tree building. Argument should be file path'
	),
	required=False
)

#command line option for uncompressed files
#parser.add_argument(
#	 '-gz',
#	 action='store_true',
#	 help=(
#		 'Are your fastq files gzipped?'
#	 ),
#	 default=False
#)

args = parser.parse_args()

outPut_Folder = args.o
metadata_file_name = args.m
input_folder = args.f
included_isolates = args.i
included_isolates_file = args.il
included_hosts = args.hf
included_hosts_file = args.hfl
threads = multiprocessing.cpu_count()
if threads > 8:
	threads =8
#gzipped = args.gz

def auto_bowtie2(outPut, fileNum,threads):
	print(fileNum, " is entering the pipeline")
	#histat 2 + samtools sort call
	command = ['bowtie2 --no-unal',
			'-p',
			str(threads),
			'-x',
			'index/70-15_small_index',
			'-U',
			file,
			'--local --very-sensitive-local',
			'|',
			'samtools',
			'sort',
			'-',
			'-@',
			'2',
			'-O',
			'bam',
			'-o',
			f'{outPut}/{fileNum}hits.bam']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)

def auto_tabix(outPut, fileNum,file_ex,run):
	command = ['tabix',
			f'{outPut}/{fileNum}{file_ex}']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	if run == 1:
		with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("tabix1 done")

	if run == 2:
		with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("tabix2 done")

def auto_mpileup(outPut,fileNum,threads):
	command = ['bcftools',
			'mpileup',
			'--threads',
			str(threads),
			'-d',
			'100000',
			#'-R',
			#'MonsterPlexRegionsFileSuperCont_small_index.txt',
			'--annotate',
			'FORMAT/AD',
			'-f',
			'index/70-15_small.fasta',
			f'{outPut}/{fileNum}hits.bam',
			'>>',
			f'{outPut}/{fileNum}.vcf']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("mpileup done")
	
def auto_call(outPut,fileNum):
	command = ['bcftools',
			'call',
			'-c',
			'--ploidy',
			'1',
			f'{outPut}/{fileNum}.vcf',
			'-o',
			f'{outPut}/seperateCall/{fileNum}call.vcf']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("call done")
	
def auto_bedtools(outPut,fileNum):
	command = ['bedtools',
			'genomecov',
			'-ibam',
			f'{outPut}/{fileNum}hits.bam',
			'-bg',
			'>',
			f'{outPut}/Coverage/{fileNum}cover.bed']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("bedtools done")

def auto_bgzip(outPut, fileNum,file_ex):
	command =['bgzip',
			f'{outPut}/{fileNum}{file_ex}']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("bgzip done")

def autoMerge(outPut, file, fileNum):
	#bg zip the bcftools call result file
	command = ['bgzip',
			f'{outPut}/seperateCall/{fileNum}call.vcf']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	#tabix the call results
	command = ['tabix',
			f'{outPut}/seperateCall/{fileNum}call.vcf.gz']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	with open(f'{outPut}/fastqListCall.txt', 'a') as append:
		append.write(f'{outPut}/seperateCall/' + file.split('/')[1].split('.')[0] + 'call.vcf.gz\n')
	with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("merge done")
		
def sampleBuilder(outPut):
	command = ['bcftools',
			'merge',
			'-l',
			f'{outPut}/fastqListCall.txt',
			'-o',
			f'{outPut}/seperateCall/{outPut}MergedCallAll.vcf']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	
	sites =[]
	sitesUsed =[]
	#reads a list of sites you want and only looks at data from there
	with open('MonsterPlexSitesList_small_index.txt', 'r') as read:
		for line in read:
			sites.append(line.strip('\n'))

	with open(f'{outPut}/seperateCall/{outPut}MergedCallAll.vcf', 'r') as read:
		seqs = list()
		check = False
		for line in read:
			if line.split('\t')[0] == '#CHROM':
				print("header past seqs made")
				fqList = line.split('\t')
				for n in range(9,len(fqList)):
					seqs.append([fqList[n], ''])
					check = True
			#elif check:
			elif check and line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1] in sites:
				#this creates a horizontal split of the line
				lineList = line.strip('\n').split('\t')
				for n in range(9,len(lineList)):
					fields = lineList[n].split(':')
					if fields[0] == '.':
						seqs[n - 9][1] += "N"
					elif len(fields[2].split(',')) == 1:
						if fields[0] == '0':
							if int(fields[2]) > 5:
								seqs[n - 9][1] += lineList[3]
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							else:
								seqs[n - 9][1] += "N"
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							#this checks alt
						elif fields[0] == '1':
							if int(fields[2]) > 5:
								seqs[n - 9][1] += lineList[4]
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							else:
								seqs[n - 9][1] += "N"
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						else:
							print("something went wrong with " + seqs[n][0] + '\n' + lineList[n])
							print("at site " + line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							quit()
					#this checks cases were both ref and alt are registered
					elif len(fields[2].split(',')) >= 2:
						#this creates a list out of the AD field
						AD = fields[2].split(',')
						#this checks if ref is blank
						if AD[0] == '.':
							if int(AD[1]) > 5:
								seqs[n - 9][1] += lineList[4].split(',')[0]
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							else:
								seqs[n - 9][1] += "N"
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						#this checks if alt is blank
						elif AD[1] == '.':
							if int(AD[0]) > 5:
								seqs[n - 9][1] += lineList[3]
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							else:
								seqs[n - 9][1] += "N"
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						#checks if ref is greater then alt
						elif int(AD[0]) > int(AD[1]):
							if int(AD[0]) > (int(AD[1]) * 20):
								seqs[n - 9][1] += lineList[3]
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							else:
								seqs[n - 9][1] += "N"
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						#checks if alt is greater than ref
						elif int(AD[1]) > int(AD[0]):
							if int(AD[1]) > (int(AD[0]) * 20):
								seqs[n - 9][1] += lineList[4].split(',')[0]
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							else:
								seqs[n - 9][1] += "N"
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						elif int(AD[1]) == int(AD[0]):
							seqs[n - 9][1] += lineList[3]
							sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						else:
							print("something went wrong with " + seqs[n][0] + '\n' + lineList[n])
							print("at site " + line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							quit()
					else:
						print("something went wrong with " + seqs[n][0] + '\n' + lineList[n])
						print("at site " + line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						quit()
	sample_metadata = metaDataBuilder(metadata_file_name)
	
	print(sample_metadata)
	
	os.mkdir(f'{outPut}/built_fasta')
	
	with open(f'{outPut}/built_fasta/{outPut}builtSeqMeta.fasta', 'a') as writeSeq:
		for read in seqs:
			seqID = read[0].split('/')[1].split('.')[0].split('hits')[0]
			if len(seqID.split("_")) > 1:
				seqID = f'{"-".join(seqID.split("_"))}'
			print(seqID)
			if (seqID) in sample_metadata:
				seqSpecies = sample_metadata[seqID][0]
				seqHost = sample_metadata[seqID][1]
				seqLineage = sample_metadata[seqID][2]
				seqCountry = sample_metadata[seqID][3]
				writeSeq.write(f'>{seqID}_{seqSpecies}_{seqHost}_{seqLineage}_{seqCountry}\n{read[1]}\n')
			else:
				writeSeq.write('>' + seqID
							   + '_._._._.' + '\n' + read[1] + '\n')
	with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("samplebuilder done")
						
def metaDataBuilder(metadata_file):
	metaData = {}
	with open(metadata_file, 'r') as read:
		for line in read:
			ID = line.split(',')[0].strip('\n')
			species = line.split(',')[1].strip('\n')
			host = line.split(',')[2].strip('\n')
			lineage = line.split(',')[3].strip('\n')
			country = line.split(',')[4].strip('\n')
			metaData[ID] = [species, host, lineage, country]
		return metaData

#filters the built seq meta by isolate
def fasta_filter(outPut,included_isolates):
	to_write = []
	with open(f'{outPut}/built_fasta/{outPut}builtSeqMeta.fasta','r') as read:
		lines = read.readlines()
		for i in range(0,len(lines)):
			if lines[i][0] == '>':
				if lines[i].split('_')[0].split(">")[1].strip() in included_isolates:
					to_write.append([lines[i],lines[i+1]])
	with open(f'{outPut}/built_fasta/{outPut}builtSeqFiltered.fasta','a') as write:
		for isolate in to_write:
			write.write(f'{isolate[0]}{isolate[1]}')
	with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("filtering done")
        
#filters the built seq meta by host
def fasta_filter_hosts(outPut,included_hosts,filtered):
	#if this has been pre filtered by isolate it adds a second layer of filtering
	if filtered:
		to_write = []
		with open(f'{outPut}/built_fasta/{outPut}builtSeqFiltered.fasta','r') as read:
			lines = read.readlines()
			for i in range(0,len(lines)):
				if lines[i][0] == '>':
					if len(lines[i].split('_')) > 2 and lines[i].split('_')[3].strip() in included_hosts:
						to_write.append([lines[i],lines[i+1]])
		with open(f'{outPut}/built_fasta/{outPut}builtSeqFiltered2.fasta','a') as write:
			for isolate in to_write:
				write.write(f'{isolate[0]}{isolate[1]}')
		os.remove(f'{outPut}/built_fasta/{outPut}builtSeqFiltered.fasta')
		os.rename(f'{outPut}/built_fasta/{outPut}builtSeqFiltered2.fasta', f'{outPut}/built_fasta/{outPut}builtSeqFiltered.fasta')
		with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("filtering done")
	#otherwise it acts the same as fasta_filter but with hosts
	else:
		to_write = []
		with open(f'{outPut}/built_fasta/{outPut}builtSeqMeta.fasta','r') as read:
			lines = read.readlines()
			for i in range(0,len(lines)):
				if lines[i][0] == '>':
					if len(lines[i].split('_')) > 2 and lines[i].split('_')[2].strip() in included_hosts:
						to_write.append([lines[i],lines[i+1]])
		with open(f'{outPut}/built_fasta/{outPut}builtSeqFiltered.fasta','a') as write:
			for isolate in to_write:
				write.write(f'{isolate[0]}{isolate[1]}')
		with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("filtering done")


def autoRAxML(outPut,filtered):
	os.mkdir(f'{outPut}/RAXML_results')
	if filtered==False:
	#command for running RAXML
		command = ['raxmlHPC',
			 '-p',
			 '1234',
			 '-f',
			 'a',
			 '-x',
			 '1234',
			 '-s',
			 f'{outPut}/built_fasta/{outPut}builtSeqMeta.fasta',
			 '-n',
			 'miniMonsterPlex.raxml',
			 '-m',
			 'GTRGAMMA',
			 '-#',
			 '1000']
		subprocess.run(' '.join(command),
				 shell=True,
				 check=True)
		subprocess.run(f'mv *.raxml {outPut}/RAXML_results/',
				 shell=True,
				 check=True)
		with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("autoraxml done")
	else:
		command = ['raxmlHPC',
			  '-p',
			   '1234',
			   '-f',
			   'a',
			   '-x',
			   '1234',
			   '-s',
			   f'{outPut}/built_fasta/{outPut}builtSeqFiltered.fasta',
			   '-n',
			   'miniMonsterPlex.raxml',
			   '-m',
			   'GTRGAMMA',
			   '-#',
			   '1000']
		subprocess.run(' '.join(command),
				  shell=True,
				  check=True)
		subprocess.run(f'mv *.raxml {outPut}/RAXML_results/',
				 shell=True,
				 check=True)
		with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("autoraxml done")

def raxmlGate(outPut_Folder,filtered):
	if filtered:
		with open(f'{outPut_Folder}/built_fasta/{outPut_Folder}builtSeqFiltered.fasta', 'r') as read:
			cnt = 0
			for line in read:
				if line[0] == ">":
					cnt += 1
		if cnt < 4:
			user_input = ''
			while user_input != 'y':
				print("Tree building requires a minium of 4 isolates and you have less then 4 in your filtering")
				user_input= input("Would you like to build a tree with all isolates(y) or quit now(n)? (y/n)")
				if user_input == 'n':
					quit()
			filtered = False
			return(filtered)
		else:
			filtered = True
			return(filtered)
	else:
		return(False)

def mlTree(outPut_Folder):
	command = ['Rscript',
		   '--vanilla',
		   'MLtree.R',
		   f'{outPut_Folder}/RAXML_results/RAxML_bestTree.miniMonsterPlex.raxml']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	command = ['mv',
	        'NA.pdf',
	        f'{outPut_Folder}/RAXML_results/{outPut_Folder}_tree.pdf']
	subprocess.run(' '.join(command),
	                        shell=True,
	                        check=True)
	with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("mltree done")

def cleanup(outPut):
	command =['mv', 
			f'{input_folder}/*.gz',
			'completed_fastq/']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	command = ['rm',
			f'{outPut_Folder}/*.bam']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	
	command = ['cat',
			f'{outPut}/seperateCall/{outPut}MergedCallAll.vcf',
			'>>',
			'totalMergedCall.vcf']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	command = ['bgzip',
			f'{outPut}/seperateCall/{outPut}MergedCallAll.vcf']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	command = ['mv',
			f'{outPut}/seperateCall/{outPut}MergedCallAll.vcf.gz',
			'processed_vcf/']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	command = ['rm',
			f'{outPut}/*.*',]
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	command = ['rm',
			'-r',
			f'{outPut}/seperateCall/',
			f'{outPut}/Coverage/']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
	command = ['cat',
			f'{outPut}/built_fasta/{outPut}builtSeqMeta.fasta',
			'>>',
			'totalFasta.mfa']
	subprocess.run(' '.join(command),
				shell=True,
				check=True)
#	os.remove(f'{outPut_Folder}/log.txt')
		
#if gzipped:
#	 fileList = glob.glob('fastq/*.gz')
#elif gzipped == False:
#	 fileListTemp = glob.glob('fastq/*.fastq')
#	 fileList2 = glob.glob('fastq/*.fq')
#	 for file in fileList2:
#		 fileListTemp.append(file)
#	 fileList =[]
#	 for file in fileList:
#		 command = ['bgzip',
#					f'fastq/{file}']
#		 subprocess.run(' '.join(command),
#						shell=True,
#						check=True)
#		 fileList.append(file + '.gz')

#this makes it so you can use the -i and -il commands at the same time
if included_isolates == None:
	included_isolates = []
if included_isolates_file != None:
	with open(included_isolates_file, 'r') as read:
		for line in read:
			if line.strip() not in included_isolates:
				included_isolates.append(line.strip())
#this makes it so you can use the -hf and -hfl commands at the same time
if included_hosts == None:
	included_hosts = []
if included_hosts_file != None:
	with open(included_hosts_file, 'r') as read:
		for line in read:
			if line.strip() not in included_hosts:
				included_hosts.append(line.strip())

threads = multiprocessing.cpu_count()
filtered = False
fileList = glob.glob(f'{input_folder}/*.gz')

try:
	with open(f'{outPut_Folder}/log.txt','r') as read:
		step = read.readlines()
		print(step)
		match step[0]:
			case "start":
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_bowtie2(outPut_Folder, fileNum, threads,1)
					auto_tabix(outPut_Folder, fileNum,'hits.bam')
					auto_mpileup(outPut_Folder, fileNum, threads)
					auto_call(outPut_Folder, fileNum)
					auto_bedtools(outPut_Folder, fileNum)
					auto_bgzip(outPut_Folder, fileNum, ".vcf")
					auto_tabix(outPut_Folder, fileNum, ".vcf.gz")
					autoMerge(outPut_Folder, file, fileNum)
					
				sampleBuilder(outPut_Folder)
				#this starts the filtering process if more then seq id is given
				if len(included_isolates) >= 1:
					fasta_filter(outPut_Folder, included_isolates)
					filtered = True
					
				if len(included_hosts) >= 1:
					fasta_filter_hosts(outPut_Folder, included_hosts,filtered)
					filtered = True
						
				filtered = raxmlGate(outPut_Folder,filtered)


				autoRAxML(outPut_Folder,filtered)

				mlTree(outPut_Folder)

				cleanup(outPut_Folder)
				
			case "autohisat2 done":
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_tabix(outPut_Folder, fileNum,'hits.bam',1)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("tabix1 done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_mpileup(outPut_Folder, fileNum, threads)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("mpileup done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_call(outPut_Folder, fileNum)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("call done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_bedtools(outPut_Folder, fileNum)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("bedtools done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_bgzip(outPut_Folder, fileNum, ".vcf")
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("bgzip done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_tabix(outPut_Folder, fileNum, ".vcf.gz",2)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("tabix2 done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					autoMerge(outPut_Folder, file, fileNum)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("merge done")
				
				sampleBuilder(outPut_Folder)
				#this starts the filtering process if more then seq id is given
				if len(included_isolates) >= 1:
					fasta_filter(outPut_Folder, included_isolates)
					filtered = True
					
				if len(included_hosts) >= 1:
					fasta_filter_hosts(outPut_Folder, included_hosts,filtered)
					filtered = True
						
				filtered = raxmlGate(outPut_Folder,filtered)


				autoRAxML(outPut_Folder,filtered)

				mlTree(outPut_Folder)

				cleanup(outPut_Folder)
			
			case "tabix1 done":
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_mpileup(outPut_Folder, fileNum, threads)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("mpileup done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_call(outPut_Folder, fileNum)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("call done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_bedtools(outPut_Folder, fileNum)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("bedtools done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_bgzip(outPut_Folder, fileNum, ".vcf")
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("bgzip done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_tabix(outPut_Folder, fileNum, ".vcf.gz",2)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("tabix2 done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					autoMerge(outPut_Folder, file, fileNum)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("merge done")
				
				sampleBuilder(outPut_Folder)
				#this starts the filtering process if more then seq id is given
				if len(included_isolates) >= 1:
					fasta_filter(outPut_Folder, included_isolates)
					filtered = True
					
				if len(included_hosts) >= 1:
					fasta_filter_hosts(outPut_Folder, included_hosts,filtered)
					filtered = True
						
				filtered = raxmlGate(outPut_Folder,filtered)


				autoRAxML(outPut_Folder,filtered)

				mlTree(outPut_Folder)

				cleanup(outPut_Folder)
				
			case "mpileup done":
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_call(outPut_Folder, fileNum)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("call done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_bedtools(outPut_Folder, fileNum)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("bedtools done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_bgzip(outPut_Folder, fileNum, ".vcf")
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("bgzip done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_tabix(outPut_Folder, fileNum, ".vcf.gz",2)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("tabix2 done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					autoMerge(outPut_Folder, file, fileNum)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("merge done")
				
				sampleBuilder(outPut_Folder)
				#this starts the filtering process if more then seq id is given
				if len(included_isolates) >= 1:
					fasta_filter(outPut_Folder, included_isolates)
					filtered = True
					
				if len(included_hosts) >= 1:
					fasta_filter_hosts(outPut_Folder, included_hosts,filtered)
					filtered = True
						
				filtered = raxmlGate(outPut_Folder,filtered)


				autoRAxML(outPut_Folder,filtered)

				mlTree(outPut_Folder)

				cleanup(outPut_Folder)
				
			case "call done":
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_bedtools(outPut_Folder, fileNum)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("bedtools done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_bgzip(outPut_Folder, fileNum, ".vcf")
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("bgzip done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_tabix(outPut_Folder, fileNum, ".vcf.gz",2)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("tabix2 done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					autoMerge(outPut_Folder, file, fileNum)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("merge done")
				
				sampleBuilder(outPut_Folder)
				#this starts the filtering process if more then seq id is given
				if len(included_isolates) >= 1:
					fasta_filter(outPut_Folder, included_isolates)
					filtered = True
					
				if len(included_hosts) >= 1:
					fasta_filter_hosts(outPut_Folder, included_hosts,filtered)
					filtered = True
						
				filtered = raxmlGate(outPut_Folder,filtered)


				autoRAxML(outPut_Folder,filtered)

				mlTree(outPut_Folder)

				cleanup(outPut_Folder)
				
			case "bedtools done":
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_bgzip(outPut_Folder, fileNum, ".vcf")
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("bgzip done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_tabix(outPut_Folder, fileNum, ".vcf.gz",2)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("tabix2 done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					autoMerge(outPut_Folder, file, fileNum)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("merge done")
				
				sampleBuilder(outPut_Folder)
				#this starts the filtering process if more then seq id is given
				if len(included_isolates) >= 1:
					fasta_filter(outPut_Folder, included_isolates)
					filtered = True
					
				if len(included_hosts) >= 1:
					fasta_filter_hosts(outPut_Folder, included_hosts,filtered)
					filtered = True
						
				filtered = raxmlGate(outPut_Folder,filtered)


				autoRAxML(outPut_Folder,filtered)

				mlTree(outPut_Folder)

				cleanup(outPut_Folder)
				
			case "bgzip done":
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					auto_tabix(outPut_Folder, fileNum, ".vcf.gz",2)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("tabix2 done")
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					autoMerge(outPut_Folder, file, fileNum)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("merge done")
				
				sampleBuilder(outPut_Folder)
				#this starts the filtering process if more then seq id is given
				if len(included_isolates) >= 1:
					fasta_filter(outPut_Folder, included_isolates)
					filtered = True
					
				if len(included_hosts) >= 1:
					fasta_filter_hosts(outPut_Folder, included_hosts,filtered)
					filtered = True
						
				filtered = raxmlGate(outPut_Folder,filtered)


				autoRAxML(outPut_Folder,filtered)

				mlTree(outPut_Folder)

				cleanup(outPut_Folder)
				
			case "tabix2 done":
				for file in fileList:
					fileNum = file.split('/')[1].split('.')[0]
					autoMerge(outPut_Folder, file, fileNum)
					with open(f'{outPut_Folder}/log.txt','w') as create:
						create.write("merge done")
				
				sampleBuilder(outPut_Folder)
				#this starts the filtering process if more then seq id is given
				if len(included_isolates) >= 1:
					fasta_filter(outPut_Folder, included_isolates)
					filtered = True
					
				if len(included_hosts) >= 1:
					fasta_filter_hosts(outPut_Folder, included_hosts,filtered)
					filtered = True
						
				filtered = raxmlGate(outPut_Folder,filtered)


				autoRAxML(outPut_Folder,filtered)

				mlTree(outPut_Folder)

				cleanup(outPut_Folder)
				
			case "merge done":
				sampleBuilder(outPut_Folder)
				#this starts the filtering process if more then seq id is given
				if len(included_isolates) >= 1:
					fasta_filter(outPut_Folder, included_isolates)
					filtered = True
					
				if len(included_hosts) >= 1:
					fasta_filter_hosts(outPut_Folder, included_hosts,filtered)
					filtered = True
						
				filtered = raxmlGate(outPut_Folder,filtered)


				autoRAxML(outPut_Folder,filtered)

				mlTree(outPut_Folder)

				cleanup(outPut_Folder)
				
			case "samplebuilder done":
				if len(included_isolates) >= 1:
					fasta_filter(outPut_Folder, included_isolates)
					filtered = True
					
				if len(included_hosts) >= 1:
					fasta_filter_hosts(outPut_Folder, included_hosts,filtered)
					filtered = True
						
				filtered = raxmlGate(outPut_Folder,filtered)


				autoRAxML(outPut_Folder,filtered)

				mlTree(outPut_Folder)

				cleanup(outPut_Folder)
				
			case "filtering done":
				autoRAxML(outPut_Folder,filtered)

				mlTree(outPut_Folder)

				cleanup(outPut_Folder)
				
			case "autoraxml done":
				mlTree(outPut_Folder)

				cleanup(outPut_Folder)
				
			case "mltree done":
				cleanup(outPut_Folder)

except:
	print("no log file")
	os.makedirs(f'{outPut_Folder}/seperateCall/')
	os.makedirs(f"{outPut_Folder}/Coverage/")
	with open(f'{outPut_Folder}/log.txt','a') as create:
		create.write("start")
	for file in fileList:
		fileNum = file.split('/')[1].split('.')[0]
		auto_bowtie2(outPut_Folder, fileNum, threads)
		with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("autohisat2 done")

	for file in fileList:
		fileNum = file.split('/')[1].split('.')[0]
		auto_tabix(outPut_Folder, fileNum,'hits.bam',1)
		with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("tabix1 done")
	for file in fileList:
		fileNum = file.split('/')[1].split('.')[0]
		auto_mpileup(outPut_Folder, fileNum, threads)
		with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("mpileup done")

	for file in fileList:
		fileNum = file.split('/')[1].split('.')[0]
		auto_call(outPut_Folder, fileNum)
		with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("call done")
	for file in fileList:
		fileNum = file.split('/')[1].split('.')[0]
		auto_bedtools(outPut_Folder, fileNum)
		with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("bedtools done")
	for file in fileList:
		fileNum = file.split('/')[1].split('.')[0]
		auto_bgzip(outPut_Folder, fileNum, ".vcf")
		with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("bgzip done")
	for file in fileList:
		fileNum = file.split('/')[1].split('.')[0]
		auto_tabix(outPut_Folder, fileNum, ".vcf.gz",2)
		with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("tabix2 done")
	for file in fileList:
		fileNum = file.split('/')[1].split('.')[0]
		autoMerge(outPut_Folder, file, fileNum)
		with open(f'{outPut_Folder}/log.txt','w') as create:
			create.write("merge done")

	sampleBuilder(outPut_Folder)
	#this starts the filtering process if more then seq id is given
	print(filtered)
	print(included_isolates)
	print(included_hosts)
	if len(included_isolates) >= 1:
		fasta_filter(outPut_Folder, included_isolates)
		filtered = True
		
	if len(included_hosts) >= 1:
		fasta_filter_hosts(outPut_Folder, included_hosts,filtered)
		filtered = True
			
	filtered = raxmlGate(outPut_Folder,filtered)

	autoRAxML(outPut_Folder,filtered)

	mlTree(outPut_Folder)

	cleanup(outPut_Folder)
