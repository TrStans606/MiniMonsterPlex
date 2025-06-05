# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 13:49:59 2023

@author: treyd
"""
#importing packages
import os
import glob
import re
import subprocess 
import multiprocessing
import shutil
import gzip

def auto_bowtie2(outPut_Folder, input_file, fileNum,threads):
	try:
		print(fileNum, " is entering the pipeline")
		#histat 2 + samtools sort call
		command = ['bowtie2 --no-unal',
			'-p',
			str(threads),
			'-x',
			'index/70-15_small_index',
			'-U',
			input_file,
			'--local --very-sensitive-local',
			'2>', f'{outPut_Folder}/bowtie_out/{fileNum}_alignment_summary.txt'
			'|',
			'samtools', 
			'sort',
			'-',
			'-@',
			'2',
			'-O',
			'bam',
			'-o',
			f'{outPut_Folder}/bowtie_out/{fileNum}hits.bam']
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(f'{outPut_Folder}/bowtie_out/{fileNum}hits.bam'):
			os.remove(f'{outPut_Folder}/bowtie_out/{fileNum}hits.bam')
		if os.path.isfile(f'{outPut_Folder}/bowtie_out/{fileNum}_alignment_summary.txt'):
			os.remove(f'{outPut_Folder}/bowtie_out/{fileNum}_alignment_summary.txt')
		print(f'Issue with file: {fileNum}')
		print(f'Errror running command: {e}')
		print(f"Something Went wrong with bowtie2 :{e.stderr}")
		quit()

def parse_alignment_summary(fileNum, outPut_Folder):
	with open(fileNum, 'r') as f:
		content = f.read()

	# Extract total reads
	total_match = re.search(r'(\d[\d,]*) reads; of these:', content)
	if not total_match:
		raise ValueError(f"Total reads not found in {fileNum}")
	total_reads = int(total_match.group(1).replace(',', ''))

	# Extract aligned 1 time
	aligned_once_match = re.search(r'(\d[\d,]*) \([\d.]+%\) aligned exactly 1 time', content)
	if not aligned_once_match:
		raise ValueError(f"Aligned exactly 1 time not found in {fileNum}")
	aligned_once = int(aligned_once_match.group(1).replace(',', ''))

	# Extract aligned >1 times
	aligned_multiple_match = re.search(r'(\d[\d,]*) \([\d.]+%\) aligned >1 times', content)
	if not aligned_multiple_match:
		raise ValueError(f"Aligned >1 times not found in {fileNum}")
	aligned_multiple = int(aligned_multiple_match.group(1).replace(',', ''))

	# Compute total aligned
	aligned_total = aligned_once + aligned_multiple
	aligned_fraction = aligned_total / total_reads

	# Return CSV string
	return f"{fileNum},{total_reads},{aligned_total},{aligned_fraction:.4f}"


def auto_mpileup(outPut, fileNum, threads):
	try:
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
				os.path.join(outPut,'bowtie_out',f'{fileNum}hits.bam'),
				'>>',
				f'{outPut}/mpileup_out/{fileNum}.vcf']
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(f'{outPut}/mpileup_out/{fileNum}.vcf'):
			os.remove(f'{outPut}/mpileup_out/{fileNum}.vcf')
		print(f'Issue with file: {fileNum}')
		print(f'Errror running command: {e}')
		print(f"Something Went wrong with mpileup :{e.stderr}")
		quit()
	
def auto_call(outPut, fileNum):
	try:
		command = ['bcftools',
				'call',
				'-c',
				'--ploidy',
				'1',
				os.path.join(outPut, 'mpileup_out', f'{fileNum}.vcf'),
				'-o',
				f'{outPut}/call_out/{fileNum}call.vcf']
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(f'{outPut}/call_out/{fileNum}call.vcf'):
			os.remove(f'{outPut}/call_out/{fileNum}call.vcf')
		print(f'Issue with file: {fileNum}')
		print(f'Errror running command: {e}')
		print(f"Something Went wrong with bctools call:{e.stderr}")
		quit()

	
def auto_bedtools(outPut, fileNum):
	try:
		command = ['bedtools',
				'genomecov',
				'-ibam',
				os.path.join(outPut, 'bowtie_out', f'{fileNum}hits.bam'),
				'-bg',
				'>',
				f'{outPut}/coverage_out/{fileNum}cover.bed']
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(f'{outPut}/coverage_out/{fileNum}cover.bed'):
			os.remove(f'{outPut}/coverage_out/{fileNum}cover.bed')
		print(f'Issue with file: {fileNum}')
		print(f'Errror running command: {e}')
		print(f"Something Went wrong with bedtools genomecov:{e.stderr}")
		quit()

def auto_bgzip(outPut, fileNum,file_ex):
	try:
		command =['bgzip',
				os.path.join(outPut,'mpileup_out',f'{fileNum}{file_ex}')]
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(os.path.join(outPut,'mpileup_out',f'{fileNum}{file_ex}.gz')):
			os.remove(os.path.join(outPut,'mpileup_out',f'{fileNum}{file_ex}.gz'))
		print(f'Issue with file: {fileNum}')
		print(f'Errror running command: {e}')
		print(f"Something Went wrong with bgzip:{e.stderr}")
		quit()

	try:
		command = ['tabix',
				os.path.join(outPut,'mpileup_out',f'{fileNum}{file_ex}.gz')]
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(os.path.join(outPut,'mpileup_out',f'{fileNum}{file_ex}.gz.tbi')):
			os.remove(os.path.join(outPut,'mpileup_out',f'{fileNum}{file_ex}.gz.tbi'))
		print(f'Issue with file: {fileNum}')
		print(f'Errror running command: {e}')
		print(f"Something Went wrong with tabix:{e.stderr}")
		quit()		

def autoVCFZip(outPut, file, fileNum):
	#bg zip the bcftools call result file
	try:
		command = ['bgzip',
				os.path.join(outPut,'call_out',f'{fileNum}call.vcf')]
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(os.path.join(outPut,'call_out',f'{fileNum}call.vcf.gz')):
			os.remove(os.path.join(outPut,'call_out',f'{fileNum}call.vcf.gz'))
		print(f'Issue with file: {fileNum}')
		print(f'Errror running command: {e}')
		print(f"Something Went wrong with bgzip:{e.stderr}")
		quit()

	#tabix the call results
	try:
		command = ['tabix',
				os.path.join(outPut,'call_out',f'{fileNum}call.vcf.gz')]
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(os.path.join(outPut,'call_out',f'{fileNum}call.vcf.gz.tbi')):
			os.remove(os.path.join(outPut,'call_out',f'{fileNum}call.vcf.gz.tbi'))
		print(f'Issue with file: {fileNum}')
		print(f'Errror running command: {e}')
		print(f"Something Went wrong with tabix:{e.stderr}")
		quit()		
	with open(f'{outPut}/fastqListCall.txt', 'a') as append:
		append.write(f'{outPut}/call_out/' + file.split('/')[1].split('.')[0] + 'call.vcf.gz\n')

def autoMerge(outPut):
	try:
		command = ['bcftools',
				'merge',
				'-l',
				f'{outPut}/fastqListCall.txt',
				'-o',
				f'{outPut}/merge_out/{outPut}MergedCallAll.vcf']
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(os.path.join(outPut,'merge_out',f'{outPut}MergedCallAll.vcf')):
			os.remove(os.path.join(outPut,'merge_out',f'{outPut}MergedCallAll.vcf'))
		print(f'Errror running command: {e}')
		print(f"Something Went wrong with bcftools merge:{e.stderr}")
		quit()		

		
def sampleBuilder(outPut,metadata_file):
	sites =[]
	sitesUsed =[]
	#reads a list of sites you want and only looks at data from there
	with open('MonsterPlexSitesList_small_index.txt', 'r') as read:
		for line in read:
			sites.append(line.strip('\n'))

	with open(f'{outPut}/merge_out/{outPut}MergedCallAll.vcf', 'r') as read:
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
		sample_metadata = metaDataBuilder(metadata_file)
		
		print(sample_metadata)
		
		os.mkdir(f'{outPut}/built_fasta')
		
		with open(f'{outPut}/built_fasta/{outPut}builtSeqMeta.fasta', 'a') as writeSeq:
			for read in seqs:
				seqID = read[0].split('/')[2].split('.')[0].split('hits')[0]
				if len(seqID.split("_")) > 1:
					seqID = f'{"-".join(seqID.split("_"))}'
				if (seqID) in sample_metadata:
					seqSpecies = sample_metadata[seqID][0]
					seqHost = sample_metadata[seqID][1]
					seqLineage = sample_metadata[seqID][2]
					seqCountry = sample_metadata[seqID][3]
					writeSeq.write(f'>{seqID}_{seqSpecies}_{seqHost}_{seqLineage}_{seqCountry}\n{read[1]}\n')
				else:
					writeSeq.write('>' + seqID
								+ '_._._._.' + '\n' + read[1] + '\n')
						
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


def autoRAxML(outPut,filtered):
	wd = os.getcwd()
	if filtered==False:
	#command for running RAXML
		try:
			command = ['raxmlHPC',
				'-w',
				os.path.join(wd,outPut),
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
					check=True,
					capture_output=True)
		except subprocess.CalledProcessError as e:
			raxmls = glob.glob(os.path.join(outPut,'*.raxml'))
			for file in raxmls:
				os.remove(file)
			print(f'Errror running command: {e}')
			print(f"Something Went wrong with raxml:{e.stderr}")
			quit()		
		os.mkdir(os.path.join(outPut,'raxml_out'))
		raxmls = glob.glob(os.path.join(outPut,'*.raxml'))
		for file in raxmls:
			shutil.move(file,os.path.join(outPut,'raxml_out'))
	else:
		try:
			command = ['raxmlHPC',
				'-w',
				os.path.join(wd,outPut),
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
					check=True,
					capture_output=True)
		except subprocess.CalledProcessError as e:
			raxmls = glob.glob(os.path.join(outPut,'*.raxml'))
			for file in raxmls:
				os.remove(file)
			print(f'Errror running command: {e}')
			print(f"Something Went wrong with raxml:{e.stderr}")
			quit()		
		os.mkdir(os.path.join(outPut,'raxml_out'))
		raxmls = glob.glob(os.path.join(outPut,'*.raxml'))
		for file in raxmls:
			shutil.move(file,os.path.join(outPut,'raxml_out'))
		

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
	try:
		command = ['Rscript',
			'--vanilla',
			'MLtree.R',
			f'{outPut_Folder}/raxml_out/RAxML_bestTree.miniMonsterPlex.raxml']
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile('NA.pdf'):
			os.remove('NA.pdf')
		print(f'Errror running command: {e}')
		print(f"Something Went wrong with mlTree:{e.stderr}")
		quit()		
	os.mkdir(os.path.join(outPut_Folder,'tree_out'))
	shutil.move('NA.pdf',f'{outPut_Folder}/tree_out/{outPut_Folder}_tree.pdf')

#series of lines for cleaing up left over temp data
def cleanup(outPut,input_folder,complete,project_name):
	files_to_move = glob.glob(os.path.join(input_folder,'*.gz'))
	for file in files_to_move:
		shutil.move(file,os.path.join(project_name,'completed_fastq/'))

	with open('totalMergedCall.vcf', 'a') as f:
		with open(f'{outPut}/merge_out/{outPut}MergedCallAll.vcf','r') as read:
			for line in read:
				f.write(line)

	with open(f'{outPut}/merge_out/{outPut}MergedCallAll.vcf', 'rb') as f_in:
		with gzip.open(f'{outPut}/merge_out/{outPut}MergedCallAll.vcf.gz', 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
	shutil.move(f'{outPut}/merge_out/{outPut}MergedCallAll.vcf.gz',os.path.join(project_name,'processed_vcf/'))

	shutil.rmtree(f'{outPut}/bowtie_out/')
	shutil.rmtree(f'{outPut}/coverage_out/')
	shutil.rmtree(f'{outPut}/call_out/')
	shutil.rmtree(f'{outPut}/mpileup_out/')
	shutil.rmtree(f'{outPut}/merge_out/')
	os.remove(f'{outPut}/fastqListCall.txt')

	
	with open('totalFasta.mfa','a') as f:
		with open(f'{outPut}/built_fasta/{outPut}builtSeqMeta.fasta','r') as read:
			for line in read:
				f.write(line)
	#does the combined anyalysis of all the data
	if complete:
		command = ['raxmlHPC',
			 '-p',
			 '1234',
			 '-f',
			 'a',
			 '-x',
			 '1234',
			 '-s',
			 'totalFasta.mfa',
			 '-n',
			 'miniMonsterPlex_full.raxml',
			 '-m',
			 'GTRGAMMA',
			 '-#',
			 '1000']
		subprocess.run(' '.join(command),
				 shell=True,
				 check=True)
		command = ['Rscript',
		   '--vanilla',
		   'MLtree.R',
		   'RAxML_bestTree.miniMonsterPlex_full.raxml']
		subprocess.run(' '.join(command),
				shell=True,
				check=True)
		files_to_delete = glob.glob(os.path.join(os.getcwd(),'*.raxml'))
		for file in files_to_delete:
			os.remove(file)
		files_to_delete = glob.glob(os.path.join(os.getcwd(),'*.raxml.support'))
		for file in files_to_delete:
			os.remove(file)

def main(project, metadata_file, complete=False):
	outPut_Folder = os.path.join(project,"output")
	metadata_file_name = os.path.join(project,"metadata",metadata_file)
	input_folder = os.path.join(project,"newFastq")
	included_isolates = None
	included_isolates_file = None
	included_hosts = None
	included_hosts_file = None
	threads = multiprocessing.cpu_count()
	if threads > 8:
		threads =8
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

	# Print a tabular summary of key variables (project, metadata_file, outPut_Folder, inPut_Folder, complete)
	print(f"{'Variable':<20} {'Value'}")
	print(f"{'-'*20} {'-'*40}")
	print(f"{'Project':<20} {project}")
	print(f"{'Metadata File':<20} {metadata_file_name}")
	print(f"{'Output Folder':<20} {outPut_Folder}")
	print(f"{'Input Folder':<20} {input_folder}")
	print(f"{'Complete mode':<20} {complete}")

	if os.path.isdir(outPut_Folder):
		pass
	else:
		os.mkdir(outPut_Folder)
		os.mkdir(os.path.join(project,'completed_fastq'))
		os.mkdir(os.path.join(project,'processed_vcf/'))
		os.mkdir(os.path.join(outPut_Folder,'bowtie_out'))
		os.mkdir(os.path.join(outPut_Folder,'mpileup_out'))
		os.mkdir(os.path.join(outPut_Folder, "processedAlignSumm"))
		os.mkdir(os.path.join(outPut_Folder,'call_out'))
		os.mkdir(os.path.join(outPut_Folder,'coverage_out'))
		os.mkdir(os.path.join(outPut_Folder, 'merge_out'))

	#everything now runs on a file by file basis skipping steps if an output file already exists
	for file in fileList:
		fileNum = file.split('/')[1].split('.')[0]
		if os.path.isfile(os.path.join(outPut_Folder,'bowtie_out',f'{fileNum}hits.bam')):
			print(f'File {fileNum} has already been processed by bowtie')
		else:
			auto_bowtie2(outPut_Folder, file,fileNum, threads)

		summary_string = parse_alignment_summary(f'{outPut_Folder}/bowtie_out/{fileNum}_alignment_summary.txt', outPut_Folder)
		# Write the result to a CSV file
		with open(f'{outPut_Folder}/processedAlignSumm/alignment_summary.csv', "a") as out_file:  # "a" = append mode
			out_file.write(summary_string + "\n")

		#capture more error info about subprocess commands
		if os.path.isfile(os.path.join(outPut_Folder,'bowtie_out',f'{fileNum}hits.bam.csi')):
			print(f"{fileNum} has already been tabixed")
		else:
			try:
				command = ['tabix',
				os.path.join(outPut_Folder,'bowtie_out',f'{fileNum}hits.bam')]
				subprocess.run(' '.join(command),
						shell=True,
						check=True,
						capture_output=True)
			except subprocess.CalledProcessError as e:
				if os.path.isfile(os.path.join(outPut_Folder,'bowtie_out',f'{fileNum}hits.bam.csi')):
					os.remove(os.path.join(outPut_Folder,'bowtie_out',f'{fileNum}hits.bam.csi'))
				print(f'Issue with file: {fileNum}')
				print(f'Errror running command: {e}')
				print(f"Something Went wrong with tabix :{e.stderr}")
				quit()

		if os.path.isfile(os.path.join(outPut_Folder,'mpileup_out',f'{fileNum}.vcf')):
			print(f'File {fileNum} has already been processed by mpileup')
		else:
			auto_mpileup(outPut_Folder, fileNum, threads)
		
		if os.path.isfile(os.path.join(outPut_Folder,'call_out',f'{fileNum}call.vcf')):
			print(f'File {fileNum} has already been processed by bcftools call')
		else:
			auto_call(outPut_Folder, fileNum)

		if os.path.isfile(f'{outPut_Folder}/coverage_out/{fileNum}cover.bed'):
			print(f'File {fileNum} has already been processed by bedtools genomcov')
		else:
			auto_bedtools(outPut_Folder, fileNum)

		file_extension = ".vcf"

		if os.path.isfile(os.path.join(outPut_Folder,'mpileup_out',f'{fileNum}{file_extension}.gz.tbi')):
			print(f'File {fileNum} has already been processed by bgzip')
		else:
			auto_bgzip(outPut_Folder, fileNum, file_extension)

		if os.path.isfile(os.path.join(outPut_Folder,'call_out',f'{fileNum}call.vcf.gz.tbi')):
			print(f'File {fileNum} has already been processed by bgzip and tabix')
		else:	
			autoVCFZip(outPut_Folder, file, fileNum)

	if os.path.isfile(os.path.join(outPut_Folder,'merge_out',f'{outPut_Folder}MergedCallAll.vcf')):
		print('Merge already done: Skipping')
		pass
	else:
		autoMerge(outPut_Folder)

	if os.path.isdir(os.path.join(outPut_Folder,'built_fasta')):
		print('Fasta already built: Skipping')
		pass
	else:
		sampleBuilder(outPut_Folder,metadata_file_name)
	#this starts the filtering process if more then seq id is given
	if len(included_isolates) >= 1:
		fasta_filter(outPut_Folder, included_isolates)
		filtered = True

	if len(included_hosts) >= 1:
		fasta_filter_hosts(outPut_Folder, included_hosts,filtered)
		filtered = True
		
	filtered = raxmlGate(outPut_Folder,filtered)

	if os.path.isdir(os.path.join(outPut_Folder,'raxml_out')):
		print("Raxml has already run")
	else:
		autoRAxML(outPut_Folder,filtered)

	mlTree(outPut_Folder)

	cleanup(outPut_Folder,input_folder,complete,project)
