#importing packages
import os
import glob
import subprocess 
import argparse

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

#command line option for RAxML versions
parser.add_argument(
    '-r',
    action='store',
    help=(
        'Name of the RAxML version used'
    ),
    default='raxmlHPC'
)

args = parser.parse_args()

outPut_Folder = args.o
metadata_file_name = args.m
RAXML_version = args.r

#autoVCF function 
def autoVCF(outPut, fileNum):
    os.makedirs(f'{outPut}/seperateCall/')
    os.makedirs(f"{outPut}/Coverage/")
    print(fileNum, " is entering the pipeline")
    #histat 2 + samtools sort call
    command = ['hisat2',
               '-p',
               '2',
               '-x',
               'index/70-15index',
               '-U',
               file,
               '--max-intronlen',
               '20',
               '--summary-file',
               f'{outPut}/{fileNum}summary.txt',
               '--dta-cufflinks',
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
    command = ['tabix',
               f'{outPut}/{fileNum}hits.bam']
    subprocess.run(' '.join(command),
                   shell=True,
                   check=True)
    #mpileup call
    command = ['bcftools',
               'mpileup',
               '--threads',
               '2',
               '-d',
               '100000',
               '-R',
               'MonsterPlexRegionsFileSuperCont.txt',
               '--annotate',
               'FORMAT/AD',
               '-f',
               'index/70-15.fasta.fasta',
               f'{outPut}/{fileNum}hits.bam',
               '>>',
               f'{outPut}/{fileNum}.vcf']
    subprocess.run(' '.join(command),
                   shell=True,
                   check=True)
    #bcftools call command
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
    #bedtools genome coverage command
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
    #zip up the results for use later in the pipeline
    command =['bgzip',
              f'{outPut}/{fileNum}.vcf']
    subprocess.run(' '.join(command),
                   shell=True,
                   check=True)
    #tabix the gzipped results
    command = ['tabix',
               f'{outPut}/{fileNum}.vcf.gz']
    subprocess.run(' '.join(command),
                   shell=True,
                   check=True)
        
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
        append.write(file.split('/')[1].split('.')[0] + 'call.vcf.gz\n')
        
def sampleBuilder(outPut):
    command = ['bcftools',
               'merge',
               '-l',
               f'{outPut}/fastqListCall.txt',
               '-o',
               f'{outPut}/seperateCall/wheatBlastMergedCallAll.vcf']
    subprocess.run(' '.join(command),
                   shell=True,
                   check=True)
    
    sites =[]
    sitesUsed =[]
    #reads a list of sites you want and only looks at data from there
    with open('MonsterPlexSitesList.txt', 'r') as read:
        for line in read:
            sites.append(line.strip('\n'))
    
    with open(f'{outPut}/seperateCall/wheatBlastMergedCallAll.vcf', 'r') as read:
        seqs = list()
        check = False
        for line in read:
            if line.split('\t')[0] == '#CHROM':
                print("header past seqs made")
                fqList = line.split('\t')
                for n in range(9,len(fqList)):
                    seqs.append([fqList[n], ''])
                    check = True
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
    
    os.mkdir(f'{outPut}/built_fasta')
    
    with open(f'{outPut}/built_fasta/builtSeqMeta.fasta', 'a') as writeSeq:
        for read in seqs:
            seqID = read[0].split('/')[1].split('.')[0]
            if seqID in sample_metadata:
                seqSpecies = sample_metadata[seqID][0]
                seqHost = sample_metadata[seqID][1]
                seqCountry = sample_metadata[seqID][2]
                writeSeq.write(f'>{seqID}_{seqSpecies}_{seqHost}_{seqCountry}\n{read[1]}\n')
            else:
                writeSeq.write('>' + read[0].split('/')[1].split('.')[0] + '\n' + read[1] + '\n')
                        
def metaDataBuilder(metadata_file):
    metaData = {}
    with open(metadata_file, 'r') as read:
        for line in read:
            ID = line.split('\t')[0].strip('\n')
            species = line.split('\t')[1].strip('\n')
            host = line.split('\t')[2].strip('\n')
            country = line.split('\t')[3].strip('\n')
            metaData[ID] = [species, host, country]
            return metaData
        
def autoRAxML(outPut,version):
    os.mkdir(f'{outPut}/RAXML_results')
    #command for running RAXML
    command = [f'./{version}',
               '-p',
               '1234',
               '-f',
               'a',
               '-x',
               '1234',
               '-s',
               f'{outPut}/built_fasta/builtSeqMeta.fasta',
               '-n',
               f'{outPut}/RAXML_results/pythreads.raxml',
               '-m',
               'GTRGAMMA',
               '-#',
               '1000']
    subprocess.run(' '.join(command),
                   shell=True,
                   check=True)
               
               
        

fileList = glob.glob('fastq/*.gz')
for file in fileList:
    fileNum = file.split('/')[1].split('.')[0]
    autoVCF(outPut_Folder, fileNum)
    autoMerge(outPut_Folder, file, fileNum)
sampleBuilder(outPut_Folder)
autoRAxML(outPut_Folder,RAXML_version)
    
