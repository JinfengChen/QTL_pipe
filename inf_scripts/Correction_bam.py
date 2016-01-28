#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
import glob
import subprocess

def usage():
    test="name"
    message='''
python MultiLib_bam.py --input ../../input/fastq/RILs_ALL_bam

Find bam files in Bam which is alternative library for other bams in RILs_ALL_bam.
We will put these in RILs_ALL_bam_multi_lib

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#GN-1    103     92 
def readtrait(infile):
    data = defaultdict(str)
    r = re.compile(r'GN-(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'GN'): 
                unit = re.split(r'\t',line)
                #print unit[0]
                m = r.search(unit[0])
                ril = m.groups(0)[0] if m else 'NA'
                #print ril
                data[ril] = 1
    return data



def bam_list(bam_files, r):
    data = defaultdict(int)
    for f in bam_files:
        m = r.search(f)
        ril = m.groups(0)[0] if m else 0
        #print f, ril
        if not ril == 0:
            data[ril] = 1
    return data
    
#return dict of ril->lib_name->bam_path
def parse_bam_all(bam_list, r):
    data = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int())))
    for lib in sorted(re.split(r'\n', bam_list)):
        unit = re.split(r' |\t', lib)
        bam = os.path.split(unit[-1])[1]
        bam = re.sub(r'.recal.bam', r'', bam)
        bam = re.sub(r'.bam', r'', bam)
        #print lib, bam
        if r.search(bam):
            ril = r.search(bam).groups(0)[0]
            data[ril][bam] = unit[-1]
            #print ril
    return data   

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    ##../input/fastq/Bam/RIL168_0_ATCACG_FC0813L1.recal.bam 
    #bam_all = subprocess.check_output('ls -all ../../input/fastq/Bam/*.recal.bam', shell=True)
    #bam_275 = subprocess.check_output('ls -all %s/*.bam' %(args.input), shell=True)
    #bam_all = glob.glob('../../input/fastq/Bam/*.recal.bam')
    ##../input/fastq/RILs_ALL_bam/GN1.bam
    #bam_in  = glob.glob('%s/*.bam' %(args.input))
    #rils    = readtrait('../../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.ALL') 

    r1 = re.compile(r'RIL(\d+)\_')
    r2 = re.compile(r'GN(\d+)\.')

    target = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_correction'
    source = '/rhome/cjinfeng/BigData/00.RD/RILs/Problem_RILs/bin/RILs_genotype/genotypes/MSU_r7.corrected'
    #/rhome/cjinfeng/BigData/00.RD/RILs/Problem_RILs/bin/RILs_genotype/genotypes/MSU_r7.corrected/RIL217_0_CTTGTA_FC251L4.recal.bam
    bams_corr = glob.glob('%s/*.bam' %(source))
    for bam in bams_corr:
        ril = re.split(r'_', os.path.split(bam)[1])[0]
        ril = re.sub(r'RIL', r'GN', ril)
        bai = re.sub(r'.bam', r'.bai', bam)
        flagstat = re.sub(r'.recal.bam', r'.dedup.flagstat', bam)
        bam_new = '%s/%s.bam' %(target, ril)
        bai_new = '%s/%s.bam.bai' %(target, ril)
        flagstat_new = '%s/%s.flagstat' %(target, ril)
        #print bam
        #print bam_new
        os.system('ln -s %s %s' %(bam, bam_new))
        os.system('ln -s %s %s' %(bai, bai_new))
        os.system('ln -s %s %s' %(flagstat, flagstat_new))

if __name__ == '__main__':
    main()

