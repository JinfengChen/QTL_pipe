#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
import glob

def usage():
    test="name"
    message='''
python NewBam.py --input ../../input/fastq/RILs_ALL_bam
Find new Bam files in ../input/Bam, which we can link into ../input/RILs_ALL_bam to update

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
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    #../input/fastq/RILs_ALL_bam/GN1.bam
    bam_all = glob.glob('../../input/fastq/Bam/*.recal.bam')
    #../input/fastq/Bam/RIL168_0_ATCACG_FC0813L1.recal.bam 
    bam_in  = glob.glob('%s/*.bam' %(args.input))
    rils    = readtrait('../../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.ALL') 

    r1 = re.compile(r'RIL(\d+)\_') 
    r2 = re.compile(r'GN(\d+)\.')
    list1 = bam_list(bam_all, r1)
    list2 = bam_list(bam_in, r2)
 
    print 'New Bam to link:'
    for bam in sorted(list1.keys(), key=int):
        if not list2.has_key(bam):
            print bam

    print 'New ril to sequence:'
    for ril in sorted(rils.keys()):
        if not list1.has_key(ril) and not list2.has_key(ril):
            print ril    

if __name__ == '__main__':
    main()

