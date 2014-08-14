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
python NewBam.py --input ../input/RILs_ALL_bam
Find new Bam files in ../input/Bam, which we can link into ../input/RILs_ALL_bam to update

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


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
    bam_all = glob.glob('../input/fastq/Bam/*.recal.bam')
    #../input/fastq/Bam/RIL168_0_ATCACG_FC0813L1.recal.bam 
    bam_in  = glob.glob('%s/*.bam' %(args.input))
   
    r1 = re.compile(r'RIL(\d+)\_') 
    r2 = re.compile(r'GN(\d+)\.')
    list1 = bam_list(bam_all, r1)
    list2 = bam_list(bam_in, r2)
 
    for bam in sorted(list1.keys(), key=int):
        if not list2.has_key(bam):
            print bam
    
if __name__ == '__main__':
    main()

