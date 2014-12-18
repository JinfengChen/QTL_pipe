#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python poorSNP.py --input 277
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#V1      V2
#0100021547A     A       T
#0100031071A     A       G
def readparent(infile, snp):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'V'): 
                unit = re.split(r'\t',line)
                if not snp.has_key(unit[0]):
                    print line
    return data

#0100031071A     GN80    G
def readsnp(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = 1
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
    snpfile = '../../input/fastq/RILs_ALL_bam/GN%s.Maq.p1.map.pileup.SNP' %(args.input)
    snp = readsnp(snpfile)
    readparent('../NB.RILs.dbSNP.SNPs.parents', snp)

if __name__ == '__main__':
    main()

