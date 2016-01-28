#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import glob

def usage():
    test="name"
    message='''
python PrepareRelocaTE.py --bam RILs_ALL_bam_correct

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam')
    parser.add_argument('-f', '--fastq')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.bam) > 0
    except:
        usage()
        sys.exit(2)

    if not args.fastq:
        #args.fastq = '/shared/wesslerlab/Rice/RIL/Illumina/'
        #args.fastq = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/Illumina_fixed_link'
        args.fastq = '/rhome/cjinfeng/Rice/RIL/Illumina_correct'

    if not args.output:
        args.output = re.sub(r'bam', 'fastq', args.bam)
        
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    #/shared/wesslerlab/Rice/RIL/Illumina/RIL39_1/RIL39_1_ACTTGA_FC133L7_p1.fq
    #/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/Bam/RIL39_1_ACTTGA_FC133L7.recal.bam
    r = re.compile(r'(((RIL\d+)_\d+)_\w+_\w+)\.recal\.bam')
    bams = glob.glob('%s/*.bam' %(args.bam))
    for bam in bams:
        bam_path = os.path.realpath(bam)
        if r.search(os.path.split(bam_path)[1]):
            ril    = r.search(os.path.split(bam_path)[1]).groups(0)[2]
            strain = r.search(os.path.split(bam_path)[1]).groups(0)[1]
            prefix = r.search(os.path.split(bam_path)[1]).groups(0)[0]
            print ril, strain, prefix
            fq1    = '%s/%s/%s_p1.fq.gz' %(args.fastq, strain, prefix)
            fq2    = '%s/%s/%s_p2.fq.gz' %(args.fastq, strain, prefix)
            subdir = '%s/%s' %(args.output, ril)
            fq1n   = '%s/%s_1.fq.gz' %(subdir, ril)
            fq2n   = '%s/%s_2.fq.gz' %(subdir, ril)
            if not os.path.exists(subdir):
                os.mkdir(subdir)
            os.system('ln -s %s %s' %(fq1, fq1n)) 
            os.system('ln -s %s %s' %(fq2, fq2n))

if __name__ == '__main__':
    main()

