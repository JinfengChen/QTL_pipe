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
python BamDir.py --input ../inf_list/Bam.Core.list --project RILs_ALL_bam_core

Read rils in Bam.Core.list, link bam and SNP file from RILs_ALL_bam to a new dir "project" without these problem rils.
After this change the dir name in "step01.genotype.sh" and rerun the QTL pipeline:
qsub step01.parent.sh
qsub -q js step01.genotype.sh
qsub -q js step02.recombination_bin.sh
qsub -q js step03.QTL.sh

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#GN104   heterozygous on chr07 21Mb
def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = 1
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if args.project is None:
        args.project = 'RILs_ALL_bam_core'

    newdir = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/%s' %(args.project)
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    count = 0    
    rils = readtable(args.input)
    bams = glob.glob('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/*.bam')
    ofile = open('../BWA.sampleRIL.list', 'w')
    for bam in sorted(bams):
        prefix, ext = os.path.splitext(bam)
        ril = os.path.basename(prefix)
        if rils.has_key(ril):
            continue
        snp = '%s.Maq.p1.map.pileup.SNP' %(prefix)
        cmd1 = 'ln -s %s %s' %(bam, newdir)
        cmd2 = 'ln -s %s %s' %(snp, newdir)
        cmd3 = 'cp %s %s/Exclued.list' %(args.input, newdir)
        print >> ofile, ril
        #print cmd1
        #print cmd2
        os.system(cmd1)
        os.system(cmd2)
        os.system(cmd3)
        count += 1
    print 'Job done: Linked %s' %(count) 
 
if __name__ == '__main__':
    main()

