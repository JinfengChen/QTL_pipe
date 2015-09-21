#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir

def usage():
    test="name"
    message='''
python Fix_Bam_ID.py
This script deal with issue 1.

Fix Bam ID by the following ways:
0. All the files from Sofia's genotype results were linked to Bam_FixID
1. For FC251, some index were assigned incorrect. We replace the whole flowcell by new process files.
2. For RILs have multi libraries, we know some are wrong. We correct these by using the right libraries.

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

#fix fc251 by replace with new one
def fix_fc251(work_dir, fc251_fixed):
    files1 = os.listdir(work_dir)
    r = re.compile(r'FC251')
    removed = defaultdict(lambda : str())
    added   = defaultdict(lambda : str())
    for f1 in sorted(files1):
        if r.search(f1):
            f1_fp = '%s/%s' %(work_dir, f1)
            unit1 = re.split(r'_', f1)
            ril1  = unit1[0]
            removed[ril1] = 1
            os.system('rm %s' %(f1_fp))
    files2 = os.listdir(fc251_fixed)
    for f2 in sorted(files2):
        if r.search(f2):
            f2_fp ='%s/%s' %(fc251_fixed, f2)
            unit2 = re.split(r'_', f2)
            ril2  = unit2[0]
            added[ril2] = 1
            os.system('rm %s/%s' %(work_dir, f2))
            os.system('ln -s %s %s/' %(f2_fp, work_dir))
    print 'Removed %s RILs: %s' %(len(removed.keys()), ','.join(sorted(removed.keys())))
    print 'Added %s RILs: %s' %(len(added.keys()), ','.join(sorted(added.keys())))

def clean_landrace(work_dir):
    list1 = glob.glob('%s/%s_*' %(work_dir, 'A160'))
    list2 = glob.glob('%s/%s_*' %(work_dir, 'RILA123'))
    for f in list1:
        print f
        os.system('rm %s' %(f))
    for f in list2:
        print f
        os.system('rm %s' %(f))

def fc_date():
    flowcell = {
    "FC133":'120710',
    "FC153":'120810',
    "FC193":'130522',
    "FC197":'130624',
    "FC205":'130708',
    "FC251":'140624',
    "FC271":'141202',
    "FC279":'141217',
    "FC0813":'130901',
    "FC1213":'130902'
    }
    return flowcell
 
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    
    work_dir    = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/Bam_fixID'    
    fc251_fixed = '/rhome/cjinfeng/BigData/00.RD/RILs/Problem_RILs/bin/RILs_genotype/genotypes/MSU_r7.corrected'
    #clean files not RILs
    clean_landrace(work_dir)
    #fix flowcell 251, use new genotype data from fc251_fixed to replace the old one
    fix_fc251(work_dir, fc251_fixed) 

    #For these libraries with multi libraries we collection information (flowcell, depth, data) and picked the newest one as representive libraries.

if __name__ == '__main__':
    main()

