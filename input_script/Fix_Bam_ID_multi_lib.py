#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from collections import OrderedDict
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
python Fix_Bam_ID_multi_lib.py
This script deal with issue 2. First, we produce a detailed table of all the libraries for all the RILs. Second, for these RILs have two libraries, we keep newest one if they are different. we keep both if they are same. Third, For these have three or more libraries we create table and correction by manual inspection. 

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

def flowcell_date():
    flowcell = {
    "FC133":'120710',  #date not lable, assigned to oldest
    "FC153":'120810',
    "FC193":'130522',
    "FC197":'130624',
    "FC205":'130708',
    "FC251":'140624',
    "FC271":'141202',
    "FC279":'141217',
    "FC0813":'130901', #date not lable, assigned to between FC205 and FC251 which roughly right
    "FC1213":'130902'  #date not lable, assigned to between FC205 and FC251 which roughly right
    }
    return flowcell

#Bam_fixID/RIL100_0_ATCACG_FC251L2.recal.bam 
def multi_lib(work_dir):
    data = defaultdict(lambda : defaultdict(lambda : str()))
    bams = glob.glob('%s/*.recal.bam' %(work_dir))
    #create dict of RIL->lib->flowcell
    for bam in sorted(bams):
        lib = re.split(r'\.', os.path.split(bam)[1])[0]
        ril = re.sub(r'RIL', r'' ,re.split(r'_', lib)[0])
        flowcell = re.split(r'_', lib)[-1][:-2]
        data[ril][lib] = flowcell
        #print lib, ril, flowcell 
    #rank the lib for each RIL by sequenced date of flowcell
    fc_date = flowcell_date()
    print 'RIL\tLib:Date'
    for ril in sorted(data.keys(), key=int):
        if len(data[ril].keys()) == 1:
            lib = data[ril].keys()[0]
            print 'RIL%s\t%s:%s' %(ril, lib, fc_date[data[ril][lib]])
        elif len(data[ril].keys()) > 1:
            lib_dict = defaultdict(lambda : str())
            for lib in data[ril].keys():
                lib_date = fc_date[data[ril][lib]]
                lib_dict[lib] = lib_date
            #sort and output
            lib_dict_sorted = OrderedDict(sorted(lib_dict.items(), key=lambda x: x[1], reverse=True))
            ril_info = []
            for lib in lib_dict_sorted:
                lib_info = '%s:%s' %(lib, lib_dict_sorted[lib])
                ril_info.append(lib_info)
            print 'RIL%s\t%s' %(ril, '\t'.join(ril_info))
       
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    
    work_dir    = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/Bam_fixID'    
    fc251_fixed = '/rhome/cjinfeng/BigData/00.RD/RILs/Problem_RILs/bin/RILs_genotype/genotypes/MSU_r7.corrected'
    #clean files not RILs
    #clean_landrace(work_dir)
    #fix flowcell 251, use new genotype data from fc251_fixed to replace the old one
    #fix_fc251(work_dir, fc251_fixed) 

    #For these libraries with multi libraries we collection information (flowcell, depth, data) and picked the newest one as representive libraries.
    multi_lib(work_dir)

if __name__ == '__main__':
    main()

