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
import subprocess

def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#0100031071A     GN278   G
def read_snp(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                pos  = int(unit[0][2:-1])
                data[pos] = unit[2]
    return data


#return dict of lib -> SNP
#Maq.p1.map.pileup.SNP
def parse_bam_all(bam_list):
    data = defaultdict(lambda : str())
    for lib in sorted(re.split(r'\n', bam_list)):
        #print lib
        unit = re.split(r' ', lib)
        if len(unit) < 2:
            continue
        #print len(unit)
        #print '%s\t%s' %(unit[-3], unit[-1])
        bam = os.path.split(unit[-1])[1]
        bam = re.sub(r'.recal.bam', r'', bam)
        bam = re.sub(r'.bam', r'', bam)
        snp = unit[-3]
        #print snp
        snp = re.sub(r'\.bam', r'.Maq.p1.map.pileup.SNP', snp)
        #print snp, bam
        data[bam] = snp
    return data  

def snp_similarity(snp1, snp2):
    snp1_dict = read_snp(snp1)
    snp2_dict = read_snp(snp2)
    total = 0
    match = 0
    if len(snp1_dict.keys()) < 20000 or len(snp2_dict.keys()) < 20000:
        return 'NA'
    for pos in snp1_dict.keys():
        if snp2_dict.has_key(pos):
            total += 1
            if snp1_dict[pos] == snp2_dict[pos]:
                match += 1
    sim = float(match)/float(total)
    return sim

def function_helper(args):
    return snp_similarity(*args)

##run function with parameters using multiprocess of #cpu
def mp_pool_function(function, parameters, cpu):
    pool = mp.Pool(int(cpu))
    imap_it = pool.map(function, tuple(parameters))
    collect_list = []
    for x in imap_it:
        #print 'status: %s' %(x)
        collect_list.append(x)
    return collect_list


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
 
    #/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_multi_lib
    bam_all   = subprocess.check_output('ls -all %s/*.bam' %('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_multi_lib'), shell=True)
    bam_dupli = subprocess.check_output('ls -all %s/*.bam' %('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_multi_lib'), shell=True)

    snp_dupli = parse_bam_all(bam_dupli) 
    snp_all   = parse_bam_all(bam_all)

    print 'Lib1\tLib2\tSimilarity' 
    for lib_d in sorted(snp_dupli.keys()):
        snp_d = snp_dupli[lib_d]
        for lib_a in sorted(snp_all.keys()):
            snp_a   = snp_all[lib_a]
            print snp_d, snp_a
            snp_sim = snp_similarity(snp_d, snp_a)
            print '%s\t%s\t%s' %(lib_d, lib_a, snp_sim)


if __name__ == '__main__':
    main()

