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
import multiprocessing as mp


def usage():
    test="name"
    message='''
python SNP_similarity_pairs.py --input SNP_similarity.list
take a list of pairs of library to check, output similary and detailed list of SNPs that matched and differed.

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
                #pos  = int(unit[0][:-1])
                #ref  = unit[0][-1]
                pos   = unit[0]
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

#0100137002A     GN232   A
def snp_similarity_list(lib1, lib2, snp1, snp2, outfile):
    #print lib1, lib2
    #print snp1, snp2
    snp1_dict = read_snp(snp1)
    snp2_dict = read_snp(snp2)
    total = 0
    match = 0
    if len(snp1_dict.keys()) < 2000 or len(snp2_dict.keys()) < 2000:
        return [lib1, lib2, 'NA', 0, 0]
    ofile = open(outfile, 'w')
    print >> ofile, 'Pos\tSNP1\tGenotype1\tSNP2\tGenotype2'
    for pos in sorted(snp1_dict.keys()):
        if snp2_dict.has_key(pos):
            total += 1
            if snp1_dict[pos] == snp2_dict[pos]:
                match += 1
            #write details
            ref = pos[-1]
            genotype1 = 0 if ref == snp1_dict[pos] else 1
            genotype2 = 0 if ref == snp2_dict[pos] else 1 
            print >> ofile, '%s\t%s\t%s\t%s\t%s' %(pos, snp1_dict[pos], genotype1, snp2_dict[pos], genotype2)
    ofile.close()
    sim = float(match)/float(total)
    return [lib1, lib2, sim, total, match]

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
    parser.add_argument('-c', '--cpu')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
 
    if not args.cpu:
        args.cpu = 2
    if not args.input:
        usage()
        exit()

    pairs_dir = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_pairs_SNPs'
    #/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_multi_lib
    bam_all   = subprocess.check_output('ls -all %s/*.bam' %('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam'), shell=True)
    bam_dupli = subprocess.check_output('ls -all %s/*.bam' %('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_multi_lib'), shell=True)

    snp_dupli = parse_bam_all(bam_dupli) 
    snp_all   = parse_bam_all(bam_all)

    ##read pairs and compare SNPs
    print 'Lib1\tLib2\tSimilarity\tTotalSNP\tMatchedSNP'
    with open (args.input, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                lib1 = unit[0]
                lib2 = unit[1]
                outfile  = '%s/%s.%s.list' %(pairs_dir, lib1, lib2)
                snp1 = snp_dupli[lib1] if snp_dupli.has_key(lib1) else snp_all[lib1]
                snp2 = snp_dupli[lib2] if snp_dupli.has_key(lib2) else snp_all[lib2]
                snp_sim = snp_similarity_list(lib1, lib2, snp1, snp2, outfile)
                print '\t'.join(map(str, snp_sim)) 

'''
    print 'Lib1\tLib2\tSimilarity'
    parameters = []
    for lib_d in sorted(snp_dupli.keys()):
        snp_d = snp_dupli[lib_d]
        for lib_a in sorted(snp_all.keys()):
            snp_a   = snp_all[lib_a]
            #print snp_d, snp_a
            parameters.append([lib_d, lib_a, snp_d, snp_a])
            #snp_sim = snp_similarity(snp_d, snp_a)
            #print '%s\t%s\t%s' %(lib_d, lib_a, snp_sim)
    collect_list_list = mp_pool_function(function_helper, parameters, args.cpu)
    for result in collect_list_list:
        print '%s\t%s\t%s\t%s\t%s' %(result[0], result[1], result[2], result[3], result[4])
'''

if __name__ == '__main__':
    main()

