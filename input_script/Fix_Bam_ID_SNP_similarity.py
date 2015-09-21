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
python Fix_Bam_ID_SNP_similarity.py --input Bam_fixID --cpu 32 > Bam_fixID.SNP.similarity
or
qsub -q highmem SNP_similarity.sh

#same ril
awk '$1~/RIL58_/ && $2~/RIL58_/' RILs_SNP.dupli_all.similarity
#all dupli
awk '$3>0.9' RILs_SNP.dupli_all.similarity |grep "NA" -v| less -S
awk '$3>0.9' RILs_SNP.dupli_all_number.similarity | grep "NA" -v > RILs_SNP.dupli_all_number.duplicate.table
#all
awk '$3>0.9' RILs_SNP.similarity | grep "NA" -v | awk '$1!=$2' | less -S
awk '$3>0.9' RILs_SNP.all_all_number.similarity | grep "NA" -v | awk '$1!=$2' > RILs_SNP.all_all_number.duplicate.table
#dupli vs dupli, no duplicate
awk '$3>0.9' RILs_SNP.dupli_dupli_number.similarity |grep "NA" -v | awk '$1!=$2' | less -S

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
                #pos  = int(unit[0][2:-1])
                pos   = unit[0]
                data[pos] = unit[2]
    return data

##CHROM  POS     REF     RIL103_0_GAGTGG_FC1213L5
#Chr1    31071   A       A/A
def read_snp_tab(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'): 
                unit = re.split(r'\t',line)
                #pos  = int(unit[0][2:-1])
                chrs  = re.sub(r'Chr', r'', unit[0]) 
                pos   = '%02d%08d%s' %(int(chrs), int(unit[1]), unit[2])
                if unit[3][0] == unit[3][2] and not unit[3][0] == '.':
                    data[pos] = unit[3][0]
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

def snp_similarity(lib1, lib2, snp1, snp2):
    #print lib1, lib2
    #print snp1, snp2
    snp1_dict = read_snp_tab(snp1)
    snp2_dict = read_snp_tab(snp2)
    total = 0
    match = 0
    if len(snp1_dict.keys()) < 2000 or len(snp2_dict.keys()) < 2000:
        return [lib1, lib2, 'NA', 0, 0, len(snp1_dict.keys()), len(snp2_dict.keys())]
    for pos in snp1_dict.keys():
        if snp2_dict.has_key(pos):
            total += 1
            if snp1_dict[pos] == snp2_dict[pos]:
                match += 1
    sim = float(match)/float(total)
    return [lib1, lib2, sim, total, match, len(snp1_dict.keys()), len(snp2_dict.keys())]

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


#make dict of lib->snp_file
def lib_snp_file(snp_files):
    data = defaultdict(lambda : str())
    for f in sorted(snp_files):
        lib = re.split(r'\.', os.path.split(f)[1])[0]
        data[lib] = os.path.abspath(f)
        #print lib, f
    return data

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
        exit(2)

    #/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_multi_lib
    #bam_all   = subprocess.check_output('ls -all %s/*.bam' %('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_multi_lib'), shell=True)
    #bam_dupli = subprocess.check_output('ls -all %s/*.bam' %('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_correction'), shell=True)

    #snp_dupli = parse_bam_all(bam_dupli) 
    #snp_all   = parse_bam_all(bam_all)

    #Bam_fixID/RIL100_0_ATCACG_FC251L2.genotype.tab
    snp_files = glob.glob('%s/*.genotype.tab' %(args.input))
    snp_dict  = lib_snp_file(snp_files)

    print 'Lib1\tLib2\tSimilarity\tTotal_Shared_SNP_Site\tTotal_Identical_SNP_Sites\tLib1_SNP\tLib2_SNP'
    parameters = []
    for lib_1 in sorted(snp_dict.keys()):
        snp_1 = snp_dict[lib_1]
        for lib_2 in sorted(snp_dict.keys()):
            snp_2   = snp_dict[lib_2]
            if not snp_1 == snp_2:
                #print snp_1, snp_2
                parameters.append([lib_1, lib_2, snp_1, snp_2])
            else:
                #print snp_1, snp_2
                pass
    collect_list_list = mp_pool_function(function_helper, parameters, args.cpu)
    for result in collect_list_list:
        print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(result[0], result[1], result[2], result[3], result[4], result[5], result[6])


if __name__ == '__main__':
    main()

