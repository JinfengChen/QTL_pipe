#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python Needreseq.py --input 1

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

'''
GN1     GN10    GN100   GN101   GN102   GN103   GN104   GN105   GN106   GN107   GN108   GN109   GN11    GN110
0100021547A     NA      NA      NA      A       A       A       NA      NA      NA      NA      A       A
0100031071A     NA      G       A       A       A       A       G       G       G       G       NA      A
'''
def snpmatrix(infile):
    data = defaultdict(int)
    data1 = defaultdict(lambda: float)
    rils = []
    total = 0
    s = re.compile(r'GN(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            m = s.search(line)
            if line.startswith('GN') and m: 
                unit = re.split(r'\t',line)
                rils = [a.replace('GN','') for a in unit]
            elif len(line) > 2:
                total = total + 1
                unit = re.split(r'\t',line)
                unit = unit[1:]
                count = 0
                for snp in unit: 
                    count = count + 1
                    if snp == 'NA':
                        data[rils[count-1]] = int(data[rils[count-1]]) + 1
    
    for ril in data.keys():
        narate = float(data[ril])/float(total)
        data1[ril] = narate
        #print ril, narate
    return data1



'''
Sample  Heading Days    Plant Height (cm) in Field      Biomass Number of Tillers       Single Plant (Grain Yield) (g)
GN-1    103     92      120.2   22      42.0
'''
def trait(infile):
    data = defaultdict(str)
    s = re.compile(r'GN-(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            m = s.search(line)
            if line.startswith('GN') and m: 
                unit = re.split(r'\t',line)
                ril = m.groups(0)[0]
                if not data.has_key(ril):
                    data[ril] = 1
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    if args.input is None:
        args.input = '1'

    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
   
    rils = trait('../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt')
    snps = snpmatrix('NB.RILs.dbSNP.SNPs.RILs')
    for ril in rils:
        if snps.has_key(ril):
            if snps[ril] > 0:
                print ril, snps[ril]
        else:
            print ril, 'No sequence'
if __name__ == '__main__':
    main()

