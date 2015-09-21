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
python Fix_ln_Illumina.py

Link the fastq files of RILs into "Illumina_fixed_link" on pigeon. The old one Illumina failed in all links.
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#lrwxrwxrwx 1 robb wesslerlab  98 Jul 30  2013 ./Illumina/A160_0_0/A160_0_ATCACG_FC197L7_p1.fq -> /shared/wesslerlab/Rice/RIL/FC197_RIL_75_140_179_234_259_A160/flowcell197_lane7_pair1_ATCACG.fastq
def fix_link(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r' ',line)
                fq_source = unit[-1]
                fq_target = unit[-3]
                #print '%s\t%s' %(fq_source, fq_target)
                fq_source = re.sub(r'/shared/wesslerlab/Rice', r'/bigdata/wesslerlab/shared/Rice', fq_source)
                fq_source = re.sub(r'/rhome/robb/Wessler-Rice', r'/bigdata/wesslerlab/shared/Rice', fq_source)
                fq_target = re.sub(r'Illumina', r'Illumina_fixed_link', fq_target)
                #print '%s\t%s' %(fq_source, fq_target)
                fq_targets_dir = os.path.dirname(fq_target)
                #print targets[-2], fq_targets_dir
                if not os.path.exists(fq_targets_dir):
                    os.mkdir(fq_targets_dir)
                os.system('ln -s %s %s' %(fq_source, fq_target))
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    outdir = 'Illumina_fixed_link'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    os.system('ls -all %s/*/*.fq > temp.fq.list' %('./Illumina'))
    fix_link('temp.fq.list')

if __name__ == '__main__':
    main()

