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
python statcore.py > RIL.bam.unique.Core.stat

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#GN115
def readcore(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'GN'): 
                unit = re.split(r'\t',line)
                data[unit[0]] = 0
                #print unit[0]
    return data


def main():
    corelist = readcore('Bam.Core.list')
    
    with open ('RIL.bam.unique.stat', 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith(r'GN'): 
                unit = re.split(r'\t',line)
                unit[0] = re.sub(r'\_\?', '', unit[0])
                #print unit[0]
                if not corelist.has_key(unit[0]):
                    print line
            else:
                print line

if __name__ == '__main__':
    main()

