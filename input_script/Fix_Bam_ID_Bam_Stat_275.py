#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
import glob
from os.path import basename

def usage():
    test="name"
    message='''
python Fix_Bam_ID_Bam_Stat.py --input Bam_fixID

Generate coverage table for all the bam from Bam_fixID
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#GN-1    103     92 
def readtrait(infile):
    data = defaultdict(str)
    r = re.compile(r'GN-(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'GN'): 
                unit = re.split(r'\t',line)
                #print unit[0]
                m = r.search(unit[0])
                ril = m.groups(0)[0] if m else 'NA'
                #print ril
                data[ril] = 1
    return data



def bam_list(bam_files, r):
    data = defaultdict(list)
    for f in bam_files:
        m = r.search(f)
        ril = m.groups(0)[0] if m else 0
        #print f, ril
        if not ril == 0:
            data[ril].append(f)
    return data

'''
14658218 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
11271499 + 0 mapped (76.90%:-nan%)
14658218 + 0 paired in sequencing
7329109 + 0 read1
7329109 + 0 read2
11120120 + 0 properly paired (75.86%:-nan%)
11229836 + 0 with itself and mate mapped
41663 + 0 singletons (0.28%:-nan%)
72911 + 0 with mate mapped to a different chr
46198 + 0 with mate mapped to a different chr (mapQ>=5)
'''
def parsestat(infile):
    data = defaultdict(lambda : int())
    count = 0
    r_total = re.compile(r'total')
    r_mapped= re.compile(r'mapped')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                count += 1
                unit = re.split(r' ', line)
                #print unit[0]
                if r_total.search(line):
                #if unit[4] == 'total':
                    data['Total'] = int(unit[0])
                #elif unit[3] == 'mapped':
                elif r_mapped.search(line) and len(unit) < 7:
                    data['Mapped'] = int(unit[0])
    return data
 

def bamstat(list1, prefix):
    r = re.compile(r'(.*)\.bam')
    r1 = re.compile(r'GN(\d+).*')
    ofile = open('%s.bam.stat' %(prefix), 'w')
    print >> ofile, 'Sample\t#Read\tAverage\tTotal\tDepth\tMapped_Depth\tMapped_rate\t#Library\tFileName'
    for riln in sorted(list1.keys(), key=int):
        bams = list1[riln]
        n = len(bams)
        for bam in bams:
            m = r.search(bam)
            if not m:
                continue
            prefix = m.groups(0)[0]
            #flagstat = '%s.flagstat' %(prefix)
            flagstat_dd = '%s.dedup.flagstat' %(prefix)
            #stat = parsestat(flagstat) if os.path.isfile(flagstat) else [0, 0]
            stat_dd = parsestat(flagstat_dd) if os.path.isfile(flagstat_dd) else [0, 0]
            name = os.path.basename(bam)
            m1 = r1.search(name)
            #print prefix, name
            if m1:
                ril = 'RIL%s' %(str(m1.groups(0)[0]))
                readn = stat_dd['Total']
                readl = 250 if ril == 'RIL242' else 100
                basen = readn*readl
                depth = float(basen)/372000000
                mapr  = float(stat_dd['Mapped'])/float(stat_dd['Total'])
                depth1 = depth*mapr
                print >> ofile, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (ril, str(readn), str(readl), str(basen), str(depth), str(depth1), str(mapr), str(n), bam)
    ofile.close()
     
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
    
    #../input/fastq/Bam/RIL168_0_ATCACG_FC0813L1.recal.bam
    bam_all = glob.glob('%s/*.bam' %(args.input))

    r1 = re.compile(r'GN(\d+)') 
    list1 = bam_list(bam_all, r1)
 
    bamstat(list1, args.input)

if __name__ == '__main__':
    main()

