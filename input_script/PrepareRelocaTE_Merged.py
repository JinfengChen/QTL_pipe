#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import glob
import fnmatch

def usage():
    test="name"
    message='''
python PrepareRelocaTE_Merged.py --bam RILs_ALL_bam_correct

Prepare fastq directory for each RIL. We can run RelocaTE or other analysis from there.
For RILs that have more than one library with same barcode, we add all these with same barcode in the directory.
For RILs that have only one barcode for selected high coverage library, we link the directory from RILs_ALL_fastq_correct.
    '''
    print message


#ril -> barcode -> abspath of fq.gz
def get_fastq_barcode(fastq_dir):
    data = defaultdict(lambda : defaultdict(lambda : list()))
    r = re.compile(r'(((RIL\d+)_\d+)_(\w+)_\w+)_p1\.fq\.gz')
    for ril_fq_dir in sorted(os.listdir(fastq_dir)):
        if fnmatch.fnmatch(ril_fq_dir, 'RIL*_*'):
            for fq in sorted(os.listdir('%s/%s' %(fastq_dir, ril_fq_dir))):
                if r.search(fq):
                    barcode= r.search(fq).groups(0)[3]
                    ril    = r.search(fq).groups(0)[2]
                    strain = r.search(fq).groups(0)[1]
                    prefix = r.search(fq).groups(0)[0]
                    data[ril][barcode].append('%s/%s/%s' %(fastq_dir, ril_fq_dir, fq))
    return data

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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam')
    parser.add_argument('-f', '--fastq')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.bam) > 0
    except:
        usage()
        sys.exit(2)

    if not args.fastq:
        #args.fastq = '/shared/wesslerlab/Rice/RIL/Illumina/'
        #args.fastq = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/Illumina_fixed_link'
        args.fastq = '/rhome/cjinfeng/Rice/RIL/Illumina_correct'

    if not args.output:
        args.output = '%s_merged' %(re.sub(r'bam', 'fastq', args.bam))
        
        
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    fastq_barcode = get_fastq_barcode(args.fastq)

    #/shared/wesslerlab/Rice/RIL/Illumina/RIL39_1/RIL39_1_ACTTGA_FC133L7_p1.fq
    #/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/Bam/RIL39_1_ACTTGA_FC133L7.recal.bam
    r = re.compile(r'(((RIL\d+)_\d+)_(\w+)_\w+)\.recal\.bam')
    bams = glob.glob('%s/*.bam' %(args.bam))
    ofile = open('%s.fastq.list' %(args.output), 'w')
    for bam in sorted(bams):
        bam_path = os.path.realpath(bam)
        if r.search(os.path.split(bam_path)[1]):
            barcode= r.search(os.path.split(bam_path)[1]).groups(0)[3]
            ril    = r.search(os.path.split(bam_path)[1]).groups(0)[2]
            strain = r.search(os.path.split(bam_path)[1]).groups(0)[1]
            prefix = r.search(os.path.split(bam_path)[1]).groups(0)[0]
            #print ril, strain, prefix, barcode
            if fastq_barcode.has_key(ril):
                if fastq_barcode[ril].has_key(barcode):
                    print >> ofile, '%s\t%s\t%s' %(ril, len(fastq_barcode[ril][barcode]), '\t'.join(sorted(fastq_barcode[ril][barcode])))
                    subdir = '%s/%s' %(args.output, ril)
                    if not os.path.exists(subdir):
                        os.mkdir(subdir)
                    for p1 in sorted(fastq_barcode[ril][barcode]):
                        #print p1
                        #continue
                        p2 = re.sub(r'p1.fq.gz', r'p2.fq.gz', p1)
                        os.system('ln -s %s %s/' %(p1, subdir))
                        os.system('ln -s %s %s/' %(p2, subdir))
                else:
                    print >> ofile, 'no fastq found for this barcode'
            else:
                print >> ofile, 'no fastq found for this ril'
                
            #fq1    = '%s/%s/%s_p1.fq.gz' %(args.fastq, strain, prefix)
            #fq2    = '%s/%s/%s_p2.fq.gz' %(args.fastq, strain, prefix)
            #subdir = '%s/%s' %(args.output, ril)
            #fq1n   = '%s/%s_1.fq.gz' %(subdir, ril)
            #fq2n   = '%s/%s_2.fq.gz' %(subdir, ril)
            #if not os.path.exists(subdir):
            #    os.mkdir(subdir)
            #os.system('ln -s %s %s' %(fq1, fq1n)) 
            #os.system('ln -s %s %s' %(fq2, fq2n))
    ofile.close()

if __name__ == '__main__':
    main()

