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
python PrepareRelocaTE_Merged_BAM.py --bam RILs_ALL_bam_correct

Prepare bam for each RIL.
For RILs that have more than one library with same barcode, we merged all these with same barcode.
For RILs that have only one barcode for selected high coverage library, we link the directory from Bam_correct.
    '''
    print message

def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource walltime=100:00:00,nodes=1:ppn=1,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

#RIL150_0_GTAGAG_FC1213L3.recal.bam
#ril -> barcode -> abspath of recal.bam
def get_bam_barcode(bam_dir):
    data = defaultdict(lambda : defaultdict(lambda : list()))
    r = re.compile(r'(((RIL\d+)_\d+)_(\w+)_\w+)\.recal\.bam')
    for file_in_dir in sorted(os.listdir(bam_dir)):
        if fnmatch.fnmatch(file_in_dir, '*.recal.bam'):
            if r.search(file_in_dir):
                barcode= r.search(file_in_dir).groups(0)[3]
                ril    = r.search(file_in_dir).groups(0)[2]
                strain = r.search(file_in_dir).groups(0)[1]
                prefix = r.search(file_in_dir).groups(0)[0]
                data[ril][barcode].append('%s/%s' %(bam_dir, file_in_dir))
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

    #Raw bam directory
    if not args.fastq:
        #args.fastq = '/shared/wesslerlab/Rice/RIL/Illumina/'
        #args.fastq = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/Illumina_fixed_link'
        args.fastq = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/Bam_correct'

    if not args.output:
        args.output = '%s_merged' %(os.path.abspath(args.bam))
        
        
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    bam_barcode = get_bam_barcode(args.fastq)

    #/shared/wesslerlab/Rice/RIL/Illumina/RIL39_1/RIL39_1_ACTTGA_FC133L7_p1.fq
    #/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/Bam/RIL39_1_ACTTGA_FC133L7.recal.bam
    r = re.compile(r'(((RIL\d+)_\d+)_(\w+)_\w+)\.recal\.bam')
    bams = glob.glob('%s/*.bam' %(args.bam))
    ofile = open('%s.bam.list' %(args.output), 'w')
    merge_shell = open('%s.merge_bam.sh' %(args.output), 'w')
    for bam in sorted(bams):
        bam_path = os.path.realpath(bam)
        bam_abspath = os.path.abspath(bam)
        if r.search(os.path.split(bam_path)[1]):
            barcode= r.search(os.path.split(bam_path)[1]).groups(0)[3]
            ril    = r.search(os.path.split(bam_path)[1]).groups(0)[2]
            strain = r.search(os.path.split(bam_path)[1]).groups(0)[1]
            prefix = r.search(os.path.split(bam_path)[1]).groups(0)[0]
            
            print ril, strain, prefix, barcode
            if bam_barcode.has_key(ril):
                if bam_barcode[ril].has_key(barcode):
                    print >> ofile, '%s\t%s\t%s' %(ril, len(bam_barcode[ril][barcode]), '\t'.join(sorted(bam_barcode[ril][barcode])))
                    #subdir = '%s/%s' %(args.output, ril)
                    #if not os.path.exists(subdir):
                    #    os.mkdir(subdir)
                    input_bam_bai  = '%s.bai' %(bam_abspath)
                    input_bam_tab  = re.sub(r'\.bam', r'.genotype.tab.gz', bam_abspath)
                    input_bam_vcf  = re.sub(r'\.bam', r'.genotype.vcf.gz', bam_abspath)
                    input_bam_stat = re.sub(r'\.bam', r'_stats', bam_abspath)
                    input_bam_flag = re.sub(r'\.bam', r'.dedup.flagstats', bam_abspath)
                    input_bam_SNP  = re.sub(r'\.bam', r'.Maq.p1.map.pileup.SNP', bam_abspath)
                    output_bam_tab = '%s/%s' %(args.output, os.path.split(bam_abspath)[1])
                    output_bam_tab = re.sub(r'\.bam', r'.genotype.tab.gz', output_bam_tab)
                    print input_bam_tab, input_bam_vcf, input_bam_SNP
                    if len(bam_barcode[ril][barcode]) > 1:
                        #merge bam files
                        input_bam_files = []
                        output_bam_file = '%s/%s' %(args.output, os.path.split(bam_abspath)[1])
                        output_bam_flagstat = re.sub(r'\.bam', r'.dedup.flagstats', output_bam_file)
                        print output_bam_file, output_bam_flagstat
                        for bam_file in sorted(bam_barcode[ril][barcode]):
                            input_bam_files.append('INPUT=%s' %(bam_file))
                        if not os.path.exists(output_bam_file):
                        #if 1:
                            print >> merge_shell, '/opt/linux/centos/7.x/x86_64/pkgs/java/jdk1.8.0_25/bin/java -Xmx10g -jar /opt/picard/1.81/MergeSamFiles.jar VALIDATION_STRINGENCY=LENIENT %s USE_THREADING=true OUTPUT=%s' %(' '.join(input_bam_files), output_bam_file)
                            print >> merge_shell, '/usr/bin/samtools index %s' %(output_bam_file)
                            print >> merge_shell, '/usr/bin/samtools flagstat %s > %s' %(output_bam_file, output_bam_flagstat)
                        if not os.path.exists(output_bam_tab):
                            print 'link for 2 libraries'
                            #print input_bam_tab, input_bam_vcf, input_bam_SNP
                            os.system('ln -s %s %s/' %(input_bam_tab, args.output))
                            os.system('ln -s %s %s/' %(input_bam_vcf, args.output))
                            os.system('ln -s %s %s/' %(input_bam_SNP, args.output))
                    else:
                        #link bam files in
                        print 'single library'
                        if not os.path.exists(output_bam_tab):
                            print 'link for 1 libraries'
                            #print input_bam_tab, input_bam_vcf, input_bam_SNP
                            os.system('ln -s %s %s/' %(bam_abspath, args.output))
                            os.system('ln -s %s %s/' %(input_bam_bai, args.output))
                            os.system('ln -s %s %s/' %(input_bam_tab, args.output))
                            os.system('ln -s %s %s/' %(input_bam_vcf, args.output))
                            os.system('ln -s %s %s/' %(input_bam_stat, args.output))
                            os.system('ln -s %s %s/' %(input_bam_flag, args.output))
                            os.system('ln -s %s %s/' %(input_bam_SNP, args.output)) 
                else:
                    print >> ofile, 'no bam found for this barcode'
            else:
                print >> ofile, 'no bam found for this ril'
                
    ofile.close()
    merge_shell.close()
    #runjob('%s.merge_bam.sh' %(args.output), 3)

if __name__ == '__main__':
    main()

