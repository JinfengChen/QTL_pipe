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

def usage():
    test="name"
    message='''
python Run_fastqc.py --fastq ../../input/fastq/Illumina_correct

fastq quality control of fq using fastqc
    '''
    print message

def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=16,walltime=2:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-fq', '--fastq')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.fastq) > 0
    except:
        usage()
        sys.exit(2)


    if not args.project:
        args.project = 'fastqc'

    adapter = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/bin/inf_scripts/adapter.list'
    sh_script = '%s.sh' %(args.project)
    ofile = open (sh_script, 'w')
    fastqs = glob.glob('%s/*/*.fq' %(os.path.abspath(args.fastq)))
    for fq in sorted(fastqs):
        #print fq
        fq_dir = os.path.split(fq)[0]
        out    = '%s_fastqc.zip' %(fq)
        if not os.path.exists(out):
            cmd = '/opt/linux/centos/7.x/x86_64/pkgs/fastqc/0.11.3/fastqc -a %s -o %s -f fastq %s' %(adapter, fq_dir, fq)
            print >> ofile, cmd
    ofile.close()
    
    #runjob(sh_script, 1)

if __name__ == '__main__':
    main()

