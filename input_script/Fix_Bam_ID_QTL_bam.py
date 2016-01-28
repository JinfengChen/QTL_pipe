#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import glob
import subprocess

def usage():
    test="name"
    message='''
python Fix_Bam_ID_QTL_bam.py --input Bam_fixID

1. Remove incorrect libraries from Bam_fixID into Bam_fixID/incorrect_ID
2. Generate a RILs_ALL_bam_fixID directory, which similar to RILs_ALL_bam that used for QTL pipeline.
We will check if the bam from Bam_fixID and RILs_ALL_bam is the same libraries for the same RIL. If so we will use the same bam file and SNP file.
If not we just link the correct bam and generate SNP file latter.

After this change the dir name in "step01.genotype.sh" and rerun the QTL pipeline:
qsub step01.parent.sh
qsub -q js step01.genotype.sh
qsub -q js step02.recombination_bin.sh
qsub -q js step03.QTL.sh

    '''
    print message

def createdir(outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir) 


def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#RIL     Lib:Date:Depth
#RIL1    RIL1_0_CGTACG_FC153L5:120810:6.13                               NA
#RIL2    RIL2_0_TTAGGC_FC153L5:120810:5.13                               NA
#Use the first lib which is newest one RIls
#Can also read from other file that manual edited if need
def read_picked_lib(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                if not len(unit[0]) < 4:
                    lib = re.split(r':', unit[1])[0]
                    data[lib] = unit[0]
                   
    return data

def archive_incorrect_lib(fixID_dir):
    incorrect_id = [
    #incorrect, resequenced between 20150920
    'RIL58_0_CTTGTA_FC193L6:130522:12.41',
    'RIL60_0_CTTGTA_FC193L6:130522:12.41',
    'RIL168_0_ATCACG_FC0813L1:130901:11.00',
    'RIL187_0_TGACCA_FC0813L1:130901:15.18',
    'RIL208_0_CAGATC_FC0813L1:130901:12.90',
    'RIL206_0_TTAGGC_FC0813L1:130901:10.39',
    'RIL211_0_GCCAAT_FC0813L1:130901:11.79',
    'RIL216_0_ACAGTG_FC1213L7:130902:1.52',
    'RIL243_0_CTTGTA_FC1213L7:130902:1.32',
    'RIL253_0_ATCACG_FC1213L7:130902:1.72'
    ]
    #create dir of archive
    archive_incorrect_id = '%s/%s' %(os.path.abspath(fixID_dir), 'archive_incorrect_id')
    createdir(archive_incorrect_id)

    for lib in incorrect_id:
        lib = re.split(r':', lib)[0]
        lib_fp = '%s/%s' %(os.path.abspath(fixID_dir), lib)
        #print lib_fp
        if os.path.isfile('%s.recal.bam' %(lib_fp)):
            os.system('mv %s.* %s' %(lib_fp, archive_incorrect_id))


#return dict of ril->lib_name->bam_path
def parse_bam_all(bam_list, r):
    data = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int())))
    for lib in sorted(re.split(r'\n', bam_list)):
        unit = re.split(r' |\t', lib)
        bam = os.path.split(unit[-1])[1]
        bam = re.sub(r'.recal.bam', r'', bam)
        bam = re.sub(r'.bam', r'', bam)
        #print lib, bam
        if r.search(bam):
            ril = r.search(bam).groups(0)[0]
            data[ril][bam] = unit[-1]
            #print ril
    return data   

#return dict of ril->lib_name->bam_path
def parse_bam_all_fixID(bam_list, r, picked):
    data = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int())))
    for bam in sorted(bam_list):
        lib = os.path.split(bam)[1]
        lib = re.sub(r'.recal.bam', r'', lib)
        if picked.has_key(lib):
            #print lib, bam
            if r.search(lib):
                ril = r.search(lib).groups(0)[0]
                data[ril][lib] = os.path.abspath(bam)
                #print ril
    return data   

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if args.project is None:
        args.project = 'RILs_ALL_bam_fixID'

    #archive incorrect ID
    archive_incorrect_lib(args.input)
   
    #create new bam directory for QTL
    newdir = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/%s' %(args.project)
    createdir(newdir)
    r1 = re.compile(r'RIL(\d+)\_')
    r2 = re.compile(r'GN(\d+)')
    ##bam origial 275
    rils_all_bams  = subprocess.check_output('ls -all %s/*.bam' %('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam'), shell=True)
    rils_all_bams_lib = parse_bam_all(rils_all_bams, r1)
    
    ##bam fixID 275
    picked = read_picked_lib('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/Bam_fixID.info')
    rils_all_bams_fix_ID = glob.glob('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/Bam_fixID/*.recal.bam')
    rils_all_bams_fix_ID_lib = parse_bam_all_fixID(rils_all_bams_fix_ID, r1, picked)   
 
    #for ril in sorted(rils_all_bams_lib.keys()):
    #    for lib in rils_all_bams_lib[ril].keys():
    #        print '%s\t%s\t%s' %(ril, lib, rils_all_bams_lib[ril][lib])
    print "These are fixed or changed library: "
    for ril in sorted(rils_all_bams_fix_ID_lib.keys()):
        for lib in rils_all_bams_fix_ID_lib[ril].keys():
            #print '%s\t%s\t%s' %(ril, lib, rils_all_bams_fix_ID_lib[ril][lib])
            if rils_all_bams_lib[ril].has_key(lib):
                bam = rils_all_bams_lib[ril][lib]
                bai = re.sub(r'\.bam', r'.bai', bam)
                flag = re.sub(r'.recal.bam', r'.dedup.flagstat', bam)
                snp = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/GN%s.Maq.p1.map.pileup.SNP' %(ril)
                newbam = '%s/GN%s.bam' %(newdir, ril)
                newbai = re.sub(r'\.bam', r'.bai', newbam)
                newflag= re.sub(r'\.bam', r'.dedup.flagstat', newbam)
                cmd1 = 'ln -s %s %s' %(bam, newbam)
                cmd2 = 'ln -s %s %s' %(bai, newbai)
                cmd3 = 'ln -s %s %s' %(flag, newflag)
                cmd4 = 'ln -s %s %s' %(snp, newdir)
                #print cmd1
                #print cmd2
                #print cmd3
                #print cmd4
                os.system(cmd1)
                os.system(cmd2)
                os.system(cmd3)
                os.system(cmd4)
            else:
                bam = rils_all_bams_fix_ID_lib[ril][lib]
                bai = re.sub(r'\.bam', r'.bai', bam)
                flag= re.sub(r'.recal.bam', r'.dedup.flagstat', bam)
                newbam = '%s/GN%s.bam' %(newdir, ril)
                newbai = re.sub(r'\.bam', r'.bai', newbam)
                newflag= re.sub(r'\.bam', r'.dedup.flagstat', newbam)
                cmd1 = 'ln -s %s %s' %(bam, newbam)
                cmd2 = 'ln -s %s %s' %(bai, newbai)
                cmd3 = 'ln -s %s %s' %(flag, newflag)
                print cmd1
                print cmd2
                print cmd3
                os.system(cmd1)
                os.system(cmd2)
                os.system(cmd3)
'''
    newdir = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/%s' %(args.project)
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    count = 0    
    rils = readtable(args.input)
    bams = glob.glob('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/*.bam')
    ofile = open('../BWA.sampleRIL.list', 'w')
    for bam in sorted(bams):
        prefix, ext = os.path.splitext(bam)
        ril = os.path.basename(prefix)
        if rils.has_key(ril):
            continue
        snp = '%s.Maq.p1.map.pileup.SNP' %(prefix)
        cmd1 = 'ln -s %s %s' %(bam, newdir)
        cmd2 = 'ln -s %s %s' %(snp, newdir)
        cmd3 = 'cp %s %s/Exclued.list' %(args.input, newdir)
        print >> ofile, ril
        #print cmd1
        #print cmd2
        #os.system(cmd1)
        #os.system(cmd2)
        #os.system(cmd3)
        count += 1
    print 'Job done: Linked %s' %(count) 
''' 

if __name__ == '__main__':
    main()

