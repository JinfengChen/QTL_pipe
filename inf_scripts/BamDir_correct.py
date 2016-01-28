#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import glob

def usage():
    test="name"
    message='''
python BamDir.py --project RILs_ALL_bam_correct

link bam and SNP file from RILs_ALL_bam to a new dir "project".
After this change the dir name in "step01.genotype.sh" and rerun the QTL pipeline:
qsub step01.parent.sh
qsub -q js step01.genotype.sh
qsub -q js step02.recombination_bin.sh
qsub -q js step03.QTL.sh

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#new sequence RIL, better than old one but may be depth is less
def must_use_list(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                data[unit[0]] = unit[1]
    return data


#Sample  #Read   Average Total   Depth   Mapped_Depth    Mapped_rate     #Library        FileName
#RIL1    23383054        100     2338305400      6.2857672043    6.13189677419   0.975520819479  1       Bam_correct/RIL1_0_CGTACG_FC153L5.recal.bam
def readtable(infile):
    must_use_lib = must_use_list('BamDir_correct.list')
    data = defaultdict(lambda : defaultdict(lambda : str()))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'RIL'): 
                unit = re.split(r'\t',line)
                ril  = unit[0]
                bam  = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/%s' %(unit[-1])
                depth= float(unit[5])
                print line
                if data.has_key(ril):
                    print 'has key: %s' %(ril)
                    if must_use_lib.has_key(ril):
                        print 'must use lib has key: %s, %s' %(data[ril][0], bam)
                        #have must_use new lib, search for new lib
                        if os.path.split(data[ril][0])[1] == os.path.split(must_use_lib[ril])[1]:
                            print 'already set must_use_lib to new lib: %s' %(data[ril][0])
                            #already set the new lib, do nothing
                            continue
                        else:
                            print 'set must_use_lib to new lib: %s, %s' %(bam, depth)
                            #new lib not set, go on searching
                            data[ril][0] = bam
                            data[ril][1] = depth
                    else:
                        print 'comparing depth'
                        #do not have muse_use new lib, use higher depth
                        if depth > float(data[ril][1]):
                            data[ril][0] = bam
                            data[ril][1] = depth
                else:
                    print 'new %s' %(ril)
                    data[ril][0] = bam
                    data[ril][1] = depth
    return data


def subset_bam_stat(sublist, stat, substat):
    libs = defaultdict(lambda : str())
    with open (sublist, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                libs[unit[1]] = unit[0]
    ofile = open(substat, 'w')
    with open (stat, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'RIL'): 
                unit = re.split(r'\t',line)
                lib = re.sub(r'.recal.bam', r'', os.path.basename(unit[-1]))
                if libs.has_key(lib):
                    print >> ofile, line
            else:
                print >> ofile, line
    ofile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.project) > 0
    except:
        usage()
        sys.exit(2)

    if args.project is None:
        args.project = 'RILs_ALL_bam_correct'
 
    data = readtable('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/Bam_correct.bam.stat')
    
    newdir = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/%s' %(args.project)
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    count = 0
    ofile = open('%s.RILs.list' %(newdir), 'w') 
    #bams = glob.glob('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/Bam_correct/*.recal.bam')
    for ril in sorted(data.keys()):
        #print ril
        bam    = data[ril][0]
        prefix = re.sub(r'.recal.bam', r'', bam)
        rilid = re.split(r'_', os.path.basename(prefix))[0]
        rilid = re.sub(r'RIL', r'GN', rilid)
        snp_tab = '%s.genotype.tab.gz' %(prefix)
        snp_vcf = '%s.genotype.vcf.gz' %(prefix)
        stat = '%s.dedup.flagstat' %(prefix)
        bai  = '%s.recal.bai' %(prefix)
        cmd1 = 'ln -s %s %s/%s.bam' %(bam, newdir, rilid)
        cmd2 = 'ln -s %s %s/%s.genotype.tab.gz' %(snp_tab, newdir, rilid)
        cmd3 = 'ln -s %s %s/%s.genotype.vcf.gz' %(snp_vcf, newdir, rilid)
        cmd4 = 'ln -s %s %s/%s.dedup.flagstats' %(stat, newdir, rilid)
        cmd5 = 'ln -s %s %s/%s.bam.bai' %(bai, newdir, rilid)
        #print cmd1
        #print cmd3
        #print cmd4
        
        if not os.path.exists('%s/%s.bam' %(newdir, rilid)):
            os.system(cmd1)
            count += 1
            print 'linked: %s\t%s' %(ril, bam)
            
        else:
            old_bam = os.path.realpath('%s/%s.bam' %(newdir, rilid))
            old_bam = os.path.basename(old_bam)
            new_bam = os.path.basename(bam)
            print 'old bam: %s' %(old_bam)
            print 'new bam: %s' %(old_bam)
            if new_bam == old_bam:
                print 'same'
                print >> ofile, '%s\t%s' %(rilid, os.path.basename(prefix))
                #same file, we do nothing and skip to next file
                continue
            else:
                print 'replaced'
                #different file, new bam must have higher coverage, update to new one
                os.system('rm %s/%s.bam' %(newdir, rilid))
                os.system('rm %s/%s.genotype.tab.gz' %(newdir, rilid))
                os.system('rm %s/%s.genotype.vcf.gz' %(newdir, rilid))
                os.system('rm %s/%s.dedup.flagstats' %(newdir, rilid))
                os.system('rm %s/%s.bam.bai' %(newdir, rilid))
                os.system('rm -R %s/%s_stats' %(newdir, rilid))
                os.system(cmd1)
                count += 1
                print 'linked: %s\t%s' %(ril, bam)
        if not os.path.exists('%s/%s.genotype.tab.gz' %(newdir, rilid)):
            os.system(cmd2)
        if not os.path.exists('%s/%s.genotype.vcf.gz' %(newdir, rilid)):
            os.system(cmd3)
        if not os.path.exists('%s/%s.dedup.flagstats' %(newdir, rilid)):
            os.system(cmd4)
        if not os.path.exists('%s/%s.bam.bai' %(newdir, rilid)):
            os.system(cmd5)
        print >> ofile, '%s\t%s' %(rilid, os.path.basename(prefix))
        #count += 1
    print 'Job done: Linked %s' %(count) 
    ofile.close()

    #get subset bam for selected RIL libraries
    subset_bam_stat('%s.RILs.list' %(newdir), '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/Bam_correct.bam.stat', '%s.bam.stat' %(newdir))
 
if __name__ == '__main__':
    main()

