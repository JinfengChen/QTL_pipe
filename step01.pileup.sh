#!/bin/sh

scripts=`pwd`/scripts

#perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/RIL_SNP_MAQ_pileup.pl --ref /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/reference/MSU_r7.fa --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_0.5X

## MSU7 Chr -> chromosome
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/RIL_SNP_BWA_pileup.pl --ref /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/reference/MSU_r7.fa --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_ALL_bam

## MSU7 orignal /rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa
#perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/RIL_SNP_BWA_pileup.pl --ref /rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_ALL

echo "done"

