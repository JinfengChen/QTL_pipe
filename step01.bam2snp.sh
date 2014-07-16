#!/bin/sh

scripts=`pwd`/scripts

## MSU7 orignal /rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/RIL_SNP_BWA_pileup.pl --ref /rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_ALL_bam

perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/RIL_SNP_MAQ_snp.pl --ref /rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa --parent /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/NB.RILs.dbSNP.SNPs.parents --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_ALL_bam

echo "done"

