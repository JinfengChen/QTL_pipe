#!/bin/sh

scripts=`pwd`/scripts

perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/RIL_SNP_MAQ_snp.pl --ref /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/reference/MSU_r7.fa --parent /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/NB.RILs.dbSNP.SNPs.parents --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_ALL

echo "done"

