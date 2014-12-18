#!/bin/sh

scripts=`pwd`/scripts

## MSU7 orignal /rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa
genome=/rhome/cjinfeng/BigData/00.RD/seqlib/MSU7_samtools0_1_16/MSU_r7.fa
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/RIL_SNP_BWA_pileup.pl --ref $genome --parent /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/NB.RILs.dbSNP.SNPs.parents --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_ALL_bam


#perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/RIL_SNP_MAQ_snp.pl --ref $genome --parent /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/NB.RILs.dbSNP.SNPs.parents --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_ALL_bam

perl scripts/trait/subtrait.pl --trait ../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.ALL --maqlist BWA.sampleRIL.list > ../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.first

echo "done"

