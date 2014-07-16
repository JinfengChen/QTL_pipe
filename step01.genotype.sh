#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR
scripts=$PBS_O_WORKDIR/scripts

#perl $scripts/genotype/RIL_VCF2Parents.pl --vcf ../input/reference/HEG4_dbSNP.vcf
#perl $scripts/genotype/RIL_SNP_MAQ.pl --ref ../input/reference/MSU_r7.fa --fastq ../input/fastq/RILs_0.5X --parents NB.RILs.dbSNP.SNPs.parents
#perl $scripts/genotype/RIL_SNP_MAQ.pl --ref ../input/reference/MSU_r7.fa --fastq ../input/fastq/RILs_ALL --parents NB.RILs.dbSNP.SNPs.parents
perl $scripts/genotype/RIL_SNP_BWA.pl --ref ../input/reference/MSU_r7.fa --fastq ../input/fastq/RILs_ALL_bam --parents NB.RILs.dbSNP.SNPs.parents

echo "done"

