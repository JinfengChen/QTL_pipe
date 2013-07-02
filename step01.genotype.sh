#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR

perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scritps/genotype/RIL_VCF2Parents.pl --vcf ../input/reference/HEG4_dbSNP.vcf
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scritps/genotype/RIL_SNP_MAQ.pl --ref ../input/reference/MSU_r7.fa --fastq ../input/fastq/002 --parents NB.RILs.dbSNP.SNPs.parents

echo "done"

