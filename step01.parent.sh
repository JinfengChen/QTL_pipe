#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR
scripts=$PBS_O_WORKDIR/scripts

#perl $scripts/genotype/RIL_VCF2Parents.pl --vcf ../input/reference/HEG4_dbSNP.vcf
#perl $scripts/genotype/RIL_VCF2Parents.pl --vcf ../input/reference/HEG4_dbSNP.VQSR.vcf
#perl $scripts/genotype/RIL_VCF2Parents.pl --vcf ../input/reference/ALL.gatk.snp.pass.vcf
perl $scripts/genotype/RIL_VCF2Parents.pl --vcf ../input/reference/ALL.gatk.snp.VQSR.pass.vcf


echo "done"

