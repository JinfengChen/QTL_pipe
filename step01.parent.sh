#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=10:00:00
#SBATCH --output=step01.parent.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./

scripts=`pwd`/scripts

#perl $scripts/genotype/RIL_VCF2Parents.pl --vcf ../input/reference/HEG4_dbSNP.vcf
#SNP called by Sofia using this one, try this one, 20160617
perl $scripts/genotype/RIL_VCF2Parents.pl --vcf ../input/reference/HEG4_dbSNP.VQSR.vcf
#perl $scripts/genotype/RIL_VCF2Parents.pl --vcf ../input/reference/ALL.gatk.snp.pass.vcf
#using this one for long time, 20160617
#perl $scripts/genotype/RIL_VCF2Parents.pl --vcf ../input/reference/ALL.gatk.snp.VQSR.pass.vcf


echo "done"

