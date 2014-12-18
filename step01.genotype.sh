#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR
scripts=$PBS_O_WORKDIR/scripts

#perl $scripts/genotype/RIL_VCF2Parents.pl --vcf ../input/reference/HEG4_dbSNP.vcf
#perl $scripts/genotype/RIL_SNP_MAQ.pl --ref ../input/reference/MSU_r7.fa --fastq ../input/fastq/RILs_0.5X --parents NB.RILs.dbSNP.SNPs.parents
#perl $scripts/genotype/RIL_SNP_MAQ.pl --ref ../input/reference/MSU_r7.fa --fastq ../input/fastq/RILs_ALL --parents NB.RILs.dbSNP.SNPs.parents

#Trait
perl scripts/trait/subtrait.pl --trait ../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.ALL --maqlist BWA.sampleRIL.list > ../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.first

#ALL
#perl $scripts/genotype/RIL_SNP_BWA.pl --fastq ../input/fastq/RILs_ALL_bam --parents NB.RILs.dbSNP.SNPs.parents
#Core
#perl $scripts/genotype/RIL_SNP_BWA.pl --fastq ../input/fastq/RILs_ALL_bam_core --parents NB.RILs.dbSNP.SNPs.parents
#261
perl $scripts/genotype/RIL_SNP_BWA.pl --fastq ../input/fastq/RILs_ALL_bam_261 --parents NB.RILs.dbSNP.SNPs.parents
echo "done"

