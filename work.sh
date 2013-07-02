
echo "prepare trait"
mkdir ../input/trait
cp /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL/input/trait/May28_2013.RIL.trait.table ./
cd ../input/trait
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scritps/trait/RIL_trait.pl --trait ../input/trait/May28_2013.RIL.trait.table


echo "prepare reference"
mkdir ../input/reference
cd ../input/reference
ln -s /rhome/cjinfeng/HEG4_cjinfeng/RILs/Depth_Evaluation/input/HEG4_dbSNP.vcf ./
ln -s /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa ./
grep -v "#" HEG4_dbSNP.vcf | awk '{print $1"\t"$2"\t"$5"\t"$4}' > HEG4_dbSNP.table
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scritps/reference/formatfa.pl --fa MSU_r7.fa --project Nipponbare
rm MSU_r7.fa
ln -s MSU_r7.reform.fa MSU_r7.fa
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scritps/reference/PseudoMaker_cjinfeng.pl HEG4_dbSNP.table MSU_r7.fa HEG4


echo "prepare fastq"
mkdir ../input/fastq/
cd ../input/fastq/
ln -s /rhome/cjinfeng/Rice/RIL/Illumina/ ./
qsub runprefastq.sh
##delete duplicate
rm GN39a_*
mkdir RILs_1x
mv GN* RILs_1X

echo "00.Mapping reads and pileup"
bash step00.mapping.sh

echo "01.genotype RILs with SNP called from resequencing HEG4"
echo "convert VCF SNPs into parents file"
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scritps/genotype/RIL_VCF2Parents.pl --vcf ../input/reference/HEG4_dbSNP.vcf
echo "genotype Maq results of RILs using parents file"
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scritps/genotype/RIL_SNP_MAQ.pl --ref ../input/reference/MSU_r7.fa --fastq ../input/fastq/002 --parents NB.RILs.dbSNP.SNPs.parents
echo "Or run two step 'genotype' in qsub"
qsub -q js step01.genotype.sh


echo "02.constructe recombiantion bin and draw bin"
echo "construct recombination bin using MPR package"
cat /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scritps/recombination_map/MPR_hmmrun.R | /rhome/cjinfeng/software/tools/R-2.15.3/bin/R --slave
#draw bin map for each RILs and for each chromosome
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scritps/recombination_map/RIL_drawbin.pl --MPR ./ --chrlen ../input/reference/MSU7.chr.inf
echo "Or run two step in qsub"
qsub -q js step02.recombination_bin.sh

echo "03.QTL"
qsub -q js step03.QTL.sh


