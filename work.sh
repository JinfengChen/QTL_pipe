
echo "prepare trait"
mkdir ../input/trait
cp /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL/input/trait/May28_2013.RIL.trait.table ./
cd ../input/trait
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/trait/RIL_trait.pl --trait ../input/trait/May28_2013.RIL.trait.table


echo "prepare reference"
mkdir ../input/reference
cd ../input/reference
ln -s /rhome/cjinfeng/HEG4_cjinfeng/RILs/Depth_Evaluation/input/HEG4_dbSNP.vcf ./
ln -s /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU7.chr.inf ./
ln -s /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa ./
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/reference/formatfa.pl --fa MSU_r7.fa --project Nipponbare
rm MSU_r7.fa
ln -s MSU_r7.reform.fa MSU_r7.fa


echo "prepare fastq"
mkdir ../input/fastq/
cd ../input/fastq/
ln -s /rhome/cjinfeng/Rice/RIL/Illumina/ ./
mkdir RIL_0.5X
cd ../../bin/
bash step00.prefastq.sh

echo "00.Mapping reads and pileup"
bash step00.mapping.sh

echo "01.genotype RILs with SNP called from resequencing HEG4"
echo "convert VCF SNPs into parents file"
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/RIL_VCF2Parents.pl --vcf ../input/reference/HEG4_dbSNP.vcf
echo "Run pileup and SNP using qsub before genotype, speed up the process (pileup already build in in Mapping process, so just run snp here)"
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/RIL_SNP_MAQ_pileup.pl --ref /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/reference/MSU_r7.fa --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_0.5X 
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/RIL_SNP_MAQ_snp.pl --ref /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/reference/MSU_r7.fa --parent /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/NB.RILs.dbSNP.SNPs.parents --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_0.5X
echo "genotype Maq results of RILs using parents file"
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/genotype/RIL_SNP_MAQ.pl --ref ../input/reference/MSU_r7.fa --fastq ../input/fastq/RILs_0.5X --parents NB.RILs.dbSNP.SNPs.parents

echo "or Run in shell"
qsub -q js step01.parent.sh
bash step01.snp.sh
qsub -q js step01.genotype.sh

echo "or Run just step01.genotype.sh, then should be very slow"
qsub -q js step01.genotype.sh

echo "add mping and cold trait"
perl scripts/trait/AddmPingTrait.pl --mping ../input/trait/Oct10_2013.mPing.table
perl scripts/trait/AddmPingTrait_USDA.pl

echo "get sub trait"
perl scripts/trait/subtrait.pl --trait ../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt --maqlist MAQ.sampleRIL.list > ../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.first
perl scripts/trait/subtrait.pl --trait ../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.ALL --maqlist BWA.sampleRIL.list > ../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.first
echo "trait correlation"
python scripts/trait/TraitCorrelation.py --input ../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.ALL > trait.cr.txt
python ../../bin/scripts/trait/TraitPlot.py --trait May28_2013.RIL.trait.table.QTL.trait.txt.ALL --parent May28_2013.RIL.trait.table.QTL.parents.txt

echo "02.constructe recombiantion bin and draw bin"
echo "construct recombination bin using MPR package"
cat /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/recombination_map/MPR_hmmrun.R | /rhome/cjinfeng/software/tools/R-2.15.3/bin/R --slave
#draw bin map for each RILs and for each chromosome
perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scripts/recombination_map/RIL_drawbin.pl --MPR ./ --chrlen ../input/reference/MSU7.chr.inf
echo "Or run two step in qsub"
qsub -q js step02.recombination_bin.sh

echo "03.QTL"
qsub -q js step03.QTL.sh



echo "Use bam files from Sofia directly"
qsub step01.parent.sh
bash step01.bam2snp.sh
qsub step01.genotype.sh
qsub -q js step02.recombination_bin.sh
qsub step03.QTL.sh

