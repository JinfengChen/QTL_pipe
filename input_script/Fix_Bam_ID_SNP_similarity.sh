#!/bin/bash
#PBS -l nodes=1:ppn=32
#PBS -l mem=100gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

bam_dir=Bam_correct
python Fix_Bam_ID_SNP_similarity.py --input $bam_dir --cpu $PBS_NP > $bam_dir\.SNP.similarity

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

