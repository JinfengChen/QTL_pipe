#0. MSU_r7.corrected is the default genotype results from Sofia's pipeline, which including multi-libraries for some RILs and even incorrect due to ID error.

#1. Bam_fixID is the working directory that fixing ID errors and testin. After testing we will have formal directory for both fastq and Bam.
#1.1 Fix FC251 and remove landrace sequence from directory
cd ..
python Fix_Bam_ID.py
#1.2 Calculate similarity bewteen all the pair libraries
python Fix_Bam_ID_SNP_similarity.py --input Bam_fixID --cpu 32 > Bam_fixID.SNP.similarity
#or
qsub Fix_Bam_ID_SNP_similarity.sh
#1.3 Multi libraries clean-up. 
python Fix_Bam_ID_multi_lib.py


