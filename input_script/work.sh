#0. MSU_r7.corrected is the default genotype results from Sofia's pipeline, which including multi-libraries for some RILs and even incorrect due to ID error.

#1. Bam_fixID is the working directory that fixing ID errors and testin. After testing we will have formal directory for both fastq and Bam.
#1.1 Fix FC251 and remove landrace sequence from directory
cd ..
python Fix_Bam_ID.py
#1.2 Calculate similarity bewteen all the pair libraries
python Fix_Bam_ID_SNP_similarity.py --input Bam_fixID --cpu 32 > Bam_fixID.SNP.similarity
#or
qsub Fix_Bam_ID_SNP_similarity.sh
#1.3 Calculate coverage for all bam
python Fix_Bam_ID_Bam_Stat.py --input Bam_fixID
#1.4 Multi libraries clean-up. 
python Fix_Bam_ID_multi_lib.py > Bam_fixID.info
#1.5 generate RILs_ALL_Bam_fix_ID for QTL pipe
python Fix_Bam_ID_QTL_bam.py --input Bam_fixID > RILs_ALL_bam_fixID.changed_id.list
python Fix_Bam_ID_Bam_Stat.py --input RILs_ALL_bam_fixID

#2. Bam_correct is the dir formally correct ID.
#2.1 Calculate coverage for all bam
python Fix_Bam_ID_Bam_Stat.py --input Bam_correct

