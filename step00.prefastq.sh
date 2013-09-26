#!/bin/bash

scripts=`pwd`/scripts
##get 0.5X of reads of each RIL
#perl $scripts/fastq/prefastq_qsub.pl --RIL /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/Illumina --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_0.5X > prefastq.list

##just link fastq file
perl $scripts/fastq/prefastq_ln.pl --RIL /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/Illumina --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_ALL > prefastq.list

##Get 10X of reads for each RILs, if fewer than 10X just link
#perl $scripts/fastq/prefastq_lowmem_qsub.pl --RIL /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/Illumina --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_ALL > prefastq.2.list


echo "Done"


