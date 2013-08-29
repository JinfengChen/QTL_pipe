#!/bin/bash

scripts=`pwd`/scripts
##get 0.5X of reads of each RIL
#perl $scripts/fastq/prefastq_qsub.pl --RIL /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/Illumina --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_0.5X > prefastq.list

##just link fastq file
perl $scripts/fastq/prefastq_ln.pl --RIL /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/Illumina --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_ALL > prefastq.list

echo "Done"


