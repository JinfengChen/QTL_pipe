#!/bin/bash

scripts=`pwd`/scripts
perl $scripts/fastq/prefastq_qsub.pl --RIL /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/Illumina --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_0.5X > prefastq.list

echo "Done"


