#!/bin/bash

scripts=$PBS_O_WORKDIR/scritps

perl $scripts/mapping/RIL_MAQ_qsub.pl --ref /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/reference/MSU_r7.fa --fastq /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_0.2X

echo "Done"

