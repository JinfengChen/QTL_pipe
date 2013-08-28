#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00

export R_LIBS="/rhome/cjinfeng/software/tools/R-2.15.3/library/"
cd $PBS_O_WORKDIR

scripts=$PBS_O_WORKDIR/scripts

perl $scripts/QTL/RIL_QTL.pl --qtlcart MPR.cross.uniq

echo "Done"

