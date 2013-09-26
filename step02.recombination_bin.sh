#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00

export R_LIBS="/rhome/cjinfeng/software/tools/R-2.15.3/library/"

cd $PBS_O_WORKDIR

scripts=$PBS_O_WORKDIR/scripts
#construct recombination bin using MPR package
#cat $scripts/recombination_map/MPR_hmmrun.R | /rhome/cjinfeng/software/tools/R-2.15.3/bin/R --slave
#only run R/qtl steps in the scripts to fill and uniq recombination bin, and output to files.
#deletion the first empty element in MPR.geno.bin!!!!!!!
#cat $scripts/recombination_map/MPR_hmmrun_Rqtl.R | /rhome/cjinfeng/software/tools/R-2.15.3/bin/R --slave
#draw bin map for each RILs and for each chromosome
perl $scripts/recombination_map/RIL_drawbin.pl --MPR ./ --chrlen ../input/reference/MSU7.chr.inf

echo "Done"
