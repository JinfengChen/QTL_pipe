#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --time=40:00:00
#SBATCH --output=step03.QTL.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./

scripts=`pwd`/scripts

export R_LIBS="/bigdata/stajichlab/cjinfeng/software/R_package_MPR/"
#export R_LIBS="/rhome/cjinfeng/software/tools/R-2.15.3/library/"

perl $scripts/QTL/RIL_QTL.pl --qtlcart MPR.cross.uniq

echo "Done"

