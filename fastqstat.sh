#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR

#find /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL/GN*_1.fq | sed 's/_1/_?/' | sort > RIL.fastq.list
perl /rhome/cjinfeng/BigData/00.RD/fastq/errorcorrection/trim/fastq_stat.pl --list RIL.fastq.list --project RIL.fastq

echo "Done"
