/opt/samtools-0.1.16/samtools mpileup -q 29 -Q 15 -d 8000 -f /rhome/cjinfeng/BigData/00.RD/seqlib/MSU7_samtools0_1_16/MSU_r7.fa /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/GN231.bam | awk '$4 > 0' > /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/GN231.Maq.p1.map.pileup
