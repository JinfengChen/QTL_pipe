QTL_pipe
========

The pipeline do the job of mapping QTL in RILs.

1. Preparation

1.1 trait

create "trait dir" and put trait file in "trait dir"

     mkdir ../input/trait
     cd ../input/trait
     cp /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL/input/trait/May28_2013.RIL.trait.table ./

parse trait file, generate trait matrix for parents and RILs, draw distribution fig for each trait.

     perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scritps/trait/RIL_trait.pl --trait ../input/trait/May28_2013.RIL.trait.table

The results files will be in "trait dir":

May28_2013.RIL.trait.table.QTL.parents.txt

May28_2013.RIL.trait.table.QTL.trait.txt

DrawQTLtrait.R

DrawQTLtrait.pdf

1.2 reference

Create "reference dir" and put reference sequence, dbSNPs in "reference dir"

     mkdir ../input/reference
     cd ../input/reference
     ln -s /rhome/cjinfeng/HEG4_cjinfeng/RILs/Depth_Evaluation/input/HEG4_dbSNP.vcf ./
     ln -s /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa ./

Format dbSNPs and reference sequence

     grep -v "#" HEG4_dbSNP.vcf | awk '{print $1"\t"$2"\t"$5"\t"$4}' > HEG4_dbSNP.table
     perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scritps/reference/formatfa.pl --fa MSU_r7.fa --project Nipponbare
     rm MSU_r7.fa
     ln -s MSU_r7.reform.fa MSU_r7.fa

Generate parent2 reference sequence (HEG4) 

     perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scritps/reference/PseudoMaker_cjinfeng.pl HEG4_dbSNP.table MSU_r7.fa HEG4


1.3 fastq

create "fastq dir" and put fastq of RILs in  "fastq dir"

     mkdir ../input/fastq/
     cd ../input/fastq/
     ln -s /rhome/cjinfeng/Rice/RIL/Illumina/ ./

Get sample of fastq of RILs (0.1 X)

     perl /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/scritps/fastq/prefastq.pl --RIL ./Illumina > prefastq.list

Or run using qsub

     qsub runprefastq.sh

Put all sample into a project dir

     mkdir RIL_1X
     mv GN* RIL_1X

2. Mapping Reads (Maq, replace with bwa?)

3. Recombination Map (MPR package, Xie et al, 2010 PNAS)

4. QTL using R (Rqtl, Broman et al, 2003 Bioinfromatics)
