#!/usr/bin/perl
=header
The scripts is designed to generate SNPs for RILs using Maq pileup command. The results will be input for MPR package developed by Qifa Zhang's lab.
--ref:   reference sequence
--fastq: dir of fastq
=cut

use Getopt::Long;
use FindBin qw($Bin $Script);

my %opt;
GetOptions(\%opt,"ref:s","fastq:s","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 --ref ../input/reference/Pseudomolecules_Nipponbare.Build4.0.fasta --fastq ../input/testfastq\n";
   exit();
}

my $maq="/opt/tyler/bin/maq";
my @map=glob("$opt{fastq}/GN*.Maq.p1.map");

pileup(\@map);


sub pileup
{
my ($map)=@_;
my @cmd;
my $maq="/opt/tyler/bin/maq";
for(my $i=0; $i< @$map; $i++){
   push @cmd, "$maq pileup -vP -q 40 $opt{ref}.bfa $map->[$i] | awk '\$4 > 0' > $map->[$i].pileup" unless (-e "$map->[$i].pileup");
   open OUT, ">pileup.sh" or die "$!";
   for(my $i=0; $i<@cmd; $i++){
      print OUT "$cmd[$i]\n";
   }
   close OUT;
}# for loop
`perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --convert no --maxjob 10 --resource nodes=1:ppn=1,mem=5G,walltime=100:00:00 pileup.sh`;
}# end of sub function


