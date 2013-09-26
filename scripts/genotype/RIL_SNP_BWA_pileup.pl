#!/usr/bin/perl
=header
This script is design to run pileup and SNP calling using qsub, which speed up the process in RIL_SNP_MAQ.pl.
--ref:   reference sequence
--parents: parents file
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

my @map=glob("$opt{fastq}/GN*.bam");

pileup(\@map);


sub pileup
{
my ($map)=@_;
my @cmd;
my $maq="/opt/tyler/bin/maq";
my $samtools="/usr/local/bin/samtools";
for(my $i=0; $i< @$map; $i++){
   my $prefix=$1 if ($map->[$i]=~/(.*)\.bam$/);
   push @cmd, "$samtools mpileup -q 40 -Q 15 -f $opt{ref} $map->[$i] | awk '\$4 > 0' > $prefix.Maq.p1.map.pileup" unless (-e "$prefix.Maq.p1.map.pileup.SNP");
   open OUT, ">pileup.sh" or die "$!";
   for(my $i=0; $i<@cmd; $i++){
      print OUT "$cmd[$i]\n";
   }
   close OUT;
}# for loop
`perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --convert no --maxjob 10 --resource nodes=1:ppn=1,mem=5G,walltime=100:00:00 pileup.sh`;
}# end of sub function


