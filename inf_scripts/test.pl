my $b = ',-1t,^O,..+2AG.^],$.CcTTt';
my $q = "???CB??!53..3";

print "$b\n";
if ($b =~ /([-|+]+[0-9]+[ACGTNacgtn]+)/){
   $b =~ s/([-|+]+[0-9]+[ACGTNacgtn]+)//g;
   print "yes: $1, $b\n";
}

if ($b =~ /([^.,ACGTNacgtn]+)/){
   $b =~ s/([^.,ACGTNacgtn]+)//g;
   print "yes: $1, $b\n";
}

#while($b=~/\d/){
#
#}

my @base=split("",$b);
my @qual=split("",$q);
for(my $i=1;$i<@base;$i++){
   my $score = ord($qual[$i]) - 33;
   print "$base[$i]\t$score\n";
}

