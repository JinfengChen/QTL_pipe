#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"csv:s","qtlcart:s","project:s","help");


my $help=<<USAGE;
perl $0 --csv test.csv 
perl $0 --qtlcart test
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

if ($opt{qtlcart}){
   $opt{project} ||= $opt{qtlcart};
   Rqtl($opt{qtlcart},"qtlcart");
}else{
   $opt{project} ||= $1.".Rqtl" if ($opt{csv}=~/(.*)\.csv/);
   Rqtl($opt{csv},"csv");
}


sub Rqtl
{
my ($qtl,$type)=@_;
my $cross="cross";
my $tranperant=90;
my $flag= $type eq "qtlcart" ? 1 : 0;
my $cmd =<<R;
library(qtl)
pdf ("$opt{project}.QTL.pdf")
## step0. read and write the data of cross
if ($flag){
   read.cross("qtlcart",dir="./",file="$qtl.cro",mapfile="$qtl.map") -> $cross
}else{
   read.cross("csv",dir="./",file="$qtl") -> $cross
}

myMaq <- pull.map($cross)
plot.map(myMaq)

if (0){
write.cross($cross,"qtlcart","$opt{project}.originalmap")
convert2riself($cross) -> $cross.riself
write.cross($cross.riself,"qtlcart","$opt{project}.riself")
est.map($cross) -> reestimate.map
replace.map($cross,reestimate.map) -> $cross.reestimate
write.cross($cross.reestimate,"qtlcart","$opt{project}.reestimatemap")
}

## step1. do single QTL analysis
## 1.1 marker regression
QTL.mr <- scanone($cross,pheno.col=c(1,2,3,4,5),method="mr")

## 1.2 do interval method stardand invertal mapping(em), Haley-Knott regression (hk), Extended HK regression (ehk), and Multiple imputation (imp)
$cross <- calc.genoprob($cross,step=1,error.prob=0.001)

## 1.2.1 estimate the threshold of LOD using permutation test
operm <- scanone($cross,pheno.col=c(1,2,3,4,5),n.perm=100, verbose=FALSE)
LOD <- summary(operm,alpha=0.05) ## get LOD threshold
QTL.mr.test <- summary(QTL.mr,perms=operm,alpha=0.05,pvalues=TRUE)
#if (FALSE){
## 1.2.2 interval mapping
QTL.em <- scanone($cross, pheno.col=c(1,2,3,4,5),method="em")
if (1){
QTL.hk <- scanone($cross, pheno.col=c(1,2,3,4,5),method="hk")
QTL.ehk <- scanone($cross, pheno.col=c(1,2,3,4,5),method="ehk")
$cross.imp <- sim.geno($cross,step=1,n.draws=64,error.prob=0.001)
QTL.imp <- scanone($cross.imp,pheno.col=c(1,2,3,4,5), method="imp")
}
## 1.2.3 compostion interval mapping
#QTL.cim.20 <- cim($cross,pheno.col=c(1,2,3,4,5),n.marcovar=3,window=20)
#QTL.cim.20.test <- summary(QTL.cim.20,perms=operm,pvalues=TRUE)

if (FALSE){

## 1.2.4 estimate QTL effect
QTL.cim.20.final <- summary(QTL.cim.20,perms=operm,alpha=0.2,pvalues=TRUE)
qtl <- makeqtl($cross.imp,chr=QTL.cim.20.final\$chr,pos=QTL.cim.20.final\$pos)
qtlname <- paste(qtl\$altname,collapse="+")
formula <- paste("y~",qtlname)
qtlfit <- fitqtl($cross.imp,qtl=qtl,formula=formula,get.ests=TRUE)
qtlvar <- qtlfit\$result.drop[,4]
qtleffect <- c(qtlfit\$ests\$ests[2],qtlfit\$ests\$ests[4])
qtlsummary <- cbind(QTL.cim.20.final, qtlvar, qtleffect)

}

## 2 write result into table and pdf
## 2.1 write QTL loci and LOD into table
write.table(QTL.mr,"$opt{project}.QTL.mr.table",sep="\t",col.names=NA)
write.table(LOD,"$opt{project}.QTL.mr.table.LOD_threshold",sep="\t",col.names=NA)
write.table(QTL.mr.test,"$opt{project}.QTL.mr.table.test",sep="\t",col.names=NA)
write.table(QTL.em,"$opt{project}.QTL.em.table",sep="\t",col.names=NA)
if(0){
write.table(QTL.hk,"$opt{project}.QTL.hk.table",sep="\t",col.names=NA)
write.table(QTL.ehk,"$opt{project}.QTL.ehk.table",sep="\t",col.names=NA)
write.table(QTL.imp,"$opt{project}.QTL.imp.table",sep="\t",col.names=NA)
}
#write.table(QTL.cim.20,"$opt{project}.QTL.cim.20",sep="\t",col.names=NA)
#write.table(QTL.cim.20.test,"$opt{project}.QTL.cim.20.test",sep="\t",col.names=NA)
#write.table(qtlsummary,"$opt{project}.QTL.cim.20.summary",sep="\t",col.names=NA)

## 2.2 write QTL curves into pdf
#pdf ("$opt{project}.QTL.pdf")
for (trait in 1:5){
     traitname <- colnames(cross\$pheno)[trait]
     ## draw all chr together
     chr <- 1:12 
        plot(QTL.mr,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (MR)", col="#000000$tranperant")
        plot(QTL.em,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (EM)", col="#FF0000$tranperant")
if(0){
        plot(QTL.hk,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (HK)", col="#808000$tranperant")
        plot(QTL.ehk,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (EHK)", col="#008000$tranperant")
        plot(QTL.imp,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (IMP)", col="#0000FF$tranperant")
        #plot(QTL.cim.20,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (CIM)", col="#800000$tranperant")
}        
        ## max of three QTL in one plot, so we use add=TRUE to plot more QTL
        plot(QTL.mr,chr=c(chr),lodcolumn=trait,main=traitname,lwd=0.2,ylab="LOD score",col="#000000$tranperant") ## black
        plot(QTL.em,QTL.hk,QTL.ehk,chr=c(chr),lodcolumn=trait,lwd=0.2,col=c("#FF0000$tranperant","#808000$tranperant","#008000$tranperant"),add=TRUE) ## red, olive, green
        plot(QTL.imp,chr=c(chr),lodcolumn=trait,lwd=0.2,col="#0000FF$tranperant",add=TRUE) ## blue
        #plot(QTL.cim.20,chr=c(chr),lodcolumn=trait,lwd=0.2,col="#800000$tranperant",add=TRUE) ## Maroon
        ## plot final QTL that parsed threshod
        #plot(qtl)
    
    ## draw chr by chr
    for (chr in 1:12){
        xtitle <- paste("Chr",chr," Map position (cM)",sep="")
        plot(QTL.mr,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (MR)", col="#000000$tranperant")
        plot(QTL.em,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (EM)", col="#FF0000$tranperant")
if(0){
        plot(QTL.hk,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (HK)", col="#808000$tranperant")
        plot(QTL.ehk,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (EHK)", col="#008000$tranperant")
        plot(QTL.imp,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (IMP)", col="#0000FF$tranperant")
        #plot(QTL.cim.20,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (CIM)", col="#800000$tranperant")
}        
        ## max of three QTL in one plot, so we use add=TRUE to plot more QTL
        plot(QTL.mr,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,lwd=0.2,ylab="LOD score",col="#000000$tranperant") ## black
        plot(QTL.em,QTL.hk,QTL.ehk,chr=c(chr),lodcolumn=trait,lwd=0.2,col=c("#FF0000$tranperant","#808000$tranperant","#008000$tranperant"),add=TRUE) ## red, olive, green
        plot(QTL.imp,chr=c(chr),lodcolumn=trait,lwd=0.2,col="#0000FF$tranperant",add=TRUE) ## blue
        #plot(QTL.cim.20,chr=c(chr),lodcolumn=trait,lwd=0.2,col="#800000$tranperant",add=TRUE) ## Maroon
        ## plot final QTL that parsed threshod
        #plot(qtl)}
    }
}

dev.off()


R

open OUT, ">$opt{project}.R" or die "$!";
     print OUT "$cmd\n";
close OUT;
`cat $opt{project}.R | /rhome/cjinfeng/software/tools/R-2.15.3/bin/R --vanilla --slave`;

}
 
