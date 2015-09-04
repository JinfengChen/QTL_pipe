library(qtl)
pdf ("MPR.cross.fill.QTL.pdf")
## step0. read and write the data of cross
if (1){
   read.cross("qtlcart",dir="./",file="MPR.cross.fill.cro",mapfile="MPR.cross.fill.map") -> cross
}else{
   read.cross("csv",dir="./",file="MPR.cross.fill") -> cross
}

myMaq <- pull.map(cross)
plot.map(myMaq)

if (0){
write.cross(cross,"qtlcart","MPR.cross.fill.originalmap")
convert2riself(cross) -> cross.riself
write.cross(cross.riself,"qtlcart","MPR.cross.fill.riself")
est.map(cross) -> reestimate.map
replace.map(cross,reestimate.map) -> cross.reestimate
write.cross(cross.reestimate,"qtlcart","MPR.cross.fill.reestimatemap")
}

## step1. do single QTL analysis
## 1.1 marker regression
QTL.mr <- scanone(cross,pheno.col=c(1,2,3,4,5,6,7,8,9,10),method="mr")

## 1.2 do interval method stardand invertal mapping(em), Haley-Knott regression (hk), Extended HK regression (ehk), and Multiple imputation (imp)
cross <- calc.genoprob(cross,step=1,error.prob=0.001)

## 1.2.1 estimate the threshold of LOD using permutation test
operm <- scanone(cross,pheno.col=c(1,2,3,4,5,6,7,8,9,10),n.perm=100, verbose=FALSE)
LOD <- summary(operm,alpha=0.05) ## get LOD threshold
QTL.mr.test <- summary(QTL.mr,perms=operm,alpha=0.05,pvalues=TRUE)
#if (FALSE){
## 1.2.2 interval mapping
QTL.em <- scanone(cross, pheno.col=c(1,2,3,4,5,6,7,8,9,10),method="em")
if (1){
QTL.hk <- scanone(cross, pheno.col=c(1,2,3,4,5,6,7,8,9,10),method="hk")
QTL.ehk <- scanone(cross, pheno.col=c(1,2,3,4,5,6,7,8,9,10),method="ehk")
cross.imp <- sim.geno(cross,step=1,n.draws=64,error.prob=0.001)
QTL.imp <- scanone(cross.imp,pheno.col=c(1,2,3,4,5,6,7,8,9,10), method="imp")
}
## 1.2.3 compostion interval mapping
#QTL.cim.20 <- cim(cross,pheno.col=c(1,2,3,4,5,6,7,8,9,10),n.marcovar=3,window=20)
#QTL.cim.20.test <- summary(QTL.cim.20,perms=operm,pvalues=TRUE)

if (FALSE){

## 1.2.4 estimate QTL effect
QTL.cim.20.final <- summary(QTL.cim.20,perms=operm,alpha=0.2,pvalues=TRUE)
qtl <- makeqtl(cross.imp,chr=QTL.cim.20.final$chr,pos=QTL.cim.20.final$pos)
qtlname <- paste(qtl$altname,collapse="+")
formula <- paste("y~",qtlname)
qtlfit <- fitqtl(cross.imp,qtl=qtl,formula=formula,get.ests=TRUE)
qtlvar <- qtlfit$result.drop[,4]
qtleffect <- c(qtlfit$ests$ests[2],qtlfit$ests$ests[4])
qtlsummary <- cbind(QTL.cim.20.final, qtlvar, qtleffect)

}

## 2 write result into table and pdf
## 2.1 write QTL loci and LOD into table
write.table(QTL.mr,"MPR.cross.fill.QTL.mr.table",sep="	",col.names=NA)
write.table(LOD,"MPR.cross.fill.QTL.mr.table.LOD_threshold",sep="	",col.names=NA)
write.table(QTL.mr.test,"MPR.cross.fill.QTL.mr.table.test",sep="	",col.names=NA)
write.table(QTL.em,"MPR.cross.fill.QTL.em.table",sep="	",col.names=NA)
if(0){
write.table(QTL.hk,"MPR.cross.fill.QTL.hk.table",sep="	",col.names=NA)
write.table(QTL.ehk,"MPR.cross.fill.QTL.ehk.table",sep="	",col.names=NA)
write.table(QTL.imp,"MPR.cross.fill.QTL.imp.table",sep="	",col.names=NA)
}
#write.table(QTL.cim.20,"MPR.cross.fill.QTL.cim.20",sep="	",col.names=NA)
#write.table(QTL.cim.20.test,"MPR.cross.fill.QTL.cim.20.test",sep="	",col.names=NA)
#write.table(qtlsummary,"MPR.cross.fill.QTL.cim.20.summary",sep="	",col.names=NA)

## 2.2 write QTL curves into pdf
#pdf ("MPR.cross.fill.QTL.pdf")
for (trait in 1:10){
     traitname <- colnames(cross$pheno)[trait]
     ## draw all chr together
     chr <- 1:12 
        plot(QTL.mr,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (MR)", col="#00000090")
        plot(QTL.em,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (EM)", col="#FF000090")
if(0){
        plot(QTL.hk,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (HK)", col="#80800090")
        plot(QTL.ehk,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (EHK)", col="#00800090")
        plot(QTL.imp,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (IMP)", col="#0000FF90")
        #plot(QTL.cim.20,chr=c(chr),lodcolumn=trait,main=traitname,ylab="LOD score (CIM)", col="#80000090")
}        
        ## max of three QTL in one plot, so we use add=TRUE to plot more QTL
        plot(QTL.mr,chr=c(chr),lodcolumn=trait,main=traitname,lwd=0.2,ylab="LOD score",col="#00000090") ## black
        plot(QTL.em,QTL.hk,QTL.ehk,chr=c(chr),lodcolumn=trait,lwd=0.2,col=c("#FF000090","#80800090","#00800090"),add=TRUE) ## red, olive, green
        plot(QTL.imp,chr=c(chr),lodcolumn=trait,lwd=0.2,col="#0000FF90",add=TRUE) ## blue
        #plot(QTL.cim.20,chr=c(chr),lodcolumn=trait,lwd=0.2,col="#80000090",add=TRUE) ## Maroon
        ## plot final QTL that parsed threshod
        #plot(qtl)
    
    ## draw chr by chr
    for (chr in 1:12){
        xtitle <- paste("Chr",chr," Map position (cM)",sep="")
        plot(QTL.mr,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (MR)", col="#00000090")
        plot(QTL.em,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (EM)", col="#FF000090")
if(0){
        plot(QTL.hk,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (HK)", col="#80800090")
        plot(QTL.ehk,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (EHK)", col="#00800090")
        plot(QTL.imp,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (IMP)", col="#0000FF90")
        #plot(QTL.cim.20,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,ylab="LOD score (CIM)", col="#80000090")
}        
        ## max of three QTL in one plot, so we use add=TRUE to plot more QTL
        plot(QTL.mr,chr=c(chr),lodcolumn=trait,main=traitname,xlab=xtitle,lwd=0.2,ylab="LOD score",col="#00000090") ## black
        plot(QTL.em,QTL.hk,QTL.ehk,chr=c(chr),lodcolumn=trait,lwd=0.2,col=c("#FF000090","#80800090","#00800090"),add=TRUE) ## red, olive, green
        plot(QTL.imp,chr=c(chr),lodcolumn=trait,lwd=0.2,col="#0000FF90",add=TRUE) ## blue
        #plot(QTL.cim.20,chr=c(chr),lodcolumn=trait,lwd=0.2,col="#80000090",add=TRUE) ## Maroon
        ## plot final QTL that parsed threshod
        #plot(qtl)}
    }
}

dev.off()



