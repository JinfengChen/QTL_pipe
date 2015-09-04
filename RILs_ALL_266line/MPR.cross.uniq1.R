library(qtl)
pdf ("MPR.cross.uniq.QTL.pdf")
nph <- 18 ## number of phenotype
## step0. read and write the data of cross
if (1){
   read.cross("qtlcart",dir="./",file="MPR.cross.uniq.cro",mapfile="MPR.cross.uniq.map") -> cross
}else{
   read.cross("csv",dir="./",file="MPR.cross.uniq") -> cross
}

myMaq <- pull.map(cross)
plot.map(myMaq)

if (0){
write.cross(cross,"qtlcart","MPR.cross.uniq.originalmap")
convert2riself(cross) -> cross.riself
write.cross(cross.riself,"qtlcart","MPR.cross.uniq.riself")
est.map(cross) -> reestimate.map
replace.map(cross,reestimate.map) -> cross.reestimate
write.cross(cross.reestimate,"qtlcart","MPR.cross.uniq.reestimatemap")
}

## step1. do single QTL analysis
## 1.1 marker regression
QTL.mr <- scanone(cross,pheno.col=1:nph,method="mr")

## 1.2 do interval method stardand invertal mapping(em), Haley-Knott regression (hk), Extended HK regression (ehk), and Multiple imputation (imp)
cross <- calc.genoprob(cross,step=1,error.prob=0.001)

## 1.2.1 estimate the threshold of LOD using permutation test
operm <- scanone(cross,pheno.col=1:nph,n.perm=100, verbose=FALSE)
LOD <- summary(operm,alpha=0.05) ## get LOD threshold
QTL.mr.test <- summary(QTL.mr,perms=operm,alpha=0.05,format="allpeak",pvalues=TRUE)

## 1.2.2 interval mapping
QTL.em <- scanone(cross, pheno.col=1:nph,method="em")
if (1){
QTL.hk <- scanone(cross, pheno.col=1:nph,method="hk")
QTL.ehk <- scanone(cross, pheno.col=1:nph,method="ehk")
cross.imp <- sim.geno(cross,step=1,n.draws=64,error.prob=0.001)
QTL.imp <- scanone(cross.imp,pheno.col=1:nph, method="imp")
}
## 1.2.3 compostion interval mapping
#QTL.cim.20 <- cim(cross,pheno.col=1:nph,n.marcovar=3,window=20)
#QTL.cim.20.test <- summary(QTL.cim.20,perms=operm,pvalues=TRUE)

if (TRUE){ ###QTL effect switch
## 1.2.4 estimate QTL effect

cross.imp <- sim.geno(cross,step=1,n.draws=64,error.prob=0.001)
cross.imp.jit <- cross.imp
#cross.imp.jit <- jittermap(cross.imp)
#cross.imp.jit <- sim.geno(cross.imp.jit,step=1,n.draws=64,error.prob=0.001)

write.table("QTL summary by fitqtl","MPR.cross.uniq.QTL.fit.summary",row.names=FALSE, col.names=FALSE, sep="	")
for (i in seq(4,length(QTL.mr.test),by=3)){
   print (names(QTL.mr.test)[i-1])
   write.table(names(QTL.mr.test)[i-1],"MPR.cross.uniq.QTL.fit.summary",append=TRUE,row.names=FALSE, col.names=FALSE ,sep="	")
   select <- (QTL.mr.test[,i] <= 0.05) 
   chr <- QTL.mr.test[,1][select]
   pos <- QTL.mr.test[,i-2][select]
   pheno <- floor(i/3)
   ###confident interval
   int_mk1 <- c()
   int_pos1 <- c()
   int_mk2 <- c()
   int_pos2 <- c()
   n <- (i-4)/3+1
   for (c in chr){ 
      interval <- lodint(QTL.mr, c, 1, drop=1.5, lodcolumn=n, expandtomarkers=TRUE)
      int_mk1 <- append(int_mk1, paste(rownames(interval)[1],rownames(interval)[2], sep="-"))
      int_pos1 <- append(int_pos1, paste(interval$pos[1],interval$pos[2],sep="-"))
  
      interval_bayes <- bayesint(QTL.mr, c, 1, prob=0.95, lodcolumn=n, expandtomarkers=TRUE)
      int_mk2 <- append(int_mk2, paste(rownames(interval_bayes)[1],rownames(interval_bayes)[2], sep="-"))
      int_pos2 <- append(int_pos2, paste(interval_bayes$pos[1],interval_bayes$pos[2],sep="-"))

   }

   if (length (chr) < 1){
      next
   }else if(length (chr) == 1){
      qtl <- makeqtl(cross.imp.jit,chr=chr,pos=pos)
      fit <- summary(fitqtl(cross.imp.jit,pheno.col=pheno,qtl=qtl,get.ests=TRUE))
      fitest <- fit$ests[2:(length(chr)+1)]
      fitse  <- fit$ests[(length(chr)+1+2):((length(chr)+1)*2)]
      fitlod <- fit$result.full[10]
      fitvar <- fit$result.full[13]
      fitpval <- fit$result.full[16]
      
   }else if (length (chr) > 1){
      qtl <- makeqtl(cross.imp.jit,chr=chr,pos=pos)
      fit <- summary(fitqtl(cross.imp.jit,pheno.col=pheno,qtl=qtl,get.ests=TRUE))
      fitest <- fit$ests[2:(length(chr)+1)]
      fitse  <- fit$ests[(length(chr)+1+2):((length(chr)+1)*2)]
      fitlod <- fit$result.drop[(length(chr)*2+1):(length(chr)*2+length(chr))]
      fitvar <- fit$result.drop[(length(chr)*3+1):(length(chr)*3+length(chr))]
      fitpval <- fit$result.drop[(length(chr)*5+1):(length(chr)*5+length(chr))]
   }
   fitsum <- cbind(chr,pos,int_pos1,int_mk1,int_pos2,int_mk2,fitlod,fitvar,fitpval,fitest,fitse)
   print (fitsum)
   write.table(fitsum,"MPR.cross.uniq.QTL.fit.summary",append=TRUE,row.names=FALSE, col.names=TRUE,sep="	") 
}
}### QTL effect switch

## 2 write result into table and pdf
## 2.1 write QTL loci and LOD into table
write.table(QTL.mr,"MPR.cross.uniq.QTL.mr.table",sep="	",col.names=NA)
write.table(LOD,"MPR.cross.uniq.QTL.mr.table.LOD_threshold",sep="	",col.names=NA)
write.table(QTL.mr.test,"MPR.cross.uniq.QTL.mr.table.test",sep="	",col.names=NA)
write.table(QTL.em,"MPR.cross.uniq.QTL.em.table",sep="	",col.names=NA)
if(0){
write.table(QTL.hk,"MPR.cross.uniq.QTL.hk.table",sep="	",col.names=NA)
write.table(QTL.ehk,"MPR.cross.uniq.QTL.ehk.table",sep="	",col.names=NA)
write.table(QTL.imp,"MPR.cross.uniq.QTL.imp.table",sep="	",col.names=NA)
}
#write.table(QTL.cim.20,"MPR.cross.uniq.QTL.cim.20",sep="	",col.names=NA)
#write.table(QTL.cim.20.test,"MPR.cross.uniq.QTL.cim.20.test",sep="	",col.names=NA)
#write.table(qtlsummary,"MPR.cross.uniq.QTL.cim.20.summary",sep="	",col.names=NA)

## 2.2 write QTL curves into pdf
#pdf ("MPR.cross.uniq.QTL.pdf")
for (trait in 1:nph){
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



