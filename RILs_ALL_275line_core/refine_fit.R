#/rhome/cjinfeng/software/tools/R-2.15.3/bin/R
library(qtl)
nph <- 18 ## number of phenotype
read.cross("qtlcart",dir="./",file="MPR.cross.uniq.cro",mapfile="MPR.cross.uniq.map") -> cross
myMaq <- pull.map(cross)
QTL.mr <- scanone(cross,pheno.col=1:nph,method="mr")
cross <- calc.genoprob(cross,step=1,error.prob=0.001)

#fit grain yield per tiller
cross.imp <- sim.geno(cross,step=1,n.draws=64,error.prob=0.001)
cross.imp.jit <- cross.imp
qtl <- makeqtl(cross.imp.jit, chr=c(1,5,6,10), pos=c(92.32, 126.38, 6.07, 4.33))
out.fq <- fitqtl(cross.imp.jit, qtl=qtl, pheno.col=15, get.ests=TRUE, formula=y~Q1+Q2+Q3+Q4)
summary(out.fq)
#interval
rqtl <- refineqtl(cross.imp.jit, qtl=qtl, pheno.col=15, formula=y~Q1+Q2+Q3+Q4)
lodint(rqtl, qtl.index=3, drop=1.5, expandtomarkers=TRUE)


#refine plant height
qtl <- makeqtl(cross.imp.jit, chr=c(1,1,1,3,6,10), pos=c(2.95, 93.85, 143, 50.90, 53.43 ,16.99))
rqtl <- refineqtl(cross.imp.jit, qtl=qtl, pheno.col=2, formula=y~Q1+Q2+Q3+Q4+Q5+Q6)
out.fq <- fitqtl(cross.imp.jit, qtl=rqtl, pheno.col=2, get.ests=TRUE, formula=y~Q1+Q2+Q3+Q4+Q5+Q6)
summary(out.fq)
#interval
lodint(rqtl, qtl.index=1, drop=1.5, expandtomarkers=TRUE)
lodint(rqtl, qtl.index=2, drop=1.5, expandtomarkers=TRUE)
lodint(rqtl, qtl.index=3, drop=1.5, expandtomarkers=TRUE)

