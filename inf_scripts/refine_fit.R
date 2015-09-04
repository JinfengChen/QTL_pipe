library(qtl)
nph <- 18 ## number of phenotype
read.cross("qtlcart",dir="./",file="MPR.cross.uniq.cro",mapfile="MPR.cross.uniq.map") -> cross
myMaq <- pull.map(cross)
QTL.mr <- scanone(cross,pheno.col=1:nph,method="mr")
cross <- calc.genoprob(cross,step=1,error.prob=0.001)

#fit grain yield per tiller
cross.imp <- sim.geno(cross,step=1,n.draws=64,error.prob=0.001)
cross.imp.jit <- cross.imp
qtl <- makeqtl(cross.imp.jit, chr=c(1,5,6,10), pos=c(91.354973, 114.584362, 5.761322, 17.780508))
out.fq <- fitqtl(cross.imp.jit, qtl=qtl, pheno.col=15, get.ests=TRUE, formula=y~Q1+Q2+Q3+Q4)
summary(out.fq)

#refine plant height
qtl <- makeqtl(cross.imp.jit, chr=c(1,1,3,10), pos=c(2.925417, 145.324312, 52.171131, 17.780508))
rqtl <- refineqtl(cross.imp.jit, qtl=qtl, pheno.col=2, formula=y~Q1+Q2+Q3+Q4)
out.fq <- fitqtl(cross.imp.jit, qtl=rqtl, pheno.col=2, get.ests=TRUE, formula=y~Q1+Q2+Q3+Q4)
summary(out.fq)


