echo "275"
cp /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_275line/MPR.cross.uniq.cro ./
cp /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_275line/MPR.cross.uniq.map ./
mydata <- read.cross("qtlcart", '', 'MPR.cross.uniq.cro', 'MPR.cross.uniq.map')
write.cross(mydata, 'csv', 'test')
python remove_phenotype.py --input test.csv > RIL275.csv
mydata <- read.cross("csv", '', 'RIL275.csv',genotypes=c("AA", "AB", "BB"), crosstype='riself')

echo "230"
cp /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_275line_core/MPR.cross.uniq.cro ./
cp /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_275line_core/MPR.cross.uniq.map ./
mydata <- read.cross("qtlcart", '', 'MPR.cross.uniq.cro', 'MPR.cross.uniq.map')
write.cross(mydata, 'csv', 'test')
python remove_phenotype.py --input test.csv > RIL230.csv
mydata <- read.cross("csv", '', 'RIL230.csv',genotypes=c("AA", "AB", "BB"), crosstype='riself')

