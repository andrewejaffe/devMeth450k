##### 
# preprocessing for ECD talk

source("00_devMeth450k_functions.R")

library(minfi)

##### PUBLIC DATA
## stem cell
## 
GEOquery::getGEOSuppFiles("GSE38216")
x = read.delim("GSE38216/GSE38216_signal_intensities.txt.gz",
	as.is=TRUE,skip=4,header=TRUE,	row.names=1)
M1 = x[,grep("Methylated", colnames(x))]
U1 = x[,grep("Unmethylated", colnames(x))]

## keep all cell types
pheno1 = c("ES","ES","NPC","NPC","DA_NEURON","DA_NEURON")
colnames(M1) = colnames(U1) = paste0(pheno1, "_", seq(along=pheno1))

### add public DLPFC
library(FlowSorted.DLPFC.450k)
mset = mapToGenome(FlowSorted.DLPFC.450k, merge=TRUE)
M2 = getMeth(mset)
U2 = getUnmeth(mset)

# put public data in order
M = cbind(M1[rownames(M2),], M2) 
U = cbind(U1[rownames(U2),], U2)

# make methylset
pheno = data.frame(CellType=c(pheno1, pData(FlowSorted.DLPFC.450k)$CellType),
	stringsAsFactors=FALSE)
rownames(pheno) = colnames(M)
Mset = GenomicMethylSet(gr=granges(mset),pData=DataFrame(pheno),
	Meth = as.matrix(M), Unmeth = as.matrix(U),
	annotation=annotation(mset),
	preprocessMethod = preprocessMethod(mset))
MsetNorm = preprocessQuantile(Mset, merge=TRUE)
compData = minfi:::pickCompProbes(MsetNorm)
coefs = compData$coefEsts
coefs = coefs[!duplicated(rownames(coefs)),] # some are dups
save(coefs, file="rdas/cellComp_estimates_cellLines_NeuNs.rda")

######### PREPROCESS LIBD SAMPLES##########
# read in phenotype data from GEO objects for reproducibility
load("/dcs01/ajaffe/Brain/DNAm/ECD2014/devMeth450k_Mset_RGset.rda")
pd = as.data.frame(pData(Mset))

path = "/dcs01/ajaffe/Brain/DNAm/ECD2014/idats/" # change for GEO download
pd$BasePath=paste0(path, pd$Chip)

# add age groups
pd$ageGroup = cut(pd$Age, breaks = c(-0.5,0,0.6,10,20,50,100))
levels(pd$ageGroup) = c("Fetal","Infant","Child","Teens","Adult","50+")
	
# read in data	
RGset = read.450k(pd$BasePath)

# add control genes
controlProbes = minfi:::.extractFromRGSet450k(RGset)
negControlPCs = prcomp(t(log2(rbind(controlProbes$greenControls$NEGATIVE,
	controlProbes$redControls$NEGATIVE)+1)))$x[,1:4]
colnames(negControlPCs) = paste0("negControl_", 	
	colnames(negControlPCs))
# pd = cbind(pd, negControlPCs)

## qc, to drop replicates
mset = mapToGenome(RGset)
qc = minfiQC(mset)
qcMat = as.data.frame(qc$qc)
# pd = cbind(pd,qcMat)

### KEEP DUPLICATES IN FOR NOW BUT FLAG
# distance to centroid to flag duplicates
bIndexes = split0(pd$BrNum)
means = colMeans(pd[,c("mMed","uMed")])
bestQC= sapply(bIndexes, function(x) {
	if(length(x) > 1) {
		tmp = rowSums((pd[x,c("mMed","uMed")] - means)^2) 
	}	else tmp = 1
	x[which.min(tmp)]
})
pd$bestQC = rep(FALSE)
pd$bestQC[bestQC] = TRUE

# via FlowSorted.DLPFC.450k
Mset = preprocessQuantile(RGset, merge=TRUE)
counts = minfi:::projectCellType(getBeta(Mset[rownames(coefs), ]), coefs)
# pd = cbind(pd,counts)

# add phenotype data to Mset
pData(Mset) = DataFrame(pd)

save(Mset, RGset, 
	file="/dcs01/ajaffe/Brain/DNAm/ECD2014/devMeth450k_Mset_RGset.rda")
