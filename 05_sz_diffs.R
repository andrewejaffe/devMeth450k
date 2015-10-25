###
# R-devel
## qsub -V -pe local 8 -l mf=40G,h_vmem=70G,h_stack=256M -cwd -b y R-devel CMD BATCH --no-save schizo_diffs.R

library(minfi)
library(limma)
library(RColorBrewer)

source("00_devMeth450k_functions.R")

# load data
load("/dcs01/ajaffe/Brain/DNAm/ECD2014/devMeth450k_Mset_RGset.rda")

# filter probes on chrX,Y and containing snps at SBE and target CpG
Mset = addSnpInfo(Mset)
Mset = Mset[is.na(getSnpInfo(Mset)$CpG_rs) & 
	is.na(getSnpInfo(Mset)$SBE_rs),]
Mset = Mset[!seqnames(rowData(Mset)) %in% c("chrX","chrY")]

## filter people
keepIndex = which(pData(Mset)$bestQC & 
	pData(Mset)$Gender == pData(Mset)$predictedSex & 
	pData(Mset)$Age > 16)
Mset = Mset[,keepIndex]
	
## extract data on remaining people
pd = pData(Mset)
p = getBeta(Mset)
map = as.data.frame(rowData(Mset))
colnames(map)[1] = "chr"

## pca
oo = order(matrixStats::rowSds(p),decreasing=TRUE)[1:100000]
pca = prcomp(t(p[oo,]))
pcaVars = getPcaVars(pca)

pdf("plots/jaffe_suppFigure9b_SzCompBatch.pdf")
palette(c("black","red"))
par(mar=c(5,6,2,2))
plot(pca$x[,1], pd$NeuN_neg,
	xlab=paste0("PC1: ", pcaVars[1], "% of Var Explain"),
	ylab="NeuN- Composition", pch = 21, 
	bg = as.numeric(factor(pd$Dx)),
	cex.axis=2,cex.lab=2)
legend("topright", levels(factor(pd$Dx)),
	col = 1:2, pch = 15, cex=1.8)
palette(brewer.pal(8,"Dark2"))
boxplot(pca$x[,2] ~ pd$Plate, outline=FALSE, 
	ylab=paste0("PC2: ", pcaVars[2], "% of Var Explain"))
points(pca$x[,2] ~ jitter(as.numeric(factor(pd$Plate)),
	amount=0.15), pch = 21, 	bg = as.numeric(factor(pd$Dx)))
tt = table(pd$Sentrix_ID, pd$Plate)
tt = apply(tt, 1, which.max)
boxplot(pca$x[,2] ~ pd$Sentrix_ID, outline=FALSE, col=tt, 
	ylab=paste0("PC2: ", pcaVars[2], "% of Var Explain"),
	cex.axis=2,cex.lab=2,xlab="Slide",xaxt="n")
dev.off()

## check for composition differences
pdf("plots/jaffe_suppFigure9a_SzCompBatch.pdf",h=6,w=8)
palette(brewer.pal(8,"Dark2"))
par(mar=c(9.5,6.5,2,2))
boxplot(NeuN_neg ~ paste0("Plate", 
	as.numeric(factor(Plate)),":",Dx), varwidth=TRUE,
		ylab="NeuN- Composition",cex.axis=1.6,cex.lab=2,
		outline=FALSE, data=as.data.frame(pd), 
		las=3,	col = c(1,2,2,3,3,4,5,5))
boxplot(NeuN_pos ~ paste0("Plate", 
	as.numeric(factor(Plate)),":",Dx), varwidth=TRUE,
		ylab="NeuN+ Composition",cex.axis=1.6,cex.lab=2,
		outline=FALSE, data=as.data.frame(pd), 
		las=3,	col = c(1,2,2,3,3,4,5,5))
boxplot(NPC ~ paste0("Plate", 
	as.numeric(factor(Plate)),":",Dx), varwidth=TRUE,
		ylab="NPC Composition",cex.axis=1.6,cex.lab=2,
		outline=FALSE, data=as.data.frame(pd), 
		las=3,	col = c(1,2,2,3,3,4,5,5))
boxplot(ES ~ paste0("Plate", 
	as.numeric(factor(Plate)),":",Dx), varwidth=TRUE,
		ylab="ES Composition",cex.axis=1.6,cex.lab=2,
		outline=FALSE, data=as.data.frame(pd), 
		las=3,	col = c(1,2,2,3,3,4,5,5))
boxplot(DA_NEURON ~ paste0("Plate", 
	as.numeric(factor(Plate)),":",Dx), varwidth=TRUE,
		ylab="Neuron Composition",cex.axis=1.6,cex.lab=2,
		outline=FALSE, data=as.data.frame(pd), 
		las=3,	col = c(1,2,2,3,3,4,5,5))
dev.off()

## counts
counts = as.data.frame(pd[,23:27])
plateFit = apply(counts, 2, function(x) summary(lm(x ~ pd$Dx *pd$Plate)))
lapply(plateFit, function(x) signif(coef(x),3))

################
# single CpGs, aka DMPs

## 
mod = model.matrix(~ Dx + Age + Race + negControl_PC1 + negControl_PC2 + 
	negControl_PC3 + negControl_PC4,data=as.data.frame(pd))

fit = lmFit(p, mod)
eb = ebayes(fit)
dmpTableSz = data.frame(meanDiff = fit$coef[,2], 
	tstat = eb$t[,2], pval = eb$p[,2], row.names=rownames(eb$t))
dmpTableSz$qval = p.adjust(dmpTableSz$pval, "fdr")
dmpTableSz$bonf = p.adjust(dmpTableSz$pval, "bonferroni")
dmpTableSz$meanControl = rowMeans(p[,pd$Dx == "Control"])

## DMPs
sum(dmpTableSz$bonf < 0.05)
sigIndex=which(dmpTableSz$bonf < 0.05)

## stats by plate
pIndexes=splitit(pd$Plate)
pIndexes = pIndexes[2:3] # only w both SZ and control
plateStats = lapply(pIndexes, function(ii) {
	cat(".")
	f = lmFit(p[,ii], mod[ii,])
	e = ebayes(f)
	out = data.frame(meanDiff=f$coef[,2], 
		tstat = e$t[,2], pval = e$p[,2],
		row.names=rownames(e$t))
})

pdf("plots/jaffe_suppFigure11_plateDiff.pdf")
par(mar=c(5,6,2,2))
plot(plateStats[[1]]$tstat ~ 
	plateStats[[2]]$tstat, subset=sigIndex,
	ylab="Plate2", xlab="Plate3",pch=21,bg="grey",
	cex.axis=1.8,cex.lab=1.8,cex.main=1.4,
	main = "T-statistics for SZ vs Control DMPs",
	xlim=c(-4,4), ylim = c(-30,30))
abline(0,1,col="red")
dev.off()

#########
## drop a plate, final analysis!
keepIndex=which(pd$Plate != "Lieber_244")
fit2 = lmFit(p[,keepIndex], mod[keepIndex,])
eb2 = ebayes(fit2)
dmpTableSz2 = data.frame(meanDiff = fit2$coef[,2], 
	tstat = eb2$t[,2], pval = eb2$p[,2], 
	row.names=rownames(eb2$t))
dmpTableSz2$qval = p.adjust(dmpTableSz2$pval, "fdr")
dmpTableSz2$bonf = p.adjust(dmpTableSz2$pval, "bonferroni")
dmpTableSz2$meanControl = rowMeans(p[,pd$Dx == "Control"])

save(dmpTableSz2, file="rdas/dropped_plate_SZ_diffs.rda")

## metrics
sigIndex2= which(dmpTableSz2$bonf < 0.05)

quantile(abs(dmpTableSz2$meanDiff[sigIndex2]))
sum(abs(dmpTableSz2$meanDiff[sigIndex2]) > 0.03)
mean(dmpTableSz2$meanDiff[sigIndex2] < 0)

## check vs mill
dmpTableSz2[c("cg26173173","cg24803255",
	"cg00903099","cg08171022"),]
	
## annotate to genes
library(bumphunter)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
theTranscripts = annotateTranscripts(
	TxDb.Hsapiens.UCSC.hg19.knownGene,codingOnly=TRUE)

an = annotateNearest(rowData(Mset), theTranscripts)
dmpTableSz2$nearestGene = as.character(theTranscripts$Gene)[an$subjectHits]
dmpTableSz2$nearestRefseq = ss(as.character(theTranscripts$Refseq), " ")[an$subjectHits]
dmpTableSz2$nearestGeneDist = an$dist	
	
sigDmpTable = dmpTableSz2[dmpTableSz2$bonf < 0.05,]
sigDmpTable = sigDmpTable[order(sigDmpTable$pval),]
gSig = matchGenes(rowRanges(Mset)[rownames(sigDmpTable)], theTranscripts)
gSig$annotation = ss(gSig$annotation, " ")
sigDmpTable = cbind(sigDmpTable, gSig)
save(sigDmpTable, file="rdas/dmp_SZ_list_annotated.rda")
write.csv(sigDmpTable, file="tables/jaffe_suppTable13_szDmps.csv")

	
#####
# pgc overlap

## read in 108 loci
pgc = read.delim("pgc2_loci.txt",as.is=TRUE)
pgc$chr = ss(pgc$Position..hg19.,":")
tmp = ss(pgc$Position..hg19.,":",2)
pgc$start = as.numeric(ss(tmp, "-"))
pgc$end = as.numeric(ss(tmp, "-",2))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))

## find overlap
ooPgc = findOverlaps(pgcGR, rowRanges(Mset))
inPgc = rep(FALSE, nrow(dmpTableSz2))
inPgc[subjectHits(ooPgc)] = TRUE
sum(inPgc)
tab = table(inPgc, dmpTableSz2$bonf < 0.05)
tab[1,1]/tab[2,1]/tab[1,2]*tab[2,2]
rowSums(tab)
chisq.test(tab)


load("/dcs01/ajaffe/Brain/DNAm/ECD2014/rdas/meQTLs_adults_cis_annotated.rda")
sum(rownames(dmpTableSz2_inPGC) %in% sig$cpg)

## explore some hits for biology
sigDmpTableNear = sigDmpTable[abs(sigDmpTable$nearestGeneDist) < 500,]
sigDmpTableNearBig = sigDmpTableNear[abs(sigDmpTableNear$meanDiff) > 0.03,]


# do GO
load("rdas/dmp_SZ_list_annotated.rda")
nullgenes =  read.delim("ref_gene_hg19.txt",header=TRUE,as.is=TRUE)
g = sigDmpTable$annotation
g = g[!(sigDmpTable$description %in% c("upstream","downstream") &
	abs(sigDmpTable$distance) > 5000)]
g = g[!is.na(g)]
goSz = dogo(g,nullgenes[,2], goP=1)
goSz$bonf = p.adjust(goSz$Pvalue, "bonf")
write.csv(goSz[,-8], file="tables/jaffe_suppTable14_szGO.csv",
	row.names=FALSE)

## check fetal
load("rdas/ageChanges_output.rda")		
pdf("plots/jaffe_suppFigure14_fetalSzEffect.pdf",useDingbat=FALSE)
par(mar=c(5,6,1,1))
plot(dmpTableSz2$meanDiff ~ dmpTable$meanDiff,
	subset=sigIndex2, xlab="Fetal versus Adult: % in Prop Meth",
	ylab="SZ versus Adult Control: % in Prop Meth",
	pch=21,bg="grey",xlim = c(-0.3,0.3), ylim = c(-0.05,0.05),
	cex.axis=2,cex.lab=2,cex=1.5)
abline(v=0,h=0,col="black",lty=2,lwd=1.6)
dev.off()

ct = cor.test(dmpTableSz2$meanDiff[sigIndex2],
	dmpTable$meanDiff[sigIndex2])
ct$p.value

sum(abs(dmpTable$meanDiff) > 0.1)
mean(abs(dmpTable$meanDiff) > 0.1)

	genes$name[repIndex]

mean(sign(dmpTableSz2$meanDiff[sigIndex2]) != 
	sign(dmpTable$meanDiff[sigIndex2]))

table(dmpTableSz2$bonf[sigIndex2] < 0.05)
	
### check composition sensitivity model
modComp = model.matrix(~ Dx + Age + Race + negControl_PC1 + negControl_PC2 + 
	negControl_PC3 + negControl_PC4 + NeuN_neg + NeuN_pos + NPC + ES + DA_NEURON,
	data=as.data.frame(pd))
fitComp = lmFit(p[,keepIndex], modComp[keepIndex,])
ebComp = ebayes(fitComp)

pdf("plots/suppFigure_compAdj_SzDiff.pdf")
par(mar=c(5,6,4,2))
plot(fitComp$coef[sigIndex2,2], fit2$coef[sigIndex2,2],
	pch = 21, bg="grey", xlab="Model + Composition",
	ylab="Model", main ="Proportion Change in DNAm\nSZ vs Control, FDR < 5% CpGs",
	xlim = c(-0.04, 0.04), ylim = c(-0.04, 0.04),
	cex.axis=2, cex.lab=2, cex.main=1.6)
abline(h=0,v=0,lty=2)
abline(0,1,lty=2, col="blue")

dev.off()

#### check smoking ####
pheno = read.csv("DNAMeth_tox_052013.csv", as.is=TRUE)
colnames(pheno)[1] = "BrNum"
pheno$BrNum = paste0("Br", pheno$BrNum)
pheno = pheno[match(pd$BrNum, pheno$BrNum),]

pd$Smoking = ifelse(pheno$Nicotine == "Positive" | 
	pheno$Cotinine == "Positive", "Yes", "No")
sIndex=which(!is.na(pd$Smoking) & pd$Plate != "Lieber_244")
modSmoking =  model.matrix(~Smoking + Dx + Age + Race + 
	negControl_PC1 + negControl_PC2 + negControl_PC3 + 
	negControl_PC4,data=as.data.frame(pd)[sIndex,])
fitSmoking = lmFit(p[,sIndex], modSmoking)
ebSmoking = ebayes(fitSmoking)
sum(p.adjust(ebSmoking$p[,3], "fdr") < 0.05)

pdf("plots/suppFigure_smokingAdj_SzDiff.pdf")
par(mar=c(5,6,4,2))
plot(fitSmoking$coef[sigIndex2,3], fit2$coef[sigIndex2,2],
	pch = 21, bg="grey", xlab="Model + Smoking",
	ylab="Model", main ="Proportion Change in DNAm\nSZ vs Control, FDR < 5% CpGs",
	xlim = c(-0.04, 0.04), ylim = c(-0.04, 0.04),
	cex.axis=2, cex.lab=2, cex.main=1.6)
abline(h=0,v=0,lty=2)
abline(0,1,lty=2, col="blue")
dev.off()

#### check treatment

pd$Antipsychotics = ifelse(pheno$Antipsychotics == "Positive", "Yes", 
	ifelse(pheno$Antipsychotics == "Negative", "No", NA))
antiIndex=which(!is.na(pd$Antipsychotics) & pd$Plate != "Lieber_244")
modAnti =  model.matrix(~Antipsychotics + Dx + Age + Race + 
	negControl_PC1 + negControl_PC2 + negControl_PC3 + 
	negControl_PC4,data=as.data.frame(pd)[antiIndex,])
fitAnti = lmFit(p[,antiIndex], modAnti)
ebAnti = ebayes(fitAnti)
sum(p.adjust(ebAnti$p[,3], "fdr") < 0.05)

pdf("plots/suppFigure_antiAdj_SzDiff.pdf")
par(mar=c(5,6,4,2))
plot(fitAnti$coef[sigIndex2,3], fit2$coef[sigIndex2,2],
	pch = 21, bg="grey", xlab="Model + Antipsychotics",
	ylab="Model", main ="Proportion Change in DNAm\nSZ vs Control, FDR < 5% CpGs",
	xlim = c(-0.04, 0.04), ylim = c(-0.04, 0.04),
	cex.axis=2, cex.lab=2, cex.main=1.6)
abline(h=0,v=0,lty=2)
abline(0,1,lty=2, col="blue")
dev.off()

# save for carolina
antiAdjStats = data.frame(meanDiff=fitAnti$coef[,2],
	tstat = ebAnti$t[,2], pval = ebAnti$p[,2])
smokingAdjStats = data.frame(meanDiff=fitSmoking$coef[,2],
	tstat = ebSmoking$t[,2], pval = ebSmoking$p[,2])
save(antiAdjStats, smokingAdjStats,
	file="/nexsan2/disk3/Feinberg_Jaffe/Carolina/antipsychAndSmokingEffects.rda")

### age of onset DMPs?
n1 = length(intersect(which(dmpTable$bonf < 0.05),
	sigIndex2))
tab1 = matrix(c(n1, length(sigIndex2)-n1, 
	sum(dmpTable$bonf < 0.05) -n1, nrow(dmpTable)-
		length(sigIndex2) - sum(dmpTable$bonf < 0.05) + n1),
			nc = 2, byrow=TRUE)
tab1[1,1]/tab1[2,1]/tab1[1,2]*tab1[2,2]	
chisq.test(tab1)$p.value

load("rdas/sensitivityAoo_controls.rda")
aooIndex=which(outAoo$bonf < 0.05)
length(aooIndex)
n2 = length(intersect(aooIndex, sigIndex2))
tab2 = matrix(c(n2, length(sigIndex2)-n2, 
	sum(outAoo$bonf < 0.05) - n2, nrow(outAoo)-
		length(sigIndex2) - sum(outAoo$bonf < 0.05) + n2),
			nc = 2, byrow=TRUE)
tab2[1,1]/tab2[2,1]/tab2[1,2]*tab2[2,2]	
chisq.test(tab2)$p.value

#### carolina model
pi1 = pd$NeuN_neg
X = as.numeric(factor(pd$Dx))-1
t2 = (1-pi1)*X ; t3 = pi1*X

modComp2 = model.matrix(~pi1[keepIndex] + 
	t2[keepIndex] + t3[keepIndex])
fitComp2 = lmFit(p[,keepIndex], modComp2)

pdf("plots/suppFigure_carolinaAdj_SzDiff.pdf")
par(mar=c(5,6,4,2))
plot(fitComp2$coef[sigIndex2,4], fit2$coef[sigIndex2,2],
	pch = 21, bg="grey", xlab="Composition Model",
	ylab="Full Model", main ="Proportion Change in DNAm\nSZ vs Control, FDR < 5% CpGs",
	xlim = c(-0.15, 0.15), ylim = c(-0.04, 0.04),
	cex.axis=2, cex.lab=2, cex.main=1.6)
dev.off()

ebComp2 = ebayes(fitComp2)

