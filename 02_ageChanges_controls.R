###
# R-devel
## qsub -V -pe local 8 -l mf=40G,h_vmem=70G,h_stack=256M -cwd -b y R-devel CMD BATCH --no-save ageChanges_controls.R

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
	pData(Mset)$Dx == "Control")
Mset = Mset[,keepIndex]
	
## extract data on remaining people
pd = as.data.frame(pData(Mset))
p = getBeta(Mset)
map = as.data.frame(rowData(Mset))
colnames(map)[1] = "chr"

## composition
counts = pd[,c("ES", "NPC","DA_NEURON","NeuN_pos","NeuN_neg")]
pdf("plots/Figure2abcd.pdf",w=5,h=3.5,
	useDingbats=FALSE)
for(i in 1:ncol(counts)) {
	par(bty="n")
	agePlotter(counts[,i], pd$Age, mainText=colnames(counts)[i],
		smoothIt = FALSE, 	ylab="Cell type proportion", ylim = c(0,1),ageLabel="top",
		ageBreaks = c(-1, -0.25, 0.5, 10, 100),jitter=FALSE)
}
dev.off()

## corr to age
# corList= apply(counts, 2, function(x) cor.test(x, pd$Age, method="spearman"))
tList= apply(counts, 2, function(x) t.test(x ~ pd$Age > 0))
signif(sapply(tList, function(x) x$p.value),3)
       # ES       NPC DA_NEURON  NeuN_pos  NeuN_neg
 # 1.21e-24  5.74e-25  9.90e-01  3.57e-24  8.09e-86

## negative control probes
pdf("plots/jaffe_suppFigure10a_RUV.pdf")
par(mar=c(5,6,2,2))
palette(brewer.pal(8, "Dark2"))
boxplot(negControl_PC1 ~ Plate, data=pd, col = 1:4,
	ylab="NegControl PC1", names = paste0("Plate", 1:4),
	cex.lab=1.8,cex.axis=1.8)
boxplot(negControl_PC2 ~ Plate, data=pd, col = 1:4,
	ylab="NegControl PC2", names = paste0("Plate", 1:4),
	cex.lab=1.8,cex.axis=1.8)
dev.off()

pdf("plots/jaffe_suppFigure10b_RUV.pdf",w=12)
par(mar=c(7,3,2,2))
tt = table(pd$Sentrix_ID, pd$Plate)
tt = apply(tt, 1, which.max)
boxplot(negControl_PC1 ~ Sentrix_ID, data=pd, 
	col = tt,	ylab="",las=2, 	cex.axis=1)
boxplot(negControl_PC2 ~ Sentrix_ID, data=pd, 
	col =tt, ylab="",las=2, 	cex.axis=1)
dev.off()

############
## differences associated with birth
############
pd$Fetal = ifelse(pd$Age < 0, 1 ,0)
mod = model.matrix(~ Fetal + negControl_PC1 + negControl_PC2 + 
	negControl_PC3 + negControl_PC4,data= pd)

## pca
oo = order(matrixStats::rowSds(p),decreasing=TRUE)[1:100000]
pca = prcomp(t(p[oo,]))
pcaVars = getPcaVars(pca)

pdf("plots/jaffe_suppFigure1_PCA.pdf")
par(mar=c(5,6,2,2))
palette(brewer.pal(8,"Set1"))
fetal3 = cut(pd$Age, c(-1,0,10,100), label=c("Fetal","Child","Adult"))
fetal3 = relevel(fetal3,ref = "Adult")
plot(pca$x, cex=1.3, pch = 21, bg = as.numeric(fetal3),
	xlab=paste0("PC1: ", pcaVars[1], "% of Var Explained"),
	ylab=paste0("PC2: ", pcaVars[2], "% of Var Explained"),
	cex.axis=1.8, cex.lab=1.8)
legend("bottomleft", levels(fetal3), col = 1:3,
	pch = 15, cex=1.8)
dev.off()

################
# single CpGs, aka DMPs
# fit = lmFit(p, mod)
# eb = ebayes(fit)
# dmpTable = data.frame(meanDiff = fit$coef[,2], 
	# tstat = eb$t[,2], pval = eb$p[,2], row.names=rownames(eb$t))
# dmpTable$qval = p.adjust(dmpTable$pval, "fdr")
# dmpTable$bonf = p.adjust(dmpTable$pval, "bonferroni")
# dmpTable$meanPost = rowMeans(p[,pd$Fetal == 0])

## DMRs
# require(doParallel)
# registerDoParallel(cores = 8)
# dmrs = bumphunterEngine(p, mod, chr = map$chr, 
	# pos = map$start, cutoff=  0.1, nullMethod = "bootstrap",
	# smooth=TRUE, B=1000)

## blocks
cobj=cpgCollapse(Mset, what="Beta")
# blocks = blockFinder(cobj$object, mod,  cutoff=  0.1, 
	# nullMethod = "bootstrap",	smooth=TRUE, B=1000) 
# save(dmpTable, dmrs, blocks, file="rdas/ageChanges_output.rda")


###############################
### plots and summaries #######

load("rdas/ageChanges_output.rda")

library(bumphunter)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
theTranscripts = annotateTranscripts(
	TxDb.Hsapiens.UCSC.hg19.knownGene,codingOnly=TRUE)

## DMPs
### plot a few examples
oo = c(order(dmpTable$meanDiff, decreasing=TRUE)[1:5],
	order(dmpTable$meanDiff, decreasing=FALSE)[1:5])
pdf("plots/Figure1_insets.pdf",h=4,w=4)
par(mar=c(3,5,2,1))
palette(brewer.pal(8,"Set1"))
for(i in oo) {
	boxplot(p[i,] ~ pd$Fetal,ylim=c(0,1),
		outline=FALSE,cex.axis=1.4,cex.lab=2,cex.main=1.5,
		names = c("Postnatal", "Fetal"), yaxt="n",
		ylab="DNAm Level",main=rownames(p)[i])

	axis(2, at=seq(0,1,0.25), c("0","0.25", 
		"0.5", "0.75","1"),cex.axis=2)
	points(p[i,] ~ jitter(pd$Fetal+1, amount=0.15), 
		pch = 21, bg=pd$Fetal+1)
}
dev.off()
	
sum(dmpTable$bonf < 0.05)
signif(mean(dmpTable$bonf < 0.05)*100,4)

an = annotateNearest(rowData(Mset), theTranscripts)
dmpTable$nearestGene = as.character(theTranscripts$Gene)[an$subjectHits]
dmpTable$nearestGeneDist = an$dist

length(unique(dmpTable$nearestGene[
	abs(dmpTable$nearestGeneDist) < 5000 & dmpTable$bonf < 0.05]))
length(unique(as.character(theTranscripts$Gene)))

sum(dmpTable$pval < 0.05)
signif(mean(dmpTable$pval < 0.05)*100,4)

pdf("plots/Figure1a.pdf")
par(mar=c(5,6,2,2))
hist(dmpTable$meanDiff, col = "grey", breaks=50,
	xlab="", cex.axis=2, cex.lab=2,main="")
dev.off()

# composition R2
modComp = model.matrix(~NeuN_pos + NeuN_neg + NPC + ES + DA_NEURON,data=pd)
fitComp = lmFit(p, modComp)
R2_Comp = getR2(p, modComp)

pdf("plots/Figure2E.pdf")
par(mar=c(5,6,2,1))
hist(R2_Comp$Adjusted_R2, breaks=50,col="grey",main="",
	xlab="",cex.axis=1.8,cex.lab=1.8, ylim = c(0,22000))
hist(R2_Comp$Adjusted_R2[dmpTable$bonf < 0.05], 
	breaks=50,col="red",main="",add=TRUE,
	xlab="",cex.axis=1.8,cex.lab=1.8, ylim = c(0,22000))
dev.off()

#### DMR table
# 10% cutoff
dmrTable = dmrs$table
dmrTable = dmrTable[dmrTable$fwer < 0.05,] # FWER filter
table(sign(dmrTable$value))

cl = clusterMaker(as.character(map$chr), map$start,
	maxGap=500)
fetal2 = factor(ifelse(pd$Age < 0, "Fetal", "Postnatal"),
	levels = c("Postnatal", "Fetal"))
genes = matchGenes(dmrTable[1:50,], theTranscripts)

### make GenomicState object from derfinder package
load("rdas/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome

pdf("plots/Figure1B_SuppFigure2.pdf",w=9,useDingbats=FALSE)
par(bty="n")
dmrPlot(dmrTable[1:10,], p, map$chr, map$start, cl, genes = genes,
	genomicState=gs, coi = fetal2, cols = "Set1", Jitter=TRUE)
dev.off()

#####################
## blocks ###

## plot
blockPlot(cset=cobj$object, blocks450 = blocks, 
	coi =fetal2,blockname = "Fetal", scale=100, N=50,
	filename="plots/Figure1C_pg21.pdf",showCancerPanel=FALSE, 
	showDiffPanel=FALSE, bty="n")

### analysis on blocks
blockTable= blocks$table	
sum(blockTable$fwer < 0.05)
fetalBlockGR = with(blockTable, GRanges(chr, IRanges(start,end),
	fwer = fwer, value = value))
fetalBlockGR = fetalBlockGR[fetalBlockGR$fwer < 0.05]
quantile(width(fetalBlockGR))

ooBlocks = findOverlaps(fetalBlockGR, theTranscripts)


### overlap to cancer?
blockFn = "rdas/cancer_blocks_hansen.rda"
cancerBlocks = sapply(blockFn, function(x) {
	xx = load(x)
	get(xx)
})[[1]]
genome(cancerBlocks)="hg19"
cancerBlocks$sign = ifelse(cancerBlocks$direction =="hypo", -1,1)
mm = annotateNearest(fetalBlockGR, cancerBlocks)

## p-value
clusterMap = split(rowRanges(Mset), cobj$blockInfo$pns)
clusterMap = unlist(range(clusterMap))

cBrain = countOverlaps(clusterMap, fetalBlockGR)
cCancer = countOverlaps(clusterMap, cancerBlocks)
tab = table(cBrain > 0, cCancer > 0, dnn = c("Brain", "Cancer"))
tab[1,1]/tab[2,1]/tab[1,2]*tab[2,2]
chisq.test(tab)$p.value

mean(mm$dist==0)
prop.table(table(mm$dist==0 , 
	cancerBlocks$sign[mm$subjectHits] == sign(fetalBlockGR$value),
	dnn = c("Overlap", "Sign")), 1)

### check differences b/w 10-25 and 25-85, 
####	e.g. before and after age of onset
Index2 = which(pd$Age > 10)
aoo = ifelse(pd$Age[Index2] > 25, 1, 0)
mod2 = model.matrix(~ aoo + negControl_PC1 + negControl_PC2 + 
	negControl_PC3 + negControl_PC4,data=pd[Index2,])

fit2 = lmFit(p[,Index2], mod2)
eb2 = ebayes(fit2)
outAoo = data.frame(meanDiff = fit2$coef[,2], 
	pval = eb2$p[,2])
outAoo$bonf = p.adjust(outAoo$pval, "bonf")

## adj for cell
fit2adj = lmFit(p[,Index2], cbind(mod2, 
	pd$NeuN_neg[Index2], pd$NeuN_pos[Index2]))
eb2adj = ebayes(fit2adj)
outAoo$meanDiffComp = fit2adj$coef[,2]
outAoo$pvalComp = eb2adj$p[,2]
outAoo$bonfComp = p.adjust(outAoo$pvalComp, "bonf")

save(outAoo, file="rdas/sensitivityAoo_controls.rda")
