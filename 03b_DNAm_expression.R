######
# R-devel

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
pdM = pData(Mset)
pM = getBeta(Mset)
mapM = rowData(Mset)

## load expression
xx=load("rdas/brain_data_from_geo.rda")
pE = p ; pdE = pdGeo ; mapE = map
rm(p,pdGeo,map)

## match up
mm = match(pdM$BrNum, pdE$FileName)
pdM = pdM[!is.na(mm),] ; pM = pM[,!is.na(mm)]
pdE = pdE[mm[!is.na(mm)],] ; pE = pE[,mm[!is.na(mm)]]

## differential expression
modE = model.matrix(~ifelse(pdM$Age < 0, 1,0))
colnames(modE)[2] = "Fetal"
fitE = lmFit(pE, modE)
ebE = ebayes(fitE)
expTable = data.frame(logFC = fitE$coef[,2], 
	tstat = ebE$t[,2], pval = ebE$p[,2],
	bonf = p.adjust(ebE$p[,2],"bonferroni"),
	Symbol = mapE$Symbol)
length(unique(expTable$Symbol[expTable$bonf < 0.05]))
mean(expTable$bonf < 0.05)

# composition R2
modComp = model.matrix(~NeuN_pos + NeuN_neg +
	NPC + ES + DA_NEURON,data=as.data.frame(pdM))
fitComp = lmFit(pE, modComp)
R2_Comp = getR2(pE, modComp)

pdf("plots/suppFigure5.pdf")
par(mar=c(5,6,2,1))
hist(R2_Comp$Adjusted_R2, col="grey",main="",
	xlab="",cex.axis=1.8,cex.lab=1.8)
dev.off()

#### load results ### 
load("rdas/ageChanges_output.rda")
mapE_gr = with(mapE, GRanges(geneChr, 
	IRanges(geneStart,geneEnd)))

## DMPs
anDmp = distanceToNearest(mapM, mapE_gr)
dmpTable$nearestGene = mapE$Symbol[subjectHits(anDmp)]
dmpTable$distToNearestGene = mcols(anDmp)$distance

pEE = pE[subjectHits(anDmp),]
dmpCor = dmpCorAdult = rep(NA,nrow(dmpTable))
adultIndex=which(pdM$Age > 13)
for(i in seq(along=dmpCor)) {
	if(i %% 10000 == 0) cat(".")
	dmpCor[i] = cor(pM[i,], pEE[i,])
	dmpCorAdult[i] = cor(pM[i,adultIndex], pEE[i,adultIndex])
}
dmpTable$corrToExprs = dmpCor
dmpTable$corrToExprsTstat = dmpTable$corrToExprs/sqrt(
	(1-dmpTable$corrToExprs^2)/(ncol(pM)-2))
dmpTable$corrToExprsPvalue = 2*pt(-abs(dmpTable$corrToExprsTstat),
	df = ncol(pM)-1)
dmpTable$corrToExprsPvalueBonf = p.adjust(dmpTable$corrToExprsPvalue,
	"bonferroni")

dmpTable$corrToExprsAdult = dmpCorAdult
dmpTable$corrToExprsAdultTstat = dmpTable$corrToExprsAdult/sqrt(
	(1-dmpTable$corrToExprsAdult^2)/(ncol(pM)-2))
dmpTable$corrToExprsAdultPvalue = 2*pt(-abs(dmpTable$corrToExprsAdultTstat),
	df = ncol(pM)-1)
save(dmpTable, file="rdas/annotated_dmp_table.rda")

# how many probes are near genes and signif?
table(dmpTable$distToNearestGene < 5000,
	dmpTable$bonf < 0.05,dnn=c("near","sig"))
prop.table(table(dmpTable$distToNearestGene < 5000,
	dmpTable$bonf < 0.05,dnn=c("near","sig")),2)

# and how many have expression changes?
mean(dmpTable$corrToExprsPvalue[
	dmpTable$distToNearestGene < 5000 & dmpTable$bonf < 0.05] 
		< 0.05)	# marginal
mean(dmpTable$corrToExprsPvalueBonf[
	dmpTable$distToNearestGene < 5000 & dmpTable$bonf < 0.05] 
		< 0.05)	# bonf

## plots
expTableToDmp = expTable[subjectHits(anDmp),]
sigDmpIndex=which(dmpTable$bonf < 0.05 & 
	dmpTable$distToNearestGene < 5000)

pdf("plots/suppFigure3_dmpVsExpression.pdf")
par(mar=c(5,6,4,2))
hist(dmpTable$corrToExprs[sigDmpIndex],
	breaks = 50, col = "grey",xlab="Correlation to Expression",
	cex.axis=2,cex.lab=2,main = "", xlim = c(-1,1))
dev.off()

cor.test(dmpTable$tstat[sigDmpIndex], 
	expTableToDmp$tstat[sigDmpIndex])

write.csv(dmpTable, file=gzfile("tables/suppTable2_DMPsWithExprs.csv.gz"))

### DMRs
library(bumphunter)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
theTranscripts = annotateTranscripts(
	TxDb.Hsapiens.UCSC.hg19.knownGene,codingOnly=TRUE)

dmrTable= dmrs$table
dmrGenes = matchGenes(dmrTable, theTranscripts)
dmrGenes$annotation = ss(dmrGenes$annotation, " ")	
dmrTable = cbind(dmrTable, dmrGenes)

## match to expression
geneMatch = lapply(dmrTable$name, grep, x=mapE$Symbol)

oo = data.frame(queryHits=rep(seq(along=geneMatch),
	elementLengths(geneMatch)), subjectHits = unlist(geneMatch))
oo = oo[!is.na(oo$subjectHits),]

meanMeth = matrix(NA, nrow = nrow(dmrTable),ncol=ncol(pM))
for(i in 1:nrow(dmrTable)) {
	if(i %% 500 == 0) cat(".")
	ii = dmrTable$indexStart[i]:dmrTable$indexEnd[i]
	if(length(ii) == 1) {
		meanMeth[i,] = pM[ii,] 
	} else meanMeth[i,] = colMeans(pM[ii,])
}
# put in order
meanMeth = meanMeth[oo$queryHits,]
pEE = pE[oo$subjectHits,]

dmrCor = rep(NA,nrow(oo))
for(i in seq(along=dmrCor)) {
	if(i %% 500 == 0) cat(".")
	dmrCor[i] = cor(meanMeth[i,], pEE[i,], use="comp")
}
dmrCorList = split(dmrCor, factor(oo$queryHits,levels=1:nrow(dmrTable)))

dmrTable$corrToExprs = as.numeric(sapply(dmrCorList, 
	function(x) x[which.max(abs(x))]))
dmrTable$corrToExprsTstat = dmrTable$corrToExprs/sqrt(
	(1-dmrTable$corrToExprs^2)/(ncol(pM)-2))
dmrTable$corrToExprsPvalue = 2*pt(-abs(dmrTable$corrToExprsTstat),
	df = ncol(pM)-1)
dmrTable$corrToExprsPvalueBonf = p.adjust(dmrTable$corrToExprsPvalue,
	"bonferroni")
save(dmrTable, file="rdas/annotated_dmr_table.rda")

write.csv(dmrTable, file="tables/jaffe_suppTable3_DMRsWithExprs.csv",
	row.names=FALSE)

pdf("plots/suppFigure4_dmrsVsExpression.pdf",h=5,w=7)
par(mar=c(7,6,2,2))
boxplot(corrToExprs ~ region, data=dmrTable, 
	outline=FALSE, subset=fwer < 0.05,xaxt="n",cex.lab=1.6,
	ylim = c(-1,1),	ylab="Corr: Meth to Exprs",cex.axis=1.6)
abline(h=0,lty=2)
abline(h=c(-1,1)*0.273, col="red")
text(seq(along=levels(dmrTable$region))-0.8, -1.35, 
	labels=levels(dmrTable$region), srt=45, pos=1, xpd=TRUE,
	cex=1.8)
points(corrToExprs~jitter(as.numeric(region), amount=0.15),
	data=dmrTable,subset=fwer < 0.05,pch=21,bg="grey")
dev.off()

# how many dmrs are near genes and signif?
nearIndex=which(!(dmrTable$description %in% c("upstream","downstream") &
	abs(dmrTable$distance) > 5000) & 
	dmrTable$fwer < 0.05)
length(nearIndex)/nrow(dmrTable)	
	
mean(dmrTable$corrToExprsPvalue[nearIndex] < 0.05,na.rm=TRUE)
sum(dmrTable$corrToExprsPvalue[nearIndex] < 0.05,na.rm=TRUE)
mean(dmrTable$corrToExprsPvalueBonf[nearIndex] < 0.05,na.rm=TRUE)
sum(dmrTable$corrToExprsPvalueBonf[nearIndex] < 0.05,na.rm=TRUE)
table(dmrTable$distance < 5000 | !dmrTable$region %in% c("upstream","downstream") & 
	dmrTable$fwer < 0.05,dnn=c("near","sig"))
prop.table(table(dmrTable$distToNearestGene < 5000,
	dmrTable$fwer < 0.05,dnn=c("near","sig")),2)

# and the correlations to expression	
mean(dmrTable$corrToExprsPvalue[
	dmrTable$fwer < 0.05 & dmrTable$distToNearestGene < 5000] < 0.05,
	na.rm=TRUE)
mean(dmrTable$corrToExprsPvalueBonf[
	dmrTable$fwer < 0.05 & dmrTable$distToNearestGene < 5000] < 0.05,
	na.rm=TRUE)
	
### blocks ####
blockTable = blocks$table
blockTable = blockTable[blockTable$fwer < 0.05,]
blockGR = with(blockTable, GRanges(chr, IRanges(start,end)))
anBlock = findOverlaps(blockGR, mapE_gr,select="all")
diffExpByBlock = split(expTable[subjectHits(anBlock),],
	factor(queryHits(anBlock), levels=seq(along=blockGR)))
blockTable$numGenesContained = elementLengths(diffExpByBlock)
diffExpByBlock = lapply(diffExpByBlock, function(x) {
	x = x[order(x$pval),]
	x[!duplicated(x$Symbol),]
})
blockTable$propSignifDeGenes = sapply(diffExpByBlock, function(x) mean(x$bonf < 0.05))
blockDiffExpT = sapply(diffExpByBlock, function(x) mean(x$tstat))

## add cancer overlap
blockFn = "rdas/cancer_blocks_hansen.rda"
cancerBlocks = sapply(blockFn, function(x) {
	xx = load(x)
	get(xx)
})[[1]]
genome(cancerBlocks)="hg19"
cancerBlocks$sign = ifelse(cancerBlocks$direction =="hypo", -1,1)
mm = annotateNearest(blockGR, cancerBlocks)
blockTable$overlapCancerBlock = mm$dist==0
blockTable$cancerBlockSign = cancerBlocks$sign[mm$subjectHits]
blockTable$cancerBlockSign[!blockTable$overlapCancerBlock ] = NA
write.csv(blockTable, file="tables/jaffe_suppTable5_blocksWithExprs.csv",row.names=FALSE)
