###
source("00_devMeth450k_functions.R")

## load composition profiles
load("rdas/cellComp_estimates_cellLines_NeuNs.rda")

# #### read in brainspan data
# p1 = read.delim("/nexsan2/disk3/ajaffe/BrainSpan/DNAm/1109_methylation_beta_values.txt",
	# as.is=TRUE, skip=14, header=TRUE,row.names=1)
# p1 = p1[cpgs,]
# pd1 = read.delim("/nexsan2/disk3/ajaffe/BrainSpan/DNAm/1109_methylation_beta_values.txt",
	# as.is=TRUE, header=FALSE,nrow = 14)

# p2 = read.delim("/nexsan2/disk3/ajaffe/BrainSpan/DNAm/1110_methylation_beta_values.txt",
	# as.is=TRUE, skip=14, header=TRUE,row.names=1)
# p2 = p2[cpgs,]
# pd2 = read.delim("/nexsan2/disk3/ajaffe/BrainSpan/DNAm/1110_methylation_beta_values.txt",
	# as.is=TRUE, header=FALSE,nrow = 14)
	
# p= cbind(p1,p2)
# colnames(p) = ss(colnames(p), "\\.")

# pd = cbind(pd1, pd2[,-1])
# pd = pd[-4,] # duplicate well
# rownames(pd) = pd[,1]
# pd = pd[,-1]
# pd = t(pd)
# pd = as.data.frame(pd, stringsAsFactors=FALSE)
# rownames(pd) = NULL
# pd$Plate = ss(pd$Plate, "\\.")
# names(pd) = gsub(" ", "", names(pd))
# pd$SampleID = gsub(" ", "_", pd$SampleID)

# identical(gsub("X","", colnames(p)), pd$CompleteBarcode) # TRUE
# colnames(p) = rownames(pd) = pd$SampleID
# save(p,pd,file="rdas/brainspan_dnam_data.rda",compress=TRUE)

# pd$Regions = ss(rownames(pd),"_",2)
# pd$Regions[pd$Regions=="AIC"] = "A1C"
# pd$Regions[pd$Regions=="MIC"] = "M1C"
# pd$Regions[pd$Regions=="VIC"] = "V1C"
# pd$Regions[pd$Regions=="SIC"] = "S1C"
# pd$Regions[pd$Regions=="STS"] = "STR"
# pd$Regions = factor(pd$Regions, levels = c("DFC","VFC","MFC",
	# "OFC","M1C","S1C", "IPC", "A1C", "STC", "ITC", "V1C", "HIP",
	# "AMY", "STR", "MD", "CBC"))
# pd$Specimen.Code = ss(pd$SampleID, "_")

# pheno = read.csv("brainspan_pheno_meth.csv",as.is=TRUE)
# tmp = as.numeric(ss(pheno$Age," ", 1))
# tmp[ss(pheno$Age," ", 2)=="M"] = tmp[ss(pheno$Age," ", 2)=="M"]/12
# pheno$Age= tmp

# pd$Age = pheno$Age[match(pd$Specimen.Code, pheno$Specimen.Code)]
# pd$Specimen.ID = pheno$Specimen.ID[match(pd$Specimen.Code, pheno$Specimen.Code)]

# countsBrainspan = minfi:::projectCellType(p[rownames(coefs), ], coefs)
# save(pd, countsBrainspan, file = "rdas/phenotype_data_brainspan.rda")
load("rdas/phenotype_data_brainspan.rda")

## plots
pdf("plots/boxplots_brainspan_byCellType.pdf",h=5,w=7)
par(mar = c(5,6,1,1))
palette(colorRampPalette(c("white", "blue"))(14))
n = c("ESCs", "NPCs", "DA Neurons", "NeuN+", "NeuN-")
for(i in 1:ncol(countsBrainspan)) {
	boxplot(countsBrainspan[,i]~Regions, data=pd, las=3,
		outline=FALSE, ylim = range(countsBrainspan[,i]),
		cex.axis=2,cex.lab=2, ylab = paste(n[i], "Proportion"))
	points(countsBrainspan[,i]~jitter(as.numeric(factor(Regions)), amount=0.1),
		data=pd, pch = 21, bg = as.numeric(factor(pd$Age)))
	abline(v = 11.5,lty=2)
	points(x=seq(2,5,len=length(unique(pd$Age))), y = rep(0.12,14),
		pch=15,col=1:14,cex=1.5)
	text(c(1.3, 6), c(0.12,0.12), c("4M", "37Y"),cex=1.5)
	text(3.5, 0.15, "Age", font=2,cex=1.4)
	rect(0.5,0.1, 6.8, 0.165)
}
dev.off()

pd$NCX = ifelse(pd$Regions %in% names(table(pd$Regions))[1:11], 0,1)
mod = model.matrix(~pd$NCX + pd$Age)
pv = t(apply(countsBrainspan[!is.na(pd$Age),], 2,
	function(x) summary(lm(x~mod-1))$coef[2:3,4]))
pv
slope = t(apply(countsBrainspan[!is.na(pd$Age),], 2,
	function(x) summary(lm(x~mod-1))$coef[2:3,1]))
slope

##################
# ## Mill data ######
# library(GEOquery)
# theData = getGEO("GSE58885")[[1]]

# pFetal = exprs(theData)
# pdFetal = pData(theData)
# pdFetal$Age = as.numeric(ss(as.character(
	# pdFetal$characteristics_ch1), ": ",2))
# map = fData(theData)
# countsFetal = minfi:::projectCellType(pFetal[rownames(coefs), ], coefs)
# save(pdFetal, countsFetal,  file="rdas/phenotype_millFetal.rda")
# save(pdFetal, countsFetal, map, pFetal, file="rdas/allData_millFetal.rda")

load("rdas/phenotype_millFetal.rda")
pdf("plots/fetalBrain_compByPcd.pdf")
par(mar=c(5,6,1,1), bty="n") 
for(i in 1:ncol(countsFetal)) {
	plot(pdFetal$Age, countsFetal[,i], 
		ylab = paste(n[i], "proportion"),
		xlab = "Days post-conception", 
		pch=21,bg="grey",cex=2,
		cex.lab=2,cex.axis=2)
	legend("topright", paste0("r = ", 
		signif(cor(pdFetal$Age, countsFetal[,i]),3)),
			cex=1.9, bty="n")
}
dev.off()

corList = apply(countsFetal, 2, cor.test, pdFetal$Age)
sapply(corList, function(x) x$p.value)

## variance explained
library(limma)
modComp = model.matrix(~NeuN_pos + NeuN_neg + NPC + ES + DA_NEURON,
	data=as.data.frame(countsFetal))
fitComp = lmFit(pFetal, modComp)
R2_Comp = getR2(pFetal, modComp)

modBoth = model.matrix(~pdFetal$Age + NeuN_pos + NeuN_neg + 
	NPC + ES + DA_NEURON,data=as.data.frame(countsFetal))
fitBoth = lmFit(pFetal, modBoth)
ebBoth = ebayes(fitBoth)
sum(p.adjust(ebBoth$p[,2],"bonf") < 0.05, na.rm=TRUE)
ii = which(p.adjust(ebBoth$p[,2],"bonf") < 0.05)
ebBoth$p[ii,]

ffBoth = getF(fitBoth, fitComp, pFetal)
ffBoth$bonf = p.adjust(ffBoth$f_pval, "bonf")

modAge = model.matrix(~pdFetal$Age)
fitAge = lmFit(pFetal, modAge)
ebAge = ebayes(fitAge)
bonfAge=  p.adjust(ebAge$p[,2], "bonf")
. 

dmps = read.csv("tables/dDMPs_spiersEtAl.csv", as.is=TRUE)
ind = which(rownames(R2_Comp) %in% dmps$Probe)
pdf("plots/comp_varExpl_DNAm_spiersEtAl.pdf")
par(mar=c(5,6,2,1))
hist(R2_Comp$Adjusted_R2, breaks=50,col="grey",main="",
	cex.axis=1.8,cex.lab=1.8, ylim = c(0,44000),
	xlab="Adjusted R2")
hist(R2_Comp$Adjusted_R2[ind], 
	breaks=50,col="red",main="",add=TRUE,
	xlab="",cex.axis=1.8,cex.lab=1.8, ylim = c(0,44000))
dev.off()
