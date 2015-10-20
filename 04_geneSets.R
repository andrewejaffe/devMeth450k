load("rdas/ageChanges_output.rda")

library(bumphunter)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
theTranscripts = annotateTranscripts(
	TxDb.Hsapiens.UCSC.hg19.knownGene,codingOnly=TRUE)

### go
nullgenes = ss(as.character(theTranscripts$Refseq), " ")[an$subjectHits]
nullgenes = unique(nullgenes[abs(an$dist) < 5000])


xx=load("rdas/annotated_dmr_table.rda")
dmrTable = dmrTable[dmrTable$fwer < 0.05,] # FWER filter

g = unique(dmrTable$annotation[!(dmrTable$description %in%
	c("upstream","downstream") & 	dmrTable$distance > 5000 )])
g = g[!is.na(g)]

goDmrs = dogo(g, nullgenes[,2],goP=1)

dmrsUp = dmrTable[sign(dmrTable$value) < 0,]
dmrsDown = dmrTable[sign(dmrTable$value) > 0,]
gUp = unique(dmrsUp$annotation[!(dmrsUp$description %in%
	c("upstream","downstream") & 	dmrsUp$distance > 5000)]) 
gDown = unique(dmrsDown$annotation[!(dmrsDown$description %in%
	c("upstream","downstream") & 	dmrsDown$distance > 5000)]) 
gUp = gUp[!is.na(gUp)] ; gDown = gDown[!is.na(gDown)]
goDmrsUp = dogo(gUp, nullgenes[,2],goP=1)
goDmrsDown = dogo(gDown, nullgenes[,2],goP=1)

goDmrs$bonf = p.adjust(goDmrs$Pvalue, "bonf")
goDmrsUp$bonf = p.adjust(goDmrsUp$Pvalue, "bonf")
goDmrsDown$bonf = p.adjust(goDmrsDown$Pvalue, "bonf")

save(goDmrs, goDmrsUp, goDmrsDown, file="rdas/GO_ageDmrs_output.rda")

## GO results ##
goOut = merge(goDmrsUp[,-c(4,8)], goDmrsDown[,-c(4,8)], all=TRUE,
	by = c("GOBPID", "Term", "Size"),suffixes=c("_up", "_down"))
goOut = goOut[order(goOut$Pvalue_up),]
rownames(goOut) = NULL
goOut$Pvalue_down[is.na(goOut$Pvalue_down)] = 1
goOut$Pvalue_up[is.na(goOut$Pvalue_up)] = 1
write.table(goOut, file="tables/GO_ageChanges.txt", 
	row.names=FALSE, sep="\t")

goOutUp = goOut[which(goOut$bonf_up < 0.05),]
goOutDown = goOut[which(goOut$bonf_down < 0.05),]
goOutDown = goOutDown[order(goOutDown$bonf_down),]

gBlocks = unique(ss(as.character(theTranscripts$Refseq), 
	" ")[subjectHits(ooBlocks)])

goBlocks = dogo(gBlocks, nullgenes[,2],goP=1)
goBlocks$bonf = p.adjust(goBlocks$Pvalue, "bonf")
sum(goBlocks$bonf < 0.05)
save(goBlocks, file="rdas/GO_block_genes.rda")
write.csv(goBlocks[,-8], file= "tables/suppTable_GO_blocks.csv",
	row.names=FALSE)

###########

##################
# pgc overlap
pgc = read.delim("pgc2_loci.txt",as.is=TRUE)
pgc$chr = ss(pgc$Position..hg19.,":")
tmp = ss(pgc$Position..hg19.,":",2)
pgc$start = as.numeric(ss(tmp, "-"))
pgc$end = as.numeric(ss(tmp, "-",2))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))

## dmps
ooPgc = findOverlaps(pgcGR, rowRanges(Mset))
inPgc = rep(FALSE, nrow(dmpTable))
inPgc[subjectHits(ooPgc)] = TRUE
sum(inPgc)

###################################
# enrichment of location of DMPs ##
tab = table(inPgc, dmpTable$bonf < 0.05)
tab[1,1]/tab[2,1]/tab[1,2]*tab[2,2]
rowSums(tab)
chisq.test(tab)


## plots
pdf("plots/dmp_by_pgc.pdf",h=5,w=7)
par(mar=c(5,6,1,1))
plot(density(dmpTable$meanDiff[!inPgc & dmpTable$bonf < 0.05]),
	lwd= 3, cex.axis=1.8, xlim=c(-.5,.5), cex.lab=1.8,
	xlab="Mean Meth Diff (Fetal-Postnatal)",main="")
lines(density(dmpTable$meanDiff[inPgc & dmpTable$bonf < 0.05]),
	col="red",lwd=3)
legend("topleft", c("PGC","NotPGC"),
	col=c("red","black"), pch=15,cex=1.6)
legend("topright", paste0("p=", signif(pv, 3)),cex=1.6)
dev.off()

## later AoO? ####
tab2 = table(inPgc, outAoo$bonf < 0.05)
tab2[1,1]/tab2[2,1]/tab2[1,2]*tab2[2,2]
chisq.test(tab2)


t.test(outAoo$meanDiffComp ~ inPgc, 
	subset=outAoo$bonfComp < 0.05)
pv2 = t.test(outAoo$meanDiff ~ inPgc, 
	subset=outAoo$bonf < 0.05)$p.value
t.test(abs(outAoo$meanDiff) ~ inPgc, 
	subset=outAoo$bonf < 0.05)
	
pdf("plots/dmp_by_pgc_aooPeriod.pdf",h=5,w=7)
par(mar=c(5,6,1,1))
plot(density(outAoo$meanDiff[!inPgc & outAoo$bonf < 0.05]),
	lwd= 3, cex.axis=1.8, xlim=c(-.2,.2), cex.lab=1.8,
	xlab="Mean Meth Diff (After vs Before AOO, 25yr)",main="")
lines(density(outAoo$meanDiff[inPgc & outAoo$bonf < 0.05]),
	col="red",lwd=3)
legend("topleft", c("PGC","NotPGC"),
	col=c("red","black"), pch=15,cex=1.6)
legend("topright", paste0("p=", signif(pv2, 3)),cex=1.6)
dev.off()
	

## dmrs
dmrsGR = makeGRangesFromDataFrame(dmrTable, keep=TRUE)
sig=dmrsGR[dmrsGR$fwer < 0.05]
ooPgc2 = findOverlaps(pgcGR, sig)
inPgc2 = rep(FALSE, nrow(dmrTable))
inPgc2[subjectHits(ooPgc2)] = TRUE

# location
quantile(width(sig))

## cluster group
tmp = split(rowRanges(Mset)[!is.na(dmrs$fitted)],
	cl[!is.na(dmrs$fitted)])
clusterGR = unlist(range(tmp))
clusterGR$numProbes = elementLengths(tmp)

table(countOverlaps(pgcGR, sig) > 0, 
	countOverlaps(pgcGR, clusterGR[! 
		(names(clusterGR) %in% sig$cluster)]) > 0,
	dnn = c("DMR", "Cluster"))

tabDmr = matrix(c(sum(inPgc2), sum(countOverlaps(clusterGR, 
	pgcGR) > 0) - sum(inPgc2), sum(!inPgc2), 
		length(clusterGR) - sum(countOverlaps(clusterGR, 
			pgcGR) > 0) -  sum(!inPgc2)), nr=2,byrow=TRUE)
tabDmr[1,1]/tabDmr[2,1]/tabDmr[1,2]*tabDmr[2,2]
chisq.test(tabDmr)

			
## permutation-based
overlap = countOverlaps(pgcGR, sig) > 0
numOverlap = sum(overlap)
numOverlap
pgc$numDMRs = countOverlaps(pgcGR, sig)
write.csv(pgc[pgc$numDMRs > 0, c(1:3,12)],
	"tables/suppTable_pgcRegionsNumDmrs.csv",
	row.names=FALSE)

## check confounding, other way
overlap = countOverlaps(sig, pgcGR) > 0	
t.test(width(sig) ~ overlap)
t.test(dmrTable$clusterL[dmrTable$fwer < 0.05] ~ overlap)

# null 
xx = load("/home/epi/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/nullRegions_otherDx_B100000.rda")
nullCounts = countOverlaps(nullRegions, sig) > 0
nullByPerm = split(nullCounts, nullRegions$permutation)
nullOverlap = sapply(nullByPerm, sum)
mean(nullOverlap > numOverlap)

pdf("plots/dmr_overlap_pgc.pdf", h=5)
par(mar=c(5,6,2,2))
hist(nullOverlap, col="grey",main="",
	xlab="# PGC Regions With 1+ DMRs",
	xlim = c(0,35),cex.axis=2,cex.lab=2)
abline(v=numOverlap, lwd=3, col="black")
legend("top", "p<1e-5",cex=1.6, bty="n")
dev.off()

## check confounding
overlapNull = countOverlaps(nullRegions, rowRanges(Mset)[!is.na(dmrs$fitted)])
overlapNull = split(overlapNull, nullRegions$permutation)
overlapNull = sapply(overlapNull, sum)
N = sum(countOverlaps(pgcGR, rowRanges(Mset)[!is.na(dmrs$fitted)]))

long = overlapNull > N
t.test(nullOverlap ~ long)

##### blocks
overlapBlock = countOverlaps(pgcGR, fetalBlockGR) > 0
numOverlapBlock = sum(overlapBlock)

nullCountsBlock = countOverlaps(nullRegions, fetalBlockGR) > 0
nullByPermBlock = split(nullCountsBlock, nullRegions$permutation)
nullOverlapBlock = sapply(nullByPermBlock, sum)
mean(nullOverlapBlock > numOverlapBlock)
