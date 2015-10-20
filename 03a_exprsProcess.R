#
source("00_devMeth450k_functions.R")

library(GEOquery)
library(limma)
library(GenomicRanges)

# to get phenotype and processed data
tmp = getGEOSuppFiles("GSE30272")
system("tar xvf GSE30272/GSE30272_RAW.tar")
system("gunzip *.gz")
system("mv *.txt GSE30272/")
fn = list.files("GSE30272", pattern="GSM",full.names=TRUE)

map = Table(getGEO("GPL4611"))

cols = list(R="SR_Mean", G = "SG_Mean", Rb = "SR_bkMean", Gb = "SG_bkMean")
RGset = read.maimages(fn,columns=cols,
	annotation =c("CloneID","Flag","CloneTitle","PlatePos",	"HKFlag",
		"GelBand",	"Cytoband","geneMap"))

# RGbg = backgroundCorrect(RGset, method="normexp", offset=50)
MA = normalizeWithinArrays(RGset,method="loess", 
	bc.method = "normexp", offset=50)

p = MA$M
rownames(p) = MA$genes$CloneID

## get phenotype
theData = getGEO("GSE30272")

pdList = lapply(theData,pData)
pd = do.call("rbind",pdList)
pd = pd[,c(1,2,grep("characteristics", names(pd)))]
pd = pd[,1:10]
pd = parsePheno(pd)

colnames(p) = rownames(pd) = pd$title
p = p[map$OligoID,] # match ordering

##### use existing data from GEMMA database
anno = read.delim("GPL4611_noParents.an.txt", 
	skip=7,as.is=TRUE,row.names=1,na.string="")
map$Symbol = anno$GeneSymbols[match(map$ID, rownames(anno))]
map$Entrez = anno$NCBIids[match(map$ID, rownames(anno))]
map$GemmaID = anno$GemmaIDs[match(map$ID, rownames(anno))]

## and align
align = read.delim("GPL4611.allalignments.txt", as.is=TRUE) 
map$aligned = ifelse(map$ID %in% align$probe, 1, 0)

## drop those that don't align
keepIndex=which(map$aligned  == 1)
p = p[keepIndex,] ; map = map[keepIndex,]

## drop those that don't match to genes
dropIndex= which(is.na(map$Symbol) | grepl("\\|", map$Symbol))
p = p[-dropIndex,] ; map = map[-dropIndex,]

## order by age
oo = order(pd$Age)
p = p[,oo] ; pd = pd[oo,]

### add additional annotation
pheno = read.csv("pheno_methToExprs.csv")
pd$FileName = pheno$FileName[match(pd$title, pheno$ArrayID)]

pdGeo = pd
colnames(p) = rownames(pdGeo) = pdGeo$FileName

## add biomaRt data
library(biomaRt)
ensembl = useMart("ENSEMBL_MART_ENSEMBL", # VERSION 75, hg19
	dataset="hsapiens_gene_ensembl",
	host="feb2014.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene",
	"chromosome_name", "start_position", "end_position"), 
	mart=ensembl)
sym = sym[sym$chromosome_name %in% c(1:22,"X","Y","MT"),]
mm = match(map$Entrez, sym$entrezgene)

map$geneChr = paste0("chr", sym$chromosome_name[mm])	
map$geneStart = sym$start_position[mm]
map$geneEnd = sym$end_position[mm]

dropIndex= which(is.na(map$geneStart))
p = p[-dropIndex,] ; map = map[-dropIndex,]

save(p,map,pdGeo, file="rdas/brain_data_from_geo.rda")
