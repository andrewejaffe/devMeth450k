###
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

# splitting variables into list
splitit = function(x) split(seq(along=x),x) # splits into list
split0 = function(x) splitit(factor(x, levels = unique(x)))

getPcaVars = function(pca)  signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100


# get R2 from matrix via limma
# mod0 must be a nested model, within mod, a la anova
getR2 = function(p,mod,mod0=NULL) {
	require(limma)
	
	fit1 = lmFit(p,mod)
	rss1 = rowSums((p-fitted(fit1))^2)
	n = ncol(p)
	k = ncol(mod)-1
	
	if(is.null(mod0)) {
		rss0 = rowSums((p-rowMeans(p))^2)
	} else {
		fit0 = lmFit(p,mod0)
		rss0 = rowSums((p-fitted(fit0))^2)
	}
	
	r2 = 1 - (rss1/rss0)
	r2adj = 1 - ((1-r2)*(n-1))/(n-k-1)
	out = data.frame(R2 = r2,Adjusted_R2 = r2adj)
	return(out)	
}


## cset: result of cpgCollapse()$cobj
## block450: results of blockFinder
## coi = covariate of interest
## N: the number of blocks to plot, default=10
## blockname = name of block track, default='coi'
## filename = where to save plots
## scale, in kb. default = 100
blockPlot = function(cset, blocks450, coi, N=10,
	blockname = "coi", filename=paste0(blockname,"_blocks.pdf"),scale=100,
	showMethPanel = TRUE, showGenePanel=TRUE, showDiffPanel=TRUE, 
	showCancerPanel = TRUE, bty= "o") {
	panels = c(showMethPanel,showGenePanel, showDiffPanel, showCancerPanel)
	
	require(GenomicRanges)

	blocksTable=with(blocks450$table, GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))
	colIds=match(c("chr","start","end"),names(blocks450$table))
	mcols(blocksTable)=blocks450$table[-colIds]

	plotRegion = blocksTable[1:N]

	## annotation based on ensembl
	cat("Loading Annotation.\n")
	load("/home/epi/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
	gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
	oo = findOverlaps(blocksTable, gs)
	anno = split(gs[subjectHits(oo)], queryHits(oo))

	# cancer blocks
	if(showCancerPanel) {
		load("/home/epi/ajaffe/Lieber/Projects/450k/devPaper/cancer_blocks_hansen.rda")
		genome(blocks)="hg19"
		cancerBlocks = blocks
	}

	cat("Ploting.")
	pdf(filename,height=5,width=10)
	par(bty=bty)

	for(i in seq(along=plotRegion)) {
		cat(".")
		r = plotRegion[i]
		tmp=subsetByOverlaps(cset,r)
		tmp450=sort(subsetByOverlaps(blocksTable,r))
		if(showCancerPanel) tmpBsmooth=subsetByOverlaps(cancerBlocks,r)
											
		beta=getBeta(tmp)
		x=start(tmp)

		ii=cset %over% r
		d=blocks450$coef[ii]
		sd=blocks450$fitted[ii]
		## which rows
		Index=split(seq_along(coi),coi)
		mns=sapply(Index,function(ind) rowMeans(beta[,ind]))  
		smns=apply(mns,2,function(y) limma::loessFit(y,x,span=.2)$fitted)

		## paneling, from hector
		mypar(1,1, brewer.name = "Set1")
		par(mar=par()$mar+c(0,3,0,0))
		omar=par()$mar
		cmar=omar
		cmar[1]=.5
		par(mar=cmar)
		
		layout(cbind(1:sum(panels)),height=c(1+2/3,1, 0.75,0.75)[1:sum(panels)])
		if(showMethPanel) {
			matplot(x,beta,col=as.numeric(factor(coi)),type="p",pch=".",
			xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,1))
			matplot(x,smns,col=1:2,type="l",lwd=2.5,add=TRUE,lty=1)
			legend("bottomright",col=seq(along=levels(factor(coi))),
				lty=1,lwd=2,legend=levels(factor(coi)),cex=.8, bty=bty)
				
			segments(min(x), .2, min(x)+scale*1000, .2)
			text(min(x),.05,labels=sprintf("%dkb",scale),pos=4,offset=0)
			axis(side=2,at=c(.2,.5,.8), cex.axis=1.8)
			mtext(sprintf("%s:%d-%d", seqnames(r), start(r), end(r)), side=3)
			mtext("Methylation",side=2, line = 2.5,cex=1.5)

			legend("topright", paste0("fwer = ",
				signif(blocks450$tab$fwer[i],3)),bty=bty)
		}	  
		# annotation
		if(showGenePanel) {
			cmar=omar
			if(!is.na(showDiffPanel)) {
				cmar[3]=0.5
			} else {
				cmar[c(1,3)]=c(0.5,0.5)
			}
			par(mar=cmar)

			plot(x,rep(0,length(x)), type="n",ylim=c(-1.5,1.5),yaxt="n",ylab="",
				 xlab="",cex.axis = 1.5, cex.lab =1.5,xaxt="n")
			a = as.data.frame(anno[[i]])
			Strand = ifelse(a$strand == "+", 1, ifelse(a$strand=="-", -1, 0))
			Col = ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue","white"))
			Lwd = ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron",0.5,0))
			axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1,cex.axis = 3)
			abline(h=0,lty=3)
			for(k in 1:nrow(a)) {
				polygon(c(a$start[k],a$end[k],a$end[k],a$start[k]),
					Strand[k]/2+c(-0.3,-0.3,0.3,0.3)*Lwd[k],col=Col[k])
			}

			## by gene
			g = split(a, sapply(a$symbol,"[", 1))
			# g = split(a, a$Symbol)
			s2 = ifelse(sapply(g, function(x) unique(x$strand))=="+",1,-1)
			g = sapply(g, function(x) (max(x$end) - min(x$start))/2 + min(x$start) )
			
			if(length(g) > 0) text(g, y=s2, names(g),font=1,pos=s2+2,cex=0.8)
					
			mtext("Genes",side=2, line = 2.5,cex=1.5)
			if(!showDiffPanel) {
				xtick=pretty(x)
				axis(side=1,at=xtick,labels=sprintf("%.1fMb",xtick / 1e6))
			}
		}

		## mean diff
		if(showDiffPanel) {
			cmar=omar
			cmar[c(1,3)]=c(.5,.5)
			par(mar=cmar)

			zz=granges(tmp)

			matplot(x,sd,xaxt="n",ylab="",xlab="",type="n",lty=1,ylim=c(-.6,.6),yaxt="n",pch=21)
			axis(side=2,at=c(-.3,0,.3),labels=c("-.3","0",".3"))
			xtick=pretty(x)
			axis(side=1,at=xtick,labels=sprintf("%.1fMb",xtick / 1e6))
			mtext("Diff",side=2, line = 2.5,cex=1.5)

			ii=which(zz$type=="OpenSea")
			blockgroup=zz$blockgroup[ii]

			blockIndexes=split(seq(along=blockgroup),blockgroup)
			for (ind in blockIndexes) {
				ind=ii[ind]
				lines(x[ind], sd[ind], lwd=2.5,col="black")
			}


			points(x[ii],d[ii],pch=21,cex=1.4,bg="black")
			axis(side=2,at=c(-2,0,2))
			abline(h=0,lty=2,col="black")

			cmar=omar
			cmar[3]=.5
			par(mar=cmar)
			matplot(x,beta,type="n",xaxt="n",yaxt="n",xlab="",
				ylab="",ylim=c(0,2),bty="n")
		}
		
		#  browser()
		if(showCancerPanel) {
			col=ifelse(tmp450$value<0 & tmp450$p.value<.05,"blue",ifelse(tmp450$value>0 & tmp450$p.value<.05,"red","black"))
			rect(start(tmp450),1+1/3,end(tmp450),1+2/3,col=col)
			if(length(tmpBsmooth) > 0)  rect(start(tmpBsmooth),1/3,end(tmpBsmooth),2/3,col=ifelse(tmpBsmooth$direction=="hypo","blue","red"))
			axis(side=2,at=c(.5,1.5),labels=c("Hansen et al.",blockname),las=1,lwd=0)
			legend("bottomleft",pt.bg=c("blue","red"),legend=c("hypo","hyper"),pch=22,cex=.8)
		}
	}
	dev.off()
}



parsePheno = function(pd) {
	require(Hmisc)
	phenoIndex=grep("characteristics",colnames(pd))
	
	#  column names
	colNames = apply(pd[,phenoIndex], 2, function(x) {
		ss(as.character(x),": ", 1)[1]
	})
	colNames = capitalize(colNames) # in Hmisc
	colnames(pd)[phenoIndex] = colNames
	
	 #variables themselves	
	pd[,phenoIndex] = apply(pd[,phenoIndex], 2, function(x) {
		ss(as.character(x),": ", 2)
	})
	pd$title = as.character(pd$title)
	pd$geo_accession = as.character(pd$geo_accession)
	
	tmp = colSums(apply(pd, 2, as.numeric))
	for(i in which(!is.na(tmp))) pd[,i] = as.numeric(pd[,i])
	
	rownames(pd) = pd$geo_accession
	return(pd)
}

# y: measurement
# age: sample age
# mod: model to fit
# mainText: what to put in the title
agePlotter = function(y, age, mod = matrix(rep(1,length(y)),ncol=1),
	mainText, smoothIt=TRUE, jitter=TRUE, ageLabel = "bottom",
	orderByAge=TRUE,ylim=NULL,ageBreaks = c(-1, 0, 1, 10, 100), 
	ylab="Adjusted Expression",pointColor = 2, lineColor= 1, alreadyFitted=NULL,
		...) {
	
	if(orderByAge) {
		oo = order(age, decreasing = FALSE)
		y = y[oo] ; age = age[oo] ; mod = mod[oo,]
		if(!is.null(alreadyFitted)) alreadyFitted = alreadyFitted[oo]
	}
	
	if(is.null(alreadyFitted)) {
		fit = fitted(lm(y~mod-1))
	} else fit = alreadyFitted
	
	fetal = cut(age, breaks =ageBreaks ,lab=FALSE)
	fIndex = splitit(fetal)	
	
	layout(matrix(c(1,1,1,2,2,3,3,4,4,4,4,4),nr = 1,byrow = TRUE))
	palette(brewer.pal(8,"Set1"))
	
	par(mar = c(4,5,3,0.45))
	if(is.null(ylim)) ylims = range(y,na.rm=TRUE) else ylims = ylim

	if(jitter) xx = jitter(age,amount=0.005) else xx=age
	plot(y ~ xx,
		subset=fIndex[[1]],
		main = "",	ylab=ylab,xlab="",
		ylim = ylims,cex.axis = 1.5, cex.lab=1.75,
		pch = 21, cex = 1.4,xaxt="n",bg = pointColor,
		xlim=c(range(age[fIndex[[1]]])+c(-0.01,0.07)),...)
	
	if(smoothIt) lines(age[fIndex[[1]]],fit[fIndex[[1]]],col=lineColor,lwd=6)
	
	axis(1,at=unique(age[fIndex[[1]]]), 
		labels = 40+52*signif(unique(age[fIndex[[1]]]),1), cex.axis=1.5)
	
	if(ageLabel == "bottom") {
		text(x = quantile(age[fIndex[[1]]],0.33), y= min(ylims), "PCW", cex=1.5)
	} else if(ageLabel == "top") {
		text(x = quantile(age[fIndex[[1]]],0.33), y= max(ylims), "PCW", cex=1.5)
	} 
	
	# infant + child
	par(mar = c(4, 0.25,3,0.25))
	for(j in 2:3) {
		plot(y ~ age,	subset=fIndex[[j]],
			main = "",ylab="",xlab="",yaxt = "n", cex=1.4,
			xlim = range(age[fIndex[[j]]])+c(-0.03,0.03),
			ylim = ylims, cex.axis = 1.5,pch = 21,  bg=pointColor)

		if(ageBreaks[2] == 0 & smoothIt) lines(age[fIndex[[j]]],fit[fIndex[[j]]],col=lineColor,lwd=6)
		if(ageBreaks[2] < 0 & smoothIt) lines(age[fIndex[[j]]][-1],fit[fIndex[[j]]][-1],col=lineColor,lwd=6)
	}
	
	# adults
	par(mar = c(4, 0.25,3,1))
	plot(y ~ age,	subset=fIndex[[4]],
			main = "",ylab="",xlab="",yaxt = "n", cex=1.4,
			xlim = range(age[fIndex[[4]]])+c(-0.01,0.01),
			ylim = ylims, cex.axis = 1.5,pch = 21, bg=pointColor)

	if(smoothIt) lines(age[fIndex[[4]]],fit[fIndex[[4]]],col=lineColor,lwd=6)

	mtext(mainText,	outer=T, line=-2.5,cex=1.35)	
	
	mtext("Age", side=1, outer=T, line=-1.5,cex=1.35)
}

require(RColorBrewer)


mypar = function(a=1,b=1,brewer.n=8,brewer.name="Dark2",...){
 par(mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
 par(mfrow=c(a,b),...)
 palette(brewer.pal(brewer.n,brewer.name))
}

makingCOI = function(coi, n , returnMean = TRUE) {
	groups = cut(coi, quantile(coi,
		prob=seq(0,1,1/(n-1)), na.rm=T),include.lowest=T,
		labels = F)
	
	if(returnMean) {
		gIndexes = splitit(groups)
		groupLabels = rep(0,length(groups))
		for(m in seq(along=gIndexes)) {
			groupLabels[gIndexes[[m]]] = signif(mean(coi[gIndexes[[m]]]),2)
		}
		groupLabels = factor(groupLabels)
		return(groupLabels)
	} else return(groups)
}

cleaningP = function(y, mod, P=ncol(mod)) {
	X=mod
	Hat=solve(t(X)%*%X)%*%t(X)
	beta=(Hat%*%t(y))
	cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
	return(cleany)
}

## from Rafa
dogo <- function(names,universe,species="human", goP = 0.01, 
	cond=FALSE, ontology = "BP"){
    if(species=="human"){
		golib="org.Hs.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Hs.egREFSEQ2EG
  } else  if (species == "mouse") {
		golib="org.Mm.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Mm.egREFSEQ2EG
  } else if (species == "rat") {
		golib="org.Rn.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Rn.egREFSEQ2EG
	}
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA))
  x=x[!is.na(x)]
  Universe=unlist(mget(as.character(universe),gomap,ifnotfound = NA))
 Universe=unique(c(Universe[!is.na(Universe)],unique(x)))

  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = ontology, pvalueCutoff = goP, conditional = cond,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)

}	

#### to-do: replace with derfinderPlot:plotCoverage
dmrPlot = function(regions, p, chr, pos, cluster, genes, coi,
	genomicState, build="hg19", species = "human", Jitter = FALSE,
	number=100,cols=NULL, lines = FALSE, linesSmooth = TRUE,
	title = TRUE, Legend = TRUE, colorRamp = FALSE, 
	meanSmooth=TRUE, plotCpG = TRUE, geneAnno = "gene") {
	
	require(bumphunter)
	require(derfinder)
	gr = GRanges(regions$chr, IRanges(regions$start, regions$end))
	
	if(build == "hg18") {
		cpg.cur = read.table("http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg18.txt",
			header = TRUE, as.is=TRUE)
		library("BSgenome.Hsapiens.UCSC.hg18")
	}
	
	if(build == "hg19") {
		cpg.cur = read.table("http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg19.txt",
			header = TRUE, as.is=TRUE)
		library("BSgenome.Hsapiens.UCSC.hg19")
	}
	
	if(build == "mm9") {
		cpg.cur = read.table("http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-mm9.txt",
			header = TRUE, as.is=TRUE)
		library("BSgenome.Mmusculus.UCSC.mm9")
	}
		
	ocpgi=data.frame(chr=I(cpg.cur[,1]), 
		start=as.numeric(cpg.cur[,2]), 
		end=as.numeric(cpg.cur[,3]))
	ADD1 = 1500; PAD = 10
	
	gr2 = GRanges(regions$chr, IRanges(regions$start - ADD1, regions$end + ADD1))
	anno = annotateRegions(gr2, gs)$annotationList
	
	# check regions
	if(is.numeric(coi)) {
		groups=coi
		gNames= sort(unique(coi))
	}

	if(is.character(coi) | is.factor(coi)) {
		groups = factor(coi)
		gNames= levels(groups)
	}
	gIndexes=split(1:length(groups),groups)
	
	brewer.n=max(3,min(length(unique(coi)),11))

	
	if(is.null(cols)) {
		mypar(brewer.n=brewer.n)
	} else if(length(cols) == 1 & cols[1] %in% rownames(brewer.pal.info)) {
		mypar(brewer.name=cols,	brewer.n=brewer.n)
	} else {
		mypar()
		palette(cols)
	} 
		
	if(colorRamp) {
		pal = colorRampPalette(palette()[3:brewer.n])
		palette(pal(brewer.n))
	}
	
	cat("Plotting.\n")
	for(j in 1:(min(nrow(regions), number))) {
	
		layout(matrix(1:2,ncol=1),heights=c(0.7,0.3))

		# first plot, region
		par(mar=c(0,4.5,0.25,1.1),oma=c(0,0,2,0))
		
		Index=(regions[j,7]-PAD):(regions[j,8]+PAD)
		Index = Index[Index %in% seq(along=cluster)]
		Index=Index[cluster[Index]==regions[j,6]]
		
		# make DMR plots
		if(Jitter) {
			posx = matrix(pos[Index], nc = ncol(p), nr = length(Index),
				byrow=  FALSE)
			posx = t(apply(posx, 1, function(x) jitter(x,amount = 12)))
			
		} else posx = pos[Index]
					
		matplot(posx, p[Index,], ylim = c(0,1),
			ylab = "", xlab = "",xaxt = "n",cex=0.7,
			cex.axis = 1.7, cex.lab = 1.7, pch=21,
			bg = as.numeric(factor(groups)),col="black",
			xlim = range(pos[Index]), yaxt="n")
		axis(2, at = c(0.2, 0.5, 0.8), cex.axis=1.7)

		xx=pos[Index]
		for(k in seq(along=gIndexes)){
			if(length(gIndexes[[k]]) == 1) yy=p[Index,gIndexes[[k]]]
			if(length(gIndexes[[k]]) > 1)	  yy=rowMeans(p[Index,gIndexes[[k]]])
			if(meanSmooth) { 
				fit1=loess(yy~xx,degree=1,span=300/(36*length(xx)),
					family="symmetric")
				lines(xx,fit1$fitted,col=k,lwd=2)
			} else 	lines(xx,yy,col=k,lwd=2)

		}
		
		mtext("Methylation",side=2, line = 2.5,cex=1.8)

		if(Legend) {
			if(length(unique(coi)) < 4) {
				legend("topleft",legend=gNames,col=1:length(gNames),
				lty=1, lwd = 4,cex=1,bty="n")
			} else {
				legend("topleft",legend=gNames,col=1:length(gNames),
					pch=19, pt.cex = 2,cex=1.1, nc = 6,bty="n")
			}	
		}
				
		abline(v=(regions$start[j]-15),lty=1)
		abline(v=(regions$end[j]+15),lty=1)

		if(title) mtext(paste0(genes$name[j],"; ", 
			genes$distance[j],"bp from tss:",genes$description[j]), outer=T,cex=1.3)
			
		# add cpgs
		if(plotCpG) {
			thechr=as.character(regions$chr[j])
			chrName = strsplit(thechr, "r")[[1]][2]
			chrName = paste("Chromosome",chrName)
			
			start = pos[Index[1]]
			end = pos[Index[length(Index)]]
			ocpgi2=ocpgi[ocpgi$chr%in%unique(as.character(thechr)),]
			
			##PLOT CPG ISLANDS
			if(species=="human") seq<-Hsapiens[[as.character(thechr) ]]
			if(species=="mouse") seq<-Mmusculus[[as.character(thechr) ]]
			
			subseq<-subseq(seq,start=start,end=end)
			cpgs=start(matchPattern("CG",subseq))+start-1

			if(plotCpG) Rug(cpgs,col="black")

			Index1 = which(ocpgi2[,1]==as.character(thechr) &
					 ((ocpgi2[,2] > start & ocpgi2[,2]< end) |
					  (ocpgi2[,3] > start & ocpgi2[,3]< end)))
			if(length(Index1)>0) sapply(Index1,function(j) Rug(unlist(ocpgi2[j,2:3]),
						 col="darkgreen",lwd=3,side=1))
		}
		
		# plot 3
		##PLOT GENES
		par(mar=c(3.5,4.5,0.25,1.1))

		plot(0,0, type="n", xlim=range(xx),ylim=c(-1.5,1.5),yaxt="n",ylab="",
			 xlab="",cex.axis = 1.5, cex.lab =1.5)
		a = as.data.frame(anno[[j]])
		Strand = ifelse(a$strand == "+", 1, ifelse(a$strand=="-", -1, 0))
		Col = ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue","white"))
		Lwd = ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron",0.5,0))
		axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1,cex.axis = 3)
		abline(h=0,lty=3)
		for(k in 1:nrow(a)) {
			polygon(c(a$start[k],a$end[k],a$end[k],a$start[k]),
				Strand[k]/2+c(-0.3,-0.3,0.3,0.3)*Lwd[k],col=Col[k])
		}
		
		if(sum(a$theRegion=="exon") > 0) {
			e = a[a$theRegion=="exon",]
			s2 = Strand[a$theRegion=="exon"]
			g = unlist(e$symbol)
			g[is.na(g)] = ""
			if(length(g) > 0) text(x=e$start + e$width/2,
				y=s2*0.75, g,font=1,pos=s2+2,cex=1.2)
		}
				
		mtext("Genes",side=2, line = 2.5,cex=1.5)
		mtext(chrName,side=1, line = 2,cex=1.4)

		abline(v=(pos[regions[j,7]]-15),lty=1)
		abline(v=(pos[regions[j,8]]+15),lty=1)
		
	}
}

Rug <- function (x, ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"), 
    quiet = getOption("warn") < 0, ...) {
    x <- as.vector(x)
    ok <- is.finite(x)
    x <- x[ok]
    if (!quiet) {
        u <- par("usr")
        u <- if (side%%2 == 1) {
            if (par("xlog")) 
                10^u[1:2]
            else u[1:2]
        }
        else {
            if (par("ylog")) 
                10^u[3:4]
            else u[3:4]
        }
        if (any(x < u[1] | x > u[2])) 
            warning("some values will be clipped")
    }
    Axis(side = side, at = x, labels = FALSE, lwd = lwd, col = col, 
        tck = ticksize, ...)
}

