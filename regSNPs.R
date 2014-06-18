#############################################################################################################
### 	regSNPs: a strategy for prioritizing regulatory single nucleotide substitutions.
### 	By: Mingxiang Teng
#############################################################################################################

#############################################################################################################
### 	Main Procedure of Promoter SNP Prioritization.
#############################################################################################################
regSNPs <- function(expt='test',SNP='SNPs.txt',EndeavourResult='endeavour_results.tsv',dir='./',GenomeDir='hg18/',HapMapSNP='HapMapSNPsCEURel27.RData',HapMapSNPsProSeq='HapMapSNPsCEURel27_ProSeq.RData',TFBS='TFBS.RData',PromoterSNPsRatio='PromoterSNPsRatio.RData',RandomSNPsRatio='RandomSNPsRatio.RData',rs=FALSE,Permutation=FALSE,RandomTimes=1000,RanStep=200){
	library(GeneR)
	data <- InputInitiation(dir=dir,rs=rs,SNP=SNP,EndeavourResult=EndeavourResult,HapMapSNP=HapMapSNP,TFBS=TFBS,PromoterSNPsRatio=PromoterSNPsRatio,HapMapSNPsProSeq=HapMapSNPsProSeq,RandomSNPsRatio=RandomSNPsRatio)
	data$Tfbs$RK <- Endeavour(EndeaResult=data$Endea$RT,Tfbs=data$Tfbs)
	data$Snps$IN <- SequenceExtractor(SNPs=data$Snps$IN,GenomeDir=GenomeDir,ExtendBP=30)
	data$Snps$SC <- SnpTfbsMaxScoreCalculation(SNPsSeq=data$Snps$IN,PSSM=data$Tfbs$PSSM)
	data$Snps$RA <- SnpTfbsRatioCalculation(data=data,TfbsTrueScore=data$Tfbs$TS,TfbsRandomScore=data$Tfbs$RS)
	data$Snps$RK <- RankSNPs(data=data,PromoterSNPsRatio=data$Ran$RA,SampleSize=nrow(data$Snps$IN))
	if(Permutation==TRUE)
	{
		if(RandomTimes<100)
		{
			stop('The random times are too low!',call.=FALSE)
		}
		data$Ran$RK <- RandomROI(RandomTimes=RandomTimes,ROIData=data,step=RanStep)
		data$Snps$PV <- array(1,dim=nrow(data$Snps$IN))
		for(i in 1:nrow(data$Snps$IN))
		{
			data$Snps$PV[i] <- sum(data$Ran$RK <= as.numeric(data$Snps$RK[i,1]))/length(data$Ran$RK)
		}
		out <- data.frame(data$Snps$IN[,1:5],data$Snps$RK[,1:2],data$Snps$PV)
		cat('Outport final results!\n')
		write.table(out,paste(dir,'PrioritizedResults',RandomTimes,'.txt',sep=''),sep='\t',row.names=FALSE,quote=FALSE)
		save(data,file=paste(dir,expt,'_result',RandomTimes,'.RData',sep=''))
	} else 	{
		out <- data.frame(data$Snps$IN[,1:5],data$Snps$RK[,1:2])####add TF information form data$Tfbs$RK
		cat('Outport final results!\n')
		write.table(out,paste(dir,'PrioritizedResults.txt',sep=''),sep='\t',row.names=FALSE,quote=FALSE)
		save(data,file=paste(dir,expt,'_result.RData',sep=''))
	}	
}

#############################################################################################################
### 	Initial The Input.
#############################################################################################################
InputInitiation <- function(dir=NULL,rs=NULL,SNP=NULL,EndeavourResult=NULL,HapMapSNP=NULL,TFBS=NULL,PromoterSNPsRatio=NULL,HapMapSNPsProSeq=NULL,RandomSNPsRatio=NULL)
{
	data <- list()
	data$Endea$RT <- read.table(paste(dir,EndeavourResult,sep=''),sep='\t',header=TRUE)
	load(paste(dir,TFBS,sep=''))
	load(paste(dir,HapMapSNP,sep=''))
	load(paste(dir,PromoterSNPsRatio,sep=''))
	load(paste(dir,HapMapSNPsProSeq,sep=''))
	load(paste(dir,RandomSNPsRatio,sep=''))
	data$Tfbs <- Tfbs
	if(rs==TRUE){
		SNPs <- read.table(paste(dir,SNP,sep=''),sep='\t',header=FALSE)
		SNPs <- HapMapSNPs[match(as.character(SNPs[,1]),as.character(HapMapSNPs[,1])),]	
	} else{
		SNPs <- read.table(paste(dir,SNP,sep=''),sep='\t',header=TRUE)
	}
	SNPs <- SNPs[sort.list(SNPs[,'pos']),]
	SNPs <- SNPs[sort.list(SNPs[,'chrom']),]
	data$Snps$IN <- SNPs
	data$ProSnps$IN <- HapMapSNPsProSeqs
	data$ProSnps$RA <- PromoterSNPsRatios
	data$Ran$RA <- RandomSNPsRatios
	cat('Initial successful.\n')
	data;
}

#############################################################################################################
### 	Rank TFBSs By Endeavour.
#############################################################################################################
Endeavour <- function(EndeaResult=NULL,Tfbs=NULL){
	TFBSPreRank <- cbind(Tfbs$TF[,1:6],EndeaResult[match(as.character(Tfbs$TF[,7]),as.character(EndeaResult[,1])),])
	TFBSPreRank[is.na(TFBSPreRank[,'rank']),'rank'] <- 1
	TFBSFinalRank <- NULL
	for(i in 1:length(Tfbs$AC))
	{
		temp <- TFBSPreRank[as.character(TFBSPreRank[,1])==Tfbs$AC[i],]
		rank <- min(as.numeric(as.character(temp[,'rank'])))
		TFBSFinalRank <- rbind(TFBSFinalRank,temp[as.numeric(as.character(temp[,'rank']))==rank,])
	}
	cat('Endeavour rank snalysis successful.\n')
	TFBSFinalRank;
}

#############################################################################################################
### 	Extract Sequence Of SNP Around Region.
#############################################################################################################
SequenceExtractor <- function(SNPs=NULL,GenomeDir=NULL,ExtendBP=NULL)
{
	chr=''
	Sequence <- NULL
	SeqLength <- NULL
	SequenceUp <- NULL
	SequenceDown <- NULL
	RefAllele <- NULL
	OthAllele <- NULL
	for(i in 1:nrow(SNPs))
	{
		if(chr!=as.character(SNPs[i,'chrom'])) 
		{
			freeAllSeq()
			chr <- as.character(SNPs[i,'chrom'])
			readFasta(paste(GenomeDir, chr, ".fa", sep=""))
			chr.size <- sizeSeq(seqno=0)
		}
		UpAndDownSeq <- getSeq(from = max(1, as.numeric(SNPs[i,'pos'])-ExtendBP), to = min(as.numeric(SNPs[i,'pos'])+ExtendBP, chr.size))
		UpSeq <- getSeq(from = max(1, as.numeric(SNPs[i,'pos'])-ExtendBP), to = min(as.numeric(SNPs[i,'pos'])-1, chr.size))
		DownSeq <- getSeq(from = max(1, as.numeric(SNPs[i,'pos'])+1), to = min(as.numeric(SNPs[i,'pos'])+ExtendBP, chr.size))
		RefA <- getSeq(from = max(1, as.numeric(SNPs[i,'pos'])), to = min(as.numeric(SNPs[i,'pos']), chr.size))
		Sequence <- c(Sequence,toupper(UpAndDownSeq))
		SeqLength <- c(SeqLength,nchar(UpAndDownSeq))
		SequenceUp <- c(SequenceUp,toupper(UpSeq))
		RefAllele <- c(RefAllele,toupper(RefA))
		Alleles <- unlist(strsplit(as.character(SNPs[i,'alleles']),split='/'))
		OthAllele <- c(OthAllele,toupper(Alleles[is.na(match(Alleles,toupper(RefA)))]))
		SequenceDown <- c(SequenceDown,toupper(DownSeq))
		#cat(i,'\n')
	}
	SNPsSeq <- cbind(SNPs,Sequence,SeqLength,RefAllele,OthAllele,SequenceUp,SequenceDown)
	cat('Sequences extracttion successful for totally ',nrow(SNPs), ' SNPs.\n')
	SNPsSeq;
}

#############################################################################################################
### 	Calculate Max PSSM Score For Each SNP-TFBS Pair.
#############################################################################################################
SnpTfbsMaxScoreCalculation <- function(SNPsSeq=NULL,PSSM=NULL)
{
	RefMaxScore <- NULL
	OthMaxScore <- NULL
	MaxScore <- list()
	for(i in 1:length(PSSM))
	{
		ScoreRefCurStr <- matrix(0,nrow(SNPsSeq),ncol(PSSM[[i]]))
		ScoreOthCurStr <- matrix(0,nrow(SNPsSeq),ncol(PSSM[[i]]))
		ScoreRefComStr <- matrix(0,nrow(SNPsSeq),ncol(PSSM[[i]]))
		ScoreOthComStr <- matrix(0,nrow(SNPsSeq),ncol(PSSM[[i]]))
		RefCurStr <- paste(substr(as.character(SNPsSeq[,'SequenceUp']),nchar(as.character(SNPsSeq[,'SequenceUp']))-ncol(PSSM[[i]])+2,nchar(as.character(SNPsSeq[,'SequenceUp']))),as.character(SNPsSeq[,'RefAllele']),substr(as.character(SNPsSeq[,'SequenceDown']),1,ncol(PSSM[[i]])-1),sep="")
		OthCurStr <- paste(substr(as.character(SNPsSeq[,'SequenceUp']),nchar(as.character(SNPsSeq[,'SequenceUp']))-ncol(PSSM[[i]])+2,nchar(as.character(SNPsSeq[,'SequenceUp']))),as.character(SNPsSeq[,'OthAllele']),substr(as.character(SNPsSeq[,'SequenceDown']),1,ncol(PSSM[[i]])-1),sep="")
		for(j in 1:nrow(SNPsSeq))
		{
			RefComStr <- strComp(RefCurStr[j])
			OthComStr <- strComp(OthCurStr[j])
			RefCurStrMotif <- MotifBinaryMatrix(RefCurStr[j])
			OthCurStrMotif <- MotifBinaryMatrix(OthCurStr[j])
			RefComStrMotif <- MotifBinaryMatrix(RefComStr)
			OthComStrMotif <- MotifBinaryMatrix(OthComStr)
			for(k in 1:ncol(PSSM[[i]]))
			{
				ScoreRefCurStr[j,k] <- round(sum(RefCurStrMotif[,k:(k+ncol(PSSM[[i]])-1)] * PSSM[[i]]),3)
				ScoreOthCurStr[j,k] <- round(sum(OthCurStrMotif[,k:(k+ncol(PSSM[[i]])-1)] * PSSM[[i]]),3)
				ScoreRefComStr[j,k] <- round(sum(RefComStrMotif[,k:(k+ncol(PSSM[[i]])-1)] * PSSM[[i]]),3)
				ScoreOthComStr[j,k] <- round(sum(OthComStrMotif[,k:(k+ncol(PSSM[[i]])-1)] * PSSM[[i]]),3)
			}
		}
		RefMaxScore <- cbind(RefMaxScore,apply(cbind(ScoreRefCurStr,ScoreRefComStr),1,max))
		OthMaxScore <- cbind(OthMaxScore,apply(cbind(ScoreOthCurStr,ScoreOthComStr),1,max))
		#cat(i,'\n')
	}
	cat('PSSM score calculation successful for ',nrow(SNPsSeq), ' SNPs and ',length(PSSM),' TFBS.\n')
	MaxScore$Ref <- RefMaxScore
	MaxScore$Oth <- OthMaxScore
	MaxScore;
}

#############################################################################################################
### 	Generate Binary Matrix For DNA Sequence Motifs.
#############################################################################################################
MotifBinaryMatrix <- function (Motif = NULL)
{
	Matrix <- matrix(0, nrow = 4, ncol = nchar(Motif))
	if (gregexpr("[aA]", Motif)[[1]][1]>0)
	{
		Matrix[1,gregexpr("[aA]", Motif)[[1]]] = 1
	}
	if (gregexpr("[cC]", Motif)[[1]][1]>0)
	{
		Matrix[2,gregexpr("[cC]", Motif)[[1]]] = 1
	}
	if (gregexpr("[gG]", Motif)[[1]][1]>0)
	{
		Matrix[3,gregexpr("[gG]", Motif)[[1]]] = 1
	}
	if (gregexpr("[tT]", Motif)[[1]][1]>0)
	{
		Matrix[4,gregexpr("[tT]", Motif)[[1]]] = 1
	}
	Matrix;
}

#############################################################################################################
### 	Calculate Ratio For Each SNP-TFBS Pair.
#############################################################################################################
SnpTfbsRatioCalculation <- function(data=NULL,TfbsTrueScore=NULL,TfbsRandomScore=NULL)
{
	Ratio <- NULL
	for(i in 1:length(data$Tfbs$AC))
	{
		RefT <- 1-pnorm(data$Snps$SC$Ref[,i],TfbsTrueScore[[i]]$MN,TfbsTrueScore[[i]]$SD)
		RefF <- 1-pnorm(data$Snps$SC$Ref[,i],mean(TfbsRandomScore[,i]),sd(TfbsRandomScore[,i]))
		OthT <- 1-pnorm(data$Snps$SC$Oth[,i],TfbsTrueScore[[i]]$MN,TfbsTrueScore[[i]]$SD)
		OthF <- 1-pnorm(data$Snps$SC$Oth[,i],mean(TfbsRandomScore[,i]),sd(TfbsRandomScore[,i]))
		Ratio <- cbind(Ratio,log2((OthT*RefF)/(RefT*OthF)))
	}
	colnames(Ratio) <- data$Tfbs$AC
	rownames(Ratio) <- as.character(data$Snps$IN[,1])
	cat('Ratio calculation successful for ',nrow(data$Snps$IN), ' SNPs and ',length(data$Tfbs$AC),' TFBS.\n')
	Ratio;
}

#############################################################################################################
### 	Rank The Candidate SNPs (both TFBS individually and overall).
#############################################################################################################
RankSNPs <- function(data=NULL,PromoterSNPsRatio=NULL,RandomTimes=1,SampleSize=NULL)
{
	Ratio <- data$Snps$RA
	RankSnpTfbs <- matrix(1,nrow(Ratio),ncol(Ratio))
	for(i in 1:ncol(RankSnpTfbs))
	{
		RankSnpTfbs[,i] <- rowSums(abs(t(PromoterSNPsRatio[,rep(i,nrow(Ratio))]))-abs(Ratio[,i])>0)/nrow(PromoterSNPsRatio)
		cat('.')
	}
	cat('\n')
	TfbsRank <- as.numeric(as.character(data$Tfbs$RK[,'rank']))
	TfbsRank <- matrix(TfbsRank,RandomTimes,length(data$Tfbs$AC),byrow=TRUE)
	TfbsRank <- TfbsRank[sort(rep(1:RandomTimes,SampleSize)),]
	OverAllScores <- matrix(1,nrow(Ratio),ncol(Ratio))
	colnames(OverAllScores) <- colnames(Ratio)
	OverAllScores <- 1-(1-RankSnpTfbs)*(1-TfbsRank)
	FinalScore <- apply(OverAllScores,1,min)
	FinalRank <- FinalScore
	if(RandomTimes==1)
	{
		FinalTFBSIndex <- NULL
		for(i in 1:length(FinalScore))
		{
			FinalTFBSIndex <- c(FinalTFBSIndex,which(OverAllScores[i,]==FinalScore[i])[1])
		}
		FinalTFBS <- colnames(Ratio)[FinalTFBSIndex]
		FinalRank <- cbind(FinalScore,FinalTFBS,OverAllScores)
		rownames(FinalRank) <- as.character(data$Snps$IN[,1])
		cat('Rank calculation successful for ',nrow(data$Snps$IN), ' SNPs and their corresponding TFBSs.\n')
	}
	FinalRank;
}

#############################################################################################################
### 	Prepare Random ROI for Final P-value (not nessesary if you just need to rank SNP candidates)
#############################################################################################################
RandomROI <- function(RandomTimes=NULL,ROIData=NULL,step=NULL){
	SNP_NUM <- nrow(ROIData$ProSnps$IN)
	ROI_NUM <- nrow(ROIData$Snps$IN)
	Scores <- NULL
	curTime <- 1
	while((curTime+step-1)<=RandomTimes)
	{
		SNP_Index <- NULL
		Rank_Index <- NULL
		for(i in 1:step)
		{
			SNP_Index <- c(SNP_Index,sample(1:SNP_NUM,ROI_NUM))
			Rank_Index <- c(Rank_Index,sample(1:nrow(ROIData$Tfbs$RK),nrow(ROIData$Tfbs$RK)))
		}
		RanData <- list()
		RanData$Snps$IN <- ROIData$ProSnps$IN[SNP_Index,]
		RanData$Snps$RA <- ROIData$ProSnps$RA[SNP_Index,]
		RanData$Tfbs$AC <- ROIData$Tfbs$AC
		RanData$Tfbs$RK <- ROIData$Tfbs$RK
		RanData$Tfbs$RK <- RanData$Tfbs$RK[rep(1:length(RanData$Tfbs$AC),step),]
		RanData$Tfbs$RK[,'rank'] <- RanData$Tfbs$RK[Rank_Index,'rank']
		Scores <- c(Scores,RankSNPs(data=RanData,PromoterSNPsRatio=ROIData$Ran$RA,RandomTimes=step,SampleSize=ROI_NUM))
		curTime <- curTime + step
		cat(curTime -1,'random times are complete.\n')
	}
	if(curTime<=RandomTimes &(curTime+step-1)>RandomTimes)
	{
		SNP_Index <- NULL
		Rank_Index <- NULL
		for(i in curTime:RandomTimes)
		{
			SNP_Index <- c(SNP_Index,sample(1:SNP_NUM,ROI_NUM))
			Rank_Index <- c(Rank_Index,sample(1:nrow(ROIData$Tfbs$RK),nrow(ROIData$Tfbs$RK)))
		}
		RanData <- list()
		RanData$Snps$IN <- ROIData$ProSnps$IN[SNP_Index,]
		RanData$Snps$RA <- ROIData$ProSnps$RA[SNP_Index,]
		RanData$Tfbs$AC <- ROIData$Tfbs$AC
		RanData$Tfbs$RK <- ROIData$Tfbs$RK
		RanData$Tfbs$RK <- RanData$Tfbs$RK[rep(1:length(RanData$Tfbs$AC),(RandomTimes-curTime+1)),]
		RanData$Tfbs$RK[,'rank'] <- RanData$Tfbs$RK[Rank_Index,'rank']
		Scores <- c(Scores,RankSNPs(data=RanData,PromoterSNPsRatio=ROIData$Ran$RA,RandomTimes=(RandomTimes-curTime+1),SampleSize=ROI_NUM))
		curTime <- RandomTimes + 1
		cat(curTime-1,'random times are complete.\n')
	}
	Scores;
}

