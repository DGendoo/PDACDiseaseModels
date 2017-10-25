# Deena M.A. Gendoo
# March 2, 2017
# Calculate Concordance across the genomic intervals between PDX-PDO pairs, for each chromosome
# FINAL_CalculateConcordance_Trios_XvsO.R

##################################################
# SETUP
##################################################
library(pheatmap)
chromosomes<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
               "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
colscheme<-c('#ffffff','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#0c2c84')  # YellowToBlue
colscheme2<-c('#ffffff','#bababa','#ffff99','#d53e4f','#fdae61','#66c2a5','#3288bd','#542788') # Multicolor
colscheme3<-c('#ffffff','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026') # YellowToRed
colscheme4<-c('#ffffff','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177') # PurpleHues

CNVs<-sort(list.files(pattern = "_XvsO_CN_Final.bed"))
CNVListTemp<-lapply(CNVs,function(x){read.table(x,header=FALSE)})

CNVs<-sapply(CNVs,function(x) gsub(x,pattern = "segments_", replacement = ""))
names(CNVListTemp)<-CNVs
names(CNVListTemp)<-gsub(x = names(CNVListTemp),pattern = "_XvsO_CN_Final.bed",replacement = "" )

Params<-read.csv("CNV_Params_PDAC_Trios.csv")
colnames(Params)[1]<-"Sample"
rownames(Params)<-Params$Sample
Params$Ploidy

for(count in 1:length(CNVListTemp))
{
  message("Processing Pair #: ",count)
  
  colnames(CNVListTemp[[count]])<-c("PDXChr","PDXStart","PDXEnd","PDXSample","PDXCN",
                                    "PDOChr","PDOStart","PDOEnd","PDOSample","PDOCN","GenomicOverlap")

  CNVListTemp[[count]]$PDXSample<-paste(names(CNVListTemp)[[count]],"_X",sep="")
  CNVListTemp[[count]]$PDOSample<-paste(names(CNVListTemp)[[count]],"_O",sep="")
  
  for (tally in 1:nrow(CNVListTemp[[count]]))
  {
    ifelse((CNVListTemp[[count]]$GenomicOverlap[tally] == 0 & CNVListTemp[[count]]$PDXCN[tally] == -1),
           CNVListTemp[[count]]$PDXCN[tally] <- NA,CNVListTemp[[count]]$PDXCN[tally] <- CNVListTemp[[count]]$PDXCN[tally])
    ifelse((CNVListTemp[[count]]$GenomicOverlap[tally] == 0 & CNVListTemp[[count]]$PDOCN[tally] == -1),
           CNVListTemp[[count]]$PDOCN[tally] <- NA,CNVListTemp[[count]]$PDOCN[tally] <- CNVListTemp[[count]]$PDOCN[tally])
  }
  
  for (rownum in 1:nrow(CNVListTemp[[count]]))
  {
    if(CNVListTemp[[count]]$GenomicOverlap[rownum] == 0)
    {
      if(is.na(CNVListTemp[[count]]$PDXCN[rownum]))
      { CNVListTemp[[count]]$GenomicOverlap[rownum] <- (CNVListTemp[[count]]$PDOEnd[rownum]-CNVListTemp[[count]]$PDOStart[rownum])}
      if(is.na(CNVListTemp[[count]]$PDOCN[rownum]))
      { CNVListTemp[[count]]$GenomicOverlap[rownum] <- (CNVListTemp[[count]]$PDXEnd[rownum]-CNVListTemp[[count]]$PDXStart[rownum])}
    } 
  }  
  
  Ploidy_X<-Params[unique(CNVListTemp[[count]]$PDXSample),]$Ploidy
  Ploidy_O<-Params[unique(CNVListTemp[[count]]$PDOSample),]$Ploidy
  CNVListTemp[[count]]$PDXCN<-CNVListTemp[[count]]$PDXCN/Ploidy_X
  CNVListTemp[[count]]$PDOCN<-CNVListTemp[[count]]$PDOCN/Ploidy_O
  
  CNVListTemp[[count]]$Jump<-abs(CNVListTemp[[count]]$PDXCN - CNVListTemp[[count]]$PDOCN)
  
}

CNVList<-CNVListTemp
rm(CNVListTemp)

# Now Calculate the Jaccard 
ConcordanceRates<-NULL
for (count in 1:length(CNVList))
{
  Sample<-CNVList[[count]] #Sample
  Sample$PDXCN[is.na(Sample$PDXCN)]<-"MISS"
  Sample$PDOCN[is.na(Sample$PDOCN)]<-"MISS"
  
  PDXSub<-Sample[,c("PDXChr","PDXStart","PDXEnd","PDXSample","PDXCN")] 
  PDXSub<-PDXSub[PDXSub$PDXCN!="MISS",] 
  PDXSub<-PDXSub[!duplicated(PDXSub), ] 
  PDXSub$IntervalLength<-PDXSub$PDXEnd-PDXSub$PDXStart
  All_PDX_Bases<-sum(as.numeric(PDXSub$IntervalLength))
  PDOSub<-Sample[,c("PDOChr","PDOStart","PDOEnd","PDOSample","PDOCN")] 
  PDOSub<-PDOSub[PDOSub$PDOCN!="MISS",] 
  PDOSub<-PDOSub[!duplicated(PDOSub), ] 
  PDOSub$IntervalLength<-PDOSub$PDOEnd-PDOSub$PDOStart
  All_PDO_Bases<-sum(as.numeric(PDOSub$IntervalLength))
  CommonSub<-Sample
  CommonSub<-CommonSub[CommonSub$PDXCN!="MISS",] 
  CommonSub<-CommonSub[CommonSub$PDOCN!="MISS",] 
  All_Common_Bases<-sum(as.numeric(CommonSub$GenomicOverlap))
  IntersectingSub<-Sample[!is.na(Sample$Jump),]
  Intersecting_Bases<-sum(as.numeric(IntersectingSub$GenomicOverlap[which(IntersectingSub$Jump<= 0.25)]))
  ConcordanceRate<-Intersecting_Bases/(All_PDX_Bases + All_PDO_Bases - All_Common_Bases)
  
  message("ConcordanceRate: ",ConcordanceRate) 
  message("") 
  ConcordanceRates[count]<-ConcordanceRate
} 
ConcordanceRates<-data.frame("Sample"=names(CNVList),"Concordance Rate"=ConcordanceRates)


write.csv(ConcordanceRates,"CNV_Jaccard_PDAC_Trios_XenograftVsPDO.csv")
