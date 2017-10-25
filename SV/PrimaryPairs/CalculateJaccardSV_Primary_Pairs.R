# Deena M.A. Gendoo
# November 1, 2016, and final update June 23, 2016
# Parse Structural Variation (SV) data (in TSV format) for Paired Tumour-PDX samples
# CalculateJaccardSV.R
# Calculate Jaccard Index 

###################################################################################################################
###################################################################################################################

library(reshape2)
library(pheatmap)


chromList<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
             "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")

######################################################################################################
# FIND ALL EVENTS PERTINENT TO TUMOUR, PDX, AS WELL AS COMMON EVENTS FOR A TUMOUR-PDX PAIR
######################################################################################################

# Read all TSVs containing structural variation events for TUMOUR
TSVs<-sort(list.files(pattern = "_Pa_P_526.annotatedSV.tsv"))
TSVList<-lapply(TSVs,function(x){read.csv(x,header = T)})
TSVs<-sapply(TSVs,function(x) gsub(x,pattern = "_Pa_P_526.annotatedSV.tsv", replacement = ""))
names(TSVList)<-TSVs
TUMOUR<-TSVList
SampleNames<-sort(unique(TSVs))
unlist(lapply(TSVList,nrow))
TumourSV<-melt(lapply(TSVList,nrow))$value

# Read all TSVs containing structural variation events for XENOGRAFT
TSVs<-sort(list.files(pattern = "_Pa_X_526.annotatedSV.tsv"))
TSVList<-lapply(TSVs,function(x){read.csv(x,header = T)})
TSVs<-sapply(TSVs,function(x) gsub(x,pattern = "_Pa_X_526.annotatedSV.tsv", replacement = ""))
names(TSVList)<-TSVs
SampleNames<-sort(unique(TSVs))
XENOGRAFT<-TSVList
unlist(lapply(TSVList,nrow))
XenoSV<-melt(lapply(TSVList,nrow))$value

# Find all Common Events in a given tumour-pdx pair
CommonSV<-CommonCalls<-COMMON<-NULL
for(count in 1:length(TSVs))
{
  Tumour<-TUMOUR[[count]]
  PDX<-XENOGRAFT[[count]]
  
  #Calculate Common SV
  PDXsv<-paste(PDX$chr1,":",PDX$pos1,"-",PDX$chr2,":",PDX$pos2,sep="")
  Tumoursv<-paste(Tumour$chr1,":",Tumour$pos1,"-",Tumour$chr2,":",Tumour$pos2,sep="")
  CommonCalls[[count]]<-intersect(Tumoursv,PDXsv)
  
  #Show differences
  Tumoursv[!(Tumoursv%in%CommonCalls[[count]])]
  PDXsv[!(PDXsv%in%CommonCalls[[count]])]
  
  #Save Events that are common (taken from the Tumour file)
  COMMON[[count]]<-Tumour[(Tumour$chr1%in%PDX$chr1 & Tumour$chr2%in%PDX$chr2 & Tumour$pos1%in%PDX$pos1 & Tumour$pos2%in%PDX$pos2),]
}
names(COMMON)<-names(TUMOUR)
CommonSV<-unlist(lapply(CommonCalls,length))

######################################################################################################
# GENERAL JACCARD SCORE CALCULATION ACROSS THE ENTIRE GENOME
######################################################################################################
# Read all TSVs containing structural variation events for TUMOUR
Total_TumourSV<-TumourSV
# Read all TSVs containing structural variation events for XENOGRAFT
Total_XenoSV<-XenoSV
# All TSVs commont to a tumour-xenograft pair
Total_CommonSV<-CommonSV

#Generate Jaccard Index calculation for total SV events
Jaccard_Total<-NULL
for(count in 1:length(TSVs))
{
  Jaccard_Total[count]<-Total_CommonSV[count]/((Total_TumourSV[count]+Total_XenoSV[count])-Total_CommonSV[count])
}
names(Jaccard_Total)<-SampleNames
sort(Jaccard_Total)

######################################################################################################
# CALCULATE EVENTS BY CHROMOSOME, AS OPPOSED TO ONLY GENERAL SCORE
######################################################################################################
Tumour_TSVList_ByChr<-NULL
for(count in 1:length(TUMOUR))
{
  Sample<-TUMOUR[[count]]
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    ChrCounts[counter]<-(length(which(Sample$chr1 %in% chromList[counter]))+length(which(Sample$chr2 %in% chromList[counter])))/2
  }
  Tumour_TSVList_ByChr[[count]]<-ChrCounts
}
Tumour_TSVList_ByChr<-data.frame(matrix(unlist(Tumour_TSVList_ByChr), nrow=10, byrow=T))
rownames(Tumour_TSVList_ByChr)<-names(TUMOUR)
colnames(Tumour_TSVList_ByChr)<-chromList

#Calculate Number of Events by Chromosome for XENOGRAFT
Xeno_TSVList_ByChr<-NULL
for(count in 1:length(XENOGRAFT))
{
  Sample<-XENOGRAFT[[count]]
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    ChrCounts[counter]<-(length(which(Sample$chr1 %in% chromList[counter]))+length(which(Sample$chr2 %in% chromList[counter])))/2
  }
  Xeno_TSVList_ByChr[[count]]<-ChrCounts
}
Xeno_TSVList_ByChr<-data.frame(matrix(unlist(Xeno_TSVList_ByChr), nrow=10, byrow=T))
rownames(Xeno_TSVList_ByChr)<-names(XENOGRAFT)
colnames(Xeno_TSVList_ByChr)<-chromList

#Calculate Number of Events by Chromosome for COMMON events
Common_TSVList_ByChr<-NULL
for(count in 1:length(COMMON))
{
  Sample<-COMMON[[count]]
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    ChrCounts[counter]<-(length(which(Sample$chr1 %in% chromList[counter]))+length(which(Sample$chr2 %in% chromList[counter])))/2
  }
  Common_TSVList_ByChr[[count]]<-ChrCounts
}
Common_TSVList_ByChr<-data.frame(matrix(unlist(Common_TSVList_ByChr), nrow=10, byrow=T))
rownames(Common_TSVList_ByChr)<-names(COMMON)
colnames(Common_TSVList_ByChr)<-chromList

# Now calculate Jaccard 
Jaccard_By_Chr<-NULL
for(sample in 1:length(COMMON))
{
  Xeno<-Xeno_TSVList_ByChr[sample,]
  Tumour<-Tumour_TSVList_ByChr[sample,]
  Common<-Common_TSVList_ByChr[sample,]
  ChrJaccard<-NULL
  for(count in 1:length(chromList))
  {
    ChrJaccard[count]<-Common[count]/((Tumour[count]+Xeno[count])-Common[count])
  }
  Jaccard_By_Chr[[sample]]<-ChrJaccard
}

Jaccard_By_Chr<-do.call(cbind,Jaccard_By_Chr)
Jaccard_By_Chr2<- data.frame(matrix(unlist(Jaccard_By_Chr), ncol=10, byrow=F))
Jaccard_By_Chr2[is.na(Jaccard_By_Chr2)]<-(-0.01)
Jaccard_By_Chr<-data.frame(Jaccard_By_Chr)

rownames(Jaccard_By_Chr2)<-chromList
colnames(Jaccard_By_Chr2)<-names(COMMON)

Jaccard_By_Chr2[Jaccard_By_Chr2<0]<-NA
colnames(Jaccard_By_Chr2)<-gsub(colnames(Jaccard_By_Chr2),pattern = "PCSI_",replacement = "")

write.csv(round(Jaccard_By_Chr2,3),file="Jaccard_By_Chr_PDAC.csv")

# Get an overall Score across 24 chromosomes
SampleScore<-NULL
for(sample in 1:ncol(Jaccard_By_Chr2))
{
  SampleTested<-Jaccard_By_Chr2[,sample]
  Empty<-length(which(is.na(SampleTested)))
  TotalCategories<-(24-Empty)
  Positive<-length(which(SampleTested>=0.6))
  SampleScore[sample]<-Positive/TotalCategories
}
names(SampleScore)<-colnames(Jaccard_By_Chr2)

message("Overall Concordance (Sc) Score: \n")
data.frame(SampleScore)

library(gplots)
pdf("Sc_Score_Plotting.pdf")
heatmap.2(as.matrix(t(rbind(SampleScore,SampleScore,Jaccard_By_Chr2))),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'))
dev.off()


#See which chromosomes have >=5 SV events in both tumour and PDX
test<-Tumour_TSVList_ByChr>=5 & Xeno_TSVList_ByChr>=5
test[test=="FALSE"]<-""
test

#Discordant Chromosomes
tester<-abs(Tumour_TSVList_ByChr - Xeno_TSVList_ByChr)>=10
tester[tester=="FALSE"]<-""
tester

library(gplots)
pdf("Jaccard_By_Chr_PDAC_Heatmap2.pdf",height = 10,width = 5)
heatmap.2(as.matrix((Jaccard_By_Chr2)),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),sepcolor = "black",tracecol = "white")
dev.off()
