# Deena M.A. Gendoo
# November 1, 2016, and final update June 23, 2017
# Parse Structural Variation (SV) data (in TSV format) for Paired Tumour-PDX samples
# CalculateJaccardSV.R
# Calculate Jaccard Index 
###################################################################################################################
###################################################################################################################

library(reshape2)

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

# Read all TSVs containing structural variation events for ORGANOID
TSVs<-sort(list.files(pattern = "_Pa_O_526.annotatedSV.tsv"))
TSVList<-lapply(TSVs,function(x){read.csv(x,header = T)})
TSVs<-sapply(TSVs,function(x) gsub(x,pattern = "_Pa_O_526.annotatedSV.tsv", replacement = ""))
names(TSVList)<-TSVs
SampleNames<-sort(unique(TSVs))
ORGANOID<-TSVList
unlist(lapply(TSVList,nrow))
OrganoidSV<-melt(lapply(TSVList,nrow))$value

CC_TumourVsPDX<-CC_TumourVsPDO<-CC_PDXvsPDO<-NULL
COMMON_TvsX<-COMMON_TvsO<-COMMON_XvsO<-NULL

for(count in 1:length(TSVs))
{
  Tumour<-TUMOUR[[count]]
  PDX<-XENOGRAFT[[count]]
  PDO<-ORGANOID[[count]]
  
  PDXsv<-paste(PDX$chr1,":",PDX$pos1,"-",PDX$chr2,":",PDX$pos2,sep="")
  Tumoursv<-paste(Tumour$chr1,":",Tumour$pos1,"-",Tumour$chr2,":",Tumour$pos2,sep="")
  PDOsv<-paste(PDO$chr1,":",PDO$pos1,"-",PDO$chr2,":",PDO$pos2,sep="")
  
  CC_TumourVsPDX[[count]]<-intersect(Tumoursv,PDXsv)
  CC_TumourVsPDO[[count]]<-intersect(Tumoursv,PDOsv)
  CC_PDXvsPDO[[count]]<-intersect(PDXsv,PDOsv)
  
  COMMON_TvsX[[count]]<-Tumour[(Tumour$chr1%in%PDX$chr1 & Tumour$chr2%in%PDX$chr2 & Tumour$pos1%in%PDX$pos1 & Tumour$pos2%in%PDX$pos2),]
  COMMON_TvsO[[count]]<-Tumour[(Tumour$chr1%in%PDO$chr1 & Tumour$chr2%in%PDO$chr2 & Tumour$pos1%in%PDO$pos1 & Tumour$pos2%in%PDO$pos2),]
  COMMON_XvsO[[count]]<-PDX[(PDX$chr1%in%PDO$chr1 & PDX$chr2%in%PDO$chr2 & PDX$pos1%in%PDO$pos1 & PDX$pos2%in%PDO$pos2),]
  
}
names(COMMON_TvsX)<-names(TSVList)
names(COMMON_TvsO)<-names(TSVList)
names(COMMON_XvsO)<-names(TSVList)

CommonSV_TvsX<-unlist(lapply(CC_TumourVsPDX,length))
CommonSV_TvsO<-unlist(lapply(CC_TumourVsPDO,length))
CommonSV_XvsO<-unlist(lapply(CC_PDXvsPDO,length))

rm(CC_TumourVsPDX)
rm(CC_TumourVsPDO)
rm(CC_PDXvsPDO)

SVTable<-cbind(SampleNames,TumourSV,XenoSV,OrganoidSV,CommonSV_TvsX,CommonSV_TvsO,CommonSV_XvsO)
#write.csv(SVTable,"Table_SV_Counts_Trios.csv")

######################################################################################################
# CALCULATE EVENTS BY CHROMOSOME, AS OPPOSED TO ONLY GENERAL SCORE
######################################################################################################

Tumour_ByChr<-NULL
for(count in 1:5)
{
  Sample<-TUMOUR[[count]]
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    ChrCounts[counter]<-(length(which(Sample$chr1 %in% chromList[counter]))+length(which(Sample$chr2 %in% chromList[counter])))/2
  }
  Tumour_ByChr[[count]]<-ChrCounts
}
Tumour_ByChr<-data.frame(matrix(unlist(Tumour_ByChr), nrow=5, byrow=T))
rownames(Tumour_ByChr)<-names(TUMOUR)
colnames(Tumour_ByChr)<-chromList

Xeno_ByChr<-NULL
for(count in 1:5)
{
  Sample<-XENOGRAFT[[count]]
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    ChrCounts[counter]<-(length(which(Sample$chr1 %in% chromList[counter]))+length(which(Sample$chr2 %in% chromList[counter])))/2
  }
  Xeno_ByChr[[count]]<-ChrCounts
}
Xeno_ByChr<-data.frame(matrix(unlist(Xeno_ByChr), nrow=5, byrow=T))
rownames(Xeno_ByChr)<-names(XENOGRAFT)
colnames(Xeno_ByChr)<-chromList

Organoid_ByChr<-NULL
for(count in 1:5)
{
  Sample<-ORGANOID[[count]]
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    ChrCounts[counter]<-(length(which(Sample$chr1 %in% chromList[counter]))+length(which(Sample$chr2 %in% chromList[counter])))/2
  }
  Organoid_ByChr[[count]]<-ChrCounts
}
Organoid_ByChr<-data.frame(matrix(unlist(Organoid_ByChr), nrow=5, byrow=T))
rownames(Organoid_ByChr)<-names(ORGANOID)
colnames(Organoid_ByChr)<-chromList


Common_ByChr_TvsX<-Common_ByChr_TvsO<-Common_ByChr_XvsO<-NULL

for(count in 1:length(COMMON_TvsX))
{
  Sample<-COMMON_TvsX[[count]]
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    ChrCounts[counter]<-(length(which(Sample$chr1 %in% chromList[counter]))+length(which(Sample$chr2 %in% chromList[counter])))/2
  }
  Common_ByChr_TvsX[[count]]<-ChrCounts
}
Common_ByChr_TvsX<-data.frame(matrix(unlist(Common_ByChr_TvsX), nrow=5, byrow=T))
rownames(Common_ByChr_TvsX)<-names(COMMON_TvsX)
colnames(Common_ByChr_TvsX)<-chromList

for(count in 1:length(COMMON_TvsO))
{
  Sample<-COMMON_TvsO[[count]]
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    ChrCounts[counter]<-(length(which(Sample$chr1 %in% chromList[counter]))+length(which(Sample$chr2 %in% chromList[counter])))/2
  }
  Common_ByChr_TvsO[[count]]<-ChrCounts
}
Common_ByChr_TvsO<-data.frame(matrix(unlist(Common_ByChr_TvsO), nrow=5, byrow=T))
rownames(Common_ByChr_TvsO)<-names(COMMON_TvsO)
colnames(Common_ByChr_TvsO)<-chromList

for(count in 1:length(COMMON_XvsO))
{
  Sample<-COMMON_XvsO[[count]]
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    ChrCounts[counter]<-(length(which(Sample$chr1 %in% chromList[counter]))+length(which(Sample$chr2 %in% chromList[counter])))/2
  }
  Common_ByChr_XvsO[[count]]<-ChrCounts
}
Common_ByChr_XvsO<-data.frame(matrix(unlist(Common_ByChr_XvsO), nrow=5, byrow=T))
rownames(Common_ByChr_XvsO)<-names(COMMON_XvsO)
colnames(Common_ByChr_XvsO)<-chromList

#Jaccard Calculations
Jaccard_By_Chr_TvsX<-Jaccard_By_Chr_TvsO<-Jaccard_By_Chr_XvsO<-NULL

for(sample in 1:5)
{
  Xeno<-Xeno_ByChr[sample,]
  Tumour<-Tumour_ByChr[sample,]
  Common<-Common_ByChr_TvsX[sample,]
  
  ChrJaccard<-NULL
  for(count in 1:length(chromList))
  {
    ChrJaccard[count]<-Common[count]/((Tumour[count]+Xeno[count])-Common[count])
  }
  Jaccard_By_Chr_TvsX[[sample]]<-ChrJaccard
}
Jaccard_By_Chr_TvsX<-do.call(cbind,Jaccard_By_Chr_TvsX)

for(sample in 1:5)
{
  Tumour<-Tumour_ByChr[sample,]
  Organo<-Organoid_ByChr[sample,]
  Common<-Common_ByChr_TvsO[sample,]
  
  ChrJaccard<-NULL
  for(count in 1:length(chromList))
  {
    ChrJaccard[count]<-Common[count]/((Tumour[count]+Organo[count])-Common[count])
  }
  Jaccard_By_Chr_TvsO[[sample]]<-ChrJaccard
}
Jaccard_By_Chr_TvsO<-do.call(cbind,Jaccard_By_Chr_TvsO)

for(sample in 1:5)
{
  Xeno<-Xeno_ByChr[sample,]
  Organo<-Organoid_ByChr[sample,]
  Common<-Common_ByChr_XvsO[sample,]
  
  ChrJaccard<-NULL
  for(count in 1:length(chromList))
  {
    ChrJaccard[count]<-Common[count]/((Xeno[count]+Organo[count])-Common[count])
  }
  Jaccard_By_Chr_XvsO[[sample]]<-ChrJaccard
}
Jaccard_By_Chr_XvsO<-do.call(cbind,Jaccard_By_Chr_XvsO)


Jaccard_By_Chr2<- data.frame(matrix(unlist(Jaccard_By_Chr_TvsX), ncol=5, byrow=F))
Jaccard_By_Chr_TvsX<-data.frame(Jaccard_By_Chr2)

Jaccard_By_Chr2<- data.frame(matrix(unlist(Jaccard_By_Chr_TvsO), ncol=5, byrow=F))
Jaccard_By_Chr_TvsO<-data.frame(Jaccard_By_Chr2)

Jaccard_By_Chr2<- data.frame(matrix(unlist(Jaccard_By_Chr_XvsO), ncol=5, byrow=F))
Jaccard_By_Chr_XvsO<-data.frame(Jaccard_By_Chr2)

rownames(Jaccard_By_Chr_XvsO)<-rownames(Jaccard_By_Chr_TvsO)<-rownames(Jaccard_By_Chr_TvsX)<-chromList
colnames(Jaccard_By_Chr_XvsO)<-colnames(Jaccard_By_Chr_TvsO)<-colnames(Jaccard_By_Chr_TvsX)<-names(TSVList)

ChromSVTable<-NULL
for (count in 1:5)
{
 tablesv<-cbind(Jaccard_By_Chr_TvsX[,count],Jaccard_By_Chr_TvsO[,count],Jaccard_By_Chr_XvsO[,count])
 rownames(tablesv)<-chromList
 colnames(tablesv)<-c(paste(names(TSVList)[count],"_T_vs_PDX",sep=""),
                      paste(names(TSVList)[count],"_T_vs_PDO",sep=""),
                      paste(names(TSVList)[count],"_PDX_vs_PDO",sep=""))
 ChromSVTable<-cbind(ChromSVTable,tablesv)
}

write.csv(ChromSVTable,file="Table_JaccardByChr_Trios.csv")

colnames(ChromSVTable)<-gsub(colnames(ChromSVTable),pattern = "PCSI_",replacement = "")

# Get an overall Score across 24 chromosomes for each of the disease model comparisons
SampleScore<-NULL
for(sample in 1:ncol(ChromSVTable))
{
  SampleTested<-ChromSVTable[,sample]
  Empty<-length(which(is.na(SampleTested)))
  TotalCategories<-(24-Empty)
  Positive<-length(which(SampleTested>=0.6))
  SampleScore[sample]<-Positive/TotalCategories
}
names(SampleScore)<-colnames(ChromSVTable)

message("Overall Concordance (Sc) Score: \n")
data.frame(SampleScore)

library(gplots)
pdf("Jaccard_By_Chr_PDAC_Trios_Heatmap1.pdf",height = 7,width = 2)
heatmap.2(as.matrix((ChromSVTable[,1:3])),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),cexRow = 0.6,cexCol = 0.6)
dev.off()

pdf("Jaccard_By_Chr_PDAC_Trios_Heatmap2.pdf",height = 7,width = 2)
heatmap.2(as.matrix((ChromSVTable[,4:6])),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),cexRow = 0.6,cexCol = 0.6)
dev.off()

pdf("Jaccard_By_Chr_PDAC_Trios_Heatmap3.pdf",height = 7,width = 2)
heatmap.2(as.matrix((ChromSVTable[,7:9])),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),cexRow = 0.6,cexCol = 0.6)
dev.off()

pdf("Jaccard_By_Chr_PDAC_Trios_Heatmap4.pdf",height = 7,width = 2)
heatmap.2(as.matrix((ChromSVTable[,10:12])),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),cexRow = 0.6,cexCol = 0.6)
dev.off()

pdf("Jaccard_By_Chr_PDAC_Trios_Heatmap5.pdf",height = 7,width = 2)
heatmap.2(as.matrix((ChromSVTable[,13:15])),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),cexRow = 0.6,cexCol = 0.6)
dev.off()

#See which chromosomes have >=5 SV events in both tumour and PDX
test<-Tumour_ByChr>=5 & Xeno_ByChr>=5
test[test=="FALSE"]<-""
test

#See which chromosomes have >=5 SV events in both tumour and PDO
test2<-Tumour_ByChr>=5 & Organoid_ByChr>=5
test2[test2=="FALSE"]<-""
test2

#See which chromosomes have >=5 SV events in both PDX and PDO
test3<-Xeno_ByChr>=5 & Organoid_ByChr>=5
test3[test3=="FALSE"]<-""
test3

#Discordant Chromosomes between Tumour-PDX
tester<-abs(Tumour_ByChr - Xeno_ByChr)>=10
tester[tester=="FALSE"]<-""
tester

#Discordant Chromosomes between Tumour-PDO
tester2<-abs(Tumour_ByChr - Organoid_ByChr)>=10
tester2[tester2=="FALSE"]<-""
tester2

#Discordant Chromosomes between PDX-PDO
tester3<-abs(Xeno_ByChr - Organoid_ByChr)>=10
tester3[tester3=="FALSE"]<-""
tester3
