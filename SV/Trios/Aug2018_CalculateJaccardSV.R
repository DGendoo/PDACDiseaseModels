# Deena M.A. Gendoo
# August 3, 2018
# Parse Structural Variation (SV) data (in TSV format) for Paired Tumour-PDX sample
# Calculate Jaccard Index and Generate Corresponding Venn Diagrams
# Jaccard calculation here is done by looking at the NUMBER OF BASES tablulated in an event
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

# Read all TSVs containing structural variation events for ORGANOID
TSVs<-sort(list.files(pattern = "_Pa_O_526.annotatedSV.tsv"))
TSVList<-lapply(TSVs,function(x){read.csv(x,header = T)})
TSVs<-sapply(TSVs,function(x) gsub(x,pattern = "_Pa_O_526.annotatedSV.tsv", replacement = ""))
names(TSVList)<-TSVs
SampleNames<-sort(unique(TSVs))
ORGANOID<-TSVList
unlist(lapply(TSVList,nrow))
OrganoidSV<-melt(lapply(TSVList,nrow))$value

# Find all Common Events in a given tumour-pdx pair [Tumour vs PDX]
CC_TumourVsPDX<-CC_TumourVsPDO<-CC_PDXvsPDO<-NULL
COMMON_TvsX<-COMMON_TvsO<-COMMON_XvsO<-NULL

for(count in 1:length(TSVs))
{
  Tumour<-TUMOUR[[count]]
  PDX<-XENOGRAFT[[count]]
  PDO<-ORGANOID[[count]]
  
  #Calculate Common SV
  PDXsv<-paste(PDX$chr1,":",PDX$pos1,"-",PDX$chr2,":",PDX$pos2,sep="")
  Tumoursv<-paste(Tumour$chr1,":",Tumour$pos1,"-",Tumour$chr2,":",Tumour$pos2,sep="")
  PDOsv<-paste(PDO$chr1,":",PDO$pos1,"-",PDO$chr2,":",PDO$pos2,sep="")
  
  CC_TumourVsPDX[[count]]<-intersect(Tumoursv,PDXsv)
  CC_TumourVsPDO[[count]]<-intersect(Tumoursv,PDOsv)
  CC_PDXvsPDO[[count]]<-intersect(PDXsv,PDOsv)
  
  #Save Events that are common (taken from the Tumour file)
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

#Show differences
#Tumoursv[!(Tumoursv%in%CommonCalls[[count]])]
#PDXsv[!(PDXsv%in%CommonCalls[[count]])]

rm(CC_TumourVsPDX)
rm(CC_TumourVsPDO)
rm(CC_PDXvsPDO)

#Sanity Check
all(unlist(lapply(CommonSV_TvsX,nrow)) == CommonSV_TvsX)
all(unlist(lapply(CommonSV_TvsO,nrow)) == CommonSV_TvsO)
all(unlist(lapply(CommonSV_XvsO,nrow)) == CommonSV_XvsO)

# View All
SampleNames
TumourSV
XenoSV
OrganoidSV
CommonSV_TvsX
CommonSV_TvsO
CommonSV_XvsO

SVTable<-cbind(SampleNames,TumourSV,XenoSV,OrganoidSV,CommonSV_TvsX,CommonSV_TvsO,CommonSV_XvsO)


######################################################################################################
# CALCULATE EVENTS BY CHROMOSOME, AS OPPOSED TO ONLY GENERAL SCORE
######################################################################################################

#Calculate Number of Events by Chromosome for TUMOUR
Tumour_ByChr<-NULL
for(count in 1:length(TUMOUR))
{
  Sample<-TUMOUR[[count]]
  message("Sample: ",names(TUMOUR)[count])
  
  #Now by Chromosome
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    message("Chr: ",chromList[counter])
    Sample_DelInvDup<-Sample[which(Sample$chr1==chromList[counter] & Sample$chr2==chromList[counter]),]
    DelInvDupTotal<-sum(as.numeric(abs(Sample_DelInvDup$pos2 - Sample_DelInvDup$pos1)))
    
    Sample_Tra<-Sample[which(Sample$type=="TRA" & (Sample$chr1==chromList[counter] | Sample$chr2==chromList[counter])),]
    TraTotal=nrow(Sample_Tra)
    
    ChrCounts[counter]<-DelInvDupTotal+TraTotal
  }
  #Assign to a sample
  Tumour_ByChr[[count]]<-ChrCounts
}
Tumour_ByChr<-data.frame(matrix(unlist(Tumour_ByChr), nrow=5, byrow=T))
rownames(Tumour_ByChr)<-names(TUMOUR)
colnames(Tumour_ByChr)<-chromList

#Calculate Number of Events by Chromosome for XENOGRAFT
Xeno_ByChr<-NULL
for(count in 1:length(XENOGRAFT))
{
  Sample<-XENOGRAFT[[count]]
  message("Sample: ",names(XENOGRAFT)[count])
  #Now by Chromosome
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    message("Chr: ",chromList[counter])
    Sample_DelInvDup<-Sample[which(Sample$chr1==chromList[counter] & Sample$chr2==chromList[counter]),]
    DelInvDupTotal<-sum(as.numeric(abs(Sample_DelInvDup$pos2 - Sample_DelInvDup$pos1)))
    
    Sample_Tra<-Sample[which(Sample$type=="TRA" & (Sample$chr1==chromList[counter] | Sample$chr2==chromList[counter])),]
    TraTotal=nrow(Sample_Tra)
    
    ChrCounts[counter]<-DelInvDupTotal+TraTotal
  }
  
  #Assign to a sample
  Xeno_ByChr[[count]]<-ChrCounts
}
Xeno_ByChr<-data.frame(matrix(unlist(Xeno_ByChr), nrow=5, byrow=T))
rownames(Xeno_ByChr)<-names(XENOGRAFT)
colnames(Xeno_ByChr)<-chromList

#Calculate Number of Events by Chromosome for ORGANOID
Organoid_ByChr<-NULL
for(count in 1:length(ORGANOID))
{
  Sample<-ORGANOID[[count]]
  message("Sample: ",names(ORGANOID)[count])
  
  #Now by Chromosome
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    message("Chr: ",chromList[counter])
    Sample_DelInvDup<-Sample[which(Sample$chr1==chromList[counter] & Sample$chr2==chromList[counter]),]
    DelInvDupTotal<-sum(as.numeric(abs(Sample_DelInvDup$pos2 - Sample_DelInvDup$pos1)))
    
    Sample_Tra<-Sample[which(Sample$type=="TRA" & (Sample$chr1==chromList[counter] | Sample$chr2==chromList[counter])),]
    TraTotal=nrow(Sample_Tra)
    
    ChrCounts[counter]<-DelInvDupTotal+TraTotal
  }
  Organoid_ByChr[[count]]<-ChrCounts
}
Organoid_ByChr<-data.frame(matrix(unlist(Organoid_ByChr), nrow=5, byrow=T))
rownames(Organoid_ByChr)<-names(ORGANOID)
colnames(Organoid_ByChr)<-chromList


Common_ByChr_TvsX<-Common_ByChr_TvsO<-Common_ByChr_XvsO<-NULL

#Calculate Number of Events by Chromosome for COMMON events between TUMOUR-PDX
for(count in 1:length(COMMON_TvsX))
{
  Sample<-COMMON_TvsX[[count]]
  message("Sample: ",names(COMMON_TvsX)[count])
  
  #Now by Chromosome
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    message("Chr: ",chromList[counter])
    Sample_DelInvDup<-Sample[which(Sample$chr1==chromList[counter] & Sample$chr2==chromList[counter]),]
    DelInvDupTotal<-sum(as.numeric(abs(Sample_DelInvDup$pos2 - Sample_DelInvDup$pos1)))
    
    Sample_Tra<-Sample[which(Sample$type=="TRA" & (Sample$chr1==chromList[counter] | Sample$chr2==chromList[counter])),]
    TraTotal=nrow(Sample_Tra)
    
    ChrCounts[counter]<-DelInvDupTotal+TraTotal
  }
  
  Common_ByChr_TvsX[[count]]<-ChrCounts
}
Common_ByChr_TvsX<-data.frame(matrix(unlist(Common_ByChr_TvsX), nrow=5, byrow=T))
rownames(Common_ByChr_TvsX)<-names(COMMON_TvsX)
colnames(Common_ByChr_TvsX)<-chromList

#Calculate Number of Events by Chromosome for COMMON events between TUMOUR-PDO
for(count in 1:length(COMMON_TvsO))
{
  Sample<-COMMON_TvsO[[count]]
  message("Sample: ",names(COMMON_TvsO)[count])
  
  #Now by Chromosome
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    message("Chr: ",chromList[counter])
    Sample_DelInvDup<-Sample[which(Sample$chr1==chromList[counter] & Sample$chr2==chromList[counter]),]
    DelInvDupTotal<-sum(as.numeric(abs(Sample_DelInvDup$pos2 - Sample_DelInvDup$pos1)))
    
    Sample_Tra<-Sample[which(Sample$type=="TRA" & (Sample$chr1==chromList[counter] | Sample$chr2==chromList[counter])),]
    TraTotal=nrow(Sample_Tra)
    
    ChrCounts[counter]<-DelInvDupTotal+TraTotal
  }
  Common_ByChr_TvsO[[count]]<-ChrCounts
}
Common_ByChr_TvsO<-data.frame(matrix(unlist(Common_ByChr_TvsO), nrow=5, byrow=T))
rownames(Common_ByChr_TvsO)<-names(COMMON_TvsO)
colnames(Common_ByChr_TvsO)<-chromList

#Calculate Number of Events by Chromosome for COMMON events between PDX-PDO
for(count in 1:length(COMMON_XvsO))
{
  Sample<-COMMON_XvsO[[count]]
  message("Sample: ",names(COMMON_XvsO)[count])
  
  #Now by Chromosome
  ChrCounts<-NULL
  for(counter in 1:length(chromList))
  {
    message("Chr: ",chromList[counter])
    Sample_DelInvDup<-Sample[which(Sample$chr1==chromList[counter] & Sample$chr2==chromList[counter]),]
    DelInvDupTotal<-sum(as.numeric(abs(Sample_DelInvDup$pos2 - Sample_DelInvDup$pos1)))
    
    Sample_Tra<-Sample[which(Sample$type=="TRA" & (Sample$chr1==chromList[counter] | Sample$chr2==chromList[counter])),]
    TraTotal=nrow(Sample_Tra)
    
    ChrCounts[counter]<-DelInvDupTotal+TraTotal
  }
  Common_ByChr_XvsO[[count]]<-ChrCounts
}
Common_ByChr_XvsO<-data.frame(matrix(unlist(Common_ByChr_XvsO), nrow=5, byrow=T))
rownames(Common_ByChr_XvsO)<-names(COMMON_XvsO)
colnames(Common_ByChr_XvsO)<-chromList


Jaccard_By_Chr_TvsX<-Jaccard_By_Chr_TvsO<-Jaccard_By_Chr_XvsO<-NULL

# Now calculate Jaccard By Chromosome for TUMOUR-PDX
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

# Now calculate Jaccard By Chromosome for TUMOUR-PDO
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

# Now calculate Jaccard By Chromosome for PDX-PDO
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

write.csv(round(ChromSVTable,4),file="Aug2018_Table_JaccardByChr_Trios.csv")
save.image("Aug2018_Jaccard_Calculations_PDAC_Trios.RData")

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



library(gplots)
pdf("Aug2018_Jaccard_By_Chr_PDAC_Trios_Heatmap1.pdf",height = 7,width = 2)
heatmap.2(as.matrix((ChromSVTable[,1:3])),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),cexRow = 0.6,cexCol = 0.6)
dev.off()

pdf("Aug2018_Jaccard_By_Chr_PDAC_Trios_Heatmap2.pdf",height = 7,width = 2)
heatmap.2(as.matrix((ChromSVTable[,4:6])),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),cexRow = 0.6,cexCol = 0.6)
dev.off()

pdf("Aug2018_Jaccard_By_Chr_PDAC_Trios_Heatmap3.pdf",height = 7,width = 2)
heatmap.2(as.matrix((ChromSVTable[,7:9])),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),cexRow = 0.6,cexCol = 0.6)
dev.off()

pdf("Aug2018_Jaccard_By_Chr_PDAC_Trios_Heatmap4.pdf",height = 7,width = 2)
heatmap.2(as.matrix((ChromSVTable[,10:12])),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),cexRow = 0.6,cexCol = 0.6)
dev.off()

pdf("Aug2018_Jaccard_By_Chr_PDAC_Trios_Heatmap5.pdf",height = 7,width = 2)
heatmap.2(as.matrix((ChromSVTable[,13:15])),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),cexRow = 0.6,cexCol = 0.6)
dev.off()
#col = c('#ffffff','#ffffcc','#a1dab4','#41b6c4','#2c7fb8','#253494')
#col = c('#ffffff','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0')
# col = c('#ffffff','#8c510a','#dfc27d','#c7eae5','#35978f')
# col = c('#ffffff','#8c510a','#dfc27d','#80cdc1','#01665e')

save.image("Aug2018_Jaccard_Calculations_PDAC_Trios.RData")
