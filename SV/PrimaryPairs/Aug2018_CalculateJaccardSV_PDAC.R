# Deena M.A. Gendoo
# August 3, 2018
# Parse Structural Variation (SV) data (in TSV format) for Paired Tumour-PDX samples
# Calculate Jaccard Index and Generate Corresponding Venn Diagrams
# Jaccard calculation here is done by looking at the NUMBER OF BASES tablulated in an event
###################################################################################################################
###################################################################################################################

library(reshape2)
library(pheatmap)


chromList<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
             "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")

######################################################################################################
######################################################################################################
######################################################################################################
# FIND ALL EVENTS PERTINENT TO TUMOUR, PDX, AS WELL AS COMMON EVENTS FOR A TUMOUR-PDX PAIR
######################################################################################################
######################################################################################################
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

#Sanity Check
all(unlist(lapply(COMMON,nrow)) == CommonSV)


# View All
SampleNames
CommonSV
TumourSV
XenoSV

######################################################################################################
######################################################################################################
######################################################################################################
# CALCULATE EVENTS BY CHROMOSOME, AS OPPOSED TO ONLY GENERAL SCORE
######################################################################################################
######################################################################################################
######################################################################################################

# Calculate Number of bases found, for all events, by Chromosome, for TUMOUR
# For DEL/INV/DUP, sum all the bases involved in the event
# For TRA, count each instance as 1 bp (ie, only 1 bp in that chromosomeA got affected by a given TRA event)

Tumour_TSVList_ByChr<-NULL
for(count in 1:length(TUMOUR))
{
  Sample<-TUMOUR[[count]]
  
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
  Common_TSVList_ByChr[[count]]<-ChrCounts
}
Common_TSVList_ByChr<-data.frame(matrix(unlist(Common_TSVList_ByChr), nrow=10, byrow=T))
rownames(Common_TSVList_ByChr)<-names(COMMON)
colnames(Common_TSVList_ByChr)<-chromList

# View All counts for chr-specific SV events
SampleNames
Common_TSVList_ByChr
Tumour_TSVList_ByChr
Xeno_TSVList_ByChr

# Now calculate Jaccard By Chromosome!
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

pdf(file="Plot_Jaccard_By_Chr.pdf",width=25,height = 7)
barplot(as.matrix(Jaccard_By_Chr2),beside = T,col = rainbow(24))
#abline(h = 1,col="red")
dev.off()

Jaccard_By_Chr2[Jaccard_By_Chr2<0]<-NA
colnames(Jaccard_By_Chr2)<-gsub(colnames(Jaccard_By_Chr2),pattern = "PCSI_",replacement = "")

library(gplots)
pdf("Aug2018_Jaccard_By_Chr_PDAC_Heatmap.pdf",height = 10,width = 10)
heatmap.2(as.matrix((Jaccard_By_Chr2)),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'))
dev.off()
#'#543005','#8c510a','#dfc27d','#c7eae5','#35978f' #Winning Scheme
#c('#543005','#8c510a','#dfc27d','#80cdc1','#01665e')) 
#col = c('#ffffff','#ffffcc','#a1dab4','#41b6c4','#2c7fb8','#253494')
#col = c('#ffffff','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0')
# col = c('#ffffff','#8c510a','#dfc27d','#c7eae5','#35978f')
# col = c('#ffffff','#8c510a','#dfc27d','#80cdc1','#01665e')

write.csv(round(Jaccard_By_Chr2,4),file="Aug2018_Jaccard_By_Chr_PDAC.csv")

#View by chromosome!
barplot(as.matrix(t(Jaccard_By_Chr2)),beside = T,col = rainbow(10))

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


pdf("Aug2018_Sc_Score_Plotting.pdf")
heatmap.2(as.matrix(t(rbind(SampleScore,SampleScore,Jaccard_By_Chr2))),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'))
dev.off()

library(gplots)
pdf("Aug2018_Jaccard_By_Chr_PDAC_Heatmap2.pdf",height = 10,width = 5)
heatmap.2(as.matrix((Jaccard_By_Chr2)),Rowv = F,Colv = F,trace = "none",na.color = "white",
          col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),sepcolor = "black",tracecol = "white")
dev.off()

#col = c('#fff

save.image("Aug2018_Jaccard_Calculations_PDAC.RData")


######################################################################################################
######################################################################################################
######################################################################################################
# CALCULATE EVENTS BY CATEGORY OF THE STRUCTURAL VARIATION. Genome-wide score still
# ERGO, GENOME-WIDE SCORE for DELETIONS, another for Inversions, etc.. Get a total of 4 scores per sample
######################################################################################################
######################################################################################################
######################################################################################################

# For DEL/INV/DUP, sum all the bases involved in the event
# For TRA, count each instance as 1 bp (ie, only 1 bp in that chromosomeA got affected by a given TRA event)

# Calculate Number of bases found, for all events, for TUMOUR
SVtypes<-c("DEL","DUP","INV","TRA")

Tumour_TSVList_ByCategory<-NULL
for(count in 1:length(TUMOUR))
{
  Sample<-TUMOUR[[count]]
  message("Sample: ",names(TUMOUR)[count])
  
  #Now by Category
  ChrCounts<-NULL
  for(category in 1:4)#categories: DEL, DUP,INV, TRA
  {
    SV_TYPE<-SVtypes[category] #Get the type assessed
    
    message("Category: ",SVtypes[category])
    SampleCat<-Sample[which(Sample$type==SV_TYPE),]
    SampleSize<-sum(as.numeric(abs(SampleCat$pos2 - SampleCat$pos1)))

    if(SV_TYPE=="TRA"){ SampleSize<-nrow(Sample)} #Just count the number of translocations that happened

    ChrCounts[category]<-SampleSize
  }
  #Assign to a sample
  Tumour_TSVList_ByCategory[[count]]<-ChrCounts
}
Tumour_TSVList_ByCategory<-data.frame(matrix(unlist(Tumour_TSVList_ByCategory), nrow=10, byrow=T))
rownames(Tumour_TSVList_ByCategory)<-names(TUMOUR)
colnames(Tumour_TSVList_ByCategory)<-SVtypes


# Calculate Number of bases found, for all events, for XENOGRAFT
SVtypes<-c("DEL","DUP","INV","TRA")

Xeno_TSVList_ByCategory<-NULL
for(count in 1:length(TUMOUR))
{
  Sample<-XENOGRAFT[[count]]
  message("Sample: ",names(XENOGRAFT)[count])
  
  #Now by Category
  ChrCounts<-NULL
  for(category in 1:4)#categories: DEL, DUP,INV, TRA
  {
    SV_TYPE<-SVtypes[category] #Get the type assessed
    
    message("Category: ",SVtypes[category])
    SampleCat<-Sample[which(Sample$type==SV_TYPE),]
    SampleSize<-sum(as.numeric(abs(SampleCat$pos2 - SampleCat$pos1)))
    
    if(SV_TYPE=="TRA"){ SampleSize<-nrow(Sample)}
    
    ChrCounts[category]<-SampleSize
  }
  #Assign to a sample
  Xeno_TSVList_ByCategory[[count]]<-ChrCounts
}
Xeno_TSVList_ByCategory<-data.frame(matrix(unlist(Xeno_TSVList_ByCategory), nrow=10, byrow=T))
rownames(Xeno_TSVList_ByCategory)<-names(XENOGRAFT)
colnames(Xeno_TSVList_ByCategory)<-SVtypes


# Calculate Number of bases found, for all events, for COMMON events
SVtypes<-c("DEL","DUP","INV","TRA")

Common_TSVList_ByCategory<-NULL
for(count in 1:length(TUMOUR))
{
  Sample<-COMMON[[count]]
  message("Sample: ",names(COMMON)[count])
  
  #Now by Category
  ChrCounts<-NULL
  for(category in 1:4)#categories: DEL, INV, DUP, TRA
  {
    SV_TYPE<-SVtypes[category] #Get the type assessed
    
    message("Category: ",SVtypes[category])
    SampleCat<-Sample[which(Sample$type==SV_TYPE),]
    SampleSize<-sum(as.numeric(abs(SampleCat$pos2 - SampleCat$pos1)))
    
    if(SV_TYPE=="TRA"){ SampleSize<-nrow(Sample)}
    
    ChrCounts[category]<-SampleSize
  }
  #Assign to a sample
  Common_TSVList_ByCategory[[count]]<-ChrCounts
}
Common_TSVList_ByCategory<-data.frame(matrix(unlist(Common_TSVList_ByCategory), nrow=10, byrow=T))
rownames(Common_TSVList_ByCategory)<-names(COMMON)
colnames(Common_TSVList_ByCategory)<-SVtypes


# Now calculate Jaccard By Category!
Jaccard_By_Cat<-NULL
for(sample in 1:length(COMMON))
{
  Xeno<-Xeno_TSVList_ByCategory[sample,]
  Tumour<-Tumour_TSVList_ByCategory[sample,]
  Common<-Common_TSVList_ByCategory[sample,]
  
  ChrJaccard<-NULL
  for(count in 1:length(SVtypes))
  {
    ChrJaccard[count]<-Common[count]/((Tumour[count]+Xeno[count])-Common[count])
  }
  
  Jaccard_By_Cat[[sample]]<-ChrJaccard
}

Jaccard_By_Cat<-do.call(cbind,Jaccard_By_Cat)
Jaccard_By_Cat2<- data.frame(matrix(unlist(Jaccard_By_Cat), ncol=10, byrow=F))
#Jaccard_By_Cat2[is.na(Jaccard_By_Cat2)]<-(-0.01)
rownames(Jaccard_By_Cat2)<-SVtypes
colnames(Jaccard_By_Cat2)<-names(COMMON)

write.csv(Jaccard_By_Cat2,file="Aug2018_JaccardBySVCategory.csv")
#Now produce barplots for the distribution of events, genome-wide, for each sample
Jaccard_By_Cat2<-as.matrix(Jaccard_By_Cat2)
pdf("Aug2018_JaccardBySVCategory.pdf",width=12,height = 5)
barplot(as.matrix(Jaccard_By_Cat2),col = c("#d73027","#74add1","#91cf60","#542788"),beside = T,ylab = "Jaccard Score")
legend("topleft", legend = c("DEL","DUP","INV","TRA"), fill = c("#d73027","#74add1","#91cf60","#542788"))
dev.off()


######################################################################################################
######################################################################################################
######################################################################################################
# FOCUS ON CHROMOTHRYPSIS: CALCULATE EVENTS BY CATEGORY OF THE STRUCTURAL VARIATION, FOR EVERY CHROMOSOME
# Ergo, SCORE for DELETIONS, another for Inversions, etc.. Get a total of 4 scores x 24 chromosomes, per sample
######################################################################################################
######################################################################################################
######################################################################################################
# For DEL/INV/DUP, sum all the bases involved in the event
# For TRA, count each instance as 1 bp (ie, only 1 bp in that chromosomeA got affected by a given TRA event)

chromList<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
             "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
SVtypes<-c("DEL","DUP","INV","TRA")


##### DO THE ANALYSIS FOR THE TUMOUR SAMPLES ######
Tumour_TSVList<-NULL

#First Loop: By Sample
for(samplecount in 1:length(TUMOUR))
{
  Sample<-TUMOUR[[samplecount]]
  message("Sample: ",names(TUMOUR)[samplecount])
  
  #Second Loop: Now by Chromosome
  ChrCounts<-NULL
  
  for(chrcounter in 1:length(chromList))
  {
    message("Chr: ",chromList[chrcounter])
    
    #Reduce to a Given Chromosome!
    Sample_DelInvDup<-Sample[which(Sample$chr1==chromList[chrcounter] & Sample$chr2==chromList[chrcounter]),]
    Sample_Tra<-Sample[which(Sample$type=="TRA" & (Sample$chr1==chromList[chrcounter] | Sample$chr2==chromList[chrcounter])),]
    
    # Third Loop: Now Assess by SV Category
    CatCount<-NULL
    for(categoryCount in 1:4)
    {
      SV_TYPE<-SVtypes[categoryCount] #Get the SV type assessed
      message("SV: ",SV_TYPE)
      
      SampleCat<-Sample_DelInvDup[which(Sample_DelInvDup$type==SV_TYPE),]
      SampleSize<-sum(as.numeric(abs(SampleCat$pos2 - SampleCat$pos1)))
      
      if(SV_TYPE=="TRA"){ SampleSize<-nrow(Sample_Tra)} #Just count the number of translocations that happened
      
      #Save counts of SV categories
      CatCount[categoryCount]<-SampleSize
      
    } #end of Third Loop
    
    
    #Second Loop: Save counts to chromosome 
    ChrCounts[[chrcounter]]<-CatCount 
  }
  
  #Third Loop: Save all chromosome info to a given sample
  ChrCounts<-data.frame(matrix(unlist(ChrCounts), nrow=24, byrow=T))
  rownames(ChrCounts)<-chromList
  colnames(ChrCounts)<-SVtypes
  
  #Assign to a sample, all scores for all chromosomes
  Tumour_TSVList[[samplecount]]<-ChrCounts
}
names(Tumour_TSVList)<-names(TUMOUR)


##### DO THE ANALYSIS FOR THE XENOGRAFT SAMPLES ######
Xeno_TSVList<-NULL

#First Loop: By Sample
for(samplecount in 1:length(XENOGRAFT))
{
  Sample<-XENOGRAFT[[samplecount]]
  message("Sample: ",names(XENOGRAFT)[samplecount])
  
  #Second Loop: Now by Chromosome
  ChrCounts<-NULL
  
  for(chrcounter in 1:length(chromList))
  {
    message("Chr: ",chromList[chrcounter])
    
    #Reduce to all SV events for a Given Chromosome!
    Sample_DelInvDup<-Sample[which(Sample$chr1==chromList[chrcounter] & Sample$chr2==chromList[chrcounter]),]
    Sample_Tra<-Sample[which(Sample$type=="TRA" & (Sample$chr1==chromList[chrcounter] | Sample$chr2==chromList[chrcounter])),]
    
    # Third Loop: Now Assess by SV Category
    CatCount<-NULL
    for(categoryCount in 1:4)
    {
      SV_TYPE<-SVtypes[categoryCount] #Get the SV type assessed
      message("SV: ",SV_TYPE)
      
      SampleCat<-Sample_DelInvDup[which(Sample_DelInvDup$type==SV_TYPE),]
      SampleSize<-sum(as.numeric(abs(SampleCat$pos2 - SampleCat$pos1)))
      
      if(SV_TYPE=="TRA"){ SampleSize<-nrow(Sample_Tra)} #Just count the number of translocations that happened
      
      #Save counts of SV categories
      CatCount[categoryCount]<-SampleSize
      
    } #end of Third Loop
    
    
    #Second Loop: Save counts to chromosome 
    ChrCounts[[chrcounter]]<-CatCount 
  }
  
  #Third Loop: Save all chromosome info to a given sample
  ChrCounts<-data.frame(matrix(unlist(ChrCounts), nrow=24, byrow=T))
  rownames(ChrCounts)<-chromList
  colnames(ChrCounts)<-SVtypes
  
  #Assign to a sample, all scores for all chromosomes
  Xeno_TSVList[[samplecount]]<-ChrCounts
}
names(Xeno_TSVList)<-names(XENOGRAFT)


##### DO THE ANALYSIS FOR THE COMMON SAMPLES ######
Common_TSVList<-NULL

#First Loop: By Sample
for(samplecount in 1:length(COMMON))
{
  Sample<-COMMON[[samplecount]]
  message("Sample: ",names(COMMON)[samplecount])
  
  #Second Loop: Now by Chromosome
  ChrCounts<-NULL
  
  for(chrcounter in 1:length(chromList))
  {
    message("Chr: ",chromList[chrcounter])
    
    #Reduce to all SV events for a Given Chromosome!
    Sample_DelInvDup<-Sample[which(Sample$chr1==chromList[chrcounter] & Sample$chr2==chromList[chrcounter]),]
    Sample_Tra<-Sample[which(Sample$type=="TRA" & (Sample$chr1==chromList[chrcounter] | Sample$chr2==chromList[chrcounter])),]
    
    # Third Loop: Now Assess by SV Category
    CatCount<-NULL
    for(categoryCount in 1:4)
    {
      SV_TYPE<-SVtypes[categoryCount] #Get the SV type assessed
      message("SV: ",SV_TYPE)
      
      SampleCat<-Sample_DelInvDup[which(Sample_DelInvDup$type==SV_TYPE),]
      SampleSize<-sum(as.numeric(abs(SampleCat$pos2 - SampleCat$pos1)))
      
      if(SV_TYPE=="TRA"){ SampleSize<-nrow(Sample_Tra)} #Just count the number of translocations that happened
      
      #Save counts of SV categories
      CatCount[categoryCount]<-SampleSize
      
    } #end of Third Loop
    
    
    #Second Loop: Save counts to chromosome 
    ChrCounts[[chrcounter]]<-CatCount 
  }
  
  #Third Loop: Save all chromosome info to a given sample
  ChrCounts<-data.frame(matrix(unlist(ChrCounts), nrow=24, byrow=T))
  rownames(ChrCounts)<-chromList
  colnames(ChrCounts)<-SVtypes
  
  #Assign to a sample, all scores for all chromosomes
  Common_TSVList[[samplecount]]<-ChrCounts
}
names(Common_TSVList)<-names(COMMON)


# NOW CALCULATE JACCARD SCORES, by Chromosome, by Category!
Jaccard_Totals<-NULL
for(sample in 1:length(COMMON))
{
  #Get the sample
  message("SAMPLE: ",names(COMMON)[sample])
  Xeno<-Xeno_TSVList[[sample]]
  Tumour<-Tumour_TSVList[[sample]]
  Common<-Common_TSVList[[sample]]
  
  ChrJaccard<-NULL
  for(chrom in 1:length(chromList))
  {
    #Get the chromosome
    XenoChrom<-Xeno[chrom,]
    TumourChrom<-Tumour[chrom,]
    CommonChrom<-Common[chrom,]
    message("Chr: ",chromList[chrom])
    
    SVJacs<-NULL
    for(count in 1:length(SVtypes))
    {
      SVJacs[count]<-CommonChrom[count]/((TumourChrom[count]+XenoChrom[count])-CommonChrom[count])
    }
    
    #Assign back to chromosome
    ChrJaccard[[chrom]]<-unlist(SVJacs)
  }
  
  #Collect all Jaccard Scores, for all Chromosomes
  ChrJaccard<-data.frame(matrix(unlist(ChrJaccard), nrow=24, byrow=T))
  rownames(ChrJaccard)<-chromList
  colnames(ChrJaccard)<-SVtypes
  # Assign back to Sample
  Jaccard_Totals[[sample]]<-ChrJaccard
}

names(Jaccard_Totals)<-names(COMMON)

#Start Plotting!
pdf("Aug2018_Jaccard_ByChrBySV.pdf",onefile = T,width = 10,height = 12)
for(sampName in 1:length(Jaccard_Totals))
{
  CurrentSample<-Jaccard_Totals[[sampName]]
  
  par(mfrow=c(4,6))
  for(tally in 1:24)
  {
    chromosome<-chromList[tally]
    barplot(as.matrix(t(round(CurrentSample[chromosome,],4))),beside = T,axes = FALSE,
            col = c("#d73027","#74add1","#91cf60","#542788"),ylim = c(0,1),main = names(Jaccard_Totals)[sampName])
    axis(2,at=seq(0,1,0.1))
  }
}
dev.off()

# for(i in 1:length(Jaccard_Totals))
# {
#   colnames(Jaccard_Totals[[i]])<-paste(names(Jaccard_Totals)[[i]],colnames(Jaccard_Totals[[i]]),sep="_")
# }
Jaccard_Totals<-(do.call(cbind,Jaccard_Totals))

write.csv(Jaccard_Totals,file="Aug2018_JaccardBySVCat_ByChr.csv")
