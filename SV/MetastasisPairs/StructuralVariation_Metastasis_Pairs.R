# Deena M.A. Gendoo
# November 1, 2016, and final update on June 23, 2017
# Parse Structural Variation (SV) data (in TSV format) for Paired Tumour-PDX samples
# FINAL_StructuralVariation_LIVER_PairsNews.R

###################################################################################################################

# Read all TSVs containing structural variation
TSVs<-sort(list.files(pattern = ".annotatedSV.tsv"))
TSVList<-lapply(TSVs,function(x){read.csv(x)})
TSVs<-sapply(TSVs,function(x) gsub(x,pattern = ".annotatedSV.tsv", replacement = ""))
names(TSVList)<-TSVs
names(TSVList)<-gsub(x = names(TSVList),pattern = "_526",replacement = "" )

#Rename TSV list to show the unique pair elements
TSVs<-sapply(TSVs,function(x){gsub(x,pattern = "_Pa_[P,X].*", replacement = "")})
TSVs<-sapply(TSVs,function(x){gsub(x,pattern = "_Lv_[M,X].*", replacement = "")})
SampleNames<-sort(unique(TSVs))

###################################################################################################################
# Generate plot for total number of SV events per sample
#Generate a loop and assess each pair individually
###################################################################################################################
library(reshape2)

SV_Counts<-matrix(nrow=4,ncol=length(TSVs),data = NA)
rownames(SV_Counts)<-c("DEL","DUP","INV","TRA")
colnames(SV_Counts)<-gsub(names(TSVList),pattern = "_526.*",replacement = "")
colnames(SV_Counts)<-gsub(colnames(SV_Counts),pattern = "_Lv",replacement = "")
for (count in 1:length(TSVs))
{
  message("Working on Sample #:",count)
  SV_Counts["DEL",count]<-length(which(TSVList[[count]]$type == "DEL"))
  SV_Counts["DUP",count]<-length(which(TSVList[[count]]$type == "DUP"))
  SV_Counts["INV",count]<-length(which(TSVList[[count]]$type == "INV"))
  SV_Counts["TRA",count]<-length(which(TSVList[[count]]$type == "TRA"))
}


#Quick formatting for figures!
ParamsTable<-as.matrix(do.call(cbind.data.frame, list(SV_Counts[,1:2],NA,SV_Counts[,3:4],NA,SV_Counts[,5:6],NA,
                                                   SV_Counts[,7:8],NA,SV_Counts[,9:10],NA,SV_Counts[,11:12])))

# Stacked barplot of total SV counts across all samples
pdf("SV_Counts_Total_Across_Samples.pdf",onefile = F,width = 10,height = 7)
par(mar=c(10,6,2,2)+0.2)
barplot(ParamsTable,col=c("#d73027","#2166ac","#66bd63","#542788"),las=2,ylab="Number of SV Events",names.arg = colnames(ParamsTable))
legend("topleft",c("Deletion","Duplication","Inversion","Translocation"),fill=c("#d73027","#2166ac","#66bd63","#542788"))
dev.off()


write.csv(SV_Counts,file="Total_SV_CountsPerSample.csv")
###################################################################################################################
# Generate plot for total number of SV events per chromosome, for each sample (To highlight chromothripsis)
###################################################################################################################

SV_Counts_Chrom<-matrix(nrow=length(TSVs),ncol=24,data = NA)
colnames(SV_Counts_Chrom)<-c(paste("chr",rep(1:22),sep = ""),"chrX","chrY")
rownames(SV_Counts_Chrom)<-gsub(names(TSVList),pattern = "_526.*",replacement = "")
rownames(SV_Counts_Chrom)<-gsub(rownames(SV_Counts_Chrom),pattern = "_Lv",replacement = "")

for (count in 1:length(TSVs))
{
  message("Working on Sample #:",count)
  Samplee<-TSVList[[count]]
  Samplee$chr1<-factor(as.character(lapply(Samplee$chr1,function(x){strsplit(as.character(x),"_")[[1]][1]})))
  Samplee$chr2<-factor(as.character(lapply(Samplee$chr2,function(x){strsplit(as.character(x),"_")[[1]][1]})))
  
  for(chromcount in 1:24)
  {
    SV_Counts_Chrom[count,chromcount]<-(length(which(Samplee$chr1 == colnames(SV_Counts_Chrom)[chromcount]))
                                        +length(which(Samplee$chr2 == colnames(SV_Counts_Chrom)[chromcount])))/2
  }
  
}

write.csv(SV_Counts_Chrom,"SV_Counts_PerChrom_Across_Samples_LIVER.csv")
