# Deena M.A. Gendoo
# November 1, 2016 and Updated February 24, 2017, final update June 23, 2017
# Parse Structural Variation (SV) data (in TSV format) for Paired Tumour-PDX samples
# FINAL_StructuralVariation_PDAC_Trios.R
###################################################################################################################
###################################################################################################################

#FLAGS!
EventTypeList=c("DEL","INV","DUP","TRA") #DEL or INV or DUP or TRA!

########################

# Read all TSVs containing structural variation
TSVs<-sort(list.files(pattern = ".annotatedSV.tsv"))
TSVList<-lapply(TSVs,function(x){read.csv(x)})
TSVs<-sapply(TSVs,function(x) gsub(x,pattern = ".annotatedSV.tsv", replacement = ""))
names(TSVList)<-TSVs
names(TSVList)<-gsub(x = names(TSVList),pattern = "_526.annotatedSV",replacement = "" )

#Rename TSV list to show the unique pair elements
TSVs<-sapply(TSVs,function(x){gsub(x,pattern = "_Pa_[O,P,X].*", replacement = "")})
SampleNames<-sort(unique(TSVs))

sampOrder<-c("PCSI_0590_Pa_P","PCSI_0590_Pa_X","PCSI_0590_Pa_O","PCSI_0592_Pa_P","PCSI_0592_Pa_X","PCSI_0592_Pa_O",
             "PCSI_0602_Pa_P","PCSI_0602_Pa_X","PCSI_0602_Pa_O","PCSI_0624_Pa_P","PCSI_0624_Pa_X","PCSI_0624_Pa_O",
             "PCSI_0642_Pa_P","PCSI_0642_Pa_X","PCSI_0642_Pa_O")
sampOrder2<-gsub(sampOrder,pattern = "_Pa",replacement = "")

###################################################################################################################
# Generate plot for total number of SV events per sample
###################################################################################################################
library(reshape2)

SV_Counts<-matrix(nrow=4,ncol=length(TSVs),data = NA)
rownames(SV_Counts)<-c("DEL","DUP","INV","TRA")
colnames(SV_Counts)<-gsub(names(TSVList),pattern = "_526.*",replacement = "")
colnames(SV_Counts)<-gsub(colnames(SV_Counts),pattern = "_Pa",replacement = "")
for (count in 1:length(TSVs))
{
  message("Working on Sample #:",count)
  SV_Counts[,count]<-melt(table(TSVList[[count]]$type))[,2]
}

SV_Counts<-SV_Counts[,sampOrder2]
ParamsTable<-as.matrix(do.call(cbind.data.frame, list(SV_Counts[,1:3],NA,SV_Counts[,4:6],NA,
                                                      SV_Counts[,7:9],NA,SV_Counts[,10:12],NA,SV_Counts[,13:15])))

# Stacked barplot of total SV counts across all samples
pdf("SV_Counts_Total_Across_Samples.pdf",onefile = F,width = 10,height = 7)
par(mar=c(7,6,2,2)+0.2)
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
rownames(SV_Counts_Chrom)<-gsub(rownames(SV_Counts_Chrom),pattern = "_Pa",replacement = "")

for (count in 1:length(TSVs))
{
  message("Working on Sample #:",count)
  Samplee<-TSVList[[count]]

  for(chromcount in 1:24)
  {
    SV_Counts_Chrom[count,chromcount]<-(length(which(Samplee$chr1 == colnames(SV_Counts_Chrom)[chromcount]))
                                        +length(which(Samplee$chr2 == colnames(SV_Counts_Chrom)[chromcount])))/2
  }
  
}

write.csv(SV_Counts_Chrom,"SV_Counts_PerChrom_Across_Samples_PDAC_TRIOS.csv")
