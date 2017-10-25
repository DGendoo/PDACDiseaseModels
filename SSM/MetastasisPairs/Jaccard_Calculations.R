# Deena M.A. Gendoo
# Metrics look at each Prim-PDX pair
# Find all the genes that are commonly mutated for a given Prim-PDX pair
# Find all the genes that are commonly mutated for a given Prim-PDX pair, across different categories of mutations
# Find which genes are commonly mutated 

# http://rpackages.ianhowson.com/bioc/VariantAnnotation/man/VCF-class.html
library(VariantAnnotation)
library(VennDiagram)
library(scales)
library(gplots)

ColTemplateDoubles<-c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462", "#b3de69", 
                      "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f","#1f78b4", "#e31a1c","#6a3d9a")
ColorListDoubles<-rep(ColTemplateDoubles, each=2) #30 colors of 15 pairs 

ColTemplateTriples<-c("#edf8b1","#7fcdbb","#2c7fb8") 
ColorListTriples<-rep(ColTemplateTriples,15) #30 colors of 15 pairs 

CatColors<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','black')

#JaccardColors<-c('#ffffff','#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58')
#JaccardColors<-c('#ffffff','#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd')
JaccardColors<-c('#ffffff',rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837')))
stringAsFactors=FALSE

################################
# COMBINE ALL THE DATA FIRST
################################

#Load the parsed VCF files (including Annovar Annotation)
load("PrimaryPDXInfo.RData")

# Read all TSVs for Oncotator
TSVs<-sort(list.files(pattern = ".tsv"))
TSVList<-lapply(TSVs,function(x){read.delim(x,header=T,skip=2)})
TSVs<-gsub(TSVs,pattern = ".snp.maf.txt.tsv",replacement = "")
names(TSVList)<-TSVs

#COMBINE ALL THE INFO FROM THE VCFs AND ONCOTATOR
COMBO<-NULL
for(Sample in 1:length(TSVList))
{
  COMBO[[Sample]]<-cbind(AlleleCounts[[Sample]],TSVList[[Sample]])
}
names(COMBO)<-names(AlleleCounts)


################################
# Look at % Overlap of Mutations across all calls
################################
#List of pairs and the variant counts
PrimPDXCounts<-NULL
for (count in seq(1,length(AlleleCounts),2))
{
  CommonVariants <- length(intersect(COMBO[[count]]$VariantMutations,COMBO[[count+1]]$VariantMutations))
  PrimaryVariants<-nrow(COMBO[[count]])
  PDXVariants<-nrow(COMBO[[count+1]])
  PrimPDXCounts[[(count+1)/2]]<-cbind(Primary=PrimaryVariants,PDX=PDXVariants,Common=CommonVariants)
}

#Rename VCF list & Counts
VCFs<-sapply(names(COMBO),function(x) gsub(x,pattern = "_Pa_[P,X].*", replacement = ""))
VCFs<-sapply(VCFs,function(x) gsub(x,pattern = "_Lv_[M,X].*", replacement = ""))
SampleNames<-sort(unique(VCFs))
names(PrimPDXCounts)<-SampleNames

#Convert to a 'plottable' matrix
df <- data.frame(matrix(unlist(PrimPDXCounts), nrow=length(PrimPDXCounts), byrow=T))
rownames(df)<-SampleNames
colnames(df)<-c("Primary","PDX","Common")
df<-as.matrix(df)
df<-t(df)

# #PLOT EACH OF THE PAIRS ON A SEPERATE PAGE
# pdf('NumberOfVariantCallsAcrossPrimPDXPairs.pdf',height=7,width=15)
# par(mar=c(8,4,2,2)) #bottom, L, top, R
# barplot(df,beside = T,col=ColorListTriples,
#         cex.axis=0.6,las=2,font=2,ylab="Number of Variant Calls Per Sample")
# abline(h = 1000,col='red',lty=5,lwd=2)
# legend("topleft",c("Tumour","PDX","Common"),fill = ColorListTriples[1:3],horiz = T)
# 
# dev.off()

################################
# Look at % Overlap of Mutations across all calls
# Overlap by Mutation Category
################################
#Get list of all mutation categories: all types of variant classifications that have been called across all samples
lapply(COMBO,function(x){table((x)$variant_classification)})
lapply(COMBO,function(x){names(table((x)$variant_classification))})
Categories<-unique(unlist(lapply(COMBO,function(x){names(table((x)$variant_classification))})))

#Find all the mutation calls pertaining to each category, for every sample
ComboSub<-NULL
for (CategoryCount in 1:length(Categories))
{
  message("Category assessed: ",Categories[CategoryCount])
  ComboSub[[CategoryCount]]<-lapply(COMBO,function(x){x[x$variant_classification==Categories[CategoryCount],]})
}

#List of pairs and the variant counts
PrimPDXCountsMutCategory<-NULL
PrimPDXCountsMut<-NULL
MutProf<-NULL
for (CategoryCount in 1:length(Categories))
{
  message("Category assessed: ",Categories[CategoryCount])
  #Load the subbed VCFs for a given mutation category type
  MutProf<-ComboSub[[CategoryCount]]
  # Traverse the files to generate list of common, Prim-specific, PDX specific across all samples
  for (count in seq(1,length(AlleleCounts),2))
  {
    CommonVariants <- length(intersect(MutProf[[count]]$VariantMutations,MutProf[[count+1]]$VariantMutations))
    PrimaryVariants<-nrow(MutProf[[count]])
    PDXVariants<-nrow(MutProf[[count+1]])
    PrimPDXCountsMut[[(count+1)/2]]<-cbind(Primary=PrimaryVariants,PDX=PDXVariants,Common=CommonVariants)
  }
  #Rename VCF list & Counts
  VCFs<-sapply(names(MutProf),function(x) gsub(x,pattern = "_Pa_[P,X].*", replacement = ""))
  VCFs<-sapply(VCFs,function(x) gsub(x,pattern = "_Lv_[M,X].*", replacement = ""))
  SampleNames<-sort(unique(VCFs))
  names(PrimPDXCountsMut)<-SampleNames
  
  #Convert to a 'plottable' matrix
  df.temp <- data.frame(matrix(unlist(PrimPDXCountsMut), nrow=length(PrimPDXCountsMut), byrow=T))
  rownames(df.temp)<-SampleNames
  colnames(df.temp)<-c("Primary","PDX","Common")
  df.temp<-as.matrix(df.temp)
  df.temp<-t(df.temp)
  
  PrimPDXCountsMutCategory[[CategoryCount]]<-df.temp
}
names(PrimPDXCountsMutCategory)<-Categories

#Add the counts for the generic mutation calls across all categories
PrimPDXCountsMutCategory$ALL_MUTATIONS<-df

# #Select some categories to see the mutation overlaps across samples
# pdf("IdenticalMutationsOverlap.pdf",height = 7,width = 12,onefile = T)
# par(mar=c(8,4,2,4),mfrow=c(1,3)) #bottom, L, top, R
# for(count in c(13,5,6))
# {
#   MatTemp<-PrimPDXCountsMutCategory[[count]]
#   #MatTemp<-MatTemp[,!(colnames(MatTemp)=="PCSI_0489")]
#   barplot(MatTemp,col=ColorListTriples,beside = T,horiz = F,las=2,main = names(PrimPDXCountsMutCategory)[[count]])
# }
# dev.off()

#Calculate Jaccard Index - easier to see?
JaccardAcrossSamples<-NULL
for(count in 1:length(PrimPDXCountsMutCategory))
{
  TempMat<-PrimPDXCountsMutCategory[[count]]
  JaccardAcrossSamples[[count]]<-apply(TempMat,2,function(x){x["Common"]/((x["Primary"]+x["PDX"])-x["Common"])})
}
names(JaccardAcrossSamples)<-names(PrimPDXCountsMutCategory)

#Convert to a 'plottable' matrix
df.temp <- data.frame(matrix(unlist(JaccardAcrossSamples), nrow=length(JaccardAcrossSamples), byrow=T))
rownames(df.temp)<-names(JaccardAcrossSamples)
colnames(df.temp)<-colnames(PrimPDXCountsMutCategory[[count]])
JaccardAcrossSamples<-df.temp
JaccardAcrossSamples<-as.matrix(JaccardAcrossSamples)
df.temp<-t(df.temp)

#OPTIONAL: Restrict to only categories that matter? 
#JaccardAcrossSamples<-JaccardAcrossSamples[c("Missense_Mutation","Nonsense_Mutation","Splice_Site","3'UTR" ),]

#Reorder Rows???
JaccardAcrossSamples_Plot<-JaccardAcrossSamples[c(9,10,3,8,1,6,7,11,12,13,4,5,2,14,15),]
rownames(JaccardAcrossSamples_Plot)<-c("Missense","Nonsense","5'UTR","lincRNA","3'UTR","IGR","Intron",
                                       "RNA","Silent","Splite_Site","De_novo_Start_InFrame","De_novo_Start_OutOfFrame","5'Flank","Start_Codon_SNP","ALL MUTATIONS")

# Reorder Columns for samples
JaccardAcrossSamples_Plot<-JaccardAcrossSamples_Plot[,c(2:6,1)]

write.csv(JaccardAcrossSamples_Plot,file="JaccardScores.csv")

pdf("JaccardScores_MutationCategories_Reduced.pdf")
par(mar=c(8,4,2,20)) #bottom, L, top, R
par(mfrow=c(1,2))
heatmap.2( JaccardAcrossSamples_Plot,
           col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),main = "LIVER",
           trace = "none", key = NA,
           scale = c("none"),
           symbreaks = min(JaccardAcrossSamples, na.rm=TRUE),
           na.color="white",
           cexRow = 0.75, cexCol = 0.75,
           Colv = FALSE,Rowv = FALSE )

dev.off()
#c('#ffffff','#8c510a','#dfc27d','#c7eae5','#35978f')





