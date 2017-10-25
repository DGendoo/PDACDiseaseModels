# Deena M.A. Gendoo
# Metrics look at each Prim-PDX pair
# Find all the genes that are commonly mutated for a given Prim-PDX pair
# Find all the genes that are commonly mutated for a given Prim-PDX pair, across different categories of mutations
# Find which genes are commonly mutated 
# FINAL_STEP_4_MetricsScript.R

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

JaccardColors<-c('#ffffff',rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837')))
stringAsFactors=FALSE

################################
# COMBINE ALL THE DATA FIRST
################################

load("PrimaryPDXInfo.RData")

TSVs<-sort(list.files(pattern = ".tsv"))
TSVList<-lapply(TSVs,function(x){read.delim(x,header=T,skip=2)})
TSVs<-gsub(TSVs,pattern = ".snp.maf.txt.tsv",replacement = "")
names(TSVList)<-TSVs

COMBO<-NULL
for(Sample in 1:length(TSVList))
{
  COMBO[[Sample]]<-cbind(AlleleCounts[[Sample]],TSVList[[Sample]])
}
names(COMBO)<-names(AlleleCounts)


################################
# Look at % Overlap of Mutations across all calls
################################
PrimPDXCounts<-NULL
for (count in seq(1,length(AlleleCounts),2))
{
  CommonVariants <- length(intersect(COMBO[[count]]$VariantMutations,COMBO[[count+1]]$VariantMutations))
  PrimaryVariants<-nrow(COMBO[[count]])
  PDXVariants<-nrow(COMBO[[count+1]])
  PrimPDXCounts[[(count+1)/2]]<-cbind(Primary=PrimaryVariants,PDX=PDXVariants,Common=CommonVariants)
}

VCFs<-sapply(names(COMBO),function(x) gsub(x,pattern = "_Pa_[P,X].*", replacement = ""))
VCFs<-sapply(VCFs,function(x) gsub(x,pattern = "_Lv_[M,X].*", replacement = ""))
SampleNames<-sort(unique(VCFs))
names(PrimPDXCounts)<-SampleNames

df <- data.frame(matrix(unlist(PrimPDXCounts), nrow=length(PrimPDXCounts), byrow=T))
rownames(df)<-SampleNames
colnames(df)<-c("Primary","PDX","Common")
df<-as.matrix(df)
df<-t(df)

################################
# Look at % Overlap of Mutations across all calls
# Overlap by Mutation Category
################################
lapply(COMBO,function(x){table((x)$variant_classification)})
lapply(COMBO,function(x){names(table((x)$variant_classification))})
Categories<-unique(unlist(lapply(COMBO,function(x){names(table((x)$variant_classification))})))

ComboSub<-NULL
for (CategoryCount in 1:length(Categories))
{
  message("Category assessed: ",Categories[CategoryCount])
  ComboSub[[CategoryCount]]<-lapply(COMBO,function(x){x[x$variant_classification==Categories[CategoryCount],]})
}

PrimPDXCountsMutCategory<-NULL
PrimPDXCountsMut<-NULL
MutProf<-NULL
for (CategoryCount in 1:length(Categories))
{
  message("Category assessed: ",Categories[CategoryCount])
  MutProf<-ComboSub[[CategoryCount]]
  for (count in seq(1,length(AlleleCounts),2))
  {
    CommonVariants <- length(intersect(MutProf[[count]]$VariantMutations,MutProf[[count+1]]$VariantMutations))
    PrimaryVariants<-nrow(MutProf[[count]])
    PDXVariants<-nrow(MutProf[[count+1]])
    PrimPDXCountsMut[[(count+1)/2]]<-cbind(Primary=PrimaryVariants,PDX=PDXVariants,Common=CommonVariants)
  }
  VCFs<-sapply(names(MutProf),function(x) gsub(x,pattern = "_Pa_[P,X].*", replacement = ""))
  VCFs<-sapply(VCFs,function(x) gsub(x,pattern = "_Lv_[M,X].*", replacement = ""))
  SampleNames<-sort(unique(VCFs))
  names(PrimPDXCountsMut)<-SampleNames
  
  df.temp <- data.frame(matrix(unlist(PrimPDXCountsMut), nrow=length(PrimPDXCountsMut), byrow=T))
  rownames(df.temp)<-SampleNames
  colnames(df.temp)<-c("Primary","PDX","Common")
  df.temp<-as.matrix(df.temp)
  df.temp<-t(df.temp)
  
  PrimPDXCountsMutCategory[[CategoryCount]]<-df.temp
}
names(PrimPDXCountsMutCategory)<-Categories

PrimPDXCountsMutCategory$ALL_MUTATIONS<-df

JaccardAcrossSamples<-NULL
for(count in 1:length(PrimPDXCountsMutCategory))
{
  TempMat<-PrimPDXCountsMutCategory[[count]]
  JaccardAcrossSamples[[count]]<-apply(TempMat,2,function(x){x["Common"]/((x["Primary"]+x["PDX"])-x["Common"])})
}
names(JaccardAcrossSamples)<-names(PrimPDXCountsMutCategory)

df.temp <- data.frame(matrix(unlist(JaccardAcrossSamples), nrow=length(JaccardAcrossSamples), byrow=T))
rownames(df.temp)<-names(JaccardAcrossSamples)
colnames(df.temp)<-colnames(PrimPDXCountsMutCategory[[count]])
JaccardAcrossSamples<-df.temp
JaccardAcrossSamples<-as.matrix(JaccardAcrossSamples)
df.temp<-t(df.temp)

JaccardAcrossSamples_Plot<-JaccardAcrossSamples[c(6,7,2,5,1,3,4,8,9,10,11,12,13),]
rownames(JaccardAcrossSamples_Plot)<-c("Missense","Nonsense","5'UTR","lincRNA","3'UTR","IGR","Intron",
                                       "RNA","Silent","Splite_Site","De_novo_Start_InFrame","De_novo_Start_OutOfFrame","ALL MUTATIONS")

pdf("JaccardScores_MutationCategories_Reduced.pdf")
par(mar=c(6,8,6,25)) #bottom, L, top, R
par(mfrow=c(1,2))
heatmap.2( JaccardAcrossSamples_Plot[1:12,],
           col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),main = "PDAC",
           trace = "none", key = NA,
           scale = c("none"),
           symbreaks = min(JaccardAcrossSamples, na.rm=TRUE),
           na.color="white",
           cexRow = 0.75, cexCol = 0.75,
           Colv = FALSE,Rowv = FALSE )

heatmap.2( JaccardAcrossSamples_Plot,
           col = c('#543005','#8c510a','#dfc27d','#c7eae5','#35978f'),main = "PDAC",
           trace = "none", key = NA,
           scale = c("none"),
           symbreaks = min(JaccardAcrossSamples, na.rm=TRUE),
           na.color="white",
           cexRow = 0.75, cexCol = 0.75,
           Colv = FALSE,Rowv = FALSE )
dev.off()



