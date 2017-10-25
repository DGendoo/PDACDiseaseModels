############################################################################################
## Deena M.A. Gendoo
## Last Updated June 25, 2017

## Assess VCF files for Primary/PDX pairs from PDAC- Get Number of Variant Calls Per Sample
## IDEA: Use this as a metric for cutoff of Prim-PDX pairs for future analysis

############################################################################################

# libraries
# http://rpackages.ianhowson.com/bioc/VariantAnnotation/man/VCF-class.html
library(VariantAnnotation)
library(VennDiagram)
library(scales)

ColTemplateDoubles<-c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462", "#b3de69", 
                      "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f","#1f78b4", "#e31a1c","#6a3d9a")
ColorListDoubles<-rep(ColTemplateDoubles, each=2) #30 colors of 15 pairs 

ColTemplateTriples<-c("#edf8b1","#7fcdbb","#2c7fb8") 
ColorListTriples<-rep(ColTemplateTriples,15) #30 colors of 15 pairs 

stringAsFactors=FALSE

load("PrimaryPDXInfo.RData")

#List of pairs and the variant counts - Assess the trio
PrimPDXCounts<-NULL
for (count in seq(1,length(AlleleCounts),3))
{
  
  OrganoidVariants<-nrow(AlleleCounts[[count]])
  PrimaryVariants<-nrow(AlleleCounts[[count+1]])
  PDXVariants<-nrow(AlleleCounts[[count+2]])
  
  Common_PDO_Pri <- length(intersect(AlleleCounts[[count]]$VariantMutations,AlleleCounts[[count+1]]$VariantMutations))
  Common_PDO_PDX <- length(intersect(AlleleCounts[[count]]$VariantMutations,AlleleCounts[[count+2]]$VariantMutations))
  Common_Pri_PDX <- length(intersect(AlleleCounts[[count+1]]$VariantMutations,AlleleCounts[[count+2]]$VariantMutations))

  Common_Pri_PDX_2 <- intersect(AlleleCounts[[count+1]]$VariantMutations,AlleleCounts[[count+2]]$VariantMutations)
  PDOmut<-AlleleCounts[[count]]$VariantMutations
  
  Common_All <- length(unique(intersect(Common_Pri_PDX_2,PDOmut)))

  PrimPDXCounts[[(count+2)/3]]<-cbind(Primary=PrimaryVariants,PDX=PDXVariants,PDO=OrganoidVariants,
                                      Common_Pri_PDX=Common_Pri_PDX,Common_PDO_PDX=Common_PDO_PDX,
                                      Common_PDO_Pri=Common_PDO_Pri,Common_All=Common_All)
}

#Rename VCF list & Counts
VCFs<-sapply(names(AlleleCounts),function(x) gsub(x,pattern = "_Pa_[P,X,O].*", replacement = ""))
SampleNames<-sort(unique(VCFs))
names(PrimPDXCounts)<-SampleNames

#Convert to a 'plottable' matrix
df <- data.frame(matrix(unlist(PrimPDXCounts), nrow=length(PrimPDXCounts), byrow=T))
rownames(df)<-SampleNames
colnames(df)<-c("Primary","PDX","PDO","Common_Primary_PDX","Common_PDO_PDX","Common_PDO_Primary","Common_All")
df<-as.matrix(df)
df<-t(df)

write.csv(df,file="Number_SSM_Across_Trios.csv")

