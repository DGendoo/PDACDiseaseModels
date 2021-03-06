############################################################################################
## Deena M.A. Gendoo
## Last Update March 3, 2017
# FINAL_STEP_2_NumberOfVariantCalls.R
## Assess VCF files for Primary/PDX pairs from PDAC- Get Number of Variant Calls Per Sample
## IDEA: Use this as a metric for cutoff of Prim-PDX pairs for future analysis

############################################################################################
# libraries
library(VariantAnnotation)
library(scales)

ColTemplateDoubles<-c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462", "#b3de69", 
                      "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f","#1f78b4", "#e31a1c","#6a3d9a")
ColorListDoubles<-rep(ColTemplateDoubles, each=2) #30 colors of 15 pairs 

ColTemplateTriples<-c("#edf8b1","#7fcdbb","#2c7fb8") 
ColorListTriples<-rep(ColTemplateTriples,15) #30 colors of 15 pairs 

stringAsFactors=FALSE

load("PrimaryPDXInfo.RData")

#List of pairs and the variant counts
PrimPDXCounts<-NULL
for (count in seq(1,length(AlleleCounts),2))
{
  CommonVariants <- length(intersect(AlleleCounts[[count]]$VariantMutations,AlleleCounts[[count+1]]$VariantMutations))
  PrimaryVariants<-nrow(AlleleCounts[[count]])
  PDXVariants<-nrow(AlleleCounts[[count+1]])
  PrimPDXCounts[[(count+1)/2]]<-cbind(Primary=PrimaryVariants,PDX=PDXVariants,Common=CommonVariants)
}

#Rename VCF list & Counts
VCFs<-sapply(names(AlleleCounts),function(x) gsub(x,pattern = "_Pa_[P,X].*", replacement = ""))
VCFs<-sapply(VCFs,function(x) gsub(x,pattern = "_Lv_[M,X].*", replacement = ""))
SampleNames<-sort(unique(VCFs))
names(PrimPDXCounts)<-SampleNames

#Convert to a 'plottable' matrix
df <- data.frame(matrix(unlist(PrimPDXCounts), nrow=length(PrimPDXCounts), byrow=T))
rownames(df)<-SampleNames
colnames(df)<-c("Primary","PDX","Common")
df<-as.matrix(df)
df<-t(df)

write.csv(df,file="Number_SSM_Across_Pairs.csv")

#Show the common variants and the variants only specific to either PDX or specific to tumour
df[1,]<-df[1,]-df[3,]
df[2,]<-df[2,]-df[3,]
test<-rbind(df[3,],df[1,],df[2,])
rownames(test)<-c("Common","Primary","Xenograft")
pdf('NumberOfVariantCallsAcrossPrimPDXPairs3.pdf',height=7,width=7)
par(mar=c(12,4,6,2)) #bottom, L, top, R
barplot(test,beside = F,col=rev(c("#d53e4f","#ffffbf","#2b83ba")),
        cex.axis=0.6,las=2,font=2,ylab="Number of Variant Calls Per Sample",ylim=c(0,8000))
legend("topright",c("Common","Tumour","PDX"),fill = rev(c("#d53e4f","#ffffbf","#2b83ba")),horiz = F)
dev.off()
