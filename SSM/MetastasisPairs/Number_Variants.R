############################################################################################
## Deena M.A. Gendoo
## Last Update March 3, 2017

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

#Split the data based on the obvious outlier samples, versus the rest
#Outlier samples have the greatest number of variant calls!
HighOutliers<-df[,c("PCSI_0489")]
HighOutliers<-data.frame("PCSI_0489"=HighOutliers)
HighOutliers$Dummy<-c(50000,50000,50000)
HighOutliers$Dummy2<-c(50000,50000,50000)
HighOutliers$Dummy3<-c(50000,50000,50000)
HighOutliers$Dummy4<-c(50000,50000,50000)
HighOutliers<-as.matrix(HighOutliers)
LowNormals<-df[,c(!colnames(df)%in%c("PCSI_0489"))]



#PLOT EACH OF THE PAIRS ON A SEPERATE PAGE
pdf('NumberOfVariantCallsAcrossPrimPDXPairs.pdf',height=7,width=18)
par(mar=c(8,4,2,2),mfrow=c(1,2)) #bottom, L, top, R
barplot(LowNormals,beside = T,col=c("#ffffbf","#d53e4f","#2b83ba"),
        cex.axis=0.6,las=2,font=2,ylab="Number of Variant Calls Per Sample")
abline(h = 1000,col='red',lty=5,lwd=2)
legend("topleft",c("PDX","Tumour","Common"),fill = c("#d53e4f","#ffffbf","#2b83ba"),horiz = T)

barplot(HighOutliers,beside = T,col=c("#ffffbf","#d53e4f","#2b83ba"),xlab="PCSI_0489",
        cex.axis=0.6,las=2,ylab="Number of Variant Calls Per Sample",font=2)
abline(h = 1000,col='red',lty=5,lwd=2)
dev.off()



# !! Show the common variants and the variants only specific to either PDX or specific to tumour!!
df[1,]<-df[1,]-df[3,]
df[2,]<-df[2,]-df[3,]

#Split the data based on the obvious outlier samples, versus the rest
#Outlier samples have the greatest number of variant calls!
HighOutliers<-df[,c("PCSI_0489")]
HighOutliers<-data.frame("PCSI_0489"=HighOutliers)
HighOutliers$Dummy<-c(50000,50000,50000)
HighOutliers$Dummy2<-c(50000,50000,50000)
HighOutliers$Dummy3<-c(50000,50000,50000)
HighOutliers$Dummy4<-c(50000,50000,50000)
HighOutliers<-as.matrix(HighOutliers)
LowNormals<-df[,c(!colnames(df)%in%c("PCSI_0489"))]




LowNorm<-rbind(LowNormals[3,],LowNormals[1,],LowNormals[2,])
rownames(LowNorm)<-c("Common","Primary","Xenograft")

HighNorm<-rbind(HighOutliers[3,],HighOutliers[1,],HighOutliers[2,])
rownames(HighNorm)<-c("Common","Primary","Xenograft")

pdf('NumberOfVariantCallsAcrossPrimPDXPairs3.pdf',height=7,width=10)
par(mar=c(8,4,2,2),mfrow=c(1,2)) #bottom, L, top, R

barplot(LowNorm,beside = F,col=rev(c("#d53e4f","#ffffbf","#2b83ba")),
        cex.axis=0.6,las=2,font=2,ylab="Number of Variant Calls Per Sample",ylim=c(0,8000))
legend("topleft",c("Common","Tumour","PDX"),fill = rev(c("#d53e4f","#ffffbf","#2b83ba")),horiz = F)

barplot(HighNorm,beside = F,col=rev(c("#d53e4f","#ffffbf","#2b83ba")),
        cex.axis=0.6,las=2,font=2,ylab="Number of Variant Calls Per Sample",ylim=c(0,51000))

dev.off()

