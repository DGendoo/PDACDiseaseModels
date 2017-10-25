############################################################################################
## Deena M.A. Gendoo
## Last Updated March 3, 2017
## Assess VCF files for Primary/PDX pairs from PDAC - Get VAF (Variant Allele Freq) between Prim-PDX pairs
#FINAL_STEP_1_VAF_PrimPDXpairs.R

############################################################################################
library(VariantAnnotation)
library(VennDiagram)
library(scales)
library(vcfR)

ColTemplateDoubles<-c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462", "#b3de69", 
               "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f","#1f78b4", "#e31a1c","#6a3d9a")
ColorListDoubles<-rep(ColTemplateDoubles, each=2) #30 colors of 15 pairs 

ColTemplateTriples<-c("#edf8b1","#7fcdbb","#2c7fb8") 
ColorListTriples<-rep(ColTemplateTriples,15) #30 colors of 15 pairs 

stringAsFactors=FALSE

######################################################################
# Look at Variant Allele Frequency (VAF) for a given pair
######################################################################

# Read all VCFs 
VCFs<-sort(list.files(pattern = ".final.vcf"))
VcfList<-lapply(VCFs,function(x){readVcf(x,genome="hg19")})
VCFs<-sapply(VCFs,function(x) gsub(x,pattern = ".final.vcf", replacement = ""))
names(VcfList)<-VCFs
names(VcfList)<-gsub(x = names(VcfList),pattern = "_526",replacement = "" )

VCFs<-sort(list.files(pattern = ".final.vcf"))
VcfListUsingVCFR<-lapply(VCFs,function(x){read.vcfR(x)})
VCFs<-sapply(VCFs,function(x) gsub(x,pattern = ".final.vcf", replacement = ""))
names(VcfListUsingVCFR)<-VCFs

# Get alleleic counts for total number of reads (DP) plus the counts for each allele
AlleleCounts<-NULL
for (count in 1:length(VcfList))
{
  message("Processing file: ",names(VcfList[count]))
  
  REF<-VcfList[[count]]@fixed$REF
  REF<-as.character(REF)
  
  test<-(VcfList[[count]]@fixed$ALT)
  test2<-lapply(test, as.character)
  ALT<-sapply(test2, `[`,1) #split the comma string, take first element!
  
  ReadDepth<-geno(VcfList[[count]])$DP[,2]
  
  GU<-geno(VcfList[[count]])$GU[,,1][,2]
  AU<-geno(VcfList[[count]])$AU[,,1][,2]
  CU<-geno(VcfList[[count]])$CU[,,1][,2]
  TU<-geno(VcfList[[count]])$TU[,,1][,2]
  
  VariantMetadata<-data.frame(VcfList[[count]]@rowRanges)
  Chromosome<-VariantMetadata$seqnames
  PositionOfVariant<-VariantMetadata$start
  
  VariantMutations<-paste(Chromosome,PositionOfVariant,REF,ALT,sep = ":")
  
  GeneName<-as.character(as.vector(VcfList[[count]]@info$GENE))
  GeneName[GeneName=="character(0)"]<-"non_coding"
  
  tester<-list(matrix(VcfListUsingVCFR[[count]]@fix[,"INFO"]))
  tester<-tester[[1]]
  test <-sapply(X = tester,FUN = function(SNV){strsplit(strsplit(SNV,'ANNOVAR=')[[1]][2],';')[[1]][1]})
  AnnovarCall<-as.character(matrix(lapply(test,FUN = function(snv2){strsplit(snv2,split = "\\,")[[1]][1]}))) #Annovar Call
  AnnovarGene<-as.character(matrix(lapply(test,FUN = function(snv3){strsplit(snv3,split = "\\,")[[1]][2]}))) #Annovar Gene

  Collated<-data.frame(cbind(REF,ALT,ReadDepth,AU,CU,GU,TU,GeneName,
                             VariantMutations,AnnovarCall,AnnovarGene),stringsAsFactors = FALSE)
  
  Collated<-Collated[!is.na(Collated[,4]),]
  
  SNVfreqList<-NULL
  for(variant in 1:nrow(Collated))
  {
    if(Collated[variant,"ALT"] =="A"){SNVfreq<-as.numeric(Collated[variant,"AU"])/as.numeric(Collated[variant,"ReadDepth"])}
    if(Collated[variant,"ALT"] =="C"){SNVfreq<-as.numeric(Collated[variant,"CU"])/as.numeric(Collated[variant,"ReadDepth"])}
    if(Collated[variant,"ALT"] =="G"){SNVfreq<-as.numeric(Collated[variant,"GU"])/as.numeric(Collated[variant,"ReadDepth"])}
    if(Collated[variant,"ALT"] =="T"){SNVfreq<-as.numeric(Collated[variant,"TU"])/as.numeric(Collated[variant,"ReadDepth"])}
    SNVfreqList[variant]<-SNVfreq
  }
  
  AlleleCounts[[count]]<-data.frame(cbind(Collated,SNVfreq=SNVfreqList),stringsAsFactors = FALSE)
}

names(AlleleCounts)<-names(VcfList)

save(AlleleCounts,file="PrimaryPDXInfo.RData")

