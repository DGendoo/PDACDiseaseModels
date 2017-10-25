############################################################################################
## Deena M.A. Gendoo
## Last Updated June 25, 2017

## Assess VCF files for Primary/PDX pairs from PDAC - Get VAF (Variant Allele Freq) between Prim-PDX pairs
# Ergo, plot the allele frequencies for all the mutants, to show the frequencies across a given Prim-PDX pair, 
# as well as highlight mutations of interest

############################################################################################

# FIRST TIME INSTALLATION ONLY
  #source("http://bioconductor.org/biocLite.R")
  #biocLite("VariantAnnotation")

# libraries
# http://rpackages.ianhowson.com/bioc/VariantAnnotation/man/VCF-class.html
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

# NOTES:
# The format of the geno file is in the format of DP:FDP:SDP:SUBDP:AU:CU:GU:TU
# DP = read depth (total number of reads for VAF and MAF analyses)
# AU/CU/GU/TU = ref allele called. These are actually a list because they are comma-split (can have multiples). 
# For the comma split, we have only used (for now) the first instance of the reads
# The AU/CU/GU/TU readouts will produce 'NA' for non-standard geno formats (DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50)

# Read all VCFs using Variant Annotation
VCFs<-sort(list.files(pattern = ".vcf"))
VcfList<-lapply(VCFs,function(x){readVcf(x,genome="hg19")})
VCFs<-sapply(VCFs,function(x) gsub(x,pattern = ".final.vcf", replacement = ""))
names(VcfList)<-VCFs
names(VcfList)<-gsub(x = names(VcfList),pattern = "_526",replacement = "" )

# NOTES:
# Major difficulty in extracting the "ANNOVAR" annotation from the INFO column using VariantAnnotation
# Instead, use vcfR package to parse this part

# Read all VCFs using vcfR
VCFs<-sort(list.files(pattern = ".vcf"))
VcfListUsingVCFR<-lapply(VCFs,function(x){read.vcfR(x)})
VCFs<-sapply(VCFs,function(x) gsub(x,pattern = ".final.vcf", replacement = ""))
names(VcfListUsingVCFR)<-VCFs

# Get alleleic counts for total number of reads (DP) plus the counts for each allele
AlleleCounts<-NULL
for (count in 1:length(VcfList))
{
  message("Processing file: ",names(VcfList[count]))
  
  # Nucleotide of the Reference Allele
  REF<-VcfList[[count]]@fixed$REF
  REF<-as.character(REF)
  
  # Nucleotide of the ALT allele
  # Note that some ALTs are comma-seperated! Take the first one
  test<-(VcfList[[count]]@fixed$ALT)
  test2<-lapply(test, as.character)
  ALT<-sapply(test2, `[`,1) #split the comma string, take first element!
  #ALT<-unlist(VcfList[[count]]@fixed$ALT)
  #ALT<-as.character(ALT)
  
  #Total Number of Reads
  ReadDepth<-geno(VcfList[[count]])$DP[,2]
  
  # Alternative allele frequencies
  # Note [,,1] stands for 'TUMOUR' sample. [,2] is the 'Tumour' column
  GU<-geno(VcfList[[count]])$GU[,,1][,2]
  AU<-geno(VcfList[[count]])$AU[,,1][,2]
  CU<-geno(VcfList[[count]])$CU[,,1][,2]
  TU<-geno(VcfList[[count]])$TU[,,1][,2]
  
  #Chromosome Location
  VariantMetadata<-data.frame(VcfList[[count]]@rowRanges)
  Chromosome<-VariantMetadata$seqnames
  PositionOfVariant<-VariantMetadata$start
  
  #Get Annotation for Variant
  VariantMutations<-paste(Chromosome,PositionOfVariant,REF,ALT,sep = ":")
  #VariantMutations<-gsub(rownames(Collated),pattern = c("_"),replacement = ":")
  #VariantMutations<-gsub(VariantMutations,pattern = c("/"),replacement = ":")
  
  # Extract Gene that has been affected
  GeneName<-as.character(as.vector(VcfList[[count]]@info$GENE))
  GeneName[GeneName=="character(0)"]<-"non_coding"
    #GeneName<-as.character(VcfList[[count]]@info$GENE)
    #GeneName[is.na(GeneName)]<-"non_coding"
  
  # Get ANNOVAR Annotations (From INFO Column of the VCF)
  # USE vcfR package for this!!
  tester<-list(matrix(VcfListUsingVCFR[[count]]@fix[,"INFO"]))
  tester<-tester[[1]]
  test <-sapply(X = tester,FUN = function(SNV){strsplit(strsplit(SNV,'ANNOVAR=')[[1]][2],';')[[1]][1]})
  AnnovarCall<-as.character(matrix(lapply(test,FUN = function(snv2){strsplit(snv2,split = "\\,")[[1]][1]}))) #Annovar Call
  AnnovarGene<-as.character(matrix(lapply(test,FUN = function(snv3){strsplit(snv3,split = "\\,")[[1]][2]}))) #Annovar Gene

  #Collate everything so far
  Collated<-data.frame(cbind(REF,ALT,ReadDepth,AU,CU,GU,TU,GeneName,
                             VariantMutations,AnnovarCall,AnnovarGene),stringsAsFactors = FALSE)
  #Rename Annovar Genes to 'non-coding' if they are NOT EXONIC
  #Collated$AnnovarGene[(Collated$AnnovarCall!="exonic")]<-"non_coding"

  #FILTERING: Remove all the variants for which Allele frequencies are 'NA'!
  Collated<-Collated[!is.na(Collated[,4]),]
  
  # Calculate allele frequencies for somatic SNVs
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

#Rename VCF list
names(AlleleCounts)<-names(VcfList)

save(AlleleCounts,file="PrimaryPDXInfo.RData")

