# Deena M.A. Gendoo
# Sept 9, 2016
# Waterfall plot showing the mutations observed for Key mutations, across all pairs
############################################################################################
# PARSE THE ONCOTATOR DATA RESULTS
############################################################################################
library(GenVisR)

TSVs<-sort(list.files(pattern = ".tsv"))
TSVList<-lapply(TSVs,function(x){read.delim(x,header=T,skip=2)})
TSVs<-gsub(TSVs,pattern = ".snp.maf.txt.tsv",replacement = "")
TSVs<-gsub(TSVs,pattern = "_Lv",replacement="")
TSVs<-gsub(TSVs,pattern = "PCSI_",replacement="")
names(TSVList)<-TSVs

lapply(TSVList,function(x){table((x)$variant_classification)})
lapply(TSVList,function(x){names(table((x)$variant_classification))})
Categories<-unique(unlist(lapply(TSVList,function(x){names(table((x)$variant_classification))})))

CatColors<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','black')

for (Sample in 1:length(TSVList))
{
  TSVList[[Sample]]$Tumor_Sample_Barcode<-TSVs[Sample]
  colnames(TSVList[[Sample]])[81]<-"Hugo_Symbol"
  colnames(TSVList[[Sample]])[23]<-"Variant_Classification"
  colnames(TSVList[[Sample]])[89]<-"Reference_Allele"
}

ALL_MUTATIONS<-do.call(rbind,TSVList)
ALL_MUTATIONS<-ALL_MUTATIONS[!ALL_MUTATIONS$Variant_Classification%in%c("De_novo_Start_OutOfFrame","Start_Codon_SNP","De_novo_Start_InFrame","lincRNA"),]

pdf("Waterfall_Plots_Mutations_KeyGenes_ALL.pdf",onefile = T,height = 3)
waterfall(ALL_MUTATIONS, plotGenes = c("ARID1A","KRAS","TP53","KDM6A","MAP2K4","SMAD4","RNF43","CDKN2A","TGFBR2"),
          mainXlabel = T,mainDropMut = T,sampOrder =c("0489_M","0489_X","0491_M","0491_X",
                                                      "0585_M","0585_X","0604_M","0604_X",
                                                      "0605_M","0605_X","0606_M","0606_X"),geneOrder = T)
dev.off()




