# Deena M.A. Gendoo
# Sept 9, 2016
# Waterfall plot showing the mutations observed for Key mutations, across all pairs
# FINAL_STEP_5_GenerateWaterfall.R

############################################################################################
library(GenVisR)

load("TSV.RData")
#TSVs<-sort(list.files(pattern = ".tsv"))
#TSVList<-lapply(TSVs,function(x){read.delim(x,header=T,skip=2)})
TSVs<-gsub(TSVs,pattern = ".snp.maf.txt.tsv",replacement = "")
TSVs<-gsub(TSVs,pattern = "_Pa",replacement="")
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
waterfall(ALL_MUTATIONS, plotGenes = c("KRAS","TP53","CDKN2A", "SMAD4", "MAP2K4", "ARID1A", "RNF43", "TGFBR2", "KDM6A"),
          mainXlabel = T,mainDropMut = T,sampOrder =c("0169_P","0169_X","0355_P","0355_X",
                                      "0589_P","0589_X","0590_P","0590_X",
                                      "0592_P","0592_X","0602_P","0602_X",
                                      "0611_P","0611_X","0624_P","0624_X",
                                      "0633_P","0633_X","0642_P","0642_X"))
dev.off()
