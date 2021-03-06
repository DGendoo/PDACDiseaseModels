# Deena M.A. Gendoo
# November 1, 2016 and final update June 23, 2017
# FINAL_CopyNumberVariation_PDAC_Trios.R
# Parse Parameters for Copy Number Data from CELLULOID
# Each parameters file has 5 columns representing:
    # 1- value: % captured by the model
    # 2- S
    # 3- N: %normal
    # 4- T1: % tumour
    # 5- ploidy: ploidy determined for that solution in celluloid

###################################################################################################################
###################################################################################################################

ColTemplateDoubles1<-c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462")
ColorListDoubles1<-rep(ColTemplateDoubles1, each=4) #30 colors of 15 pairs 

# Read all parameters file
CNVs<-sort(list.files(pattern = "parameters"))
CNVList<-lapply(CNVs,function(x){read.table(x,header=TRUE)})
CNVs<-sapply(CNVs,function(x) gsub(x,pattern = "parameters_", replacement = ""))
names(CNVList)<-CNVs
names(CNVList)<-gsub(x = names(CNVList),pattern = "_526",replacement = "" )
names(CNVList)<-gsub(x = names(CNVList),pattern = ".txt",replacement = "" )
names(CNVList)<-gsub(x = names(CNVList),pattern = "Pa_",replacement = "" )

# Compress all into one table
ParamsTable<-do.call(rbind.data.frame, CNVList)

ParamsTable<-ParamsTable[c("PCSI_0590_P","PCSI_0590_X","PCSI_0590_O","PCSI_0624_P","PCSI_0624_X","PCSI_0624_O",
                           "PCSI_0642_P","PCSI_0642_X","PCSI_0642_O"),]
write.csv(ParamsTable,"CNV_Params_PDAC_Trios.csv")

#Quick formatting for figures!
ParamsTable<-as.matrix(do.call(rbind.data.frame, c(CNVList[1:3],NA,CNVList[4:6],NA,CNVList[7:9])))
ParamsTable<-ParamsTable[c("PCSI_0590_P","PCSI_0590_X","PCSI_0590_O","4","PCSI_0624_P","PCSI_0624_X","PCSI_0624_O","8",
                           "PCSI_0642_P","PCSI_0642_X","PCSI_0642_O"),]
ParamsTable<-data.frame(ParamsTable)

# Plot the Tumour & Normal component across the pairs
pdf("MergedInfo_PDAC_Trios.pdf",height = 7)
par(mfrow=c(2,1))
barplot(ParamsTable$Ploidy,col = ColorListDoubles1,las=2,
        names.arg = rownames(ParamsTable),ylab="Ploidy",ylim=c(0,4))
abline(a = 2,b = 0,col="red",lty=2)
barplot(ParamsTable$T1*100,col = ColorListDoubles1,las=2,
        ylim = c(0,100),ylab = "% Tumour")
dev.off()

