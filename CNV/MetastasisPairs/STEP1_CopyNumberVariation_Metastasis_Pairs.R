# Deena M.A. Gendoo
# November 1, 2016 and Modified May 11, 2017, final update June 23, 2017
# FINAL_CopyNumberVariation_LIVER.R 
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
ColorListDoubles1<-rep(ColTemplateDoubles1, each=3) #30 colors of 15 pairs 


# Read all parameters file
CNVs<-sort(list.files(pattern = "parameters"))
CNVList<-lapply(CNVs,function(x){read.table(x,header=TRUE)})
CNVs<-sapply(CNVs,function(x) gsub(x,pattern = "parameters_", replacement = ""))
names(CNVList)<-CNVs
names(CNVList)<-gsub(x = names(CNVList),pattern = "_526",replacement = "" )
names(CNVList)<-gsub(x = names(CNVList),pattern = ".txt",replacement = "" )
names(CNVList)<-gsub(x = names(CNVList),pattern = "Lv_",replacement = "" )

# Compress all into one table
ParamsTable<-do.call(rbind.data.frame, CNVList)
write.csv(ParamsTable,"CNV_Params_LIVER.csv")

#Quick formatting for figures!
ParamsTable<-as.matrix(do.call(rbind.data.frame, c(CNVList[1:2],NA,CNVList[3:4],NA,CNVList[5:6],NA,
                                                   CNVList[7:8],NA,CNVList[9:10],NA,CNVList[11:12])))
ParamsTable<-data.frame(ParamsTable)



# Plot the Tumour & Normal component across the pairs
pdf("MergedInfo_LIVER_Pairs.pdf",height = 7)
par(mfrow=c(2,1))
barplot(ParamsTable$Ploidy,col = ColorListDoubles1,las=2,
        names.arg = rownames(ParamsTable),ylab="Ploidy")
abline(a = 2,b = 0,col="red",lty=2)
barplot(ParamsTable$T1*100,col = ColorListDoubles1,las=2,
        ylim = c(0,100),ylab = "% Tumour")
dev.off()

