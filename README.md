**PDAC Disease Models**
========

## Publication

This code and analysis pipeline was designed and developed by **Deena M.A. Gendoo** for the following publication: 

**_Whole Genomes Define Concordance of Matched Primary, Xenograft, and Organoid Models of Pancreas Cancer_**

**Please cite:** 
Gendoo DMA et al. Whole genomes define concordance of matched primary, xenograft, and organoid models of pancreas cancer. PLoS Comput Biol. 2019 Jan 10;15(1):e1006596. doi: 10.1371/journal.pcbi.1006596. PMID: 30629588; PMCID: PMC6328084.
URL to manuscript: https://pubmed.ncbi.nlm.nih.gov/30629588/


**Questions or Comments:** 
Please email deena.gendoo1984@gmail.com or d.gendoo@bham.ac.uk

## Introduction to the Analysis

This repository hosts code to analyze whole-genome sequencing (WGS) data for the following data types: 

- SSM : simple somatic mutation 
- SV : structural variation 
- CNV : copy number variation 

The analysis is conducted on **Pancreatic Ductal Adenocarcinoma (PDAC)**, across three cohorts:

1. Primary resected tumours and matching Xenografts (Primary-PDX pairs) (n=10 pairs)
2. Liver metastasis and matching Xenografts (Metastasis-PDX pairs) (n=6 pairs)
3. Trios: resected primary with matched Xenograft (PDX) and matching organoid (PDO) (Primary-PDX-PDO trios) (n=5 trios)

NB: The trios are a subset of the 10 primary-PDX pairs.  

## The Analysis 

We describe how to reproduce the statistical analysis as reported in the manuscript. To do this, please proceed to:

1. Set up the software environment
2. Run the R scripts

### Set up the software environment

We developed and tested our analysis pipeline using R running on Mac OS X platforms.

To mimic our software environment the following R packages should be installed. All these packages are available on CRAN or Bioconductor.

```
R version 3.3.1 (2016-06-21)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.3 (Yosemite)

locale:
[1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

attached base packages:
 [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] GenVisR_1.0.4              plotrix_3.6-4              vcfR_1.4.0                
 [4] scales_0.4.1               VennDiagram_1.6.17         futile.logger_1.4.3       
 [7] VariantAnnotation_1.18.7   Rsamtools_1.24.0           Biostrings_2.40.2         
[10] XVector_0.12.1             SummarizedExperiment_1.2.3 Biobase_2.32.0            
[13] GenomicRanges_1.24.3       GenomeInfoDb_1.8.7         IRanges_2.6.1             
[16] S4Vectors_0.10.3           copynumber_1.12.0          BiocGenerics_0.18.0       
[19] gplots_3.0.1               RCircos_1.2.0              pheatmap_1.0.8            
[22] reshape2_1.4.2            

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.9             ape_4.0                 lattice_0.20-34         gtools_3.5.0           
 [5] rprojroot_1.2           assertthat_0.1          digest_0.6.11           plyr_1.8.4             
 [9] backports_1.0.5         futile.options_1.0.0    evaluate_0.10           RSQLite_1.1-2          
[13] ggplot2_2.2.1           zlibbioc_1.18.0         GenomicFeatures_1.24.5  lazyeval_0.2.0         
[17] gdata_2.17.0            vegan_2.4-2             Matrix_1.2-7.1          rmarkdown_1.6          
[21] pinfsc50_1.1.0          BiocParallel_1.6.6      stringr_1.2.0           RCurl_1.95-4.8         
[25] biomaRt_2.28.0          munsell_0.4.3           rtracklayer_1.32.2      mgcv_1.8-16            
[29] htmltools_0.3.5         gridExtra_2.2.1         tibble_1.2              XML_3.98-1.5           
[33] permute_0.9-4           viridisLite_0.1.3       GenomicAlignments_1.8.4 MASS_7.3-45            
[37] bitops_1.0-6            nlme_3.1-128            gtable_0.2.0            DBI_0.5-1              
[41] magrittr_1.5            KernSmooth_2.23-15      stringi_1.1.2           viridis_0.3.4          
[45] lambda.r_1.1.9          RColorBrewer_1.1-2      tools_3.3.1             FField_0.1.0           
[49] BSgenome_1.40.1         AnnotationDbi_1.34.4    colorspace_1.3-2        cluster_2.0.5          
[53] caTools_1.17.1          memoise_1.0.0           knitr_1.15.1                         

```
### Running the R Scripts

Once the packages are installed, please download this github repository. 

This repository contains three folders: 

1. **SSM** : code for analysis for simple somatic mutation (SSM) data.
2. **SV**: code for analysis of structural variation (SV) data.
3. **CNV**: code for analysis of copy number variation (CNV) data. 

Each folder is sub-divided into the following folders:

**PrimaryPairs**: for analysis conducted on 10 primary-PDX pairs

**MetastasisPairs**: for analysis conducted on 6 metastasis-PDX pairs

**Trios**: for analysis conducted on 5 primary-PDX-PDO trios (unless otherwise indicated)


To run the analysis, in each subfolder, please first run the R Script that starts with 'STEP1'. 
Once you run this script, all other scripts in that folder can be run individually, in no particular order. 

For subfolders that don't have a 'STEP1', this also means that any of the scripts can be run stand-alone in that folder. 

## List of Datasets Used in the Analysis

Please consult the attached Supplementary data file in the manuscript, which contains whole-genome sequencing (WGS) data for the following data types: 

- SSM : simple somatic mutation 
- SV : structural variation 
- CNV : copy number variation 

### Files needed to run the analysis

The format of the Supplementary Data file follows the same format as the code repository. 
This Data files are split three folders. Below is description of these folders and the types of files contained within them. 

1. **SSM** : Please consult files *.final.vcf or *snp.maf.txt.tsv. 
2. **SV**: Please consult files *.annotatedSV.tsv.
3. **CNV**: CELLULOID parameters are provided in the txt files starting with (parameters). Actual copy number segments are provided in the files (segments*.txt). The code analyzes files BED files (.bed) that show the overlapping genomic segments between tumours and matched disease models.These files end with (_CN_Final.bed)



Each folder is sub-divided into the following folders:

**PrimaryPairs**: for analysis conducted on 10 primary-PDX pairs

**MetastasisPairs**: for analysis conducted on 6 metastasis-PDX pairs

**Trios**: for analysis conducted on 5 primary-PDX-PDO trios (unless otherwise indicated)



Below is a listing of the files found in each folder and subfolder, for every data type:

##### SSM - Primary Pairs

```
PCSI_0169_Pa_P.snp.maf.txt.tsv	PCSI_0590_Pa_X.snp.maf.txt.tsv	PCSI_0624_Pa_P.snp.maf.txt.tsv
PCSI_0169_Pa_P_526.final.vcf	PCSI_0590_Pa_X_526.final.vcf	PCSI_0624_Pa_P_526.final.vcf
PCSI_0169_Pa_X.snp.maf.txt.tsv	PCSI_0592_Pa_P.snp.maf.txt.tsv	PCSI_0624_Pa_X.snp.maf.txt.tsv
PCSI_0169_Pa_X_526.final.vcf	PCSI_0592_Pa_P_526.final.vcf	PCSI_0624_Pa_X_526.final.vcf
PCSI_0355_Pa_P.snp.maf.txt.tsv	PCSI_0592_Pa_X.snp.maf.txt.tsv	PCSI_0633_Pa_P.snp.maf.txt.tsv
PCSI_0355_Pa_P_526.final.vcf	PCSI_0592_Pa_X_526.final.vcf	PCSI_0633_Pa_P_526.final.vcf
PCSI_0355_Pa_X.snp.maf.txt.tsv	PCSI_0602_Pa_P.snp.maf.txt.tsv	PCSI_0633_Pa_X.snp.maf.txt.tsv
PCSI_0355_Pa_X_526.final.vcf	PCSI_0602_Pa_P_526.final.vcf	PCSI_0633_Pa_X_526.final.vcf
PCSI_0589_Pa_P.snp.maf.txt.tsv	PCSI_0602_Pa_X.snp.maf.txt.tsv	PCSI_0642_Pa_P.snp.maf.txt.tsv
PCSI_0589_Pa_P_526.final.vcf	PCSI_0602_Pa_X_526.final.vcf	PCSI_0642_Pa_P_526.final.vcf
PCSI_0589_Pa_X.snp.maf.txt.tsv	PCSI_0611_Pa_P.snp.maf.txt.tsv	PCSI_0642_Pa_X.snp.maf.txt.tsv
PCSI_0589_Pa_X_526.final.vcf	PCSI_0611_Pa_P_526.final.vcf	PCSI_0642_Pa_X_526.final.vcf
PCSI_0590_Pa_P.snp.maf.txt.tsv	PCSI_0611_Pa_X.snp.maf.txt.tsv
PCSI_0590_Pa_P_526.final.vcf	PCSI_0611_Pa_X_526.final.vcf

```

##### SSM - Metastasis Pairs

```
PCSI_0489_Lv_M.snp.maf.txt.tsv	PCSI_0585_Lv_M.snp.maf.txt.tsv	PCSI_0605_Lv_M.snp.maf.txt.tsv
PCSI_0489_Lv_M_526.final.vcf	PCSI_0585_Lv_M_526.final.vcf	PCSI_0605_Lv_M_526.final.vcf
PCSI_0489_Lv_X.snp.maf.txt.tsv	PCSI_0585_Lv_X.snp.maf.txt.tsv	PCSI_0605_Lv_X.snp.maf.txt.tsv
PCSI_0489_Lv_X_526.final.vcf	PCSI_0585_Lv_X_526.final.vcf	PCSI_0605_Lv_X_526.final.vcf
PCSI_0491_Lv_M.snp.maf.txt.tsv	PCSI_0604_Lv_M.snp.maf.txt.tsv	PCSI_0606_Lv_M.snp.maf.txt.tsv
PCSI_0491_Lv_M_526.final.vcf	PCSI_0604_Lv_M_526.final.vcf	PCSI_0606_Lv_M_526.final.vcf
PCSI_0491_Lv_X.snp.maf.txt.tsv	PCSI_0604_Lv_X.snp.maf.txt.tsv	PCSI_0606_Lv_X.snp.maf.txt.tsv
PCSI_0491_Lv_X_526.final.vcf	PCSI_0604_Lv_X_526.final.vcf	PCSI_0606_Lv_X_526.final.vcf

```

##### SSM - Trios

```
PCSI_0590_Pa_O.snp.maf.txt.tsv	PCSI_0592_Pa_X.snp.maf.txt.tsv	PCSI_0624_Pa_P.snp.maf.txt.tsv
PCSI_0590_Pa_O_526.final.vcf	PCSI_0592_Pa_X_526.final.vcf	PCSI_0624_Pa_P_526.final.vcf
PCSI_0590_Pa_P.snp.maf.txt.tsv	PCSI_0602_Pa_O.snp.maf.txt.tsv	PCSI_0624_Pa_X.snp.maf.txt.tsv
PCSI_0590_Pa_P_526.final.vcf	PCSI_0602_Pa_O_526.final.vcf	PCSI_0624_Pa_X_526.final.vcf
PCSI_0590_Pa_X.snp.maf.txt.tsv	PCSI_0602_Pa_P.snp.maf.txt.tsv	PCSI_0642_Pa_O.snp.maf.txt.tsv
PCSI_0590_Pa_X_526.final.vcf	PCSI_0602_Pa_P_526.final.vcf	PCSI_0642_Pa_O_526.final.vcf
PCSI_0592_Pa_O.snp.maf.txt.tsv	PCSI_0602_Pa_X.snp.maf.txt.tsv	PCSI_0642_Pa_P.snp.maf.txt.tsv
PCSI_0592_Pa_O_526.final.vcf	PCSI_0602_Pa_X_526.final.vcf	PCSI_0642_Pa_P_526.final.vcf
PCSI_0592_Pa_P.snp.maf.txt.tsv	PCSI_0624_Pa_O.snp.maf.txt.tsv	PCSI_0642_Pa_X.snp.maf.txt.tsv
PCSI_0592_Pa_P_526.final.vcf	PCSI_0624_Pa_O_526.final.vcf	PCSI_0642_Pa_X_526.final.vcf

```

##### SV - Primary Pairs

```
PCSI_0169_Pa_P_526.annotatedSV.tsv	PCSI_0590_Pa_X_526.annotatedSV.tsv	PCSI_0624_Pa_P_526.annotatedSV.tsv
PCSI_0169_Pa_X_526.annotatedSV.tsv	PCSI_0592_Pa_P_526.annotatedSV.tsv	PCSI_0624_Pa_X_526.annotatedSV.tsv
PCSI_0355_Pa_P_526.annotatedSV.tsv	PCSI_0592_Pa_X_526.annotatedSV.tsv	PCSI_0633_Pa_P_526.annotatedSV.tsv
PCSI_0355_Pa_X_526.annotatedSV.tsv	PCSI_0602_Pa_P_526.annotatedSV.tsv	PCSI_0633_Pa_X_526.annotatedSV.tsv
PCSI_0589_Pa_P_526.annotatedSV.tsv	PCSI_0602_Pa_X_526.annotatedSV.tsv	PCSI_0642_Pa_P_526.annotatedSV.tsv
PCSI_0589_Pa_X_526.annotatedSV.tsv	PCSI_0611_Pa_P_526.annotatedSV.tsv	PCSI_0642_Pa_X_526.annotatedSV.tsv
PCSI_0590_Pa_P_526.annotatedSV.tsv	PCSI_0611_Pa_X_526.annotatedSV.tsv

```

##### SV - Metastasis Pairs

```
PCSI_0489_Lv_M_526.annotatedSV.tsv	PCSI_0585_Lv_M_526.annotatedSV.tsv	PCSI_0605_Lv_M_526.annotatedSV.tsv
PCSI_0489_Lv_X_526.annotatedSV.tsv	PCSI_0585_Lv_X_526.annotatedSV.tsv	PCSI_0605_Lv_X_526.annotatedSV.tsv
PCSI_0491_Lv_M_526.annotatedSV.tsv	PCSI_0604_Lv_M_526.annotatedSV.tsv	PCSI_0606_Lv_M_526.annotatedSV.tsv
PCSI_0491_Lv_X_526.annotatedSV.tsv	PCSI_0604_Lv_X_526.annotatedSV.tsv	PCSI_0606_Lv_X_526.annotatedSV.tsv

```

##### SV - Trios

```
PCSI_0590_Pa_O_526.annotatedSV.tsv	PCSI_0592_Pa_X_526.annotatedSV.tsv	PCSI_0624_Pa_P_526.annotatedSV.tsv
PCSI_0590_Pa_P_526.annotatedSV.tsv	PCSI_0602_Pa_O_526.annotatedSV.tsv	PCSI_0624_Pa_X_526.annotatedSV.tsv
PCSI_0590_Pa_X_526.annotatedSV.tsv	PCSI_0602_Pa_P_526.annotatedSV.tsv	PCSI_0642_Pa_O_526.annotatedSV.tsv
PCSI_0592_Pa_O_526.annotatedSV.tsv	PCSI_0602_Pa_X_526.annotatedSV.tsv	PCSI_0642_Pa_P_526.annotatedSV.tsv
PCSI_0592_Pa_P_526.annotatedSV.tsv	PCSI_0624_Pa_O_526.annotatedSV.tsv	PCSI_0642_Pa_X_526.annotatedSV.tsv

```

##### CNV - Primary Pairs

```
parameters_PCSI_0169_Pa_P_526.txt	parameters_PCSI_0642_Pa_P_526.txt	segments_PCSI_0611_CN_Final.bed
parameters_PCSI_0169_Pa_X_526.txt	parameters_PCSI_0642_Pa_X_526.txt	segments_PCSI_0611_Pa_P_526.txt
parameters_PCSI_0355_Pa_P_526.txt	segments_PCSI_0169_CN_Final.bed		segments_PCSI_0611_Pa_X_526.txt
parameters_PCSI_0355_Pa_X_526.txt	segments_PCSI_0169_Pa_P_526.txt		segments_PCSI_0624_CN_Final.bed
parameters_PCSI_0589_Pa_P_526.txt	segments_PCSI_0169_Pa_X_526.txt		segments_PCSI_0624_Pa_P_526.txt
parameters_PCSI_0589_Pa_X_526.txt	segments_PCSI_0355_CN_Final.bed		segments_PCSI_0624_Pa_X_526.txt
parameters_PCSI_0590_Pa_P_526.txt	segments_PCSI_0355_Pa_P_526.txt		segments_PCSI_0633_CN_Final.bed
parameters_PCSI_0590_Pa_X_526.txt	segments_PCSI_0355_Pa_X_526.txt		segments_PCSI_0633_Pa_P_526.txt
parameters_PCSI_0611_Pa_P_526.txt	segments_PCSI_0589_CN_Final.bed		segments_PCSI_0633_Pa_X_526.txt
parameters_PCSI_0611_Pa_X_526.txt	segments_PCSI_0589_Pa_P_526.txt		segments_PCSI_0642_CN_Final.bed
parameters_PCSI_0624_Pa_P_526.txt	segments_PCSI_0589_Pa_X_526.txt		segments_PCSI_0642_Pa_P_526.txt
parameters_PCSI_0624_Pa_X_526.txt	segments_PCSI_0590_CN_Final.bed		segments_PCSI_0642_Pa_X_526.txt
parameters_PCSI_0633_Pa_P_526.txt	segments_PCSI_0590_Pa_P_526.txt
parameters_PCSI_0633_Pa_X_526.txt	segments_PCSI_0590_Pa_X_526.txt

```

##### CNV - Metastasis Pairs

```
parameters_PCSI_0489_Lv_M_526.txt	parameters_PCSI_0606_Lv_M_526.txt	segments_PCSI_0585_Lv_X_526.txt
parameters_PCSI_0489_Lv_X_526.txt	parameters_PCSI_0606_Lv_X_526.txt	segments_PCSI_0604_CN_Final.bed
parameters_PCSI_0491_Lv_M_526.txt	segments_PCSI_0489_CN_Final.bed		segments_PCSI_0604_Lv_M_526.txt
parameters_PCSI_0491_Lv_X_526.txt	segments_PCSI_0489_Lv_M_526.txt		segments_PCSI_0604_Lv_X_526.txt
parameters_PCSI_0585_Lv_M_526.txt	segments_PCSI_0489_Lv_X_526.txt		segments_PCSI_0605_CN_Final.bed
parameters_PCSI_0585_Lv_X_526.txt	segments_PCSI_0491_CN_Final.bed		segments_PCSI_0605_Lv_M_526.txt
parameters_PCSI_0604_Lv_M_526.txt	segments_PCSI_0491_Lv_M_526.txt		segments_PCSI_0605_Lv_X_526.txt
parameters_PCSI_0604_Lv_X_526.txt	segments_PCSI_0491_Lv_X_526.txt		segments_PCSI_0606_CN_Final.bed
parameters_PCSI_0605_Lv_M_526.txt	segments_PCSI_0585_CN_Final.bed		segments_PCSI_0606_Lv_M_526.txt
parameters_PCSI_0605_Lv_X_526.txt	segments_PCSI_0585_Lv_M_526.txt		segments_PCSI_0606_Lv_X_526.txt


```


##### CNV - Trios

```
parameters_PCSI_0590_Pa_O_526.txt	parameters_PCSI_0642_Pa_X_526.txt	segments_PCSI_0624_Pa_X_526.txt
parameters_PCSI_0590_Pa_P_526.txt	segments_PCSI_0590_Pa_O_526.txt		segments_PCSI_0624_PvsO_CN_Final.bed
parameters_PCSI_0590_Pa_X_526.txt	segments_PCSI_0590_Pa_P_526.txt		segments_PCSI_0624_XvsO_CN_Final.bed
parameters_PCSI_0624_Pa_O_526.txt	segments_PCSI_0590_Pa_X_526.txt		segments_PCSI_0642_Pa_O_526.txt
parameters_PCSI_0624_Pa_P_526.txt	segments_PCSI_0590_PvsO_CN_Final.bed	segments_PCSI_0642_Pa_P_526.txt
parameters_PCSI_0624_Pa_X_526.txt	segments_PCSI_0590_XvsO_CN_Final.bed	segments_PCSI_0642_Pa_X_526.txt
parameters_PCSI_0642_Pa_O_526.txt	segments_PCSI_0624_Pa_O_526.txt		segments_PCSI_0642_PvsO_CN_Final.bed
parameters_PCSI_0642_Pa_P_526.txt	segments_PCSI_0624_Pa_P_526.txt		segments_PCSI_0642_XvsO_CN_Final.bed

```

### 
