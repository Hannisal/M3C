# time:20/06/2021
# Auther:ShuaiWang
library(survival)
library(plyr)
## note:the expression data and survival data was download from GDC TCGA (cancertype) project,while the phenotype files were download from the TCGA (cancertype) project.Both them was obtained from UCSC.
## BRCA
## 1.Match expression data and clinical profiles to obtain a cohort containing all patients with no clinical feature missing
exp<-read.table("1.BRCA/HiSeqV2",header = T,row.names = 1,check.names = F)
sur<-read.table("1.BRCA/TCGA-BRCA.survival.tsv",sep = "\t",header = T)
pheno<-read.table("1.BRCA/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix",sep = "\t",header = T)
colnames(pheno)
pheno<-pheno[,c("sampleID","Age_at_Initial_Pathologic_Diagnosis_nature2012","Gender_nature2012","histological_type","Tumor_nature2012","Node_nature2012","Metastasis_nature2012","ER_Status_nature2012","PR_Status_nature2012","HER2_Final_Status_nature2012","additional_pharmaceutical_therapy","additional_radiation_therapy","additional_surgery_locoregional_procedure","additional_surgery_metastatic_procedure","history_of_neoadjuvant_treatment")]

## 2.Delete patients with clinical feature missing
pheno<-pheno[!is.na(pheno$Age_at_Initial_Pathologic_Diagnosis_nature2012),]
pheno<-pheno[!pheno$Gender_nature2012=="",]
pheno<-pheno[!pheno$Tumor_nature2012=="",]
pheno<-pheno[!pheno$Node_nature2012=="",]
pheno<-pheno[!pheno$Metastasis_nature2012=="",]
pheno<-pheno[!pheno$ER_Status_nature2012=="",]
pheno<-pheno[!pheno$PR_Status_nature2012=="",]
pheno<-pheno[!pheno$HER2_Final_Status_nature2012=="",]


## 3.Delete patients with indeterminate clinical feature
table(pheno$Gender_nature2012)
table(pheno$Tumor_nature2012)
table(pheno$Node_nature2012)
table(pheno$Metastasis_nature2012)
table(pheno$ER_Status_nature2012)
table(pheno$PR_Status_nature2012)
table(pheno$HER2_Final_Status_nature2012)
### T stage
pheno<-pheno[!pheno$Tumor_nature2012=="TX",]
### ER
pheno<-pheno[!pheno$ER_Status_nature2012=="Indeterminate",]
### PR
pheno<-pheno[!pheno$PR_Status_nature2012=="Indeterminate",]
### HER2
pheno<-pheno[!pheno$HER2_Final_Status_nature2012=="Equivocal",]

sum(!duplicated(pheno$sampleID))
pheno<-pheno[!duplicated(pheno$sampleID),]

pheno<-pheno[pheno$sampleID%in%sur$sample,]

sur<-sur[sur$sample%in%pheno$sampleID,-2]
sum(!duplicated(sur$sample))

sur<-sur[!duplicated(sur$sample),]
exp<-exp[,colnames(exp)%in%sur$sample]
exp<-exp[,!duplicated(colnames(exp))]

pheno<-pheno[match(colnames(exp),pheno$sampleID),]
sur<-sur[match(colnames(exp),sur$sample),]

## 3.Finally,719 patients was chosen
save(exp,pheno,sur,file = "TCGA-UCSC-BRCA.Rdata")

## KIRC
## 1.Match expression data and clinical profiles to obtain a cohort containing all patients with no clinical feature missing
exp<-read.table("2.KIRC/HiSeqV2",header = T,row.names = 1,check.names = F)
sur<-read.table("2.KIRC/TCGA-KIRC.survival.tsv",sep = "\t",header = T)
pheno<-read.table("2.KIRC/TCGA.KIRC.sampleMap_KIRC_clinicalMatrix",sep = "\t",header = T)
colnames(pheno)

pheno<-pheno[,c("sampleID","age_at_initial_pathologic_diagnosis","gender","histological_type","neoplasm_histologic_grade","pathologic_T","pathologic_N","pathologic_M","hemoglobin_result","lactate_dehydrogenase_result","platelet_qualitative_result","white_cell_count_result","additional_pharmaceutical_therapy","additional_radiation_therapy","additional_surgery_locoregional_procedure","additional_surgery_metastatic_procedure","history_of_neoadjuvant_treatment","followup_treatment_success")]

## 2.Delete patients with clinical feature missing
pheno<-pheno[!is.na(pheno$age_at_initial_pathologic_diagnosis),]
pheno<-pheno[!pheno$gender=="",]
pheno<-pheno[!pheno$neoplasm_histologic_grade=="",]
pheno<-pheno[!pheno$pathologic_T=="",]
pheno<-pheno[!pheno$pathologic_N=="",]
pheno<-pheno[!pheno$pathologic_M=="",]


## 3.Delete patients with indeterminate clinical feature
table(pheno$gender)
table(pheno$pathologic_T)
table(pheno$pathologic_N)
table(pheno$pathologic_M)
table(pheno$neoplasm_histologic_grade)


### T stage
pheno$pathologic_T<-substr(pheno$pathologic_T,1,2)
### N stage
pheno<-pheno[!pheno$pathologic_N=="NX",]
### M stage
pheno<-pheno[!pheno$pathologic_M=="MX",]
### Grage
pheno<-pheno[!pheno$neoplasm_histologic_grade=="GX",]


sum(!duplicated(pheno$sampleID))
pheno<-pheno[!duplicated(pheno$sampleID),]


pheno<-pheno[pheno$sampleID%in%sur$sample,]

sur<-sur[sur$sample%in%pheno$sampleID,-2]
sum(!duplicated(sur$sample))
sur<-sur[!duplicated(sur$sample),]




exp<-exp[,colnames(exp)%in%sur$sample]
sum(!duplicated(colnames(exp)))
exp<-exp[,!duplicated(colnames(exp))]

pheno<-pheno[match(colnames(exp),pheno$sampleID),]
sur<-sur[match(colnames(exp),sur$sample),]

## 3.Finally,288 patients was chosen
save(exp,pheno,sur,file = "TCGA-UCSC-KIRC.Rdata")

## LIHC
## 1.Match expression data and clinical profiles to obtain a cohort containing all patients with no clinical feature missing
exp<-read.table("3.LIHC/HiSeqV2",header = T,row.names = 1,check.names = F)
sur<-read.table("3.LIHC/TCGA-LIHC.survival.tsv",sep = "\t",header = T)
pheno<-read.table("3.LIHC/TCGA.LIHC.sampleMap_LIHC_clinicalMatrix",sep = "\t",header = T)
colnames(pheno)

pheno<-pheno[,c("sampleID","age_at_initial_pathologic_diagnosis","gender","histological_type","neoplasm_histologic_grade","pathologic_T","pathologic_N","pathologic_M","additional_pharmaceutical_therapy","additional_radiation_therapy","history_of_neoadjuvant_treatment")]

## 2.Delete patients with clinical feature missing
pheno<-pheno[!is.na(pheno$age_at_initial_pathologic_diagnosis),]
pheno<-pheno[!pheno$gender=="",]
pheno<-pheno[!pheno$neoplasm_histologic_grade=="",]
pheno<-pheno[!pheno$pathologic_T=="",]
pheno<-pheno[!pheno$pathologic_N=="",]
pheno<-pheno[!pheno$pathologic_M=="",]


## 3.Delete patients with indeterminate clinical feature
table(pheno$gender)
table(pheno$pathologic_T)
table(pheno$pathologic_N)
table(pheno$pathologic_M)
table(pheno$neoplasm_histologic_grade)


### T stage
pheno<-pheno[!pheno$pathologic_T=="[Discrepancy]",]
pheno$pathologic_T<-substr(pheno$pathologic_T,1,2)
### N stage
pheno<-pheno[!pheno$pathologic_N=="NX",]
### M stage
pheno<-pheno[!pheno$pathologic_M=="MX",]
### Grage
pheno<-pheno[!pheno$neoplasm_histologic_grade=="GX",]


sum(!duplicated(pheno$sampleID))
pheno<-pheno[!duplicated(pheno$sampleID),]


pheno<-pheno[pheno$sampleID%in%sur$sample,]

sur<-sur[sur$sample%in%pheno$sampleID,-2]
sum(!duplicated(sur$sample))
sur<-sur[!duplicated(sur$sample),]


exp<-exp[,!duplicated(colnames(exp))]

exp<-exp[,colnames(exp)%in%sur$sample]
sum(!duplicated(colnames(exp)))


pheno<-pheno[match(colnames(exp),pheno$sampleID),]
sur<-sur[match(colnames(exp),sur$sample),]

## 3.Finally,263patients was chosen
save(exp,pheno,sur,file = "TCGA-UCSC-LIHC.Rdata")
## LUAD
## 1.Match expression data and clinical profiles to obtain a cohort containing all patients with no clinical feature missing
exp<-read.table("4.LUAD/HiSeqV2",header = T,row.names = 1,check.names = F)
sur<-read.table("4.LUAD/TCGA-LUAD.survival.tsv",sep = "\t",header = T)
pheno<-read.table("4.LUAD/TCGA.LUAD.sampleMap_LUAD_clinicalMatrix",sep = "\t",header = T)
colnames(pheno)

pheno<-pheno[,c("sampleID","age_at_initial_pathologic_diagnosis","gender","histological_type","pathologic_T","pathologic_N","pathologic_M","additional_pharmaceutical_therapy","additional_radiation_therapy","additional_surgery_locoregional_procedure","additional_surgery_metastatic_procedure","history_of_neoadjuvant_treatment","followup_treatment_success","egfr_mutation_result","kras_mutation_result")]

## 2.Delete patients with clinical feature missing
pheno<-pheno[!is.na(pheno$age_at_initial_pathologic_diagnosis),]
pheno<-pheno[!pheno$gender=="",]
pheno<-pheno[!pheno$pathologic_T=="",]
pheno<-pheno[!pheno$pathologic_N=="",]
pheno<-pheno[!pheno$pathologic_M=="",]


## 3.Delete patients with indeterminate clinical feature
table(pheno$gender)
table(pheno$pathologic_T)
table(pheno$pathologic_N)
table(pheno$pathologic_M)


### T stage
pheno$pathologic_T<-substr(pheno$pathologic_T,1,2)
### N stage
pheno<-pheno[!pheno$pathologic_N=="NX",]
### M stage
pheno$pathologic_M<-substr(pheno$pathologic_M,1,2)
pheno<-pheno[!pheno$pathologic_M=="MX",]



sum(!duplicated(pheno$sampleID))
pheno<-pheno[!duplicated(pheno$sampleID),]


pheno<-pheno[pheno$sampleID%in%sur$sample,]

sur<-sur[sur$sample%in%pheno$sampleID,-2]
sum(!duplicated(sur$sample))
sur<-sur[!duplicated(sur$sample),]


exp<-exp[,!duplicated(colnames(exp))]

exp<-exp[,colnames(exp)%in%sur$sample]
sum(!duplicated(colnames(exp)))


pheno<-pheno[match(colnames(exp),pheno$sampleID),]
sur<-sur[match(colnames(exp),sur$sample),]

## 3.Finally,381 patients was chosen
save(exp,pheno,sur,file = "TCGA-UCSC-LUAD.Rdata")
## STAD
## 1.Match expression data and clinical profiles to obtain a cohort containing all patients with no clinical feature missing
exp<-read.table("5.STAD/HiSeqV2",header = T,row.names = 1,check.names = F)
sur<-read.table("5.STAD/TCGA-STAD.survival.tsv",sep = "\t",header = T)
pheno<-read.table("5.STAD/TCGA.STAD.sampleMap_STAD_clinicalMatrix",sep = "\t",header = T)
colnames(pheno)

pheno<-pheno[,c("sampleID","age_at_initial_pathologic_diagnosis","gender","histological_type","neoplasm_histologic_grade","pathologic_T","pathologic_N","pathologic_M","additional_pharmaceutical_therapy","additional_radiation_therapy","additional_surgery_locoregional_procedure","additional_surgery_metastatic_procedure","history_of_neoadjuvant_treatment","followup_treatment_success")]

## 2.Delete patients with clinical feature missing
pheno<-pheno[!is.na(pheno$age_at_initial_pathologic_diagnosis),]
pheno<-pheno[!pheno$gender=="",]
pheno<-pheno[!pheno$neoplasm_histologic_grade=="",]
pheno<-pheno[!pheno$pathologic_T=="",]
pheno<-pheno[!pheno$pathologic_N=="",]
pheno<-pheno[!pheno$pathologic_M=="",]


## 3.Delete patients with indeterminate clinical feature
table(pheno$gender)
table(pheno$pathologic_T)
table(pheno$pathologic_N)
table(pheno$pathologic_M)
table(pheno$neoplasm_histologic_grade)


### T stage
pheno$pathologic_T<-substr(pheno$pathologic_T,1,2)
pheno<-pheno[!pheno$pathologic_T=="TX",]
### N stage
pheno<-pheno[!pheno$pathologic_N=="[Discrepancy]",]
pheno$pathologic_N<-substr(pheno$pathologic_N,1,2)
pheno<-pheno[!pheno$pathologic_N=="NX",]
### M stage
pheno<-pheno[!pheno$pathologic_M=="MX",]
### Grage
pheno<-pheno[!pheno$neoplasm_histologic_grade=="GX",]


sum(!duplicated(pheno$sampleID))
pheno<-pheno[!duplicated(pheno$sampleID),]


pheno<-pheno[pheno$sampleID%in%sur$sample,]

sur<-sur[sur$sample%in%pheno$sampleID,-2]
sum(!duplicated(sur$sample))
sur<-sur[!duplicated(sur$sample),]


exp<-exp[,!duplicated(colnames(exp))]

exp<-exp[,colnames(exp)%in%sur$sample]
sum(!duplicated(colnames(exp)))


pheno<-pheno[match(colnames(exp),pheno$sampleID),]
sur<-sur[match(colnames(exp),sur$sample),]

## 3.Finally,397 patients was chosen
save(exp,pheno,sur,file = "TCGA-UCSC-STAD.Rdata")
