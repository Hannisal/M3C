# time:04/06/2021
# Auther:ShuaiWang
options(digits = 4)
rm(list = ls())
library(survival)
library(plyr)
## BRCA
## Based on the expression data from the included BRCA patients with Male deleted
### METTL2A
cancertype<-"BRCA"
gene<-"METTL2A"
load("TCGA-UCSC-BRCA.Rdata")
exp_mettl2a<-exp[gene,]
## 1.univariate analysis
unicox<-cbind(pheno[,-c(4,11:15)],sur[,2:3],t(exp_mettl2a))
### transform clinical feature into numeric
table(unicox$Age_at_Initial_Pathologic_Diagnosis_nature2012)
median(unicox$Age_at_Initial_Pathologic_Diagnosis_nature2012)
unicox$Age_at_Initial_Pathologic_Diagnosis_nature2012<-ifelse(unicox$Age_at_Initial_Pathologic_Diagnosis_nature2012<=median(unicox$Age_at_Initial_Pathologic_Diagnosis_nature2012),0,1)

class(unicox)
table(unicox$Gender_nature2012)
unicox$Gender_nature2012<-ifelse(unicox$Gender_nature2012=="FEMALE",0,1)

table(unicox$Tumor_nature2012)
unicox$Tumor_nature2012<-ifelse(unicox$Tumor_nature2012=="T1",1,ifelse(unicox$Tumor_nature2012=="T2",2,ifelse(unicox$Tumor_nature2012=="T3",3,4)))

table(unicox$Node_nature2012)
unicox$Node_nature2012<-ifelse(unicox$Node_nature2012=="N0",0,ifelse(unicox$Node_nature2012=="N1",1,ifelse(unicox$Node_nature2012=="N2",2,3)))

table(unicox$Metastasis_nature2012)
unicox$Metastasis_nature2012<-ifelse(unicox$Metastasis_nature2012=="M0",0,1)

table(unicox$ER_Status_nature2012)
unicox$ER_Status_nature2012<-ifelse(unicox$ER_Status_nature2012=="Negative",0,1)

table(unicox$PR_Status_nature2012)
unicox$PR_Status_nature2012<-ifelse(unicox$PR_Status_nature2012=="Negative",0,1)

table(unicox$HER2_Final_Status_nature2012)
unicox$HER2_Final_Status_nature2012<-ifelse(unicox$HER2_Final_Status_nature2012=="Negative",0,1)

colnames(unicox)<-c("sampleID","Age","Gender","T.stage","N.stage","M.stage","ER.status","PR.status","HER2.status","OS","OS.time","METTL2A")



save(unicox,file = paste0(cancertype,"-",gene,"-","unicox30.Rdata"))

## delete the male petients and patients with OS.time<30
unicox<-unicox[unicox$Gender=="0",]
unicox<-unicox[!unicox$OS.time<30,]
unicox<-unicox[,-3]
BaSurv<-Surv(time = unicox$OS.time,event = unicox$OS)
unicox$BaSurv<-with(unicox,BaSurv)

GCox<-coxph(BaSurv~METTL2A,data = unicox)
Gsum<-summary(GCox)
HR<-round(Gsum$coefficients[,2],2)
PValue<-round(Gsum$coefficients[,5],4)
CI<-paste0(round(Gsum$conf.int[,3:4],2),collapse = "-")
uniouttalbe<-data.frame('Characteristics'='METTL2A',
                        'Hazard Ratio'=HR,
                        'CI95'=CI,
                        'P Value'=PValue)
uniouttalbe<-function(x){
  FML<-as.formula(paste0('BaSurv~',x))
  GCox<-coxph(FML,data = unicox)
  Gsum<-summary(GCox)
  HR<-round(Gsum$coefficients[,2],2)
  PValue<-round(Gsum$coefficients[,5],4)
  CI<-paste0(round(Gsum$conf.int[,3:4],2),collapse = "-")
  uniouttalbe<-data.frame('Characteristics'=x,
                          'Hazard Ratio'=HR,
                          'CI95'=CI,
                          'P Value'=PValue)
  return(uniouttalbe)
}
uniouttalbe("METTL2A")

VarNames<-colnames(unicox)[c(2:8,11)]
Univar<-lapply(VarNames, uniouttalbe)
Univar<-ldply(Univar,data.frame)
Univar$Characteristics[Univar$P.Value<0.1]

## 2.mutivariate analysis
fml<-as.formula(paste0('BaSurv~',paste0(Univar$Characteristics[Univar$P.Value<0.1],collapse = "+")))
Multicox<-coxph(fml,data = unicox)
MultiSum<-summary(Multicox)

MultiName<-as.character(Univar$Characteristics[Univar$P.Value<0.1])
MHR<-round(MultiSum$coefficients[,2],2)
MPV<-round(MultiSum$coefficients[,5],4)
MCIL<-round(MultiSum$conf.int[,3],2)
MCIU<-round(MultiSum$conf.int[,4],2)
MCI<-paste0(MCIL,"-",MCIU)
Multicox<-data.frame('Characteristics'=MultiName,
                     'Hazard Ratio'=MHR,
                     'CI95'=MCI,
                     'P Value'= MPV)
write.csv(Univar,file = paste0(cancertype,"-",gene,"-","Univar.csv"))
write.csv(Multicox,file = paste0(cancertype,"-",gene,"-","Multicox.csv"))

## KIRC or LIHC or LUAD or STAD
### METTL2A or METTL2B or METTL6 or METTL8
cancertype<-"STAD"
gene<-"METTL6"

load(paste0("TCGA-UCSC-",cancertype,".Rdata"))
exp_mettl2a<-exp[gene,]
## 1.univariate analysis
unicox<-cbind(pheno[,-c(4,9:14)],sur[,2:3],t(exp_mettl2a))
### transform clinical feature into numeric
class(unicox)
median(unicox$age_at_initial_pathologic_diagnosis)
unicox$age_at_initial_pathologic_diagnosis<-ifelse(unicox$age_at_initial_pathologic_diagnosis<=median(unicox$age_at_initial_pathologic_diagnosis),0,1)


table(unicox$gender)
unicox$gender<-ifelse(unicox$gender=="FEMALE",0,1)

table(unicox$pathologic_T)
unicox$pathologic_T<-ifelse(unicox$pathologic_T=="T1",1,ifelse(unicox$pathologic_T=="T2",2,ifelse(unicox$pathologic_T=="T3",3,4)))

table(unicox$pathologic_N)
unicox$pathologic_N<-ifelse(unicox$pathologic_N=="N0",0,ifelse(unicox$pathologic_N=="N1",1,ifelse(unicox$pathologic_N=="N2",2,3)))

table(unicox$pathologic_M)
unicox$pathologic_M<-ifelse(unicox$pathologic_M=="M0",0,1)

table(unicox$neoplasm_histologic_grade)
unicox$neoplasm_histologic_grade<-ifelse(unicox$neoplasm_histologic_grade=="G1",1,ifelse(unicox$neoplasm_histologic_grade=="G2",2,ifelse(unicox$neoplasm_histologic_grade=="G3",3,4)))


colnames(unicox)<-c("sampleID","Age","Gender","Grade","T.stage","N.stage","M.stage","OS","OS.time",gene)



save(unicox,file = paste0(cancertype,"-",gene,"-","unicox30.Rdata"))


unicox<-unicox[!unicox$OS.time<30,]

BaSurv<-Surv(time = unicox$OS.time,event = unicox$OS)
unicox$BaSurv<-with(unicox,BaSurv)

GCox<-coxph(BaSurv~Age,data = unicox)
Gsum<-summary(GCox)
HR<-round(Gsum$coefficients[,2],2)
PValue<-round(Gsum$coefficients[,5],4)
CI<-paste0(round(Gsum$conf.int[,3:4],2),collapse = "-")
uniouttalbe<-data.frame('Characteristics'='Age',
                        'Hazard Ratio'=HR,
                        'CI95'=CI,
                        'P Value'=PValue)
uniouttalbe<-function(x){
  FML<-as.formula(paste0('BaSurv~',x))
  GCox<-coxph(FML,data = unicox)
  Gsum<-summary(GCox)
  HR<-round(Gsum$coefficients[,2],2)
  PValue<-round(Gsum$coefficients[,5],4)
  CI<-paste0(round(Gsum$conf.int[,3:4],2),collapse = "-")
  uniouttalbe<-data.frame('Characteristics'=x,
                          'Hazard Ratio'=HR,
                          'CI95'=CI,
                          'P Value'=PValue)
  return(uniouttalbe)
}
uniouttalbe(gene)

VarNames<-colnames(unicox)[c(2:7,10)]
Univar<-lapply(VarNames, uniouttalbe)
Univar<-ldply(Univar,data.frame)
Univar$Characteristics[Univar$P.Value<0.1]

## 2.mutivariate analysis
fml<-as.formula(paste0('BaSurv~',paste0(Univar$Characteristics[Univar$P.Value<0.1],collapse = "+")))
Multicox<-coxph(fml,data = unicox)
MultiSum<-summary(Multicox)

MultiName<-as.character(Univar$Characteristics[Univar$P.Value<0.1])
MHR<-round(MultiSum$coefficients[,2],2)
MPV<-round(MultiSum$coefficients[,5],4)
MCIL<-round(MultiSum$conf.int[,3],2)
MCIU<-round(MultiSum$conf.int[,4],2)
MCI<-paste0(MCIL,"-",MCIU)
Multicox<-data.frame('Characteristics'=MultiName,
                     'Hazard Ratio'=MHR,
                     'CI95'=MCI,
                     'P Value'= MPV)
write.csv(Univar,file = paste0(cancertype,"-",gene,"-","Univar.csv"))
write.csv(Multicox,file = paste0(cancertype,"-",gene,"-","Multicox.csv"))
