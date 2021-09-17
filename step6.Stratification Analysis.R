# time:04/06/2021
# Auther:ShuaiWang
## Based on the expression data from the included BRCA patients with Male deleted
rm(list = ls())
## 1.preparation for stratification analysis
load("TCGA-UCSC-BRCA.Rdata")
pheno<-pheno[pheno$Gender_nature2012=="FEMALE",]
exp<-exp[,pheno$sampleID]
sur<-sur[sur$sample%in%pheno$sampleID,]
exp_mettl2a<-exp["METTL2A",]
unicox<-cbind(pheno[,-c(3,4,11:15)],sur[,2:3],t(exp_mettl2a))
colnames(unicox)<-c("sampleID","Age","T.stage","N.stage","M.stage","ER.status","PR.status","HER2.status","OS","OS.time","METTL2A")
unicox$Age<-ifelse(unicox$Age>median(unicox$Age),paste0(">",median(unicox$Age)),paste0("<=",median(unicox$Age)))

## Age
## 2.bartlett.test for each clinical feature
bartlett.test(METTL2A~Age,unicox)
### p<0.05,heterogeneity of variance
## 3.Shapiro-Wilk test for each clinical feature
### Age<=58
v1<-unicox[unicox$Age=="<=58",ncol(unicox)]
qqnorm(v1)
qqline(v1)
shapiro.test(v1)
### p<0.05,no normal distribution
### Age>58
v2<-unicox[unicox$Age==">58",ncol(unicox)]
qqnorm(v2)
qqline(v2)
shapiro.test(v2)
### p<0.05,no normal distribution

## T.stage
## 2.bartlett.test for each clinical feature
bartlett.test(METTL2A~T.stage,unicox)
### p<0.05,heterogeneity of variance
## 3.Shapiro-Wilk test for each clinical feature
### T1
v1<-unicox[unicox$T.stage=="T1",ncol(unicox)]
qqnorm(v1)
qqline(v1)
shapiro.test(v1)
### p<0.05,no normal distribution
### T2
v2<-unicox[unicox$T.stage=="T2",ncol(unicox)]
qqnorm(v2)
qqline(v2)
shapiro.test(v2)
### p<0.05,no normal distribution
### T3
v3<-unicox[unicox$T.stage=="T3",ncol(unicox)]
qqnorm(v3)
qqline(v3)
shapiro.test(v3)
### p<0.05,no normal distribution
### T4
v4<-unicox[unicox$T.stage=="T4",ncol(unicox)]
qqnorm(v4)
qqline(v4)
shapiro.test(v4)
### p>0.05,normal distribution

## N.stage
## 2.bartlett.test for each clinical feature
bartlett.test(METTL2A~N.stage,unicox)
### p<0.05,heterogeneity of variance
## 3.Shapiro-Wilk test for each clinical feature
### N0
v1<-unicox[unicox$N.stage=="N0",ncol(unicox)]
qqnorm(v1)
qqline(v1)
shapiro.test(v1)
### p<0.05,no normal distribution
### N1
v2<-unicox[unicox$N.stage=="N1",ncol(unicox)]
qqnorm(v2)
qqline(v2)
shapiro.test(v2)
### p<0.05,no normal distribution
### N2
v3<-unicox[unicox$N.stage=="N2",ncol(unicox)]
qqnorm(v3)
qqline(v3)
shapiro.test(v3)
### p<0.05,no normal distribution
### N3
v4<-unicox[unicox$N.stage=="N3",ncol(unicox)]
qqnorm(v4)
qqline(v4)
shapiro.test(v4)
### p>0.05,normal distribution

## M.stage
## 2.bartlett.test for each clinical feature
bartlett.test(METTL2A~M.stage,unicox)
### p>0.05,Homogeneity of variance
## 3.Shapiro-Wilk test for each clinical feature
### M0
v1<-unicox[unicox$M.stage=="M0",ncol(unicox)]
qqnorm(v1)
qqline(v1)
shapiro.test(v1)
### p<0.05,no normal distribution
### M1
v2<-unicox[unicox$M.stage=="M1",ncol(unicox)]
qqnorm(v2)
qqline(v2)
shapiro.test(v2)
### p>0.05,normal distribution

## ER.status
## 2.bartlett.test for each clinical feature
bartlett.test(METTL2A~ER.status,unicox)
### p>0.05,Homogeneity of variance
## 3.Shapiro-Wilk test for each clinical feature
### Negative
v1<-unicox[unicox$ER.status=="Negative",ncol(unicox)]
qqnorm(v1)
qqline(v1)
shapiro.test(v1)
### p>0.05,normal distribution
### Positive
v2<-unicox[unicox$ER.status=="Positive",ncol(unicox)]
qqnorm(v2)
qqline(v2)
shapiro.test(v2)
### p<0.05,normal distribution

## PR.status
## 2.bartlett.test for each clinical feature
bartlett.test(METTL2A~PR.status,unicox)
### p>0.05,Homogeneity of variance
## 3.Shapiro-Wilk test for each clinical feature
### Negative
v1<-unicox[unicox$ER.status=="Negative",ncol(unicox)]
qqnorm(v1)
qqline(v1)
shapiro.test(v1)
### p>0.05,normal distribution
### Positive
v2<-unicox[unicox$PR.status=="Positive",ncol(unicox)]
qqnorm(v2)
qqline(v2)
shapiro.test(v2)
### p<0.05,normal distribution

## HER2.status
## 2.bartlett.test for each clinical feature
bartlett.test(METTL2A~HER2.status,unicox)
### p<0.05,heterogeneity of variance
## 3.Shapiro-Wilk test for each clinical feature
### Negative
v1<-unicox[unicox$ER.status=="Negative",ncol(unicox)]
qqnorm(v1)
qqline(v1)
shapiro.test(v1)
### p>0.05,normal distribution
### Positive
v2<-unicox[unicox$ER.status=="Positive",ncol(unicox)]
qqnorm(v2)
qqline(v2)
shapiro.test(v2)
### p<0.05,normal distribution


## ALL above clinical feature can't fit the condition about normal distribution&Homogeneity of variance

## 2.So the Wilcoxon test, Kruskal-Wallis H test and Nemenyi method are chose
## BaseLine Table
BaseLine<-unicox[,c("sampleID","Age","T.stage","N.stage","M.stage","ER.status","PR.status","HER2.status","OS","OS.time","METTL2A")]

BaseLine$Expression_METTL2A<-ifelse(BaseLine$METTL2A>median(BaseLine$METTL2A),"High","Low")
### Proportion
Base<-subset(BaseLine,select = -c(OS,OS.time,METTL2A))
Freq<-lapply(Base[,2:9],table)
Prop<-lapply(Freq[1:8], prop.table)
### Table
Character<-c(names(Freq[1]),names(Freq[[1]]))
Noc<-c(NA,paste0(Freq[[1]],'(',Prop[[1]],")"))
Characteristics<-data.frame('Characteristics'=Character,'Number of Cases'=Noc)

Char<-NULL
for (i in 1:8) {
  Character<-c(names(Freq[i]),names(Freq[[i]]))
  Noc<-c(NA,paste0(Freq[[i]],'(',round(Prop[[i]],4)*100,")"))
  Characteristics<-data.frame('Characteristics'=Character,'Number of Cases(%)'=Noc)
  Char<-rbind(Char,Characteristics)
}
write.csv(Char,'Characteristics.csv')

## The Association Between METTL2A Expression And Clinical Features Parameters
###     Chi-square test for fourfold table
###     Warning message:In chisq.test(ftable) : Chi-squared approximation may be incorrect
###     Concerning this warning message, the theoretical frequency of male patients is less than 5 and BRCA mainly presents in famale patients. So the male patients are deleted to reduce bias.
###     reture to the obove codes
fst<-NULL
for (j in 1:6) {
  Character<-colnames(Base[,-c(1,3,4)])[j]
  data<-Base[,c(Character,"Expression_METTL2A")]
  colnames(data)<-c("Vector","Expression")
  ftable<-xtabs(~Vector+Expression,data)
  test<-chisq.test(ftable)
  Xsquared<-test$statistic
  p.value<-test$p.value
  ftable<-as.data.frame(ftable)
  
  NocH<-c(NA,ftable[1:2,3])
  NocL<-c(NA,ftable[3:4,3])
  NocX<-c(NA,Xsquared,NA)
  NocP<-c(NA,p.value,NA)
  
  fourtable<-data.frame('Characteristics'=c(Character,as.vector(ftable[1:2,1])),"High expression"=NocH,"Low expression"=NocL,"X-squared"=NocX,"P.value"=NocP)
  fst<-rbind(fst,fourtable)
}

write.csv(fst,'fourtablechisuqare.csv')
###     Chi-square test for R*C table
###     T.stage
data<-Base[,c(3,4,9)]

rctable<-xtabs(~T.stage+Expression_METTL2A,data)
test<-chisq.test(rctable)
Xsquared<-test$statistic
p.value<-test$p.value
rctable<-as.data.frame(rctable)

NocH<-c(NA,rctable[1:4,3])
NocL<-c(NA,rctable[5:8,3])
NocX<-c(NA,Xsquared,NA,NA,NA)
NocP<-c(NA,p.value,NA,NA,NA)

rctable<-data.frame('Characteristics'=c("T.stage",as.vector(rctable[1:4,1])),"High expression"=NocH,"Low expression"=NocL,"X-squared"=NocX,"P.value"=NocP)
rcst<-NULL
rcst<-rbind(rcst,rctable)
###     N.stage
rctable<-xtabs(~N.stage+Expression_METTL2A,data)
test<-chisq.test(rctable)
Xsquared<-test$statistic
p.value<-test$p.value
rctable<-as.data.frame(rctable)

NocH<-c(NA,rctable[1:4,3])
NocL<-c(NA,rctable[5:8,3])
NocX<-c(NA,Xsquared,NA,NA,NA)
NocP<-c(NA,p.value,NA,NA,NA)

rctable<-data.frame('Characteristics'=c("N.stage",as.vector(rctable[1:4,1])),"High expression"=NocH,"Low expression"=NocL,"X-squared"=NocX,"P.value"=NocP)

rcst<-rbind(rcst,rctable)
write.csv(rcst,'rctablechisuqare.csv')


## preparation for boxplot
write.csv(unicox,file = "TCGA-BRCA-Clinicopathological.csv")
