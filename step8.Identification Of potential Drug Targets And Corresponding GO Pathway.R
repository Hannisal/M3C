# time:04/06/2021
# Auther:ShuaiWang
## PRISM
setwd("")
library(stringr)
exp<-read.csv("CCLE_expression.csv",check.names = F,row.names = 1)
colnames(exp)<-do.call(rbind,str_split(colnames(exp)," "))[,1]
sinfo<-read.csv("sample_info.csv",check.names = F,row.names = 1)
table(sinfo$primary_disease)
sinfo<-sinfo[sinfo$primary_disease=="Breast Cancer",]
PRISM_AUC<-read.csv("Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_19Q4.csv",check.names = F,row.names = 1)

## KNN
outtable<-data.frame("ID"=NULL)
for (i in colnames(PRISM_AUC)) {
  outtable<-rbind(outtable,sum(is.na(PRISM_AUC[,i]))>(1440/5))
}
colnames(outtable)<-"ID"
PRISM_AUC<-PRISM_AUC[,!outtable$ID]
library(impute)
mydata<-impute.knn(as.matrix(PRISM_AUC))
mydata<-mydata$data

## Mydata corresponds to the AUC of PRISM
## Matches EXP and PRISM
exp<-exp[rownames(exp)%in%rownames(mydata),]
mydata<-mydata[rownames(mydata)%in%rownames(exp),]
mydata<-mydata[match(rownames(exp),rownames(mydata)),]


exp<-exp[rownames(exp)%in%rownames(sinfo),]
mydata<-mydata[rownames(mydata)%in%rownames(sinfo),]

## METTL2A
exp_mettl2a<-exp[,"METTL2A"]
AUC_Cor<-cbind(exp_mettl2a,mydata)
colnames(AUC_Cor)<-do.call(rbind,str_split(colnames(AUC_Cor)," "))[,1]
save(AUC_Cor,file = "CELL_METTL2Acor_AUC.Rdata")

## GDSC
rm(list = ls())
setwd("")
sinfo<-read.csv("sample_info.csv",check.names = F,row.names = 1)
table(sinfo$primary_disease)
sinfo<-sinfo[sinfo$primary_disease=="Breast Cancer",]

exp<-read.table("Cell_line_RMA_proc_basalExp.txt",header = T,check.names = F,sep = "\t")
exp_mettl2a<-exp[714,]

sum(str_detect(exp$GENE_title,"METTL2A"))
str_detect(exp$GENE_title,"METTL2A")

colnames(exp_mettl2a)<-do.call(rbind,str_split(colnames(exp_mettl2a),"DATA."))[,2]
exp_mettl2a<-exp_mettl2a[,colnames(exp_mettl2a)%in%sinfo$COSMICID]
sinfo<-sinfo[sinfo$COSMICID%in%colnames(exp_mettl2a),]
exp_mettl2a<-exp_mettl2a[,match(sinfo$COSMICID,colnames(exp_mettl2a))]

GDSC_AUC1<-read.csv("Drug_sensitivity_AUC_(Sanger_GDSC1).csv",check.names = F,row.names = 1)
GDSC_AUC2<-read.csv("Drug_sensitivity_AUC_(Sanger_GDSC2).csv",check.names = F,row.names = 1)

GDSC_AUC1<-GDSC_AUC1[rownames(GDSC_AUC1)%in%rownames(GDSC_AUC2),]
GDSC_AUC2<-GDSC_AUC2[rownames(GDSC_AUC2)%in%rownames(GDSC_AUC1),]

GDSC_AUC<-cbind(GDSC_AUC1,GDSC_AUC2)
# knn
outtable<-data.frame("ID"=NULL)
for (i in colnames(GDSC_AUC)) {
  outtable<-rbind(outtable,sum(is.na(GDSC_AUC[,i]))>(790/5))
}
colnames(outtable)<-"ID"
GDSC_AUC<-GDSC_AUC[,!outtable$ID]
library(impute)
mydata<-impute.knn(as.matrix(GDSC_AUC))
mydata<-mydata$data

mydata<-mydata[rownames(mydata)%in%rownames(sinfo),]
mydata<-mydata[match(rownames(sinfo),rownames(mydata)),]

AUC_Cor<-cbind(mydata,t(exp_mettl2a))
colnames(AUC_Cor)[ncol(AUC_Cor)]<-"exp_mettl2a"
colnames(AUC_Cor)<-do.call(rbind,str_split(colnames(AUC_Cor)," "))[,1]
save(AUC_Cor,file="GDSC_METTL2Acor_AUC.Rdata")

## PRISM
rm(list = ls())
load("CELL_METTL2Acor_AUC.Rdata")
PRISM<-as.data.frame(AUC_Cor)

write.csv(t(PRISM),file = "PRISM.AUC.csv")
## logFC of AUC
AUC<-t(PRISM[,-1])
group_list<-as.factor(ifelse(PRISM$exp_mettl2a>median(PRISM$exp_mettl2a),"High","Low"))
group_list<-factor(group_list,levels = c("Low","High"))


table(group_list)
pvalueFilter=0.05
logFCfilter=0.1


outTab=data.frame()
dimnames=list(rownames(AUC),colnames(AUC))
data=matrix(as.numeric(as.matrix(AUC)),nrow=nrow(AUC),dimnames=dimnames)

for(i in rownames(data)){
  MedicineName=i
  rt=rbind(auc=data[i,],group_list=group_list)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(auc ~ group_list, data=rt)
  conGeneMeans=mean(data[i,group_list=="Low"])
  treatGeneMeans=mean(data[i,group_list=="High"])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,group_list=="Low"])
  treatMed=median(data[i,group_list=="High"])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(MedicineName=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

write.table(outTab,file="PRISM-ALL.xls",sep="\t",row.names=F,quote=F)

outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$pValue))<pvalueFilter),]
write.table(outDiff,file="PRISM-diff.xls",sep="\t",row.names=F,quote=F)


## GDSC
rm(list = ls())
load("GDSC_METTL2Acor_AUC.Rdata")
GDSC<-as.data.frame(AUC_Cor)

write.csv(t(GDSC),file = "GDSC.AUC.csv")
## logFC of AUC
AUC<-t(GDSC[,-ncol(GDSC)])
group_list<-as.factor(ifelse(GDSC$exp_mettl2a>median(GDSC$exp_mettl2a),"High","Low"))
group_list<-factor(group_list,levels = c("Low","High"))


table(group_list)
pvalueFilter=0.05
logFCfilter=0.1


outTab=data.frame()
dimnames=list(rownames(AUC),colnames(AUC))
data=matrix(as.numeric(as.matrix(AUC)),nrow=nrow(AUC),dimnames=dimnames)

for(i in rownames(data)){
  MedicineName=i
  rt=rbind(auc=data[i,],group_list=group_list)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(auc ~ group_list, data=rt)
  conGeneMeans=mean(data[i,group_list=="Low"])
  treatGeneMeans=mean(data[i,group_list=="High"])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,group_list=="Low"])
  treatMed=median(data[i,group_list=="High"])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(MedicineName=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

write.table(outTab,file="GDSC-ALL.xls",sep="\t",row.names=F,quote=F)

outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$pValue))<pvalueFilter),]
write.table(outDiff,file="GDSC-diff.xls",sep="\t",row.names=F,quote=F)


## boxplot of AUC
## GDSC
rm(list = ls())
load("GDSC_METTL2Acor_AUC.Rdata")
medicine<-read.table("medicine.txt",header = T)

input<-AUC_Cor[,c(medicine$name,colnames(AUC_Cor)[ncol(AUC_Cor)])]
write.csv(input,file = "GDSC_input.csv",row.names = T,quote = F)

## PRISM
rm(list = ls())
load("CELL_METTL2Acor_AUC.Rdata")
medicine<-read.table("medicine.txt",header = T)

input<-AUC_Cor[,c(medicine$name,colnames(AUC_Cor)[1])]
write.csv(input,file = "PRISM_input.csv",row.names = T,quote = F)