# time:04/06/2021
# Auther:ShuaiWang
## Based on the expression data from the included BRCA patients with Male deleted
rm(list = ls())
## preparation for GSEA
load("TCGA-UCSC-BRCA.Rdata")
pheno<-pheno[pheno$Gender_nature2012=="FEMALE",]
exp<-exp[,pheno$sampleID]
sur<-sur[sur$sample%in%pheno$sampleID,]
exp_mettl2a<-as.numeric(exp["METTL2A",])
group_list<-ifelse(exp_mettl2a>median(exp_mettl2a),"High","Low")
ncol(exp)
#[1] 711
nrow(exp)
#[1] 20530
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=2^data-1
exp<-cbind(rownames(exp),Description="BRCA",exp)
colnames(exp)[1]<-"NAME"

write.table(exp,file = "GSEA_input.txt",sep = "\t",row.names = F,quote = F)


write.csv(paste0(group_list,collapse = " "),file = "group_GSEA.csv")

## correlation analysis
## preparation for GSEA
load("TCGA-UCSC-BRCA.Rdata")
pheno<-pheno[pheno$Gender_nature2012=="FEMALE",]
exp<-exp[,pheno$sampleID]
sur<-sur[sur$sample%in%pheno$sampleID,]
exp_mettl2a<-as.numeric(exp["METTL2A",])
exp<-exp[!rownames(exp)%in%"METTL2A",]

outTab=data.frame()
outTab_cor<-data.frame()


## Normality test and Pearson correlation analysis
ZPearson<-function(x,y){
  x<-as.numeric(x)
  y<-as.numeric(y)
  xt<-shapiro.test(x)
  yt<-shapiro.test(y)
  if ((xt$p.value>0.05)&(yt$p.value>0.05)) {
    if(sd(y)>0.01){
      df1=as.data.frame(cbind(x,y))
      corT=cor.test(x,y,method="pearson")
      cor=corT$estimate
      pValue=corT$p.value
      outVector=rbind(outVector,pValue)
      outCor=rbind(outCor,cor)
    }
    else{
      outVector=rbind(outVector,pValue=1)
      outCor=rbind(outCor,cor=0)
    }
  }else{
    if(sd(y)>0.01){
      df1=as.data.frame(cbind(x,y))
      corT=cor.test(x,y,method="spearman")
      cor=corT$estimate
      pValue=corT$p.value
      outVector=rbind(outVector,pValue)
      outCor=rbind(outCor,cor)
    }
    else{
      outVector=rbind(outVector,pValue=1)
      outCor=rbind(outCor,cor=0)
    }
  }
}

shapiro.test(exp_mettl2a)
## expression of METTL2A isn't a normal distribution
outVector=data.frame("METTL2A"=NULL)
outCor=data.frame("METTL2A"=NULL)
for(i in rownames(exp)){
  x=as.numeric(exp_mettl2a)
  y=as.numeric(exp[i,])
  if(sd(y)>0.01){
    df1=as.data.frame(cbind(x,y))
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pValue=corT$p.value
    outVector=rbind(outVector,pValue)
    outCor=rbind(outCor,cor)
  }
  else{
    outVector=rbind(outVector,pValue=1)
    outCor=rbind(outCor,cor=0)
  }
}
outTab=cbind(outCor,outVector)
colnames(outTab)<-c("cor","Pvalue")
rownames(outTab)<-rownames(exp)
write.csv(outTab,file="METTL2Acor.result.csv",row.names=T,quote=F)
## GO analysis
library("clusterProfiler")
library("org.Hs.eg.db")
gene<-read.table("Enrichment_table.txt",header = T,sep = "\t",check.names = F)
gene<-bitr(gene$`Gene Symbol`,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
ego<-enrichGO(gene = gene$ENTREZID,
              OrgDb = org.Hs.eg.db,
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05,
              ont = "ALL",
              readable = T)
library(ggplot2)

class(ego)
head(ego)
write.csv(ego[1:308],file = "goresult.csv")

# count排序前30个
go<-read.csv("goresult.csv",header = T,row.names = 1)
GO_term_order=factor(1:30,labels = go$Description)
COLS<-c('#66C3A5','#8DA1CB','#FD8D62')
ggplot(data = go,
       aes(x=GO_term_order,y=Count,
           fill= ONTOLOGY))+
  geom_bar(stat = "identity",
           width=0.8)+
  scale_fill_manual(values = COLS)+
  theme_bw()+
  xlab('GO term')+ylab("Num of Genes")+
  labs(title = "The Top30 Enrichment Go Terms")+
  theme(axis.text.x = element_text(face = "bold",color = "gray50",angle = 70,vjust = 1,hjust = 1))
ggsave("GO.pdf", width = 14, height = 8)
