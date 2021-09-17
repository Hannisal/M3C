# time:02/06/2021
# Auther:ShuaiWang
## Based on the expression data from the included BRCA patients with Male deleted and the somatic mutation files for this patients
rm(list = ls())
library(maftools)
load("TCGA-UCSC-BRCA.Rdata")
pheno<-pheno[pheno$Gender_nature2012=="FEMALE",]
exp<-exp[,pheno$sampleID]
sur<-sur[sur$sample%in%pheno$sampleID,]
## Preparation for Genetic alteration analysis

maf <- "1.BRCA/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf"
laml = read.maf(maf = maf,isTCGA = T)

### set group based on the expression of METTL2A in the cohort built(n=711)
exp<-t(exp)
highmettl2a<-rownames(exp)[exp[,"METTL2A"]>median(as.numeric(exp[,"METTL2A"]))]
lowmettl2a<-rownames(exp)[!exp[,"METTL2A"]>median(as.numeric(exp[,"METTL2A"]))]

### subsetMaf
laml.high <- subsetMaf(maf=laml, tsb=highmettl2a, isTCGA=T)
laml.low <- subsetMaf(maf=laml, tsb=lowmettl2a, isTCGA=T)
### mafCompare
hvsl <- mafCompare(m1=laml.high, laml.low, m1Name="METTL2A High Epression", m2Name="METTL2A Low expression", minMut=5)

### output
write.table(hvsl$results, file="high_vs_low.tsv", quote=FALSE, row.names=FALSE, sep="\t")
write.table(hvsl$results, file="high_vs_low.csv", quote=FALSE, row.names=FALSE, sep=",")
pdf("mafforest.pdf", width = 8, height = 5)
forestPlot(mafCompareRes=hvsl, pVal=0.005, color=c("maroon", "royalblue"), geneFontSize=0.8)
dev.off()


pdf("mafh.pdf", width = 14, height = 10)
oncoplot(maf=laml.high,top=30,fontSize = 0.6,showTumorSampleBarcodes = F,drawColBar = F)
dev.off()
pdf("mafl.pdf", width = 14, height = 10)
oncoplot(maf=laml.low,top=30,fontSize = 0.6,showTumorSampleBarcodes = F,drawColBar = F)
dev.off()
