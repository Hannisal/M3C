# time:20/06/2021
# Auther:ShuaiWang
options(digits = 4)
library(ggplot2)
## load the expression file of four M3C related genes in TCGA database
count_files = dir("M3CGeneExp/",pattern = "^singleGeneExp",recursive = T)

for (i in 1:4) {
  exp<-read.table(count_files[i],header = T)
  exp$Type<-as.factor(exp$Type)
  levels(exp$Type)<-c("Normal","Tumor")
  ylabname <- paste(colnames(exp)[2], "expression(log2(FPKM+1))")
  colnames(exp)[2]<-"Gene"
  # rmove tissue without normal
  exp_withoutNormal<-NULL
  exp_withNormal <- NULL
  for (j in unique(exp$CancerType)) {
    exp1<-exp[exp$CancerType==j,]
    if (length(unique(exp1$Type))==1) {
      exp_withoutNormal<-rbind(exp_withoutNormal,exp1)
    }
      else{
      exp_withNormal<-rbind(exp_withNormal,exp1)
    }
  }
  # calculate p value
  pvalues <- sapply(exp_withNormal$CancerType, function(x) {
    res <- wilcox.test(as.numeric(Gene) ~ Type, data = subset(exp_withNormal, CancerType == x)) # two group，wilcox.test or t.test； multigroups，kruskal.test or aov(one-way ANOVA test)
    res$p.value
  })
  pv <- data.frame(gene = exp_withNormal$CancerType, pvalue = pvalues)
  pv$sigcode <- cut(pv$pvalue, c(0,0.0001, 0.001, 0.01, 0.05, 1), 
                    labels=c('****','***', '**', '*', 'ns'))
  logFoldchanges<-sapply(exp_withNormal$CancerType, function(x) {
    data = subset(exp_withNormal, CancerType == x)
    lfc<-log2(mean(2^data$Gene[data$Type=="Tumor"]-1)/mean(2^data$Gene[data$Type=="Normal"]-1))
  })
  lf <- data.frame(gene = exp_withNormal$CancerType, logFoldchange = logFoldchanges)
  lf$logFoldchange<-sprintf("%0.2f",lf$logFoldchange)
  # box plot
  p.box <- ggplot(exp_withNormal, aes(x=CancerType, y=Gene, color=Type, fill=Type)) +
    geom_boxplot(alpha = .5) + 
    theme_classic() + 
    scale_fill_manual(values=c(Normal="#3082AD",Tumor="#CC3333")) + 
    scale_color_manual(values=c(Normal="#3082AD",Tumor="#CC3333")) + 
    
    theme(axis.text.x = element_text(colour="black", size = 11,
                                     angle = 45, hjust = .5, vjust = .5)) +
    geom_text(aes(x=gene, y=max(exp_withNormal$Gene) * 1.1,
                  label = pv$sigcode),
              data=pv, 
              inherit.aes=F) +
    geom_text(aes(x=gene,y=max(exp_withNormal$Gene) * 1.2,
                  label = lf$logFoldchange),
              data = lf,
              size = 4,
              inherit.aes = F)+
    ylab(ylabname)
  p.box
  p.box.dot <- p.box + geom_point(shape = 21, size=.5,
                                  position = position_jitterdodge(),
                                  alpha = .5) 
  p.box.dot
  ggsave(paste0(ylabname,"ggplotboxDot.pdf"), width = 14, height = 5)
}
