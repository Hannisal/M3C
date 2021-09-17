# time:25/06/2021
# Auther:ShuaiWang
## Based on the expression data from the included BRCA patients with Male deleted
rm(list = ls())
library(stringr)
## expression matrix of GEO

# GSE1456
# change GSE.matrix.txt and GPL filename for different ananlysis
exp<-read.table("GEO BRCA/GSE1456-GPL97_series_matrix.txt",comment.char = "!",header = T,sep = "\t")


cli<-read.csv("GEO BRCA/GSE1456-GPL97-clinical.csv",header = T,check.names = F)

anno<-data.table::fread("GEO BRCA/GPL97-17394.txt")

## METTL2A
probe<-anno$ID[str_detect(anno$`Gene Symbol`,"METTL2A")]

exp<-exp[exp$ID_REF==probe,]
rownames(exp)<-"METTL2A";exp<-exp[,-1]

exp<-exp[,colnames(exp)%in%cli$ID]
exp<-exp[,match(cli$ID,colnames(exp))]
# match
mat<-cbind(cli,t(exp))

write.csv(mat,"GSE1456-GPL97_input.csv",row.names = F)




# GSE3494
rm(list = ls())
# change GSE.matrix.txt and GPL filename for different ananlysis
exp<-read.table("GSE3494-GPL97_series_matrix.txt",comment.char = "!",header = T,sep = "\t")

cli<-read.table("clinical.txt",header = T,comment.char = "#",sep = "\t",check.names = F)
anid<-read.csv("GSE3494-GPL97-clinical.csv",header = T,check.names = F)
anid$`INDEX(ID)`<-do.call(rbind,str_split(anid$`INDEX(ID)`," "))[,1]
cli<-cli[match(anid$`INDEX(ID)`,cli$`INDEX (ID)`),]
cli$ID<-anid$ID

anno<-data.table::fread("GPL97-17394.txt")

## METTL2A
probe<-anno$ID[str_detect(anno$`Gene Symbol`,"METTL2A")]

exp<-exp[exp$ID_REF==probe,]
rownames(exp)<-"METTL2A";exp<-exp[,-1]

exp<-exp[,colnames(exp)%in%cli$ID]
exp<-exp[,match(cli$ID,colnames(exp))]
# match
mat<-cbind(cli,t(exp))

write.csv(mat,"GSE3494-GPL97_input.csv",row.names = F)



# GSE4922
rm(list = ls())
# change GSE.matrix.txt and GPL filename for different ananlysis
exp<-read.table("GSE4922-GPL97_series_matrix.txt",comment.char = "!",header = T,sep = "\t")


cli<-read.csv("GSE4922-GPL97-clinical.csv",header = T,check.names = F)

anno<-data.table::fread("GPL97-17394.txt")

## METTL2A
probe<-anno$ID[str_detect(anno$`Gene Symbol`,"METTL2A")]

exp<-exp[exp$ID_REF==probe,]
rownames(exp)<-"METTL2A";exp<-exp[,-1]

exp<-exp[,colnames(exp)%in%cli$ID]
exp<-exp[,match(cli$ID,colnames(exp))]
# match
mat<-cbind(cli,t(exp))

write.csv(mat,"GSE4922-GPL97_input.csv",row.names = F)



# GSE6532
rm(list = ls())
# change GSE.matrix.txt and GPL filename for different ananlysis
exp<-read.table("GSE6532-GPL97_series_matrix.txt",comment.char = "!",header = T,sep = "\t")


cli<-read.csv("GSE6532-GPL97-clinical.csv",header = T,check.names = F)

anno<-data.table::fread("GPL97-17394.txt")

## METTL2A
probe<-anno$ID[str_detect(anno$`Gene Symbol`,"METTL2A")]

exp<-exp[exp$ID_REF==probe,]
rownames(exp)<-"METTL2A";exp<-exp[,-1]

exp<-exp[,colnames(exp)%in%cli$ID]
exp<-exp[,match(cli$ID,colnames(exp))]
# match
mat<-cbind(cli,t(exp))

write.csv(mat,"GSE6532-GPL97_input.csv",row.names = F)


# recycle
GSE_id<-c("GSE1456","GSE4922","GSE6532")

#exp_files = dir("GSE1456/",pattern = "*series_matrix.txt$",recursive = T)
#anno_files = dir("GPL/",pattern = "*.txt$",recursive = T)
i<-GSE_id[2]
j<-exp_files[1]
for (i in GSE_id) {
  exp_files = dir(paste0(i,"/"),pattern = "*series_matrix.txt$",recursive = T)
  anno_files = dir("GPL/",pattern = "*.txt$",recursive = T)
  for (j in exp_files) {
    anno_files0<-anno_files[str_detect(j,do.call(rbind,str_split(anno_files,pattern = "-"))[,1])]
    clinical_files<-paste0(str_split(j,pattern = "_")[[1]][1],"-clinical.csv")
    ## expression matrix of GEO
    exp<-read.table(paste0(i,"/",j),comment.char = "!",header = T,sep = "\t")
    
    cli<-read.csv(paste0(i,"/",clinical_files),header = T)
    
    anno<-data.table::fread(paste0("GPL/",anno_files0))
    
    ## METTL2A
    probe<-anno$ID[str_detect(anno$`Gene Symbol`,"METTL2A")]
    
    exp<-exp[exp$ID_REF==probe,]
    rownames(exp)<-"METTL2A";exp<-exp[,-1]
    exp<-exp[,colnames(exp)%in%cli$ID]
    exp<-exp[,match(cli$ID,colnames(exp))]
    # match
    mat<-cbind(cli,t(exp))
    
    write.csv(mat,paste0(str_split(j,pattern = "_")[[1]][1],"_input.csv"),row.names = F)
  }
}
