# Install required packages
#install.packages("ggplot2")  # For data visualization
#install.packages("stringr")  # For string manipulation

# Bioconductor packages
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
 # install.packages("BiocManager")
#}


#BiocManager::install("org.Hs.eg.db")       # Annotation package for human genes
#BiocManager::install("clusterProfiler")     # For enrichment analysis
#BiocManager::install("enrichplot")          # For enrichment analysis visualizations
#BiocManager::install("DOSE")                # For Disease Ontology Semantic enrichment analysis
#BiocManager::install("ggnewscale")          # For handling multiple color scales in ggplot

library("org.Hs.eg.db")  
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("ggnewscale")
library("enrichplot")
library("DOSE")
library(stringr)

setwd("D:\\Glioblastoma2\\09GO_KEGG\\GSE116520")


pvalueFilter=0.05         
qvalueFilter=1  
showNum=7

rt=read.table("core_hub.txt",sep="\t",check.names=F,header=F)      
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID
gene=unique(gene)

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}


kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

write.table(GO,file="GO.xls",sep="\t",quote=F,row.names = F)


if(nrow(GO)<30){
  showNum=nrow(GO)
}

pdf(file="GO_barplot.pdf",width = 9,height =7)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
print(bar)
dev.off()


pdf(file="GO_bubble.pdf",width = 9,height =7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
print(bub)
dev.off()
