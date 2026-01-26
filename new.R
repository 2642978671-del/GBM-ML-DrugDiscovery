#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install(c("clusterProfiler", "enrichplot", "org.Hs.eg.db", "WGCNA", "GSEABase"))

#install.packages(c("ggplot2", "limma", "pheatmap", "ggsci", "dplyr"))

library(ggplot2)             #引用包
logFCfilter=2                #logFC过滤条件
adj.P.Val.Filter=0.05        #矫正后的p值过滤条件
inputFile="all.txt"          #输入文件
setwd("D:\\Glioblastoma\\06volcano")      #设置工作目录

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
#定义显著性
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")
rt=cbind(rt, Sig=Sig)

# **【修改部分 1: 计算数量】**
sig_counts = table(Sig)
up_count = sig_counts["Up"]
down_count = sig_counts["Down"]
not_count = sig_counts["Not"]
# 创建标题文本
plot_title = paste0("DEGs (Up:", up_count, " | Down:", down_count, ")")

# **【修改部分 2: 绘制火山图，并将数量加入标题】**
p=ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
    geom_point(aes(col=Sig))+
    scale_color_manual(values=c("green", "black", "red"))+
    xlim(-5,5)+
    labs(title = plot_title, 
         x = expression(log[2](FC)), 
         y = expression(-log[10](adj.P.Val)))+ # 添加轴标签，使用纯文本请自行调整
    geom_vline(xintercept=c(-logFCfilter,logFCfilter), col="blue", cex=1, linetype=2)+
    geom_hline(yintercept= -log10(adj.P.Val.Filter), col="blue", cex=1, linetype=2)+
    theme(plot.title=element_text(size=16, hjust=0.5, face="bold"))
p=p+theme_bw()

#输出火山图
pdf(file="volcano.pdf", width=6, height=5.1)
print(p)
dev.off()

# **【修改部分 3: 控制台输出数量】**
print("--- 差异基因计数 ---")
print(sig_counts)
print(paste0("总计筛选出显著DEGs: ", up_count + down_count, "个"))