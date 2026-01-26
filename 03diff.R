 #Install packages if not already installed
 # if (!requireNamespace("BiocManager", quietly = TRUE))
  #   install.packages("BiocManager")
 # BiocManager::install("limma")
 # install.packages("pheatmap")
  #install.packages("ggplot2")

# Load libraries
library(limma)
library(pheatmap)
library(ggplot2)

inputFile = "geneMatrix.txt"     # Expression data file
conFile = "s1.txt"               # Control group sample info file
treatFile = "s2.txt"             # Experimental group sample info file
logFCfilter = 2                  # logFC filtering condition (logFC=0.585, 1.5x change; logFC=1, 2x change; logFC=2, 4x change)
adj.P.Val.Filter = 0.05          # Adjusted p-value filtering condition

geoID = "GSE116520"                # GEO database study ID
setwd("D:\\Glioblastoma\\01download\\GSE116520")      #���ù���Ŀ¼a

# Read input file and organize data
rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[, 1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
rt = avereps(data)

# If data is not log2 transformed, automatically log2 transform
qx = as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC = ((qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0))
if (LogC) {
    rt[rt < 0] = 0
    rt = log2(rt + 1)
}
data = normalizeBetweenArrays(rt)

# Read sample information files (control group and experimental group)
sample1 = read.table(conFile, header = FALSE, sep = "\t", check.names = FALSE)
sample2 = read.table(treatFile, header = FALSE, sep = "\t", check.names = FALSE)
sampleName1 = gsub("^ | $", "", as.vector(sample1[, 1]))
sampleName2 = gsub("^ | $", "", as.vector(sample2[, 1]))
conData = data[, sampleName1]
treatData = data[, sampleName2]
data = cbind(conData, treatData)
conNum = ncol(conData)
treatNum = ncol(treatData)

# Differential analysis
Type = c(rep("con", conNum), rep("treat", treatNum))
design <- model.matrix(~0 + factor(Type))
colnames(design) <- c("con", "treat")
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(treat - con, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Output differential results for all genes
allDiff = topTable(fit2, adjust = 'fdr', number = 200000)
allDiffOut = rbind(id = colnames(allDiff), allDiff)
write.table(allDiffOut, file = paste0(geoID, ".all.txt"), sep = "\t", quote = FALSE, col.names = FALSE)

# Output normalized expression data for all genes
Type = c(rep("Control", conNum), rep("Treat", treatNum))
outData = rbind(id = paste0(colnames(data), "_", Type), data)
write.table(outData, file = paste0(geoID, ".normalize.txt"), sep = "\t", quote = FALSE, col.names = FALSE)

# Filter differential results for significant genes
diffSig = allDiff[with(allDiff, (abs(logFC) > logFCfilter & adj.P.Val < adj.P.Val.Filter)), ]
diffSigOut = rbind(id = colnames(diffSig), diffSig)
write.table(diffSigOut, file = paste0(geoID, ".diff.txt"), sep = "\t", quote = FALSE, col.names = FALSE)

# Identify up-regulated and down-regulated genes
upRegulated = diffSig[diffSig$logFC > logFCfilter, ]
downRegulated = diffSig[diffSig$logFC < -logFCfilter, ]

# Save up-regulated genes to up.txt
write.table(upRegulated, file = "up.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Number of up-regulated genes: ", nrow(upRegulated), "\n")
cat("Up-regulated genes:\n", paste(rownames(upRegulated), collapse = ", "), "\n")

# Save down-regulated genes to down.txt
write.table(downRegulated, file = "down.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Number of down-regulated genes: ", nrow(downRegulated), "\n")
cat("Down-regulated genes:\n", paste(rownames(downRegulated), collapse = ", "), "\n")

# Plot heatmap of differential genes
geneNum = 20    # Set number of genes
diffSig = diffSig[order(as.numeric(as.vector(diffSig$logFC))), ]
diffGeneName = as.vector(rownames(diffSig))
diffLength = length(diffGeneName)
hmGene = c()
if (diffLength > (2 * geneNum)) {
    hmGene = diffGeneName[c(1:geneNum, (diffLength - geneNum + 1):diffLength)]
} else {
    hmGene = diffGeneName
}
hmExp = data[hmGene, ]
# Define annotation file
Type = c(rep("Control", conNum), rep("Treat", treatNum))
names(Type) = colnames(data)
Type = as.data.frame(Type)
# Plot heatmap
pdf(file = paste0(geoID, ".heatmap.pdf"), width = 9, height = 7)
pheatmap(hmExp,
         annotation = Type,
         color = colorRampPalette(c("#d596d8", "white", "#e64b35"))(50),
         cluster_cols = F,
         show_colnames = F,
         scale = "row",
         fontsize = 7,
         fontsize_row = 5,
         fontsize_col = 7)
dev.off()

# Define significance
allDiff$logFC[allDiff$logFC > 20] = 20
allDiff$logFC[allDiff$logFC < -20] = -20
Significant = ifelse((allDiff$adj.P.Val < adj.P.Val.Filter & abs(allDiff$logFC) > logFCfilter), ifelse(allDiff$logFC > logFCfilter, "Up", "Down"), "Not")

# Count number of up and down regulated genes for annotation
upCount = sum(Significant == "Up")
downCount = sum(Significant == "Down")

# Plot volcano plot
p = ggplot(allDiff, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(col = Significant), alpha = 0.5) +
    scale_color_manual(values = c("#d596d8", "gray", "#76df76")) +
    scale_size_continuous(range = c(1, 3)) +  # Reduce bubble size range
    labs(title = "Differential Gene Expression", 
         x = "Log Fold Change (logFC)", 
         y = "-log10(Adjusted P-Value)") +
    annotate("text", x = max(allDiff$logFC), y = max(-log10(allDiff$adj.P.Val)), label = paste("Up:", upCount), hjust = 1, vjust = -0.4, color = "#76df76") +
    annotate("text", x = min(allDiff$logFC), y = max(-log10(allDiff$adj.P.Val)), label = paste("Down:", downCount), hjust = 0, vjust = -0.4, color = "#d596d8") +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.position = "none",  # Hide the legend
          plot.margin = unit(c(1, 1, 1, 1), "cm"))  # Adjust margins for a better layout
p = p + theme_bw()

# Output volcano plot
pdf(paste0(geoID, ".vol.pdf"), width = 5.5, height = 5)
print(p)
dev.off()





# Plot bubble plot (for Up and Down-regulated genes)
p = ggplot(allDiff, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(size = abs(logFC), col = Significant), alpha = 0.6) +  # Small bubbles for logFC magnitude
    scale_color_manual(values = c("#d596d8", "gray", "#76df76")) +  # Colors for up, down, and non-significant
    scale_size_continuous(range = c(1, 3)) +  # Reduce bubble size range
    labs(title = "Differential Gene Expression", 
         x = "Log Fold Change (logFC)", 
         y = "-log10(Adjusted P-Value)") +
    annotate("text", x = max(allDiff$logFC), y = max(-log10(allDiff$adj.P.Val)), 
             label = paste("Up:", upCount), hjust = 1, vjust = -0.7, color = "#76df76", size = 3) +
    annotate("text", x = min(allDiff$logFC), y = max(-log10(allDiff$adj.P.Val)), 
             label = paste("Down:", downCount), hjust = 0, vjust = -0.7, color = "#d596d8", size = 3) +
    theme_bw() +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.position = "none",  # Hide the legend
          plot.margin = unit(c(1, 1, 1, 1), "cm"))  # Adjust margins for a better layout

# Output the bubble plot
pdf(paste0(geoID, ".bubble_plot_adjusted.pdf"), width = 6, height = 6)
print(p)
dev.off()