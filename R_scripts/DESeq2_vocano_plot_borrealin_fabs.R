library(DESeq2)
library(knitr)
library(plyr)

clean_cache()

seq_len_file <- "C:/RSYNC/AP_MS_C2H2/uniprot_human_sequence_length.tab"
setwd("C:/RSYNC/LUMIER_C2H2/Borealin")

seq_len_table <- read.delim(seq_len_file, stringsAsFactors=F)
head(seq_len_table)
class(seq_len_table$Entry)


spc_table <- read.csv("borealin_fabs.csv", stringsAsFactors=F)
head(spc_table)
dim(spc_table)
class(spc_table$ProteinID)

# remove duplicated preys
dup <- spc_table$ProteinID[duplicated(spc_table$ProteinID)]
spc_table <- spc_table[!spc_table$GeneID %in% dup,]

# attach protein length to spc table
spc_length_table <- merge(spc_table, seq_len_table, by.x="ProteinID", by.y="Entry")
head(spc_length_table)
dim(spc_length_table)

# remove proteins that are not mapped to Uniprot id
notmapped <- spc_table$ProteinID[!spc_table$ProteinID %in% seq_len_table$Entry]
removed <- c(dup, notmapped)
reason <- c(rep("duplacated", length(dup)), rep("not mapped", length(notmapped)))
rem <- cbind(removed, reason)
write.table(rem, "removed_ids.txt", row.names=F, quote=F)

# normalized spectral count by protein length. The longest protein length is scaled to 1
length_normalized_spc <- round(spc_length_table[,4:21]/t(spc_length_table[,27]/max(spc_length_table[,27])))
head(length_normalized_spc)
row.names(length_normalized_spc) <- spc_length_table$GeneName
head(length_normalized_spc)

# select the first 4 samples and remove rows that spectral counts are 0 in these samples
CDCA8 <- as.matrix(length_normalized_spc[,1:4])
colnames(CDCA8) <- c("X10800.CDCA8_293_2", "X10799.CDCA8_293_1", "X10798.CDCA8_293_2", "X10797.CDCA8_293_1")
head(CDCA8)
meanCDCA8 <- apply(CDCA8, 1, mean)
meanCDCA8 <- meanCDCA8[meanCDCA8>0]
CDCA8 <- CDCA8[row.names(CDCA8) %in% names(meanCDCA8),]
dim(CDCA8)
heatmap(CDCA8)

# Differential analysis with DESeq2
conditions <- factor(c("TREAT","TREAT","CONTROL","CONTROL"))
ddsM <- DESeqDataSetFromMatrix(CDCA8, DataFrame(conditions), ~ conditions)


dds <- DESeq(ddsM)
resultsNames(dds)


res <- results(dds, contrast=c("conditions", "TREAT","CONTROL"))
res <- na.omit(res)

  
plotMA(res, ylim=c(-2,2))
  

summary(res, alpha=0.05)


colData(dds)
metadata(res)

plotDispEsts(dds)

# make vocano plot
alpha <- 0.05 # Threshold on the adjusted p-value
lfc <- 2.0 # Threshold on the log2foldchange
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))

topT <- as.data.frame(res)

#Adjusted P values (FDR Q values)
pdf("vocano_plot.pdf", height=8, width=10)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))

with(subset(topT, padj<0.05 & abs(log2FoldChange)>lfc), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=c(-lfc,lfc), col="black", lty=4, lwd=1.2)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=1.2)
text(topT$log2FoldChange[gn.selected_pos], -log10(topT$padj)[gn.selected_pos], lab=rownames(topT)[gn.selected_pos], cex=1.0)

dev.off()

# select and output significant preys
gn.selected_pos <- topT$log2FoldChange > lfc & topT$padj < alpha 
gn.selected_neg <- topT$log2FoldChange < -lfc & topT$padj < alpha 
gn.selected <- c()
topT_pos <- topT[gn.selected_pos,]
topT_neg <- topT[gn.selected_neg,]

topT_pos <- topT_pos[order(topT_pos$padj),]
topT_neg <- topT_neg[order(topT_neg$padj),]
topT_out <- rbind(topT_pos, topT_neg)
write.table(topT_out, "DESeq2_results.tab", col.names=NA, sep="\t", quote=F)
