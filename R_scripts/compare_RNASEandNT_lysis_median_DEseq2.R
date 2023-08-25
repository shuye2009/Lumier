library(gplots)
library(MASS)
library(DESeq2)
library(knitr)
library(plyr)
library(psych)

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

wdir <- "C:/RSYNC/LUMIER_C2H2/Archive_oct032018/Archive1_3/Results_oct032018_control_distribution_normalization"
outdir <- file.path(wdir,"DEseq2")
dir.create(outdir)
setwd(outdir)

for(treat in c("NT", "RN", "SNDREAD")){ #"DNASE")){ #
  control_file <- file.path(wdir, paste(treat,"control_of_all_96well_plates.tab", sep="_")) ## scaled by each 96well plate's median
  control_table <- read.table(control_file, header=T, sep="\t")
  
  pdf(paste(treat,"distribution_of_controls.pdf", sep="_"))
  par(mfrow=c(2,2))
  for(i in 1:length(colnames(control_table))){
    coll <- colnames(control_table)[i]
    print(hist(control_table[,i], xlab=coll, main=coll))
  }
  dev.off()
}


lum_file <- file.path(wdir, "table_for_zscore.tab") ## 
lum_table <- read.table(lum_file, header=T, sep="\t")
dim(lum_table)
head(lum_table)

baits <- row.names(lum_table)
baits <- baits[!baits %in% c("H2O-GFP", "Flag-GFP")]
lum_table <- lum_table[baits,]
dummy <- c("Flag", "H2O")
example_baits <- c("ZNF391","ZNF320","ZFP161","ZNF16","ZNF483")
mixed_baits <- c(dummy, example_baits)

col_medians <- apply(lum_table, 2, median)
col_medians <- na.omit(col_medians)
row_medians <- apply(lum_table, 1, median)
row_medians <- na.omit(row_medians)
median_of_medians <- median(col_medians)

lum_table[is.na(lum_table)] <- median_of_medians
lum_table_t <- t(lum_table)
#pairs(lum_table_t[,dummy])

hist(row_medians)
hist(col_medians)

col_sds <- apply(lum_table, 2, sd)
col_sds <- na.omit(col_sds)
row_sds <- apply(lum_table, 1, sd)
row_sds <- na.omit(row_sds)
sd_of_sds <- sd(col_sds)


#pairs(lum_table_t[,dummy])

hist(row_sds)
hist(col_sds)
sort(row_sds)
sort(col_sds)


pairs.panels(t(lum_table)[,mixed_baits], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)


lum_table_minus_col_median <- t(t(lum_table) - col_medians)
heatmap.2(as.matrix(lum_table_minus_col_median), trace="none", srtCol=45, cexCol=0.5)#, Rowv=F, Colv=F, dendrogram = "none")

lum_table_minus_row_median <- lum_table - row_medians
heatmap.2(as.matrix(lum_table_minus_row_median), trace="none", srtCol=45, cexCol=0.5)#, Rowv=F, Colv=F, dendrogram = "none")

lum_table_summed <- lum_table_minus_col_median + lum_table_minus_row_median
heatmap.2(as.matrix(lum_table_summed), trace="none", srtCol=45, cexCol=0.5, Rowv=F, Colv=F, dendrogram = "none")


pairs.panels(t(lum_table_summed)[,mixed_baits], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)



or_file <- file.path(wdir, "table_for_original_lum.tab")  ## original lumier score
or_table <- read.table(or_file, header=T, sep="\t")
or_table <- log(or_table)
dim(or_table)
head(or_table)

medians_or <- apply(or_table, 2, median)
medians_or <- na.omit(medians_or)
median_of_medians <- median(medians_or)
or_table[is.na(or_table)] <- median_of_medians
boxplot(or_table[,1:12])

pairs.panels(t(or_table)[,mixed_baits], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

pdf("original_lum_heatmap.pdf", width=40, height=40)
heatmap.2(as.matrix(or_table), trace="none", srtCol=45, cexCol=0.5, Rowv=F, Colv=F, dendrogram = "none")
dev.off()
pdf("normalized_lum_heatmap.pdf", width=40, height=40)
heatmap.2(as.matrix(lum_table), trace="none", srtCol=45, cexCol=0.5, Rowv=F, Colv=F, dendrogram = "none")
dev.off()

#plot(medians_or, medians, ylim=c(0,145000))

qor_table <- quantile_normalisation(or_table)  ## quantile normalized over all baits of each replicate on 4 96well plates

medians_qor <- apply(qor_table, 2, median)
medians_qor
boxplot(qor_table[,1:12])

qlum_table <- quantile_normalisation(lum_table)  ## quantile normalized over all baits of each replicate on 4 96well plates

medians_lum <- apply(qlum_table, 2, median)
medians_lum
boxplot(qlum_table[,1:12])

pdf("original_lum_quantile_heatmap.pdf", width=40, height=40)
heatmap.2(as.matrix(qor_table), trace="none", srtCol=45, cexCol=0.5)
dev.off()
pdf("normalized_lum_quantile_heatmap.pdf", width=40, height=40)
heatmap.2(as.matrix(qlum_table), trace="none", srtCol=45, cexCol=0.5)
dev.off()

################## check medians  ###############
sample <- colnames(or_table)
prey_names <- unique(unlist(lapply(sample, function(x){unlist(strsplit(x, "_", fixed=T))[1]})))
treatment <- unique(unlist(lapply(sample, function(x){unlist(strsplit(x, "_", fixed=T))[2]})))
prey_names
treatment

if(0){
table_list <- list(or_table, qor_table, lum_table, qlum_table)
names(table_list) <- c("original", "quantile", "zscore", "quantile_zscore")
for(dataName in c("original", "quantile")){
  df <- table_list[[dataName]]
  print(dataName)
 

  for(treat in treatment[2:length(treatment)]){
    base_treat <- treatment[1]
  
    pdf(paste("vocano_plot", treat, "vs", base_treat, dataName, "normalized.pdf",sep="_"))
    par(mfrow=c(3,3))
    DE_table <- data.frame()
    for(prey in sort(prey_names)){
      #prey <- "ZNF669"
      print(prey)
      select <- round(df[,sample[grep(paste(prey, "_",sep=""), sample, fixed=T)]])
      
      conditions <- factor(sort(rep(treatment, 2)))
      ddsM <- DESeqDataSetFromMatrix(select, DataFrame(conditions), ~ conditions)
      
      
      dds <- DESeq(ddsM)
      resultsNames(dds)
      
      
      res <- results(dds, contrast=c("conditions", base_treat, treat))
      res <- na.omit(res)
      
      res
      #plotMA(res, ylim=c(-6,6))
      
      
      summary(res, alpha=0.05)
      
      
      colData(dds)
      metadata(res)
      
      #plotDispEsts(dds)
      
      # make vocano plot
      alpha <- 0.2 # Threshold on the adjusted p-value
      lfc <- 2.0 # Threshold on the log2foldchange
      cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
      
      topT <- as.data.frame(res)
      gn.selected_pos <- subset(topT, padj<alpha & log2FoldChange > lfc)
      gn.selected_neg <- subset(topT, padj<alpha & log2FoldChange < -lfc)
     
      #pdf(paste(prey, "vocano_plot.pdf", sep="_"), height=8, width=10)
      print(with(topT, plot(log2FoldChange, -log10(padj), pch=20, main=prey, cex=1.0, xlim=c(-10, 10), xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value))))
      
      print(with(subset(topT, padj<alpha & abs(log2FoldChange)>lfc), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=1.2)))
      print(abline(v=0, col="black", lty=3, lwd=1.0))
      print(abline(v=c(-lfc,lfc), col="black", lty=4, lwd=1.2))
      print(abline(h=-log10(alpha), col="black", lty=4, lwd=1.2))
      #text(gn.selected_pos$log2FoldChange, -log10(gn.selected_pos$padj), lab=rownames(gn.selected_pos), cex=1.0)
      
      #dev.off()
      DE_results <- rbind(gn.selected_pos, gn.selected_neg)
      preys <- rep(prey, dim(DE_results)[1])
      baits <- row.names(DE_results)
      DE_results <- cbind(baits, preys, DE_results)
      head(DE_results)
      row.names(DE_results) <- NULL
      DE_table <- rbind(DE_table, DE_results)
    }
    dev.off()
    DE_table <- DE_table[order(DE_table$preys),]
    write.table(DE_table, paste("results_of_differential_analysis",treat, "vs", base_treat, dataName, "normalization.tab",sep="_"), row.names=F, sep="\t")
  }
}
}

if(1){
outdir <- file.path(wdir,"SAINT_INT")
dir.create(outdir)
setwd(outdir)
### format quantile normalized data for SAINT score
pl_file <- "C:/RSYNC/AP_MS_C2H2/human_name_length.txt"
protein_length <- read.delim(pl_file, header=F, sep="\t")
head(protein_length)
colnames(protein_length) <- c("protein", "seq.length")
dummy <- c("Flag", "Flag-GFP", "H2O", "H2O-GFP")
saint_table <- or_table
# create bait file
bait_names <- sort(rep(row.names(saint_table),2))
purif <- seq(1, length(bait_names), 1)
purif <- paste(bait_names, purif, sep="_")
group <- rep("T", length(bait_names))
group[bait_names %in% dummy] <- "C"
bait_df <- data.frame(purif, bait_names, group)
write.table(bait_df, "SAINT_bait_file.txt", sep="\t", row.names=F, col.names=F, quote=F)

# create prey file
prey_names <- unique(unlist(lapply(sample, function(x){unlist(strsplit(x, "_", fixed=T))[1]})))
#prey_names <- c(prey_names,toupper(row.names(saint_table)))
prey_names[!prey_names %in% protein_length$protein]
prey_length <- protein_length[protein_length$protein %in% prey_names, ]
dim(prey_length)
prey_df <- data.frame(prey_length$protein, prey_length$protein)
write.table(prey_length, "SAINT_prey_file.txt", sep="\t", row.names=F, col.names=F, quote=F)
#writeLines(prey_names, "SAINT_prey_file.txt")

# create interaction file
for(treat in treatment){
  #treat <- "NT"
  inte_df <- NULL
  
  select <- round(saint_table[,sample[grep(paste(treat, "_",sep=""), sample, fixed=T)]])
  for(bait_ind in 1:dim(bait_df)[1]){
    bait <- as.character(bait_df$bait_names[bait_ind])
    purification <- as.character(bait_df$purif[bait_ind])
    print(c(treat, bait, purification))
    for(prey_col in colnames(select)){
      prey <- unlist(strsplit(prey_col, "_", fixed=T))[1]
      rep <- as.numeric(unlist(strsplit(prey_col, "_", fixed=T))[3])
      if(bait_ind %% 2 == rep %% 2){
        lum <- select[bait, prey_col]
        inte_df <- rbind(inte_df, c(purification, bait, prey, lum))
      }
      
    }
    
  }
  write.table(inte_df, paste(treat, "SAINT_interaction_file.txt", sep="_"), sep="\t", row.names=F, col.names=F, quote=F)
 # system(paste("saint-int-ctrl", paste(treat, "SAINT_interaction_file.txt", sep="_"), "SAINT_prey_file.txt", "SAINT_bait_file.txt"), show.output.on.console=T)
#  saint <- read.table("list.txt", sep="\t", header=T)
  # saint-int-ctrl NT_SAINT_interaction_file.txt SAINT_prey_file.txt SAINT_bait_file.txt
#  saint <- saint[order(saint$SaintScore, decreasing=T),]
#  write.table(saint, paste(treat, "SAINT_score.txt", sep="_"), sep="\t", row.names=F, col.names=T, quote=F)
}

}


#system("perl C:/RSYNC/worksapce2015/MASS_PEC_SCORE/compass.pl RN_table_for_COMPASSscore.tab RN_COMPASS_score.txt", show.output.on.console=T)
#system("perl C:/RSYNC/worksapce2015/MASS_PEC_SCORE/compass.pl SNDREAD_table_for_COMPASSscore.tab SNDREAD_COMPASS_score.txt", show.output.on.console=T)

