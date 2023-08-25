library(gplots)
library(MASS)
library(DESeq2)
library(knitr)
library(plyr)
library(psych)



wdir <- "C:/RSYNC/LUMIER_C2H2/Archive_oct032018/Archive1_3/Results_oct032018_control_distribution_normalization"
#outdir <- file.path(wdir,"Median_subtraction_zscore") ## not working properly
outdir <- file.path(wdir,"Median_subtraction")
dir.create(outdir)
setwd(outdir)


#or_file <- file.path(wdir, "table_for_zscore.tab")  ## original lumier score:"table_for_original_lum.tab"
or_file <- file.path(wdir, "table_for_original_lum.tab")  
or_table <- read.table(or_file, header=T, sep="\t")
or_table <- log(or_table)
baits <- row.names(or_table)
baits <- baits[!baits %in% c("H2O-GFP", "Flag-GFP")]

dummy <- c("Flag", "H2O")
example_baits <- c("ZNF391","ZNF320","ZFP161","ZNF16","ZNF483")
mixed_baits <- c(dummy, example_baits)
or_table <- or_table[baits,]
dim(or_table)
head(or_table)

test_df <- or_table[mixed_baits,]
a <- apply(test_df, 1, median)
b <- apply(test_df, 1, min)
c <- apply(test_df, 1, max)
d <- apply(test_df, 1, sd)

c - a
log(c-a)
log(c)-log(a)

rcol <- sample(seq(1,300, 1),6)

medians_or <- apply(or_table, 2, median)
medians_or <- na.omit(medians_or)
median_of_medians <- median(medians_or)
or_table[is.na(or_table)] <- median_of_medians

pdf("original_lum_pairwise_bait_correlation.pdf", width=10, height=10)
pairs.panels(t(or_table)[,mixed_baits], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
dev.off()

pdf("original_lum_pairwise_prey_correlation.pdf", width=10, height=10)
pairs.panels(or_table[,rcol], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
dev.off()

### subtract prey median
or_table_minus_col_median <- t(t(or_table) - medians_or)

test_df <- or_table_minus_col_median[mixed_baits,]
a <- apply(test_df, 1, median)
b <- apply(test_df, 1, min)
c <- apply(test_df, 1, max)
d <- apply(test_df, 1, sd)
#or_table_minus_col_median <- apply(or_table, 2, scale)
#row.names(or_table_minus_col_median) <- row.names(or_table)
pdf("prey_median_subtracted_lum_pairwise_bait_correlation.pdf", width=10, height=10)
pairs.panels(t(or_table_minus_col_median)[,mixed_baits], 
             method = "spearman", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
dev.off()

pdf("prey_median_subtracted_lum_pairwise_prey_correlation.pdf", width=10, height=10)
pairs.panels(or_table_minus_col_median[,rcol], 
             method = "spearman", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
dev.off()

### subtract bait median after subtract prey median
row_medians <- apply(or_table_minus_col_median, 1, median)

or_table_minus_row_median <- or_table_minus_col_median - row_medians
pdf("preyandbait_median_subtracted_lum_pairwise_bait_correlation.pdf", width=10, height=10)
pairs.panels(t(or_table_minus_row_median)[,mixed_baits], 
             method = "spearman", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
dev.off()
pdf("preyandbait_median_subtracted_lum_pairwise_prey_correlation.pdf", width=10, height=10)
pairs.panels(or_table_minus_row_median[,rcol], 
             method = "spearman", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
dev.off()

pairs(t(or_table[dummy,]))
pairs(t(or_table_minus_col_median[dummy,]))
pairs(t(or_table_minus_row_median[dummy,]))

or_table[dummy, order(or_table["H2O",])]
hist(medians_or)
or_table_minus_col_median[dummy, order(or_table["H2O",])]

pdf("original_lum_heatmap.pdf", width=40, height=40)
heatmap.2(as.matrix(or_table), trace="none", srtCol=45, cexCol=0.5)#, Rowv=F, Colv=F, dendrogram = "none")
dev.off()
pdf("prey_median_subtracted_lum_heatmap.pdf", width=40, height=40)
heatmap.2(as.matrix(or_table_minus_col_median), trace="none", srtCol=45, cexCol=0.5)#, Rowv=F, Colv=F, dendrogram = "none")
dev.off()
pdf("prey_and_bait_median_subtracted_lum_heatmap.pdf", width=40, height=40)
heatmap.2(as.matrix(or_table_minus_row_median), trace="none", srtCol=45, cexCol=0.5)#, Rowv=F, Colv=F, dendrogram = "none")
dev.off()




################## convert matrix to network  ###############
sample <- colnames(or_table)
prey_names <- unique(unlist(lapply(sample, function(x){unlist(strsplit(x, "_", fixed=T))[1]})))
treatment <- unique(unlist(lapply(sample, function(x){unlist(strsplit(x, "_", fixed=T))[2]})))
prey_names
treatment

setwd(wdir)
for(treat in treatment){
  #treat <- "NT"
  inte_df <- NULL
  
  select <- or_table_minus_col_median[,sample[grep(paste(treat, "_",sep=""), sample, fixed=T)]]
  for(bait_ind in 1:length(baits)){
    bait <- as.character(baits[bait_ind])
    
    print(c(treat, bait))
    for(prey_col in colnames(select)){
      prey <- unlist(strsplit(prey_col, "_", fixed=T))[1]
      rep <- unlist(strsplit(prey_col, "_", fixed=T))[3]
      rep <- as.numeric(rep) - 1
      rep_id <- paste("rep", rep, sep="")
      lum <- select[bait, prey_col]
      inte_df <- rbind(inte_df, c(prey, rep_id, bait, lum))
    }
    
  }
  #out_file <- paste(treat, "replicates_only_subtraction_normalized_zscore.tab", sep="_")
  #write.table(inte_df, out_file, row.names=F, col.names=c("prey", "rep_id", "bait", "normalized_zscore"), quote=F, sep="\t")
  out_file <- paste(treat, "replicates_only_subtraction_normalized_lum.tab", sep="_")
  write.table(inte_df, out_file, row.names=F, col.names=c("prey", "rep_id", "bait", "normalized_lum"), quote=F, sep="\t")
}


