library(gplots)
library(MASS)
library(DESeq2)
library(knitr)
library(plyr)
library(psych)
library(GPArotation)
library(data.table)

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
outdir <- file.path(wdir,"Clustering")
dir.create(outdir)
setwd(outdir)


or_file <- file.path(wdir, "table_for_original_lum.tab")  ## original lumier score:"table_for_original_lum.tab"
or_table <- read.table(or_file, header=T, sep="\t")
or_table <- log(or_table)

baits <- row.names(or_table)
baits <- baits[!baits %in% c("Control1", "Control2", "Control3", "Control4")] #"Trim28", "ZNF777", "H2O", "Flag", "H2O-GFP", "Flag-GFP")] # "Control1", "Control2", "Control3", "Control4",

or_table <- or_table[baits, ]

#hist(or_table[, c("ZNF265_NT_1",  "ZNF265_NT_2", "ZNF768_NT_1", "ZNF768_NT_2")])
or_table_t <- data.frame(t(or_table))

pairs.panels(or_table_t[row.names(or_table_t) %like% "_NT_", c("H2O.GFP", "Flag.GFP", sample(colnames(or_table_t), 2))], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
################## convert matrix to network  ###############
sample <- colnames(or_table)
prey_names <- unique(unlist(lapply(sample, function(x){unlist(strsplit(x, "_", fixed=T))[1]})))
treatment <- unique(unlist(lapply(sample, function(x){unlist(strsplit(x, "_", fixed=T))[2]})))
sample
prey_names
treatment
baits

for(treat in treatment){
  
  inte_df <- NULL
  
  select <- or_table[,sample[grep(paste(treat, "_",sep=""), sample, fixed=T)]]
  select_t <- data.frame(t(select))
  bait_medians <- apply(select, 1, median)
  prey_medians <- apply(select, 2, median)
  bait_means <- apply(select, 1, mean)
  prey_means <- apply(select, 2, mean)
 # m <- mardia(select)
#  kurtosi(select)
#  sk_prey <- skew(select)
#  sk_bait <- skew(select_t)
#  names(sk_prey) <- colnames(select)
#  names(sk_bait) <- colnames(select_t)
  for(bait_ind in 1:length(baits)){
    bait <- as.character(baits[bait_ind])
    
    print(c(treat, bait, bait_ind))
    rep1_info <- NULL
    for(prey_col in colnames(select)){
      prey <- unlist(strsplit(prey_col, "_", fixed=T))[1]
      rep <- unlist(strsplit(prey_col, "_", fixed=T))[3]
      rep <- as.numeric(rep)
      
      lum <- select[bait, prey_col]
      relative_to_other_baits <- (lum - prey_medians[prey_col])*prey_medians[prey_col]
      relative_to_other_preys <- (lum - bait_medians[bait])*bait_medians[bait]
      md_prey <- (prey_means[prey_col] - prey_medians[prey_col])
      md_bait <- (bait_means[bait] - bait_medians[bait])
      
      #median_product <- bait_medians[bait]*prey_medians[prey_col]
      #sd_product <- bait_sds[bait]*prey_sds[prey_col]
      
      av <- c(lum, relative_to_other_baits, relative_to_other_preys)
      if(rep == 1){
        rep1_info <- av
      }else{
        inte_df <- rbind(inte_df, c(paste(bait,prey,sep="_"), rep1_info, av))
      }
      
    }
    
  }
  
  colnames(inte_df) <- c("interaction", "rep1_lum", "rep1_relativeB", "rep1_relativeP",  "rep2_lum", "rep2_relativeB", "rep2_relativeP")
  row.names(inte_df) <- inte_df[,1]
  inte_df <- inte_df[, !colnames(inte_df) %in% "interaction"]
  out_file <- paste(treat, "bait_prey_matrix_for_clustering.tab", sep="_")
  write.table(inte_df, out_file, row.names=T, col.names=T, quote=F, sep="\t")
  
  inte_matrix <- quantile_normalisation(type.convert(as.matrix(inte_df)))
  dim(inte_matrix)
  
  if(1){
  pdf(paste(treat, "feature_correlation.pdf", sep="_"))
  pairs.panels(inte_matrix, 
               method = "pearson", # correlation method
               hist.col = "#00AFBB",
               density = TRUE,  # show density plots
               ellipses = TRUE # show correlation ellipses
  )
  dev.off()
  }
  pc <- pca(inte_matrix, 2, rotate="Varimax")
  pdf(paste(treat, "PCA_plot.pdf", sep="_"))
  biplot(pc)
  dev.off()
  #heatmap.2(inte_matrix, trace="none")
  k <- kmeans(inte_matrix, 7, iter.max=10000, algorithm = "Hartigan-Wong", nstart=10)
  k
  k$size
  k$centers
  row_order <- row.names(inte_matrix)
  clusters <- k$cluster[row_order]
  sumCol <- apply(inte_matrix, 1, sum)
  
  head(inte_matrix)
  negatives <- clusters + 10
  pchs <- rep(8, dim(inte_matrix)[1])
  names(pchs) <- row.names(inte_matrix)
  pchs[grepl("GFP", names(negatives), fixed=T)] <- 20
  negatives[grepl("GFP", names(negatives), fixed=T)] <- "red"
  bait_average <- inte_matrix[,3] + inte_matrix[,6]
  prey_average <- inte_matrix[,2] + inte_matrix[,5]
  
  pdf(paste(treat, "cluster_visualization.pdf", sep="_"))
  plot(bait_average, prey_average, xlab="distance to bait median", ylab="distance to prey median", col=negatives, pch=pchs)
  abline(h=0, v=0, lty=3)
  dev.off()
  
  km <- cbind(inte_matrix, sumCol, clusters)
  km <- km[order(km[,dim(km)[2]], km[,dim(km)[2]-1]), ]
  out_file <- paste(treat, "bait_prey_matrix_kmeans_clusters.tab", sep="_")
  write.table(km, out_file, row.names=T, col.names=NA, quote=F, sep="\t")
  
  
}


