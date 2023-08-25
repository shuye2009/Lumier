# Author: Shuye Pu
# Created on: May 06, 2019
# Last modified on: Aug 24, 2023
# This is the final version for defining interaction networks

library("MASS")

find_threshold <- function(positive, negative, FDR_threshold, plot = FALSE){
  false_positive_rate <- vector()
  true_positive_rate <- vector()
  false_discovery_rate <- vector()
  p <- length(positive)
  n <- length(negative)
  cutoff <- seq(0, max(positive), 0.1)
  zscore_threshold <- max(positive)
  FPR_at_threshold <- 1
  TPR_at_threshold <- 0
  precision_at_threshold <- 0
  for(i in 1:length(cutoff)){
    FP <- length(negative[negative > cutoff[i]])
    TP <- length(positive[positive > cutoff[i]])
    
    false_positive_rate[i] <- FP/n
    true_positive_rate[i] <- TP/p
    false_discovery_rate[i] <- FP/(TP+FP)
    if(false_discovery_rate[i] < FDR_threshold && cutoff[i] < zscore_threshold){
      zscore_threshold <- cutoff[i]
      FPR_at_threshold <- FP/n
      TPR_at_threshold <- TP/p
      precision_at_threshold <- TP/(TP+FP)
    }
  }
  if(plot){
    plot(false_positive_rate, true_positive_rate, xlim=c(0,1), ylim=c(0,1), xlab="FPR", ylab="TPR", main="ROC")
   
    abline(v=FPR_at_threshold)
    abline(0,1, lty=3)
    plot(cutoff, false_discovery_rate, xlab="zscore cutoff", ylab="FDR", main="FDR")
  }
  
  return(c(zscore_threshold, FPR_at_threshold, TPR_at_threshold, precision_at_threshold))
}

#wdir <- "C:/RSYNC/LUMIER_C2H2/Archive_oct032018/Archive1_3/Results_oct032018_control_distribution_normalization"


wdir <- "C:/RSYNC/LUMIER_C2H2/Archive_oct032018/Archive1_3/Results_oct032018_control_distribution_normalization"
outdir <- file.path(wdir,"replicates")
#outdir <- file.path(wdir,"Median_subtraction")

dir.create(outdir)
setwd(outdir)
FDR <- 0.15
dummy <- c("Control1", "Control2", "Control3", "Blank", "Control4", "Flag", "H2O")
NT <- read.table(file.path(wdir,"NT_replicates_only_control_normalized_lum.tab"), header=T)
#NT <- read.table(file.path(wdir,"NT_replicates_only_subtraction_normalized_lum.tab"), header=T)

NT_network_positive <- NT[is.element(NT[,3], "Control4"),]
NT_network_background <- NT[is.element(NT[,3], "Flag"),]
NT_network <- NT[!is.element(NT[,3], dummy),]

nt_pos <- dim(NT_network_positive)[1]

nt_cutoffs <- NULL

iter <- 1000
for(i in 1:iter){
  NT_network_negative <- NT[sample(row.names(NT_network), size=nt_pos),]
  nt_cutoff <-  find_threshold(NT_network_positive[,4], NT_network_negative[,4], FDR)
  nt_cutoffs <- rbind(nt_cutoffs, nt_cutoff)
}

brks <- seq(-5.0, 5.0, 0.1)
cutoff_row <- nt_cutoffs[which(nt_cutoffs[,1] == median(nt_cutoffs[,1]))[1],]
nt_cutoff <- cutoff_row[1]
FPR_cutoff <- cutoff_row[2]
TPR_cutoff <- cutoff_row[3]
precision <- cutoff_row[4]

png("Zscore_cutoff_plot_of_NT.png", 800, 600)
hist(nt_cutoffs[,1], xlab="Zscore_cutoff", main="NT")
#text(2,200, sprintf("FPR = %.3f",FPR_cutoff), adj=c(0,1))
#text(2,180, sprintf("TPR = %.3f",TPR_cutoff), adj=c(0,1))
#text(2,160, sprintf("Precision = %.3f",precision), adj=c(0,1))
legend("topleft", legend=c(sprintf("FPR = %.3f",FPR_cutoff),sprintf("TPR = %.3f",TPR_cutoff), sprintf("Precision = %.3f",precision)))
dev.off()

nt_values <- NT_network[,4]

x <- nt_values
fit <- fitdistr(x, "normal")
class(fit)
para <- fit$estimate
png("Density_plot_of_NT_rep.png", 800, 600)
hist(x, prob = TRUE, main="Density plot of NT")
curve(dnorm(x, para[1], para[2]), col = 3, add=T)
dev.off()

png("QQ_plot_NT_rep.png", 800, 600)
qqnorm(x, main="Normal Q-Q Plot: NT rep")
qqline(x)
dev.off()

png("Density_plot_of_NT_background.png", 800, 600)
plot(NULL, NULL, xlim=c(-5,5), ylim=c(0,1), xlab="Zscore", ylab="Density", main="Density plot of NT with positive and negative controls")
lines(density(NT_network_positive[,4]), col="red", lwd=2)
lines(density(NT_network_background[,4]), col="blue", lwd=2)
lines(density(NT_network[,4]), col="black", lwd=2)
abline(v=nt_cutoff, lty=2)
legend("topleft", legend=c("Positive control", "Negative control", "Sample"), col=c("red", "blue", "black"), lwd=2)
dev.off()

#nt_cutoff <- qnorm(0.95, para[1], para[2])
#curve(dexp(x, para[1]), col = 3, add=T)
#nt_cutoff <- qexp(0.95, para[1])
#nt_cutoff


# display NT reproducibility
NT_network_rep0 <- NT_network[NT_network[,2]=="rep0",]
NT_network_rep1 <- NT_network[NT_network[,2]=="rep1",]

dim(NT_network_rep0)
dim(NT_network_rep1)

rep0_baits <- as.character(unique(NT_network_rep0[,3]))
rep1_baits <- as.character(unique(NT_network_rep1[,3]))
rep0_preys <- as.character(unique(NT_network_rep0[,1]))
rep1_preys <- as.character(unique(NT_network_rep1[,1]))

common_baits <- rep0_baits[rep0_baits %in% rep1_baits]
common_preys <- rep0_preys[rep0_preys %in% rep1_preys]

common_rows_nt0 <- matrix(NA,ncol=length(rep0_preys),nrow=length(rep0_baits), dimnames=list(rep0_baits,rep0_preys))
common_rows_nt1 <- matrix(NA,ncol=length(rep1_preys),nrow=length(rep1_baits), dimnames=list(rep1_baits,rep1_preys))
common_rows_nt <- NULL
for(i in 1:dim(NT_network_rep0)[1]){
  bait <- as.character(NT_network_rep0[i,3])
  prey <- as.character(NT_network_rep0[i,1])
  common_rows_nt0[bait,prey] <- NT_network_rep0[i,4]
  
}

for(i in 1:dim(NT_network_rep1)[1]){
  bait <- as.character(NT_network_rep1[i,3])
  prey <- as.character(NT_network_rep1[i,1])
  common_rows_nt1[bait,prey] <- NT_network_rep1[i,4]
  
}

for(i in 1:length(common_baits)){
  bait <- common_baits[i]
  for(j in 1:length(common_preys)){
    prey <- common_preys[j]
    rep0 <- common_rows_nt0[bait,prey]
    rep1 <- common_rows_nt1[bait,prey]
    if(!is.na(rep0) && !is.na(rep1)){
      common_rows_nt <- rbind(common_rows_nt, t(c(prey, bait, rep0, rep1)))
    }
  }
}
head(common_rows_nt)

png("Reproducibility_of_NT.png", 800, 600)
plot(type.convert(common_rows_nt[,3]), type.convert(common_rows_nt[,4]), xlab="NT_replicate_0", ylab="NT_replicate_1",main="Reproducibility of NT")
abline(lm(type.convert(common_rows_nt[,4]) ~ type.convert(common_rows_nt[,3])), col="red")
abline(v=nt_cutoff,h=nt_cutoff)
ct <- cor.test(type.convert(common_rows_nt[,3]), type.convert(common_rows_nt[,4]))
co <- ct$estimate
text(-4,8, sprintf("r = %.3f",co), adj=c(1,1))
text(-4,7.5, sprintf("cutoff = %.3f",nt_cutoff), adj=c(1,1))
dev.off()

NT_network_thresholded <- common_rows_nt[(as.double(common_rows_nt[, 3]) > nt_cutoff & as.double(common_rows_nt[, 4]) > nt_cutoff), ]
write.table(NT_network_thresholded, paste("NT", "thresholded_network_replicates.tab",sep="_"), col.names=c("Bait", "Prey", paste("NT>",sprintf("%.1f",nt_cutoff),sep=""), paste("NT>",sprintf("%.1f",nt_cutoff),sep="")), row.names=F, sep="\t")
bait_int_count <- as.data.frame(apply(NT_network_thresholded[,1:2], 2, table)[[2]])
prey_int_count <- as.data.frame(apply(NT_network_thresholded[,1:2], 2, table)[[1]])
bait_int_count <- bait_int_count[order(bait_int_count$Freq, decreasing=T),]
prey_int_count <- prey_int_count[order(prey_int_count$Freq, decreasing=T),]
write.table(bait_int_count, paste("NT", "thresholded_network_replicates_bait_int.tab",sep="_"), col.names=c("Bait", "interaction_count"), row.names=F, sep="\t")
write.table(prey_int_count, paste("NT", "thresholded_network_replicates_prey_int.tab",sep="_"), col.names=c("Prey", "interaction_count"), row.names=F, sep="\t")


for(treatment in c("RN", "SNDREAD")){
#for(treatment in c("DNASE")){

  DEseq2_table <- read.table(file.path(wdir, "DEseq2", paste("results_of_differential_analysis",treatment, "vs_NT_original_normalization.tab", sep="_")),header=T)  
  #RN_LYSIS <- read.table(file.path(wdir,paste(treatment,"replicates_only_control_normalized_lum.tab", sep="_")), header=T)
  RN_LYSIS <- read.table(file.path(wdir,paste(treatment,"replicates_only_subtraction_normalized_lum.tab", sep="_")), header=T)
  
  RN_LYSIS_network_positive <- RN_LYSIS[is.element(RN_LYSIS[,3], "Control4"),]
  RN_LYSIS_network_background <- RN_LYSIS[is.element(RN_LYSIS[,3], "Flag"),]
  
  RN_LYSIS_network <- RN_LYSIS[!is.element(RN_LYSIS[,3], dummy),]

  rn_pos <- dim(RN_LYSIS_network_positive)[1]

  rn_cutoffs <- NULL
  
  for(i in 1:iter){
    RN_LYSIS_network_negative <- RN_LYSIS[sample(row.names(RN_LYSIS_network), size=rn_pos),]
    RN_cutoff <-  find_threshold(RN_LYSIS_network_positive[,4], RN_LYSIS_network_negative[,4], FDR)
    rn_cutoffs <- rbind(rn_cutoffs, RN_cutoff)
  }
  
  cutoff_row <- rn_cutoffs[which(rn_cutoffs[,1] == median(rn_cutoffs[,1]))[1],]
  RN_LYSIS_cutoff <- cutoff_row[1]
  FPR_cutoff <- cutoff_row[2]
  TPR_cutoff <- cutoff_row[3]
  precision <- cutoff_row[4]
  
  png(paste("Zscore_cutoff_plot_of_",treatment,".png",sep=""), 800, 600)
  hist(rn_cutoffs[,1], xlab="Zscore_cutoff", main=treatment)
  legend("topleft", legend=c(sprintf("FPR = %.3f",FPR_cutoff),sprintf("TPR = %.3f",TPR_cutoff), sprintf("Precision = %.3f",precision)))
  
  dev.off()
  
  #define RN_LYSIS data cutoff, old code
  
  RN_LYSIS_values <- RN_LYSIS_network[,4]
  
  x <- RN_LYSIS_values
  fit <- fitdistr(x, "normal")
  class(fit)
  para1 <- fit$estimate
  png(paste("Density_plot_of",treatment,"rep.png",sep="_"), 800, 600)
  hist(x, prob = TRUE, main=paste("Density plot of",treatment,sep=" "))
  curve(dnorm(x, para1[1], para1[2]), col = 3, add=T)
  dev.off()
  
  
  png(paste("QQ_plot",treatment,"rep.png",sep="_"), 800, 600)
  qqnorm(x,main=paste("Normal Q-Q Plot:",treatment,"rep",sep=" "))
  qqline(x)
  dev.off()
  
  png(paste("Density_plot_of", treatment, "background.png", sep="_"), 800, 600)
  plot(NULL, NULL, xlim=c(-5,5), ylim=c(0,1), xlab="Zscore", ylab="Density", main="Density plot of positive and negative controls")
  lines(density(RN_LYSIS_network_positive[,4]), col="red", lwd=2)
  lines(density(RN_LYSIS_network_background[,4]), col="blue", lwd=2)
  lines(density(RN_LYSIS_network[,4]), col="black", lwd=2)
  abline(v=RN_LYSIS_cutoff, lty=2)
  legend("topleft", legend=c("Positive control", "Negative control", "Sample"), col=c("red", "blue", "black"), lwd=2)
  dev.off()
  
        # display RN_LYSIS reproducibility
  RN_LYSIS_network_rep0 <- RN_LYSIS_network[RN_LYSIS_network[,2]=="rep0",]
  RN_LYSIS_network_rep1 <- RN_LYSIS_network[RN_LYSIS_network[,2]=="rep1",]
  
  dim(RN_LYSIS_network_rep0)
  dim(RN_LYSIS_network_rep1)
  
  rep0_baits <- as.character(unique(RN_LYSIS_network_rep0[,3]))
  rep1_baits <- as.character(unique(RN_LYSIS_network_rep1[,3]))
  rep0_preys <- as.character(unique(RN_LYSIS_network_rep0[,1]))
  rep1_preys <- as.character(unique(RN_LYSIS_network_rep1[,1]))
  
  common_rows_RN_LYSIS0 <- matrix(NA,ncol=length(rep0_preys),nrow=length(rep0_baits), dimnames=list(rep0_baits,rep0_preys))
  common_rows_RN_LYSIS1 <- matrix(NA,ncol=length(rep1_preys),nrow=length(rep1_baits), dimnames=list(rep1_baits,rep1_preys))
  common_rows_RN_LYSIS <- NULL
  
  for(i in 1:dim(RN_LYSIS_network_rep0)[1]){
    bait <- as.character(RN_LYSIS_network_rep0[i,3])
    prey <- as.character(RN_LYSIS_network_rep0[i,1])
    common_rows_RN_LYSIS0[bait,prey] <- RN_LYSIS_network_rep0[i,4]
    
  }
  
  for(i in 1:dim(RN_LYSIS_network_rep1)[1]){
    bait <- as.character(RN_LYSIS_network_rep1[i,3])
    prey <- as.character(RN_LYSIS_network_rep1[i,1])
    common_rows_RN_LYSIS1[bait,prey] <- RN_LYSIS_network_rep1[i,4]
    
  }
  
  common_baits <- rep0_baits[rep0_baits %in% rep1_baits]
  common_preys <- rep0_preys[rep0_preys %in% rep1_preys]
  common_rows_nt0_RN_LYSIS0 <- NULL
  common_rows_nt1_RN_LYSIS1 <- NULL
  for(i in 1:length(common_baits)){
    bait <- common_baits[i]
    for(j in 1:length(common_preys)){
      prey <- common_preys[j]
      rep0_RN_LYSIS <- common_rows_RN_LYSIS0[bait,prey]
      rep1_RN_LYSIS <- common_rows_RN_LYSIS1[bait,prey]
      if(!is.na(rep0_RN_LYSIS) && !is.na(rep1_RN_LYSIS)){
        common_rows_RN_LYSIS <- rbind(common_rows_RN_LYSIS, t(c(prey, bait, rep0_RN_LYSIS, rep1_RN_LYSIS)))
      }
      rep0_nt <- common_rows_nt0[bait,prey]
      rep1_nt <- common_rows_nt1[bait,prey]
      if(!is.na(rep0_nt) && !is.na(rep0_RN_LYSIS)){
        common_rows_nt0_RN_LYSIS0 <- rbind(common_rows_nt0_RN_LYSIS0, t(c(prey, bait, rep0_nt, rep0_RN_LYSIS)))
      }
      if(!is.na(rep1_nt) && !is.na(rep1_RN_LYSIS)){
        common_rows_nt1_RN_LYSIS1 <- rbind(common_rows_nt1_RN_LYSIS1, t(c(prey, bait, rep1_nt, rep1_RN_LYSIS)))
      }
    }
  }
  head(common_rows_RN_LYSIS)
  
  png(paste("Reproducibility_of_",treatment,".png", sep=""),800, 600)
  plot(type.convert(common_rows_RN_LYSIS[,3]), type.convert(common_rows_RN_LYSIS[,4]), xlab=paste(treatment,"replicate_0",sep="_"), ylab=paste(treatment,"replicate_1",sep="_"), main=paste("Reproducibility of",treatment, sep=" "))
  abline(lm(type.convert(common_rows_RN_LYSIS[,4]) ~ type.convert(common_rows_RN_LYSIS[,3])), col="red")
  abline(v=RN_LYSIS_cutoff,h=RN_LYSIS_cutoff)
  ct <- cor.test(type.convert(common_rows_RN_LYSIS[,3]), type.convert(common_rows_RN_LYSIS[,4]))
  co <- ct$estimate
  text(-3,8, sprintf("r = %.3f",co), adj=c(1,1))
  text(-3,7, sprintf("cutoff = %.3f",RN_LYSIS_cutoff), adj=c(1,1))
  dev.off()
  
  png(paste("Correlation_between_NT_and",treatment,"repicate_0.png", sep="_"), 800, 600)
  plot(type.convert(common_rows_nt0_RN_LYSIS0[,3]), type.convert(common_rows_nt0_RN_LYSIS0[,4]), xlab="NT_replicate_0", ylab=paste(treatment,"replicate_0",sep="_"),main=paste("Correlation between NT and",treatment, "repicate_0", sep=" "))
  abline(lm(type.convert(common_rows_nt0_RN_LYSIS0[,4]) ~ type.convert(common_rows_nt0_RN_LYSIS0[,3])), col="red")
  abline(v=nt_cutoff,h=RN_LYSIS_cutoff)
  ct <- cor.test(type.convert(common_rows_nt0_RN_LYSIS0[,3]), type.convert(common_rows_nt0_RN_LYSIS0[,4]))
  co <- ct$estimate
  text(-3,8.0, sprintf("r = %.3f",co), adj=c(1,1))
  text(-3,7.5, sprintf("NT cutoff = %.3f",nt_cutoff), adj=c(1,1))
  text(-3,7.0, sprintf(paste(treatment,"cutoff = %.3f", sep=" "), RN_LYSIS_cutoff), adj=c(1,1))
  dev.off()
  
  png(paste("Correlation_between_NT_and",treatment,"repicate_1.png",sep="_"), 800, 600)
  plot(type.convert(common_rows_nt1_RN_LYSIS1[,3]), type.convert(common_rows_nt1_RN_LYSIS1[,4]), xlab="NT_replicate_1", ylab=paste(treatment,"replicate_1",sep="_"),main=paste("Correlation between NT and",treatment, "repicate_1", sep=" "))
  abline(lm(type.convert(common_rows_nt1_RN_LYSIS1[,4]) ~ type.convert(common_rows_nt1_RN_LYSIS1[,3])), col="red")
  abline(v=nt_cutoff,h=RN_LYSIS_cutoff)
  ct <- cor.test(type.convert(common_rows_nt1_RN_LYSIS1[,3]), type.convert(common_rows_nt1_RN_LYSIS1[,4]))
  co <- ct$estimate
  text(-3,8.0, sprintf("r = %.3f",co), adj=c(1,1))
  text(-3,7.5, sprintf("NT cutoff = %.3f",nt_cutoff), adj=c(1,1))
  text(-3,7.0, sprintf(paste(treatment,"cutoff = %.3f",sep=" "),RN_LYSIS_cutoff), adj=c(1,1))
  dev.off()
  ## define thresholded network
  
  RN_LYSIS_network_thresholded <- common_rows_RN_LYSIS[(as.double(common_rows_RN_LYSIS[, 3]) > RN_LYSIS_cutoff & as.double(common_rows_RN_LYSIS[, 4]) > RN_LYSIS_cutoff), ]

  #RN_LYSIS_network_thresholded[RN_LYSIS_network_thresholded[,2]=="YY1",]
  
  
  write.table(RN_LYSIS_network_thresholded, paste(treatment, "thresholded_network_replicates.tab",sep="_"), col.names=c("Prey", "Bait",  paste(treatment,">",sprintf("%.3f",RN_LYSIS_cutoff), sep=""), paste(treatment,">",sprintf("%.3f",RN_LYSIS_cutoff), sep="")), row.names=F, sep="\t")
  bait_int_count <- as.data.frame(apply(RN_LYSIS_network_thresholded[,1:2], 2, table)[[2]])
  prey_int_count <- as.data.frame(apply(RN_LYSIS_network_thresholded[,1:2], 2, table)[[1]])
  bait_int_count <- bait_int_count[order(bait_int_count$Freq, decreasing=T),]
  prey_int_count <- prey_int_count[order(prey_int_count$Freq, decreasing=T),]
  write.table(bait_int_count, paste(treatment, "thresholded_network_replicates_bait_int.tab",sep="_"), col.names=c("Bait", "interaction_count"), row.names=F, sep="\t")
  write.table(prey_int_count, paste(treatment, "thresholded_network_replicates_prey_int.tab",sep="_"), col.names=c("Prey", "interaction_count"), row.names=F, sep="\t")
  
  
  UNION_network <- NULL
  common_int <- 0
  
  for(b in common_baits){
    for(p in common_preys){
  
      NT_v <- min(common_rows_nt0[b,p],common_rows_nt1[b,p])
      RN_LYSIS_v <- min(common_rows_RN_LYSIS0[b,p], common_rows_RN_LYSIS1[b,p])
      NT_above <- F
      RN_LYSIS_above <- F
      
      for(x in 1:dim(NT_network_thresholded)[1]){
        if ((NT_network_thresholded[x,1] == p) && (NT_network_thresholded[x,2] == b)){
          #NT_v <- (type.convert(NT_network_thresholded[x, 3]) + type.convert(NT_network_thresholded[x, 4]))/2
          NT_above <- T
        }
      }
      for(y in 1:dim(RN_LYSIS_network_thresholded)[1]){
        if ((RN_LYSIS_network_thresholded[y,1] == p) && (RN_LYSIS_network_thresholded[y,2] == b)){
          #RN_LYSIS_v <- (type.convert(RN_LYSIS_network_thresholded[y, 3]) + type.convert(RN_LYSIS_network_thresholded[y, 4]))/2
          RN_LYSIS_above <- T
        }
      }
      
      if(NT_above || RN_LYSIS_above){
        
        if(NT_above){
          label <- "NT"
          if(RN_LYSIS_above){
            label <- "Both"
            common_int <- common_int + 1
          }
        }else{
          label <- treatment
        }
        
        line <- t(c(p, b, sprintf("%.3f",NT_v), sprintf("%.3f",RN_LYSIS_v), label))
        UNION_network <- rbind(UNION_network, line)
      }
    
    }
  }
  colnames(UNION_network) <- c("Prey", "Bait", paste("NT>",sprintf("%.1f",nt_cutoff),sep=""),  paste(treatment,">",sprintf("%.1f",RN_LYSIS_cutoff), sep=""), "Note")
  UNION_network_merge <- merge(UNION_network, DEseq2_table, by.x=c("Bait", "Prey"), by.y=c("baits", "preys"), all=T)
  write.table(UNION_network_merge, paste("NT",treatment,"joint_thresholded_network_replicates.tab",sep="_"), row.names=F, sep="\t")
  UNION_network_merge_DE <- UNION_network_merge[which(as.character(UNION_network_merge$Note) != "NA" & as.character(UNION_network_merge$baseMean) != "NA"), ]
  UNION_network_merge_DE <- UNION_network_merge_DE[order(UNION_network_merge_DE$Note),]
  write.table(UNION_network_merge_DE, paste("NT",treatment,"joint_thresholded_network_replicates_DE.tab",sep="_"), row.names=F, sep="\t")
  dim(NT_network_thresholded)
  dim(RN_LYSIS_network_thresholded)
  dim(UNION_network)
  common_int
  
  png(paste("NT",treatment,"correlation_final.png",sep="_"), 800, 600)
  plot(type.convert(UNION_network[,3]), type.convert(UNION_network[,4]), xlab="NT_zscore", ylab=paste(treatment,"zscore",sep="_"),main=paste("NT",treatment,"correlation final",sep="_"))
  abline(v=nt_cutoff,h=RN_LYSIS_cutoff)
  text(2, 8.5, sprintf("NT cutoff = %.3f",nt_cutoff), adj=c(1,1))
  text(2, 8.0, sprintf(paste(treatment,"cutoff = %.3f",sep=" "), RN_LYSIS_cutoff), adj=c(1,1))
  dev.off()
}
