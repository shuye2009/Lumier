
library("MASS")
setwd("C:/RSYNC/LUMIER_C2H2/secondProtocol")

RN_LYSIS <- read.table("RN_LYSIS_bait_prey_log_matrix_cluster.tab")
NT <- read.table("NT_bait_prey_log_matrix_cluster.tab")
head(RN_LYSIS)
head(NT)

rowNT <- row.names(NT)
rowRN_LYSIS <- row.names(RN_LYSIS)
rowdiff <- rowRN_LYSIS[!(rowRN_LYSIS %in% rowNT)]
rowdiff

DIFF <- RN_LYSIS - NT
head(DIFF)

if(0){ #add row and col mean to the margins
  
  DIFF_sq = apply(DIFF, c(1,2), function(x) x^2)
  
  bait_diff <- apply(DIFF_sq, 1, mean)
  bait_diff_m <- mean(bait_diff)
  
  DIFF <- cbind(DIFF, bait_diff)
  prey_diff <- apply(DIFF_sq, 2, mean)
  
  bait_diff <- bait_diff_m
  prey_diff <- t(c(prey_diff, bait_diff))
  row.names(prey_diff) <- "prey_diff"
  colnames(prey_diff) <- colnames(DIFF)
  DIFF <- rbind(DIFF, prey_diff)
}

write.table(DIFF, "RN_LYSIS-NT_bait_prey_matrix_difference.tab", col.names=NA, row.names=T, sep="\t")

baits <- rownames(NT)
preys <- colnames(NT)

RN_LYSIS_network <- NULL
NT_network <- NULL

dummy <- c("Blank", "Control1", "Control2", "Control3", "Control4", "Control5", "Control6", "Flag-GFP", "H2O", "H2O-GFP")

for(i in 1:length(baits)){
  bait <- baits[i]
  for(j in 1:length(preys)){
    prey <- preys[j]
    if(!is.element(bait, dummy)){
      RN_LYSIS_lum <- RN_LYSIS[bait, prey]
      nt_lum <- NT[bait, prey]
      if(!is.nan(nt_lum)){
        nt_row <- t(c(bait, prey, nt_lum))
        NT_network <- rbind(NT_network, nt_row)
      }
      if(!is.nan(RN_LYSIS_lum)){
        RN_LYSIS_row <- t(c(bait, prey, RN_LYSIS_lum))
        RN_LYSIS_network <- rbind(RN_LYSIS_network, RN_LYSIS_row)
      }
    }
  }
}

RN_LYSIS_values <- type.convert(RN_LYSIS_network[,3])
nt_values <- type.convert(NT_network[,3])



x <- nt_values
fit <- fitdistr(x, "normal")
class(fit)
para <- fit$estimate

png("Density_plot_of_NT.png", 800, 600)
hist(x, prob = TRUE, main="Density plot of NT")
curve(dnorm(x, para[1], para[2]), col = 3, add=T)
dev.off();
nt_cutoff <- qnorm(0.95, para[1], para[2])
#curve(dexp(x, para[1]), col = 3, add=T)
#nt_cutoff <- qexp(0.95, para[1])
nt_cutoff

x <- RN_LYSIS_values
fit <- fitdistr(x, "normal")
class(fit)
para1 <- fit$estimate

png("Density_plot_of_RN_LYSIS.png", 800, 600)
hist(x, prob = TRUE, main="Density plot of RN_LYSIS")
curve(dnorm(x, para1[1], para1[2]), col = 3, add=T)
dev.off()
RN_LYSIS_cutoff <- qnorm(0.95, para1[1], para1[2])
#curve(dexp(x, para1[1]), col = 3, add=T)
#RN_LYSIS_cutoff <- qexp(0.95, para1[1])
RN_LYSIS_cutoff



NT_network_thresholded <- NT_network[type.convert(NT_network[, 3]) > nt_cutoff, ]
RN_LYSIS_network_thresholded <- RN_LYSIS_network[type.convert(RN_LYSIS_network[, 3]) > RN_LYSIS_cutoff, ]

NT_network_thresholded
RN_LYSIS_network_thresholded[RN_LYSIS_network_thresholded[,2]=="YY1",]

UNION_network <- NULL
common_int <- 0

for(i in 1:dim(NT_network)[1]){
  b <- NT_network[i,1]
  p <- NT_network[i,2]
  NT_v <- 0
  RN_LYSIS_v <- 0
  
  for(x in 1:dim(NT_network_thresholded)[1]){
    if ((NT_network_thresholded[x,1] == b) && (NT_network_thresholded[x,2] == p)){
      NT_v <- type.convert(NT_network_thresholded[x, 3])
    }
  }
  for(y in 1:dim(RN_LYSIS_network_thresholded)[1]){
    if ((RN_LYSIS_network_thresholded[y,1] == b) && (RN_LYSIS_network_thresholded[y,2] == p)){
      RN_LYSIS_v <- type.convert(RN_LYSIS_network_thresholded[y, 3])
    }
  }
  
  if(NT_v + RN_LYSIS_v > 0){
    label <- "NT"
    
    if(NT_v > 0){
      if(RN_LYSIS_v > 0){
        label <- "Both"
        common_int <- common_int + 1
      }
    }else{
      label <- "RN_LYSIS"
    }
    NT_v <- sprintf("%.3f",NT_v)
    RN_LYSIS_v <- sprintf("%.3f",RN_LYSIS_v)
    line <- t(c(b,p, NT_v, RN_LYSIS_v, label))
    UNION_network <- rbind(UNION_network, line)
  }
  

}

write.table(UNION_network, "NT_RN_LYSIS_joint_thresholded_network.tab", col.names=c("Bait", "Prey", paste("NT>",sprintf("%.3f",nt_cutoff),sep=""),  paste("RN_LYSIS>",sprintf("%.3f",RN_LYSIS_cutoff), sep=""), "Note"), row.names=F, sep="\t")
dim(NT_network_thresholded)
dim(RN_LYSIS_network_thresholded)
dim(UNION_network)
common_int
