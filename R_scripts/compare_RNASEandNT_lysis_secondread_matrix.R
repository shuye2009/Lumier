
library("MASS")
wdir <- "C:/RSYNC/LUMIER_C2H2/SecondProtocol/Results_mar082018"
outdir <- file.path(wdir,"matrix")
dir.create(outdir)
setwd(outdir)
#SNDREAD and RN_LYSIS are interchangeable
for(treatment in c("RN_LYSIS", "SNDREAD")){
  RN_LYSIS <- read.table(file.path(wdir,paste(treatment,"bait_prey_log_matrix_cluster.tab", sep="_")))
  NT <- read.table(file.path(wdir,"NT_bait_prey_log_matrix_cluster.tab"))
  head(RN_LYSIS)
  head(NT)
  
  rowNT <- row.names(NT)
  rowRN_LYSIS <- row.names(RN_LYSIS)
  rowdiff <- rowRN_LYSIS[!(rowRN_LYSIS %in% rowNT)]
  rowdiff
  
  colNT <- colnames(NT)
  colRN_LYSIS <- colnames(RN_LYSIS)
  coldiff <- colNT[!(colNT %in% colRN_LYSIS)]
  coldiff
  
  NT <- NT[,colRN_LYSIS] #size of RN_LYSIS and NT are different
  
  DIFF <- RN_LYSIS - NT  
  #head(DIFF)
  
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
  
  write.table(DIFF, paste(treatment, "NT_bait_prey_matrix_difference.tab", sep="_"), col.names=NA, row.names=T, sep="\t")
  
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
  
  png("Density_plot_of_NT_matrix.png", 800, 600)
  hist(x, prob = TRUE, ylim=c(0,0.6), main="Density plot of NT")
  curve(dnorm(x, para[1], para[2]), col = 3, add=T)
  dev.off();
  
  png("QQ_plot_NT_matrix.png", 800, 600)
  qqnorm(x, main="Normal Q-Q Plot: NT matrix")
  qqline(x)
  dev.off() 
  
  nt_cutoff <- qnorm(0.95, para[1], para[2])
  #curve(dexp(x, para[1]), col = 3, add=T)
  #nt_cutoff <- qexp(0.95, para[1])
  nt_cutoff
  
  x <- RN_LYSIS_values
  fit <- fitdistr(x, "normal")
  class(fit)
  para1 <- fit$estimate
  
  png(paste("Density_plot_of",treatment,"matrix.png",sep="_"), 800, 600)
  hist(x, prob = TRUE, ylim=c(0,0.6), main=paste("Density plot of", treatment, sep=" "))
  curve(dnorm(x, para1[1], para1[2]), col = 3, add=T)
  dev.off()
  
  png(paste("QQ_plot",treatment,"matrix.png",sep="_"), 800, 600)
  qqnorm(x,main=paste("Normal Q-Q Plot:",treatment,"matrix",sep=" "))
  qqline(x)
  dev.off() 
  
  RN_LYSIS_cutoff <- qnorm(0.95, para1[1], para1[2])
  #curve(dexp(x, para1[1]), col = 3, add=T)
  #RN_LYSIS_cutoff <- qexp(0.95, para1[1])
  RN_LYSIS_cutoff
  
  
  
  NT_network_thresholded <- NT_network[type.convert(NT_network[, 3]) > nt_cutoff, ]
  RN_LYSIS_network_thresholded <- RN_LYSIS_network[type.convert(RN_LYSIS_network[, 3]) > RN_LYSIS_cutoff, ]
  
  NT_network_thresholded
  RN_LYSIS_network_thresholded[RN_LYSIS_network_thresholded[,2]=="YY1",]
  
  write.table(NT_network_thresholded, paste("NT", "thresholded_network_matrix.tab",sep="_"), col.names=c("Bait", "Prey", paste("NT>",sprintf("%.3f",nt_cutoff),sep="")), row.names=F, sep="\t")
  write.table(RN_LYSIS_network_thresholded, paste(treatment, "thresholded_network_matrix.tab",sep="_"), col.names=c("Bait", "Prey",  paste(treatment,">",sprintf("%.3f",RN_LYSIS_cutoff), sep="")), row.names=F, sep="\t")
  
  
  UNION_network <- NULL
  common_int <- 0
  
  for(i in 1:dim(NT_network)[1]){
    b <- NT_network[i,1]
    p <- NT_network[i,2]
    
    NT_v <- type.convert(NT_network[i,3])
    RN_LYSIS_v <- 0
    subset1 <- RN_LYSIS_network[RN_LYSIS_network[,1]==b & RN_LYSIS_network[,2]==p,3]
    
    if(length(subset1) > 0){
      RN_LYSIS_v <- type.convert(subset1)
    }
    
    NT_above <- F
    RN_LYSIS_above <- F
    
    for(x in 1:dim(NT_network_thresholded)[1]){
      if ((NT_network_thresholded[x,1] == b) && (NT_network_thresholded[x,2] == p)){
        #NT_v <- type.convert(NT_network_thresholded[x, 3])
        NT_above <- T
      }
    }
    for(y in 1:dim(RN_LYSIS_network_thresholded)[1]){
      if ((RN_LYSIS_network_thresholded[y,1] == b) && (RN_LYSIS_network_thresholded[y,2] == p)){
        #RN_LYSIS_v <- type.convert(RN_LYSIS_network_thresholded[y, 3])
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
      NT_v <- sprintf("%.3f",NT_v)
      RN_LYSIS_v <- sprintf("%.3f",RN_LYSIS_v)
      line <- t(c(b,p, NT_v, RN_LYSIS_v, label))
      UNION_network <- rbind(UNION_network, line)
    }
    
  
  }
  
  write.table(UNION_network, paste("NT", treatment, "joint_thresholded_network_matrix.tab",sep="_"), col.names=c("Bait", "Prey", paste("NT>",sprintf("%.3f",nt_cutoff),sep=""),  paste(treatment,">",sprintf("%.3f",RN_LYSIS_cutoff), sep=""), "Note"), row.names=F, sep="\t")
  dim(NT_network_thresholded)
  dim(RN_LYSIS_network_thresholded)
  dim(UNION_network)
  common_int
}
