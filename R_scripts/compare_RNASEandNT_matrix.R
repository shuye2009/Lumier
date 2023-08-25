library("inflection")
library("MASS")
library("outliers")
setwd("C:/RSYNC/LUMIER_C2H2/NOV232016")

RNASE <- read.table("RNASE_bait_prey_matrix_cluster.tab")
NT <- read.table("NT_bait_prey_matrix_cluster.tab")
head(RNASE)
head(NT)

rowNT <- row.names(NT)
rowRNASE <- row.names(RNASE)
rowdiff <- rowRNASE[!(rowRNASE %in% rowNT)]
rowdiff

DIFF <- NT - RNASE
head(DIFF)

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

write.table(DIFF, "RNASE-NT_bait_prey_matrix_difference.tab", col.names=NA, row.names=T, sep="\t")

baits <- rownames(NT)
preys <- colnames(NT)

RNASE_network <- NULL
NT_network <- NULL

dummy <- c("Blank", "Control1", "Control2", "Control3", "Control4", "Control5", "Control6", "Flag-GFP")

for(i in 1:length(baits)){
  bait <- baits[i]
  for(j in 1:length(preys)){
    prey <- preys[j]
    if(!is.element(bait, dummy)){
      rnase_lum <- RNASE[bait, prey]
      nt_lum <- NT[bait, prey]
      if(!is.nan(nt_lum)){
        nt_row <- t(c(bait, prey, nt_lum))
        NT_network <- rbind(NT_network, nt_row)
      }
      if(!is.nan(rnase_lum)){
        rnase_row <- t(c(bait, prey, rnase_lum))
        RNASE_network <- rbind(RNASE_network, rnase_row)
      }
    }
  }
}

rnase_values <- type.convert(RNASE_network[,3])
nt_values <- type.convert(NT_network[,3])



x <- nt_values
fit <- fitdistr(x, "exponential")
class(fit)
para <- fit$estimate
hist(x, prob = TRUE, main="Density plot of NT")
#curve(dnorm(x, para[1], para[2]), col = 3, add=T)
#nt_cutoff <- qnorm(0.95, para[1], para[2])
curve(dexp(x, para[1]), col = 3, add=T)
nt_cutoff <- qexp(0.95, para[1])
nt_cutoff

x <- rnase_values
fit <- fitdistr(x, "exponential")
class(fit)
para1 <- fit$estimate
hist(x, prob = TRUE, main="Density plot of RNASE")
#curve(dnorm(x, para1[1], para1[2]), col = 3, add=T)
#rnase_cutoff <- qnorm(0.95, para1[1], para1[2])
curve(dexp(x, para1[1]), col = 3, add=T)
rnase_cutoff <- qexp(0.95, para1[1])
rnase_cutoff

if(0){ #old way for picking cutoff which may not be valid

    quantile(nt_values, probs=c(0.88, 0.90,0.95))
    quantile(rnase_values, probs=c(0.88, 0.90,0.95))
    
    ntdist <- ecdf(nt_values) #empirical cdf
    rnasedist <- ecdf(rnase_values)
    
    
    un <- knots(ntdist)  # get unique values
    ur <- knots(rnasedist)
    
    
    inflection_nt <- edeci(un, ntdist(un), 0) # find the inflection point on cdf, edeci is the original function
    inflection_nt
    inflection_rnase <- edeci(ur, rnasedist(ur), 0)
    inflection_rnase
    
    
    
    ntdist(inflection_nt[,3]) # cdf p value for inflection point
    rnasedist(inflection_rnase[,3])
    
    pnorm(inflection_nt[,3]) # cdf p value for inflection point
    pnorm(inflection_rnase[,3])
    
    plot(ntdist, col.points = "blue")
    abline(v=inflection_nt[,3], h=ntdist(inflection_nt[,3]))
    plot(rnasedist, do.points=TRUE, col.points = "blue")
    abline(v=inflection_rnase[,3], h=rnasedist(inflection_rnase[,3]))
    
    br <- c(seq(-6, 5), by=0.5)
    hist(nt_values, prob=T)
    abline(v=inflection_nt[,3])
    hist(rnase_values, prob=T)
    abline(v=inflection_rnase[,3])
    
    nt_cutoff <- inflection_nt[,3]
    rnase_cutoff <- inflection_rnase[,3]
}


NT_network_thresholded <- NT_network[type.convert(NT_network[, 3]) > nt_cutoff, ]
RNASE_network_thresholded <- RNASE_network[type.convert(RNASE_network[, 3]) > rnase_cutoff, ]

NT_network_thresholded
RNASE_network_thresholded[RNASE_network_thresholded[,2]=="YY1",]

UNION_network <- NULL
common_int <- 0

for(i in 1:dim(NT_network)[1]){
  b <- NT_network[i,1]
  p <- NT_network[i,2]
  NT_v <- 0
  RNASE_v <- 0
  
  for(x in 1:dim(NT_network_thresholded)[1]){
    if ((NT_network_thresholded[x,1] == b) && (NT_network_thresholded[x,2] == p)){
      NT_v <- type.convert(NT_network_thresholded[x, 3])
    }
  }
  for(y in 1:dim(RNASE_network_thresholded)[1]){
    if ((RNASE_network_thresholded[y,1] == b) && (RNASE_network_thresholded[y,2] == p)){
      RNASE_v <- type.convert(RNASE_network_thresholded[y, 3])
    }
  }
  
  if(NT_v + RNASE_v > 0){
    label <- "NT"
    
    if(NT_v > 0){
      if(RNASE_v > 0){
        label <- "Both"
        common_int <- common_int + 1
      }
    }else{
      label <- "RNASE"
    }
    NT_v <- sprintf("%.3f",NT_v)
    RNASE_v <- sprintf("%.3f",RNASE_v)
    line <- t(c(b,p, NT_v, RNASE_v, label))
    UNION_network <- rbind(UNION_network, line)
  }
  

}

write.table(UNION_network, "NT_RNASE_joint_thresholded_network.tab", col.names=c("Bait", "Prey", paste("NT>",sprintf("%.3f",nt_cutoff),sep=""),  paste("RNASE>",sprintf("%.3f",rnase_cutoff), sep=""), "Note"), row.names=F, sep="\t")
dim(NT_network_thresholded)
dim(RNASE_network_thresholded)
dim(UNION_network)
common_int
