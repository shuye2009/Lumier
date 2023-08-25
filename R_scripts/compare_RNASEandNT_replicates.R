
library("MASS")
library("outliers")
setwd("C:/RSYNC/LUMIER_C2H2/NOV232016")

RNASE <- read.table("RNASE_replicates_only_log_normalized_lum.tab", header=T)
NT <- read.table("NT_replicates_only_log_normalized_lum.tab", header=T)
head(RNASE)
head(NT)

dummy <- c("Blank", "Control1", "Control2", "Control3", "Control4", "Control5", "Control6", "Flag-GFP", "H2O")

NT_network <- NT[!is.element(NT[,3], dummy),]
RNASE_network <- RNASE[!is.element(RNASE[,3], dummy),]

#define NT data cutoff
nt_values <- NT_network[,4]

x <- nt_values
fit <- fitdistr(x, "normal")
class(fit)
para <- fit$estimate
hist(x, prob = TRUE, main="Density plot of NT")
curve(dnorm(x, para[1], para[2]), col = 3, add=T)
qqnorm(x)
qqline(x)
nt_cutoff <- qnorm(0.95, para[1], para[2])
#curve(dexp(x, para[1]), col = 3, add=T)
#nt_cutoff <- qexp(0.95, para[1])
nt_cutoff


      # display NT reproducibility
NT_network_rep0 <- NT_network[NT_network[,2]=="rep0",]
NT_network_rep1 <- NT_network[NT_network[,2]=="rep1",]

dim(NT_network_rep0)
dim(NT_network_rep1)

rep0_baits <- as.character(unique(NT_network_rep0[,3]))
rep1_baits <- as.character(unique(NT_network_rep1[,3]))
rep0_preys <- as.character(unique(NT_network_rep0[,1]))
rep1_preys <- as.character(unique(NT_network_rep1[,1]))

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

for(i in 1:length(rep0_baits)){
  bait <- rep0_baits[i]
  for(j in 1:length(rep0_preys)){
    prey <- rep0_preys[j]
    rep0 <- common_rows_nt0[bait,prey]
    rep1 <- common_rows_nt1[bait,prey]
    if(!is.na(rep0) && !is.na(rep1)){
      common_rows_nt <- rbind(common_rows_nt, t(c(prey, bait, rep0, rep1)))
    }
  }
}
head(common_rows_nt)


plot(type.convert(common_rows_nt[,3]), type.convert(common_rows_nt[,4]), xlab="NT_replicate_0", ylab="NT_replicate_1",main="Reproducibility of NT")
abline(lm(type.convert(common_rows_nt[,4]) ~ type.convert(common_rows_nt[,3])), col="red")
abline(v=nt_cutoff,h=nt_cutoff)
ct <- cor.test(type.convert(common_rows_nt[,3]), type.convert(common_rows_nt[,4]))
co <- ct$estimate
text(-4,4, sprintf("r = %.3f",co), adj=c(1,1))
text(-4,3.5, sprintf("cutoff = %.3f",nt_cutoff), adj=c(1,1))

#define RNASE data cutoff

rnase_values <- RNASE_network[,4]

x <- rnase_values
fit <- fitdistr(x, "normal")
class(fit)
para1 <- fit$estimate
hist(x, prob = TRUE, main="Density plot of RNASE")
curve(dnorm(x, para1[1], para1[2]), col = 3, add=T)
qqnorm(x)
qqline(x)
rnase_cutoff <- qnorm(0.95, para1[1], para1[2])
#curve(dexp(x, para1[1]), col = 3, add=T)
#rnase_cutoff <- qexp(0.95, para1[1])
rnase_cutoff
      # display RNASE reproducibility
RNASE_network_rep0 <- RNASE_network[RNASE_network[,2]=="rep0",]
RNASE_network_rep1 <- RNASE_network[RNASE_network[,2]=="rep1",]

dim(RNASE_network_rep0)
dim(RNASE_network_rep1)

rep0_baits <- as.character(unique(RNASE_network_rep0[,3]))
rep1_baits <- as.character(unique(RNASE_network_rep1[,3]))
rep0_preys <- as.character(unique(RNASE_network_rep0[,1]))
rep1_preys <- as.character(unique(RNASE_network_rep1[,1]))

common_rows_rnase0 <- matrix(NA,ncol=length(rep0_preys),nrow=length(rep0_baits), dimnames=list(rep0_baits,rep0_preys))
common_rows_rnase1 <- matrix(NA,ncol=length(rep1_preys),nrow=length(rep1_baits), dimnames=list(rep1_baits,rep1_preys))
common_rows_rnase <- NULL
for(i in 1:dim(RNASE_network_rep0)[1]){
  bait <- as.character(RNASE_network_rep0[i,3])
  prey <- as.character(RNASE_network_rep0[i,1])
  common_rows_rnase0[bait,prey] <- RNASE_network_rep0[i,4]
  
}

for(i in 1:dim(RNASE_network_rep1)[1]){
  bait <- as.character(RNASE_network_rep1[i,3])
  prey <- as.character(RNASE_network_rep1[i,1])
  common_rows_rnase1[bait,prey] <- RNASE_network_rep1[i,4]
  
}

common_rows_nt0_rnase0 <- NULL
common_rows_nt1_rnase1 <- NULL
for(i in 1:length(rep0_baits)){
  bait <- rep0_baits[i]
  for(j in 1:length(rep0_preys)){
    prey <- rep0_preys[j]
    rep0_rnase <- common_rows_rnase0[bait,prey]
    rep1_rnase <- common_rows_rnase1[bait,prey]
    if(!is.na(rep0_rnase) && !is.na(rep1_rnase)){
      common_rows_rnase <- rbind(common_rows_rnase, t(c(prey, bait, rep0_rnase, rep1_rnase)))
    }
    rep0_nt <- common_rows_nt0[bait,prey]
    rep1_nt <- common_rows_nt1[bait,prey]
    if(!is.na(rep0_nt) && !is.na(rep0_rnase)){
      common_rows_nt0_rnase0 <- rbind(common_rows_nt0_rnase0, t(c(prey, bait, rep0_nt, rep0_rnase)))
    }
    if(!is.na(rep1_nt) && !is.na(rep1_rnase)){
      common_rows_nt1_rnase1 <- rbind(common_rows_nt1_rnase1, t(c(prey, bait, rep1_nt, rep1_rnase)))
    }
  }
}
head(common_rows_rnase)


plot(type.convert(common_rows_rnase[,3]), type.convert(common_rows_rnase[,4]), xlab="RNASE_replicate_0", ylab="RNASE_replicate_1",main="Reproducibility of RNASE")
abline(lm(type.convert(common_rows_rnase[,4]) ~ type.convert(common_rows_rnase[,3])), col="red")
abline(v=rnase_cutoff,h=rnase_cutoff)
ct <- cor.test(type.convert(common_rows_rnase[,3]), type.convert(common_rows_rnase[,4]))
co <- ct$estimate
text(-4,4, sprintf("r = %.3f",co), adj=c(1,1))
text(-4,3.5, sprintf("cutoff = %.3f",rnase_cutoff), adj=c(1,1))

plot(type.convert(common_rows_nt0_rnase0[,3]), type.convert(common_rows_nt0_rnase0[,4]), xlab="NT_replicate_0", ylab="RNASE_replicate_0",main="Correlation between NT and RNASE repicate_0")
abline(lm(type.convert(common_rows_nt0_rnase0[,4]) ~ type.convert(common_rows_nt0_rnase0[,3])), col="red")
abline(v=nt_cutoff,h=rnase_cutoff)
ct <- cor.test(type.convert(common_rows_nt0_rnase0[,3]), type.convert(common_rows_nt0_rnase0[,4]))
co <- ct$estimate
text(-4,4, sprintf("r = %.3f",co), adj=c(1,1))
text(-4,3.5, sprintf("NT cutoff = %.3f",nt_cutoff), adj=c(1,1))
text(-4,3.0, sprintf("RNASE cutoff = %.3f",rnase_cutoff), adj=c(1,1))

plot(type.convert(common_rows_nt1_rnase1[,3]), type.convert(common_rows_nt1_rnase1[,4]), xlab="NT_replicate_1", ylab="RNASE_replicate_1",main="Correlation between NT and RNASE repicate_1")
abline(lm(type.convert(common_rows_nt1_rnase1[,4]) ~ type.convert(common_rows_nt1_rnase1[,3])), col="red")
abline(v=nt_cutoff,h=rnase_cutoff)
ct <- cor.test(type.convert(common_rows_nt1_rnase1[,3]), type.convert(common_rows_nt1_rnase1[,4]))
co <- ct$estimate
text(-3,4, sprintf("r = %.3f",co), adj=c(1,1))
text(-3,3.5, sprintf("NT cutoff = %.3f",nt_cutoff), adj=c(1,1))
text(-3,3.0, sprintf("RNASE cutoff = %.3f",rnase_cutoff), adj=c(1,1))

## define thresholded network

NT_network_thresholded <- common_rows_nt[(as.double(common_rows_nt[, 3]) > nt_cutoff & as.double(common_rows_nt[, 4]) > nt_cutoff), ]
RNASE_network_thresholded <- common_rows_rnase[(as.double(common_rows_rnase[, 3]) > nt_cutoff & as.double(common_rows_rnase[, 4]) > nt_cutoff), ]

NT_network_thresholded
RNASE_network_thresholded[RNASE_network_thresholded[,2]=="YY1",]

UNION_network <- NULL
common_int <- 0

for(b in rep0_baits){
  for(p in rep0_preys){

    NT_v <- (common_rows_nt0[b,p] + common_rows_nt1[b,p])/2
    RNASE_v <- (common_rows_rnase0[b,p] + common_rows_rnase1[b,p])/2
    NT_above <- F
    RNASE_above <- F
    
    for(x in 1:dim(NT_network_thresholded)[1]){
      if ((NT_network_thresholded[x,1] == p) && (NT_network_thresholded[x,2] == b)){
        #NT_v <- (type.convert(NT_network_thresholded[x, 3]) + type.convert(NT_network_thresholded[x, 4]))/2
        NT_above <- T
      }
    }
    for(y in 1:dim(RNASE_network_thresholded)[1]){
      if ((RNASE_network_thresholded[y,1] == p) && (RNASE_network_thresholded[y,2] == b)){
        #RNASE_v <- (type.convert(RNASE_network_thresholded[y, 3]) + type.convert(RNASE_network_thresholded[y, 4]))/2
        RNASE_above <- T
      }
    }
    
    if(NT_above || RNASE_above){
      
      if(NT_above){
        label <- "NT"
        if(RNASE_above){
          label <- "Both"
          common_int <- common_int + 1
        }
      }else{
        label <- "RNASE"
      }
      
      line <- t(c(p, b, sprintf("%.3f",NT_v), sprintf("%.3f",RNASE_v), label))
      UNION_network <- rbind(UNION_network, line)
    }
  
  }
}

write.table(UNION_network, "NT_RNASE_joint_thresholded_network_replicates.tab", col.names=c("Prey", "Bait", paste("NT>",sprintf("%.3f",nt_cutoff),sep=""),  paste("RNASE>",sprintf("%.3f",rnase_cutoff), sep=""), "Note"), row.names=F, sep="\t")
dim(NT_network_thresholded)
dim(RNASE_network_thresholded)
dim(UNION_network)
common_int

plot(type.convert(UNION_network[,3]), type.convert(UNION_network[,4]), xlab="NT_zscore", ylab="RNASE_zscore",main="NT_RNASE correlation final")
abline(v=nt_cutoff,h=rnase_cutoff)
text(0,-0.2, sprintf("NT cutoff = %.3f",nt_cutoff), adj=c(1,1))
text(0,-0.4, sprintf("RNASE cutoff = %.3f",rnase_cutoff), adj=c(1,1))
