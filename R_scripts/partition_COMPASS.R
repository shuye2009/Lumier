
library("plyr")
library("ggplot2")
#library("lattice")
# this script implements Otsu's method of partition a set into two clusters by maximizing intercluster dissimilarity
# dis <- w1w2(mu1-m2)^2, where w1 and w2 are fractions of each cluster and mu1 and mu2 are means of each cluster

## not used, not working well

partition_scores <- function(aScore){
  aScore_min <- min(aScore)
  aScore_max <- floor(max(aScore))
  
  aScore_len <- length(aScore)
  
  cutoff <- 0
  max_dist <- 0
  max_w2 <- 0
  
  max_list <- list()
  stat_df <- NULL
  
  for(i in seq(aScore_min, aScore_max, 0.1)){
    set1 <- aScore[aScore <= i]
    set2 <- aScore[aScore > i]
    w1 <- length(set1)/aScore_len
    w2 <- length(set2)/aScore_len
    mu1 <- mean(set1)
    mu2 <- mean(set2)
    dist <- w1*w2*(mu2-mu1)*(mu2-mu1)
    if(dist > max_dist){
      max_dist <- dist
      cutoff <- i
      max_w2 <- length(set2)
      max_list <- list(set1, set2)
    }
    
    diff <- mu2 - mu1
    #print(paste(i, w1, w2, mu1, mu2, diff, dist, sep="   "))
    stat_df <- rbind(stat_df, c(i, w1, w2, mu1, mu2, diff, dist))
    
  }
  
  stat_df
  plot(stat_df[,1], stat_df[,7], xlab="Score", ylab="Inter-class distance")
  points(stat_df[,1], stat_df[,3], col="blue")
  abline(v=cutoff)
  print(paste("cutoff", cutoff, "max dist", max_dist, "fraction w2", max_w2, aScore_len, sep=" "))
  return(max_list)
}

wdir <- "C:/RSYNC/LUMIER_C2H2/Archive_oct032018/Archive1_3/Results_oct032018_control_distribution_normalization/HGscore_from_APMS_180913"

outdir <- file.path(wdir,"maxquant")
dir.create(outdir)
setwd(outdir)



NT_compass <- read.table("untreated_COMPASS_score.tab", header=T, sep="\t")
RN_compass <- read.table("Rnase_COMPASS_score.tab", header=T, sep="\t")

head(NT_compass)

print("NT")
aScore <- sort(as.numeric(NT_compass$wd))
summary(aScore)

plot(aScore)
seqx <- seq(1, length(aScore), 1)
score_spl <- smooth.spline(aScore, df=20)
lines(score_spl, col="red")

curvature <- predict(score_spl, x = seqx, deriv = 2)

plot(curvature, type = "l")
max_index <- curvature$x[which(curvature$y==max(curvature$y))]
cutoff <- aScore[max_index]

print("RN")
aScore <- sort(as.numeric(RN_compass$wd))
summary(aScore)

plot(aScore)
seqx <- seq(1, length(aScore), 1)
score_spl <- smooth.spline(aScore, df=20)
lines(score_spl, col="red")

curvature <- predict(score_spl, x = seqx, deriv = 2)

plot(curvature, type = "l")
max_index <- curvature$x[which(curvature$y==max(curvature$y))]
aScore[max_index]#RN_partition <- partition_scores(aScore)

