
library("plyr")
library("ggplot2")
#library("lattice")

spline_cutoff <- function(aScore, df){
  aScore <- sort(aScore)

  seqx <- seq(1, length(aScore), 1)
  score_spl <- smooth.spline(seqx, aScore, df=df)
  
  
  curvature <- predict(score_spl, x = seqx, deriv = 2)
  max_curvature_cutoff <- aScore[which(curvature$y==max(curvature$y))]
  plot(curvature, type = "l")
  plot(aScore)
  lines(score_spl, col="red")
  abline(v=curvature$x[which(curvature$y==max(curvature$y))], lty=3)
  
  return(max_curvature_cutoff)
}

wdir <- "C:/RSYNC/LUMIER_C2H2/Archive_oct032018/Archive1_3/Results_oct032018_control_distribution_normalization/HGscore_from_APMS_180913"

outdir <- file.path(wdir,"prohits")
dir.create(outdir)
setwd(outdir)


prohits <- read.delim(file.path(outdir,"7_comparison_matrix_prohits.tab"), header=T, sep="\t")

head(prohits)

prohits_matrix <- prohits[, 2:dim(prohits)[2]]
row.names(prohits_matrix) <- as.character(prohits[,1])
prohits_matrix[is.na(prohits_matrix)] <- 0
head(prohits_matrix)

prohits_matrix <- as.matrix(prohits_matrix)
pdf("Total_peptide_count_Prohits_heatmap.pdf", width=10, height=10)
heatmap(prohits_matrix, margins=c(10,10))
dev.off()

### remove bad samples

cols <- colnames(prohits_matrix)
cols
bad <- c("X12119.ZNF99_1.N.EGFP..7652", "X12130.ZNF148_RNAse_2.N.EGFP..7707", "X12129.ZNF148_RNAse_1.N.EGFP..7707")
cols <- cols[!cols %in% bad]
prohits_matrix <- prohits_matrix[,cols]
preys <- row.names(prohits_matrix)

## format for compass scoring
Rnase <- NULL
untreated <- NULL
Rc <- 0
Uc <- 0
for(i in 1:length(colnames(prohits_matrix))){
  #i <- 5
  purifs <- colnames(prohits_matrix)[i]
  purif <- unlist(strsplit(purifs, ".", fixed=T))[2]
  bait <- unlist(strsplit(purif, "_", fixed=T))[1]
  for(j in 1:length(preys)){
    prey <- as.character(preys[j])
    spc <- as.numeric(prohits_matrix[j,i])
    if(!is.na(spc)){
      arow <- c(purif, bait, prey, spc)
      if(grepl("RNAse", purif, fixed=T)){
        Rc <- Rc + 1
        Rnase <- rbind(Rnase, c(purif, bait, prey, spc))
      }else{
        Uc <- Uc + 1
        untreated <- rbind(untreated, c(purif, bait, prey, spc))
      }
    }
  }
}
Rc
Uc
colnames(Rnase) <- c("Purification", "Bait", "Prey", "Spectral counts")
colnames(untreated) <- c("Purification", "Bait", "Prey", "Spectral counts")

write.table(Rnase, "Rnase_prohits_for_COMPASS.txt", row.names=F, col.names=T, quote=F, sep="\t")
write.table(untreated, "untreated_prohits_for_COMPASS.txt", row.names=F, col.names=T, quote=F, sep="\t")

system("perl C:/RSYNC/worksapce2015/MASS_PEC_SCORE/compass.pl Rnase_prohits_for_COMPASS.txt Rnase_COMPASS_score.tab 0 0", show.output.on.console=T)
system("perl C:/RSYNC/worksapce2015/MASS_PEC_SCORE/compass.pl untreated_prohits_for_COMPASS.txt untreated_COMPASS_score.tab 0 0", show.output.on.console=T)


### compare Rnase against untreated compass score
NT_compass <- read.table("untreated_COMPASS_score.tab", header=T, sep="\t")
Rnase_compass <- read.table("Rnase_COMPASS_score.tab", header=T, sep="\t")

head(NT_compass)

df <- 40
compass_NT_cutoff <- spline_cutoff(as.numeric(NT_compass$wd), df)
compass_RN_cutoff <- spline_cutoff(as.numeric(Rnase_compass$wd), df)

merge_compass <- merge(NT_compass, Rnase_compass, by=c("bait", "prey"), all=F)
merge_compass <- merge_compass[,c("bait", "prey", "wd.x", "wd.y")]
colnames(merge_compass) <- c("bait", "prey", "untreated", "Rnase")
png("Rnase_untreated_COMPASS_correlation_prohits.png")
plot(log2(merge_compass$untreated), log2(merge_compass$Rnase))
abline(h=log2(compass_RN_cutoff), v=log2(compass_NT_cutoff))
dev.off()

selected <- merge_compass[(merge_compass$untreated > compass_NT_cutoff | merge_compass$Rnase > compass_RN_cutoff),]
mutate(selected, note="BOTH")
selected[selected$untreated > compass_NT_cutoff & selected$Rnase < compass_RN_cutoff, "note"] <- "NT"
selected[selected$untreated < compass_NT_cutoff & selected$Rnase > compass_RN_cutoff, "note"] <- "RNASE"
selected[selected$untreated > compass_NT_cutoff & selected$Rnase > compass_RN_cutoff, "note"] <- ""
selected <- selected[order(selected$bait, selected$note),]
colnames(selected) <- c("bait", "prey", paste("untreated", sprintf("%.2f",compass_NT_cutoff)), paste("Rnase", sprintf("%.2f", compass_RN_cutoff)), "Note")
write.table(selected, "Rnase_untreated_COMPASS_difference_prohits.tab", row.names=F, col.names=T, quote=F, sep="\t")


colramp = colorRampPalette(c('white', 'blue', 'green', 'yellow', 'red'))
smoothScatter(log(merge_compass$untreated), log(merge_compass$Rnase), colramp = colramp, nbin=c(500,500))

x1 <- log2(merge_compass$untreated) + log2(merge_compass$Rnase)
x2 <- log2(merge_compass$untreated) - log2(merge_compass$Rnase)

plot(x1, x2)
abline(h=c(-1, 1), v=6)

### ### check correlation betwee replicates of maxquant peptide counts

Rnase_spc_maxquant <- read.table("../maxquant/Rnase_for_HGscore.tab", header=T, sep="\t")
untreated_spc_maxquant <- read.table("../maxquant/untreated_for_HGscore.tab", header=T, sep="\t")
Rnase_spc_maxquant_purif <- as.character(unique(sort(Rnase_spc_maxquant$Purification)))
purif_index <- unlist(lapply(Rnase_spc_maxquant_purif, function(x){unlist(strsplit(x, "_", fixed=T))[12]}))

Rnase_spc_maxquant_purif_rep1 <- Rnase_spc_maxquant_purif[as.numeric(purif_index) %% 2 == 1]
Rnase_spc_maxquant_purif_rep2 <- Rnase_spc_maxquant_purif[as.numeric(purif_index) %% 2 == 0]
Rnase_spc_maxquant_rep1 <- Rnase_spc_maxquant[Rnase_spc_maxquant$Purification %in% Rnase_spc_maxquant_purif_rep1, ]
Rnase_spc_maxquant_rep2 <- Rnase_spc_maxquant[Rnase_spc_maxquant$Purification %in% Rnase_spc_maxquant_purif_rep2, ]

untreated_spc_maxquant_purif <- as.character(unique(sort(untreated_spc_maxquant$Purification)))
purif_index <- unlist(lapply(untreated_spc_maxquant_purif, function(x){unlist(strsplit(x, "_", fixed=T))[12]}))

untreated_spc_maxquant_purif_rep1 <- untreated_spc_maxquant_purif[as.numeric(purif_index) %% 2 == 1]
untreated_spc_maxquant_purif_rep2 <- untreated_spc_maxquant_purif[as.numeric(purif_index) %% 2 == 0]
untreated_spc_maxquant_rep1 <- untreated_spc_maxquant[untreated_spc_maxquant$Purification %in% untreated_spc_maxquant_purif_rep1, ]
untreated_spc_maxquant_rep2 <- untreated_spc_maxquant[untreated_spc_maxquant$Purification %in% untreated_spc_maxquant_purif_rep2, ]

untreated_mq_rep_merge <- merge(untreated_spc_maxquant_rep1, untreated_spc_maxquant_rep2, by=c("Bait", "Prey"), all=F)
untreated_rep1 <- as.numeric(untreated_mq_rep_merge[,4])
untreated_rep2 <- as.numeric(untreated_mq_rep_merge[,6])
png("Maxquant_replicates_untreated.png")
plot(untreated_rep1, untreated_rep2, col="blue", main="Maxquant")
dev.off()
Rnase_mq_rep_merge <- merge(Rnase_spc_maxquant_rep1, Rnase_spc_maxquant_rep2, by=c("Bait", "Prey"), all=F)
Rnase_rep1 <- as.numeric(Rnase_mq_rep_merge[,4])
Rnase_rep2 <- as.numeric(Rnase_mq_rep_merge[,6])
png("Maxquant_replicates_Rnase.png")
plot(Rnase_rep1, Rnase_rep2, col="blue", main="Maxquant")
dev.off()

Rnase_mq_rep_merge[Rnase_mq_rep_merge[,4] < 2 & Rnase_mq_rep_merge[,6] > 10, c(1,2,4,6)]
untreated_mq_rep_merge[untreated_mq_rep_merge[,4] > 10 & untreated_mq_rep_merge[,6] < 2, c(1,2,4,6)]

### check correlation betwee replicates of prohits peptide counts


Rnase_spc_prohits <- as.data.frame(Rnase)
untreated_spc_prohits <- as.data.frame(untreated)

Rnase_spc_prohits_purif <- unique(sort(Rnase_spc_prohits$Purification))
Rnase_spc_prohits_purif_rep1 <- Rnase_spc_prohits_purif[seq(1, length(Rnase_spc_prohits_purif), 2)]
Rnase_spc_prohits_purif_rep2 <- Rnase_spc_prohits_purif[seq(2, length(Rnase_spc_prohits_purif), 2)]
Rnase_spc_prohits_rep1 <- Rnase_spc_prohits[Rnase_spc_prohits$Purification %in% Rnase_spc_prohits_purif_rep1, ]
Rnase_spc_prohits_rep2 <- Rnase_spc_prohits[Rnase_spc_prohits$Purification %in% Rnase_spc_prohits_purif_rep2, ]

untreated_spc_prohits_purif <- unique(sort(untreated_spc_prohits$Purification))
untreated_spc_prohits_purif_rep1 <- untreated_spc_prohits_purif[seq(1, length(untreated_spc_prohits_purif), 2)]
untreated_spc_prohits_purif_rep2 <- untreated_spc_prohits_purif[seq(2, length(untreated_spc_prohits_purif), 2)]
untreated_spc_prohits_rep1 <- untreated_spc_prohits[untreated_spc_prohits$Purification %in% untreated_spc_prohits_purif_rep1, ]
untreated_spc_prohits_rep2 <- untreated_spc_prohits[untreated_spc_prohits$Purification %in% untreated_spc_prohits_purif_rep2, ]

untreated_prohits_rep_merge <- merge(untreated_spc_prohits_rep1, untreated_spc_prohits_rep2, by=c("Bait", "Prey"), all=F)
untreated_rep1 <- as.numeric(untreated_prohits_rep_merge[,4])
untreated_rep2 <- as.numeric(untreated_prohits_rep_merge[,6])
png("Prohits_replicates_untreated.png")
plot(untreated_rep1, untreated_rep2, col="blue", main="Prohits")
dev.off()

Rnase_prohits_rep_merge <- merge(Rnase_spc_prohits_rep1, Rnase_spc_prohits_rep2, by=c("Bait", "Prey"), all=F)
Rnase_rep1 <- as.numeric(Rnase_prohits_rep_merge[,4])
Rnase_rep2 <- as.numeric(Rnase_prohits_rep_merge[,6])
png("Prohits_replicates_Rnase.png")
plot(Rnase_rep1, Rnase_rep2, col="blue", main="Prohits")
dev.off()

### check correlation between prohits and maxquant


untreated_prohits_mq_rep1_merge <- merge(untreated_spc_prohits_rep1, untreated_spc_maxquant_rep1, by=c("Bait", "Prey"), all=F)
Rnase_prohits_mq_rep1_merge <- merge(Rnase_spc_prohits_rep1, Rnase_spc_maxquant_rep1, by=c("Bait", "Prey"), all=F)
colnames(untreated_prohits_mq_rep1_merge) <- c("Bait", "Prey", "Purification_Prohits", "Prohits", "Purification_Maxquant","Maxquant")
colnames(Rnase_prohits_mq_rep1_merge) <- c("Bait", "Prey", "Purification_Prohits", "Prohits", "Purification_Maxquant","Maxquant")

untreated_prohits <- as.numeric(untreated_prohits_mq_rep1_merge[,4])
untreated_maxquant <- as.numeric(untreated_prohits_mq_rep1_merge[,6])
png("Prohits_Maxquant_correlation_untreated_rep1.png")
plot(untreated_prohits, untreated_maxquant, col="blue", main="untreated")
dev.off()

Rnase_prohits <- as.numeric(Rnase_prohits_mq_rep1_merge[,4])
Rnase_maxquant <- as.numeric(Rnase_prohits_mq_rep1_merge[,6])
png("Prohits_Maxquant_correlation_Rnase_rep1.png")
plot(Rnase_prohits, Rnase_maxquant, col="blue", main="Rnase")
dev.off()

untreated_prohits_mq_rep2_merge <- merge(untreated_spc_prohits_rep2, untreated_spc_maxquant_rep2, by=c("Bait", "Prey"), all=F)
Rnase_prohits_mq_rep2_merge <- merge(Rnase_spc_prohits_rep2, Rnase_spc_maxquant_rep2, by=c("Bait", "Prey"), all=F)
colnames(untreated_prohits_mq_rep2_merge) <- c("Bait", "Prey", "Purification_Prohits", "Prohits", "Purification_Maxquant","Maxquant")
colnames(Rnase_prohits_mq_rep2_merge) <- c("Bait", "Prey", "Purification_Prohits", "Prohits", "Purification_Maxquant","Maxquant")

untreated_prohits <- as.numeric(untreated_prohits_mq_rep2_merge[,4])
untreated_maxquant <- as.numeric(untreated_prohits_mq_rep2_merge[,6])
png("Prohits_Maxquant_correlation_untreated_rep2.png")
plot(untreated_prohits, untreated_maxquant, col="blue", main="untreated")
dev.off()

Rnase_prohits <- as.numeric(Rnase_prohits_mq_rep2_merge[,4])
Rnase_maxquant <- as.numeric(Rnase_prohits_mq_rep2_merge[,6])
png("Prohits_Maxquant_correlation_Rnase_rep2.png")
plot(Rnase_prohits, Rnase_maxquant, col="blue", main="Rnase")
dev.off()





### compute compass score for maxquant 

system("perl C:/RSYNC/worksapce2015/MASS_PEC_SCORE/compass.pl ../maxquant/Rnase_for_HGscore.tab ../maxquant/Rnase_COMPASS_score.tab 0 0", show.output.on.console=T)
system("perl C:/RSYNC/worksapce2015/MASS_PEC_SCORE/compass.pl ../maxquant/untreated_for_HGscore.tab ../maxquant/untreated_COMPASS_score.tab 0 0", show.output.on.console=T)

NT_compass_mq <- read.table("../maxquant/untreated_COMPASS_score.tab", header=T, sep="\t")
Rnase_compass_mq <- read.table("../maxquant/Rnase_COMPASS_score.tab", header=T, sep="\t")

head(NT_compass_mq)
df <- 40
compass_mq_NT_cutoff <- spline_cutoff(as.numeric(NT_compass_mq$wd), df)
compass_mq_RN_cutoff <- spline_cutoff(as.numeric(Rnase_compass_mq$wd), df)
merge_compass_mq <- merge(NT_compass_mq, Rnase_compass_mq, by=c("bait", "prey"), all=F)
merge_compass_mq <- merge_compass_mq[,c("bait", "prey", "wd.x", "wd.y")]
colnames(merge_compass_mq) <- c("bait", "prey", "untreated", "Rnase")
png("../maxquant/Rnase_untreated_COMPASS_correlation_maxquant.png")
plot(log2(merge_compass_mq$untreated), log2(merge_compass_mq$Rnase))
abline(h=log2(compass_mq_RN_cutoff), v=log2(compass_mq_NT_cutoff))
dev.off()

selected <- NULL
selected <- merge_compass_mq[(merge_compass_mq$untreated > compass_mq_NT_cutoff | merge_compass_mq$Rnase > compass_mq_RN_cutoff),]
mutate(selected, note="BOTH")
selected[selected$untreated > compass_mq_NT_cutoff & selected$Rnase < compass_mq_RN_cutoff, "note"] <- "NT"
selected[selected$untreated < compass_mq_NT_cutoff & selected$Rnase > compass_mq_RN_cutoff, "note"] <- "RNASE"
selected[selected$untreated > compass_mq_NT_cutoff & selected$Rnase > compass_mq_RN_cutoff, "note"] <- ""
selected <- selected[order(selected$bait, selected$note),]
colnames(selected) <- c("bait", "prey", paste("untreated", sprintf("%.2f",compass_mq_NT_cutoff)), paste("Rnase", sprintf("%.2f", compass_mq_RN_cutoff)), "Note")

write.table(selected, "../maxquant/Rnase_untreated_COMPASS_difference_maxquant.tab", row.names=F, col.names=T, quote=F, sep="\t")


########### HGscore processing
NT_hg_mq <- read.table("../maxquant/untreated_ihgscore0p5.txt", header=T, sep="\t")
Rnase_hg_mq <- read.table("../maxquant/Rnase_ihgscore0p5.txt", header=T, sep="\t")

head(NT_hg_mq)
df <- 50
hg_mq_NT_cutoff <- spline_cutoff(as.numeric(NT_hg_mq$Score), df)
hg_mq_RN_cutoff <- spline_cutoff(as.numeric(Rnase_hg_mq$Score), df)
hg_mq_NT_cutoff
hg_mq_RN_cutoff

merge_hg_mq <- merge(NT_hg_mq, Rnase_hg_mq, by=c("Protein1", "Protein2"), all=F)
merge_hg_mq <- merge_hg_mq[,c("Protein1", "Protein2", "Score.x", "Score.y")]
colnames(merge_hg_mq) <- c("Protein1", "Protein2", "untreated", "Rnase")
png("../maxquant/Rnase_untreated_hg_correlation_maxquant.png")
plot(log2(merge_hg_mq$untreated), log2(merge_hg_mq$Rnase))
abline(h=log2(hg_mq_RN_cutoff), v=log2(hg_mq_NT_cutoff))
dev.off()

selected <- NULL
selected <- merge_hg_mq[(merge_hg_mq$untreated > hg_mq_NT_cutoff | merge_hg_mq$Rnase > hg_mq_RN_cutoff),]
mutate(selected, note="BOTH")
selected[selected$untreated > hg_mq_NT_cutoff & selected$Rnase < hg_mq_RN_cutoff, "note"] <- "NT"
selected[selected$untreated < hg_mq_NT_cutoff & selected$Rnase > hg_mq_RN_cutoff, "note"] <- "RNASE"
selected[selected$untreated > hg_mq_NT_cutoff & selected$Rnase > hg_mq_RN_cutoff, "note"] <- ""
selected <- selected[order(selected$Protein1, selected$note),]
colnames(selected) <- c("Protein1", "Protein2", paste("untreated", sprintf("%.2f",hg_mq_NT_cutoff)), paste("Rnase", sprintf("%.2f", hg_mq_RN_cutoff)), "Note")

write.table(selected, "../maxquant/Rnase_untreated_hg_difference_maxquant.tab", row.names=F, col.names=T, quote=F, sep="\t")


## SAINT score processing, SAINT not working: Segmentation fault (core dumped)
if(0){

NT_saint_mq <- read.table("../maxquant/untreated_saint_score.txt", header=T, sep="\t")
Rnase_saint_mq <- read.table("../maxquant/Rnase_saint_score.txt", header=T, sep="\t")

head(NT_saint_mq)
df <- 50
saint_mq_NT_cutoff <- spline_cutoff(as.numeric(NT_saint_mq$Score), df)
saint_mq_RN_cutoff <- spline_cutoff(as.numeric(Rnase_saint_mq$Score), df)
saint_mq_NT_cutoff
saint_mq_RN_cutoff

merge_saint_mq <- merge(NT_saint_mq, Rnase_saint_mq, by=c("Protein1", "Protein2"), all=F)
merge_saint_mq <- merge_saint_mq[,c("Protein1", "Protein2", "Score.x", "Score.y")]
colnames(merge_saint_mq) <- c("Protein1", "Protein2", "untreated", "Rnase")
png("../maxquant/Rnase_untreated_saint_correlation_maxquant.png")
plot(log2(merge_saint_mq$untreated), log2(merge_saint_mq$Rnase))
abline(h=log2(saint_mq_RN_cutoff), v=log2(saint_mq_NT_cutoff))
dev.off()

selected <- NULL
selected <- merge_saint_mq[(merge_saint_mq$untreated > saint_mq_NT_cutoff | merge_saint_mq$Rnase > saint_mq_RN_cutoff),]
mutate(selected, note="BOTH")
selected[selected$untreated > saint_mq_NT_cutoff & selected$Rnase < saint_mq_RN_cutoff, "note"] <- "NT"
selected[selected$untreated < saint_mq_NT_cutoff & selected$Rnase > saint_mq_RN_cutoff, "note"] <- "RNASE"
selected[selected$untreated > saint_mq_NT_cutoff & selected$Rnase > saint_mq_RN_cutoff, "note"] <- ""
selected <- selected[order(selected$Protein1, selected$note),]
colnames(selected) <- c("Protein1", "Protein2", paste("untreated", sprintf("%.2f",saint_mq_NT_cutoff)), paste("Rnase", sprintf("%.2f", saint_mq_RN_cutoff)), "Note")

write.table(selected, "../maxquant/Rnase_untreated_saint_difference_maxquant.tab", row.names=F, col.names=T, quote=F, sep="\t")

}





