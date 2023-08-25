
library("MASS")
library(plyr)
setwd("C:/RSYNC/LUMIER_C2H2/Nov_10_2016_EM")

#SNDREAD and RN_LYSIS are interchangeable

INT <- read.table("inter_Nov_10_2016_EM.dat", header=F, stringsAsFactors = F)
PREY <- read.table("prey_Nov_10_2016_EM.dat", header=F, stringsAsFactors = F)
colnames(INT) <- c("purification", "bait", "prey", "total_spc", "uniq_spc")
colnames(PREY) <- c("prey", "mol_weight")

for(i in 1:dim(PREY)[1]){
  if(grepl(",", PREY[i,1])){
    print(PREY[i,1])
    temp <- unlist(strsplit(as.character(PREY[i, 1]), ", ", fixed=T))
    PREY[i,1] <- temp[[1]]
    
    print(temp[[1]])
   
  }
}

PREY <- aggregate(PREY$mol_weight, by=list(PREY$prey), FUN=max)  #remove duplicated preys
colnames(PREY) <- c("prey", "mol_weight")

head(PREY)
dim(PREY)
dim(INT)

total_num_purification <- length(unique(INT[,"purification"]))
total_num_prey <- length(unique(INT[,"prey"]))
total_num_bait <- length(unique(INT[,"bait"]))

int_median_t <- aggregate(INT$total_spc, by=list(INT$prey), FUN=median)
int_median_u <- aggregate(INT$uniq_spc, by=list(INT$prey), FUN=median)
colnames(int_median_t) <- c("prey", "median_total_spc")  # median of total spectral count
colnames(int_median_u) <- c("prey", "median_uniq_spc")  # median of unique spectral count

int_prey_frequency <- count(INT, "prey")

colnames(int_prey_frequency) <- c("prey", "frequency")

int_prey_frequency <- mutate(int_prey_frequency, percentage=(frequency*100)/total_num_purification)


int_median <- merge(int_median_t, int_median_u)
int_median

int_fr_mw <- merge(PREY, int_prey_frequency)

prey_table <- merge(int_fr_mw, int_median)
prey_table <- arrange(prey_table, desc(percentage))
dim(prey_table)
head(prey_table)

write.table(prey_table, "Prey_stats.tab", sep="\t", col.names=NA)
cat(c(paste("\nNumber of baits",total_num_bait, sep="\t"),
      paste("Number of preys",total_num_prey, sep="\t"),
      paste("Number of purifications",total_num_purification, sep="\t")), 
    file="Prey_stats.tab", sep="\n", append=T)

prey_frequency <- int_prey_frequency$percentage

pdf("Density_plot_prey_freqency.pdf")
bk <- seq(0,100, 0.1)
hist(prey_frequency, breaks=bk, freq = F)
#plot(ecdf(x))
dev.off()

x=prey_frequency

# fit normal distribution
fit <- fitdistr(prey_frequency, "normal")
fit
para1 <- fit$estimate[1]
para2 <- fit$estimate[2]

plot(density(x))
curve(dnorm(x, mean=para1, sd=para2), col="magenta", add=T)
qqplot(rnorm(length(x), mean=para1, sd=para2), x)
qqline(x)

frequency_cutoff <- qnorm(0.95, mean=para1, sd=para2)
frequency_cutoff

# fit exponential distribution
fit <- fitdistr(prey_frequency, "exponential")
fit
para1 <- fit$estimate

plot(density(x))
curve(dexp(x, para1), col="magenta", add=T)
qqplot(rexp(length(x), para1), x)
qqline(x)

frequency_cutoff <- qexp(0.95, para1)
frequency_cutoff

# fit gamma distribution
fit <- fitdistr(prey_frequency, "gamma")
fit
para1 <- fit$estimate[1]
para2 <- fit$estimate[2]

plot(density(x))
curve(dgamma(x, shape=para1, rate=para2),  col="magenta", add=T)
qqplot(rgamma(length(x), shape=para1, rate=para2), x)
qqline(x)

frequency_cutoff <- qgamma(0.95, shape=para1, rate=para2)
frequency_cutoff

# fit weibull distribution
fit <- fitdistr(prey_frequency, densfun=dweibull,start=list(scale=2,shape=1))
fit
para1 <- fit$estimate[1]
para2 <- fit$estimate[2]

plot(density(x))
curve(dweibull(x, scale=para1, shape=para2),  col="magenta", add=T)
qqplot(rweibull(length(x), scale=para1, shape=para2), x)
qqline(x)

frequency_cutoff <- qweibull(0.99, scale=para1, shape=para2)
frequency_cutoff





