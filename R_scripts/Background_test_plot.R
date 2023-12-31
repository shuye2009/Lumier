
library("MASS")

wdir <- "C:/RSYNC/LUMIER_C2H2/Archive_oct032018/Archive1_3/For_Dicer"
setwd(wdir)
BG <- read.csv(file.path(wdir,"Background_test.csv"), header=T)
head(BG)
lum <- log(BG$Lum)
plate1 <- matrix(lum[1:96], ncol=12, byrow=T)
plate2 <- matrix(lum[97:192], ncol=12, byrow=T)
plate3 <- matrix(lum[193:288], ncol=12, byrow=T)
plate4 <- matrix(lum[289:384], ncol=12, byrow=T)
plate_matrix <- rbind(cbind(plate1, plate2), cbind(plate3, plate4))
head(plate_matrix)
par(mfrow=c(2,2))
image(t(plate1[c(8:1),]), col=topo.colors(96))
image(t(plate2[c(8:1),]), col=topo.colors(96))
image(t(plate3[c(8:1),]), col=topo.colors(96))
image(t(plate4[c(8:1),]), col=topo.colors(96))

filled.contour(t(plate1[c(8:1),]), col=topo.colors(30))
filled.contour(t(plate2[c(8:1),]), col=topo.colors(30))
filled.contour(t(plate3[c(8:1),]), col=topo.colors(30))
filled.contour(t(plate4[c(8:1),]), col=topo.colors(30))

par(mfrow=c(1,1))
png("plate_image_pure_background.png", 800, 600)
image(t(plate_matrix[c(16:1),]), col=topo.colors(384))
dev.off()
png("plate_contour_pure_background.png", 800, 600)
filled.contour(t(plate_matrix[c(16:1),]), col=topo.colors(20))
dev.off()
cor_matrix1 <- cor(plate1)
cor_matrix1t <- cor(t(plate1))
heatmap.2(cor_matrix1, trace="none", col=topo.colors(20))
heatmap.2(cor_matrix1t, trace="none", col=topo.colors(20))
heatmap.2(plate_matrix[c(16:1),], trace="none", col=topo.colors(20))
png("histogram_of_pure_background.png", 800, 600)
par(mfrow=c(2,2))
hist(lum[1:96], xlim=c(8,16), main="H2O", xlab="log(lum)")
hist(lum[97:192], xlim=c(8,16), main="GFP", xlab="log(lum)")
hist(lum[193:288], xlim=c(8,16), main="ILF2", xlab="log(lum)")
hist(lum[289:384], xlim=c(8,16), main="TRIM28", xlab="log(lum)")
dev.off()

png("Density_plot_of_pure_background.png", 800, 600)
par(mfrow=c(1,1))
plot(NULL, NULL, xlim=c(8,16), ylim=c(0,5), xlab="log(lum)", ylab="Density", main="Density plot of background tests")
lines(density(lum[1:96]), col="red", lwd=2)
lines(density(lum[97:192]), col="blue", lwd=2)
lines(density(lum[193:288]), col="cyan", lwd=2)
lines(density(lum[289:384]), col="green", lwd=2)

legend("topleft", legend=c("H2O", "GFP", "ILF2", "TRIM28"), col=c("red", "blue", "cyan", "green"), lwd=2)
dev.off()