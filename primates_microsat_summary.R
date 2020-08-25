setwd("/Users/jannawilloughby/GDrive/primate_heterozygosity/")
library(scales)
library(vioplot)

data = read.table("IUCN_t3.csv", header=T, sep=",")
colnames(data) = c("Publications", "Family", "Sp", "NMicrosats", "IUCN", "t_prim_old", "t_prim", "N", "Ne", "Ht_prim", "Ho_mean", "Ho_sd", "Ne", "Nlocations", "Nsamples")

length(table(data$Sp))
sum(data$Nsamples, na.rm=T)
sum(data$NMicrosats, na.rm=T)

####estimate and plot family H####
family = as.character(unique(data$Family))
family = family[order(family)]

#set up plot
pdf("../figure1.pdf", height=7, width=7, onefile=T)
par(bty='n')
plot(-100,-100, xlim=c(0.25,5.75), ylim=c(0,1.06), ylab="", xlab="", axes=F)
axis(side=1, at=seq(1,5,1), labels=family, tick=T, pos=0)
axis(side=2, at=seq(0,1,0.25), labels=seq(0,1,0.25), tick=T, pos=0.5)
segments(x0=0.5, x1=5.6, y0=0,    y1=0, lwd=1)
segments(x0=0.5, x1=0.5, y0=0,    y1=1.05, lwd=1)
segments(x0=0.5, x1=5.6, y0=1.05, y1=1.05, lwd=1)
segments(x0=5.6, x1=5.6, y0=0,    y1=1.05, lwd=1)
segments(x0=0.5, x1=5.6, y0=0,    y1=0, lwd=1)
mtext(text="heterozygosity", side=2)

#calcualte median and plot data
OUT = NULL
for(f in 1:length(family)){
  t = data[data$Family==as.character(family[f]),,drop=F]
  t$xval = f - sample(seq(0,0.2, 0.005), nrow(t))
  points(x=t$xval, y=t$Ho_mean, col=alpha("grey50", 1), pch=19, cex=0.75)
  vioplot(x=t$Ho_mean, col = alpha("grey50", 0.9), plotCentre = "line", side = "right", ylim=c(0,0.3), add=T, at=(f+0.1))
  OUT = rbind(OUT, c(as.character(family[f]), mean(t$Ho_mean, na.rm=T), median(t$Ho_mean, na.rm=T), sd(t$Ho_mean, na.rm=T), sum(t$Nlocations)))
}
familydata = as.data.frame(OUT)
colnames(familydata) = c("family", "Ho_mean", "Ho_median","Ho_sd", "npops")
write.table(familydata, "../familydata.csv", col.names = T, row.names = F, sep=",")
dev.off()

####IUCN rank mean/sd H####
IUCN = unique(data$IUCN)
OUT = NULL
for(i in 1:length(IUCN)){
  t = data[data$IUCN==as.character(IUCN[i]),,drop=F]
  OUT = rbind(OUT, c(as.character(IUCN[i]), mean(t$Ho_mean, na.rm=T), sd(t$Ho_mean, na.rm=T), sum(t$Nlocations)))
}
IUCNdata = as.data.frame(OUT)
colnames(IUCNdata) = c("IUCN", "Ho_mean", "Ho_sd", "npops")
write.table(IUCNdata, "../IUCNdata.csv", col.names = T, row.names = F, sep=",")

####IUCN rank mean/sd H####
IUCN = unique(data$IUCN)
OUT = NULL
for(i in 1:length(IUCN)){
  t = data[data$IUCN==as.character(IUCN[i]),,drop=F]
  t$t_prim[t$t_prim<0] = 0
  OUT = rbind(OUT, c(as.character(IUCN[i]), mean(t$t_prim, na.rm=T), sd(t$t_prim, na.rm=T), sum(t$Nlocations)))
}
IUCNTdata = as.data.frame(OUT)
colnames(IUCNTdata) = c("IUCN", "t_prim", "Ho_sd", "npops")
write.table(IUCNTdata, "../IUCNTdata.csv", col.names = T, row.names = F, sep=",")

####t cutoff IUCN rank####
IUCN = unique(data$IUCN)
OUT = NULL
for(i in 1:length(IUCN)){
  t = data[data$IUCN==as.character(IUCN[i]),,drop=F]
  below0 = below10 = below50 = below1000 = over100 = 0
  below0   = sum(t$Nlocations[t$t_prim<0])
  below10  = sum(t$Nlocations[t$t_prim<10])  - below0
  below50  = sum(t$Nlocations[t$t_prim<50])  - below0 - below10
  below100 = sum(t$Nlocations[t$t_prim<100]) - below0 - below10 - below50
  over100  = sum(t$Nlocations[t$t_prim>=100])

  OUT = rbind(OUT, c(as.character(IUCN[i]), below0, below10, below50, below100, over100, sum(t$Nlocations)))
}
IUCNttdata = as.data.frame(OUT)
colnames(IUCNttdata) = c("IUCN", "below0", "below10", "below50", "below100", "over100", "npops")
write.table(IUCNttdata, "../IUCNttdata.csv", col.names = T, row.names = F, sep=",")


