markerfile <- read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/ooa5k.marker.txt")
mfile <- cbind(markerfile$MarkerID, (markerfile$AlleleCount1/(markerfile$AlleleCount0+markerfile$AlleleCount1)))   #mfile store the frequency of all markers from the simulation, e.g. N3
#mfrare1 <- mfile[mfile[,2] <= 0.01, ]           #all markers with derived allele <= 1% (used when we want to consider only rare variants under a certain threshold)
#mfrare <- mfrare1[,1]              #all markers with derived allele <= 1%
mfrare <- mfile[,1] #used when all variants are considered (have to neutralize it when graphs/stats for rare variants)

a1<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_1.txt")      #true vs inferred TMRCA calculated thanks to patrick TMRCA algorithm (cf link in report)
a2<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_2.txt")
a3<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_3.txt")
a4<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_4.txt")
a5<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_5.txt")
a6<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_6.txt")
a7<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_7.txt")
a8<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_8.txt")
a9<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_9.txt")
a10<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_10.txt")
asum=rbind(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)    #all true vs inferred TMRCA from N3 gathered together
asum=unique(asum$MarkerID)         
asum=as.data.frame(asum)                              #we make a dataframe of the markers dated
names(asum) <- "asum"
e1<-read.table(header=FALSE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_real_1000_allele_mut1.txt")    #file created by my algorithm to find and store "true" allele age, considering mutation clock
e1=cbind(asum, e1)
e1=subset(e1, e1$asum %in% mfrare)     #we only consider variants present in mfrare (usefull when we want to consider only a specific set of variants among the 1000)
e1=e1$V1
e1=as.data.frame(e1)
e2<-read.table(header=FALSE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_real_1000_allele_rec1.txt")
e2=cbind(asum, e2)
e2=subset(e2, e2$asum %in% mfrare)
e2=e2$V1
e2=as.data.frame(e2)
e3<-read.table(header=FALSE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_real_1000_allele_comb1.txt")
e3=cbind(asum, e3)
e3=subset(e3, e3$asum %in% mfrare)
e3=e3$V1
e3=as.data.frame(e3)

f1<-read.table(header=FALSE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_estimate_1000_allele_mut1.txt")    #file created by my algorithm to store estimated allele age, considering mutation clock
f1=cbind(asum, f1)
f1=subset(f1, f1$asum %in% mfrare)
f1=f1$V1
f1=as.data.frame(f1)
f2<-read.table(header=FALSE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_estimate_1000_allele_rec1.txt")
f2=cbind(asum, f2)
f2=subset(f2, f2$asum %in% mfrare)
f2=f2$V1
f2=as.data.frame(f2)
f3<-read.table(header=FALSE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_estimate_1000_allele_comb1.txt")
f3=cbind(asum, f3)
f3=subset(f3, f3$asum %in% mfrare)
f3=f3$V1
f3=as.data.frame(f3)
names(e1)<-"V1"
names(e2)<-"V1"
names(e3)<-"V1"
names(f1)<-"V1"
names(f2)<-"V1"
names(f3)<-"V1"

e=rbind(e1, e2, e3)         #the "true" ages, all clock models considered
f=rbind(f1, f2, f3)         #the estimated allele ages, all clock models considered

RMSL10E=function(a,b) {       #function used to calculate RMSLE
  n=length(a)
  a=log10(a)
  b=log10(b)
  sqrt(sum((a-b)**2/n))
}

clock=vector(mode = "character", length=length(e[[1]]))
clock1=vector(mode = "character", length=3)
r2=vector(mode = "numeric", length=3)
rh=vector(mode = "numeric", length=3)
rmsle=vector(mode = "numeric", length=3)

i=1
k=1
frontiers=c(0, length(e1[[1]]), length(e2[[1]]), length(e3[[1]]))
ranges=vector(mode = "numeric", length=4) 
for (m in 1:4) {
  ranges[m]=sum(frontiers[1:m])           
}
type_clock=c("Mutation clock", "Recombination clock", "Combined clock")
while (k<=length(e[[1]])) {
  if (ranges[i]<k & k<=ranges[i+1]) {
    clock[k]=type_clock[i]
    if (k==ranges[i]+1) {
      clock1[i]=type_clock[i]
      b=ranges[i]+1
      c=ranges[i+1]
      b0=cor(log(e$V1[b:c]), log(f$V1[b:c]))**2   #R², not used in the final report
      r2[i]=as.numeric(format(b0, digits=3))
      b1=cor.test(~e$V1[b:c]+f$V1[b:c], cbind(e,f), method="spearman", continuity=FALSE, conf.level=0.95)
      b1=b1$estimate     
      rh[i]=as.numeric(format(b1, digits=3))
      b2=RMSL10E(e$V1[b:c], f$V1[b:c])
      rmsle[i]=as.numeric(format(b2, digits=3))
    } else if (k==ranges[i+1]) {
      i=i+1
    }
    
  }
  k=k+1
}


az=cbind(e, f, clock)                 #the dataset used to create density plots with ggplot2

#install.packages("ggplot2")
#install.packages("cowplot")
#install.packages("Metrics")
#install.packages("ggthemes")
library(ggthemes)
library(Metrics)
library(ggplot2)
library(cowplot)
library(broom)
library(dplyr)

clock=clock1
dataz=cbind(clock, r2, rh, rmsle)
datazz=as.data.frame(dataz)          #the dataset of the metrics corresponding to the graphs

g = ggplot(az, aes(x=az[[1]], y=az[[2]]))
g = g + facet_grid(.~clock) +
  geom_raster(aes(az[[1]], az[[2]], fill = (..count..)/tapply(..count..,..PANEL..,max)[..PANEL..]), position = "identity", stat = "bin2d", binwidth = c(0.1, 0.1)) +
  geom_tile(aes(az[[1]], az[[2]], fill = (..count..)/tapply(..count..,..PANEL..,max)[..PANEL..]), position = "identity", stat = "bin2d", binwidth =   c(0.1, 0.1)) +
  geom_abline(intercept = c(0,0), slope = 1, alpha = 1/2)
g = g + geom_text(data=datazz, inherit.aes=FALSE, size=6, face="bold", aes(x=3000, y=3, label=paste("rho =", rh, "\n RMSLE =", rmsle)))
g = g + coord_cartesian(xlim = c(0.5, 99000), ylim = c(0.5, 99000), expand = F) +
  scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
  scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
  scale_fill_gradientn(colours = c("ivory", "firebrick2", "firebrick3", "firebrick4"), na.value = "white", limits = c(-0.05,1.05), breaks = c(0,0.25,0.5,0.75,1)) +
  theme_few() +
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 9),
        panel.border = element_rect(fill = NA, colour = "grey20", size = 2/3),
        panel.background = element_rect(fill = NA, colour = "grey20", size = 2/3),
        strip.text.x = element_text(face = "bold", size = 11),
        strip.text.y = element_text(size = 13),
        legend.justification=c(0,1), legend.position=c(1 - (3*0.985)/3, (1*0.98)/1),
        legend.background = element_rect(fill = "grey90", colour = "grey70", size = 1/4),
        axis.ticks = element_line(colour = "grey20"),
        axis.text.x = element_text(hjust = 0),
        axis.title = element_text(size = 11),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.margin = margin(-1,1,0.5,0.5, "mm"),
        legend.text = element_text(size = 7),
        legend.title = element_blank()) +
  ylab("Estimated allele age (generations)") + xlab("True allele age (generations)")

