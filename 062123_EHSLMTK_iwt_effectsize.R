setwd("/Users/sdy897/Library/CloudStorage/OneDrive-UniversityofTexasatSanAntonio/Research/BRL-liver/testing-SF/Rcode")

library(readxl)
library(MASS)
source("Pini_source_edit.R")

############
# EHS and LMTK
tmp.EHS = as.data.frame(read_xlsx("EHS+LMTK.xlsx", sheet=3, skip=2, col_names = TRUE))
freq.E = as.matrix(tmp.EHS[,seq(1,ncol(tmp.EHS)-1, by=2)])
dat.E = 10*log10(tmp.EHS[,seq(2,ncol(tmp.EHS), by=2)])

tmp.LMTK = as.data.frame(read_xlsx("EHS+LMTK.xlsx", sheet=4, skip=2, col_names = TRUE))
freq.L = as.matrix(tmp.LMTK[,seq(1,ncol(tmp.LMTK)-1, by=2)])
dat.L = 10*log10(tmp.LMTK[,seq(2,ncol(tmp.LMTK), by=2)])

###########
# Figure for paper
par(mar = c(4, 4, 2, 2))
freq.all = sort(unique(c(as.numeric(freq.L), as.numeric(freq.E))))
plot(freq.E[,1], dat.E[,1], xlim=c(min(freq.all),100), ylim=c(min(dat.L), -0.7), type="n", xlab="Frequency (MHz)", ylab="BSC (dB)")
for(j in 1:ncol(dat.L)){
  lines(freq.L[,j], dat.L[,j], col="grey60", lwd=2)
}
for(i in 1:ncol(dat.E)){
  lines(freq.E[,i], dat.E[,i], lwd=2)
}
legend("bottomright", lwd=2, col=c("black", "grey60"), legend=c("EHS","LMTK"))


###########
# interpolation
freq <- seq(12, 100, by = 0.4)

Y.E=matrix(nrow=ncol(dat.E), ncol=length(freq))
Y.L=matrix(nrow=ncol(dat.L), ncol=length(freq))

for(i in 1:ncol(dat.E)){
  approx.f=approxfun(freq.E[,i], dat.E[,i])
  Y.E[i,]=approx.f(freq)
}

for(i in 1:ncol(dat.L)){
  approx.f=approxfun(freq.L[,i], dat.L[,i])
  Y.L[i,]=approx.f(freq)
}

##############
# M-est
huber.k=0.5
huber.unscale.E=c(); huber.E=c(); med.E=c()
for (t in 1:length(freq)){
  tmp.s=mad(Y.E[,t],na.rm=TRUE)
  huber.E=c(huber.E, huber(Y.E[,t], k = huber.k, tol = 1e-06)$mu)
  huber.unscale.E=c(huber.unscale.E, huber(Y.E[,t], k = huber.k/tmp.s, tol = 1e-06)$mu)
  med.E = c(med.E, median(Y.E[,t], na.rm=TRUE))
}
mean.E = colMeans(Y.E, na.rm=TRUE)

plot(huber.E, type="l");lines(huber.unscale.E, col="green")
lines(med.E, col="red"); lines(mean.E, col="blue")


huber.unscale.L=c(); huber.L=c(); med.L=c()
for (t in 1:length(freq)){
  tmp.s=mad(Y.L[,t],na.rm=TRUE)
  huber.L=c(huber.L, huber(Y.L[,t], k = huber.k, tol = 1e-06)$mu)
  huber.unscale.L=c(huber.unscale.L, huber(Y.L[,t], k = huber.k/tmp.s, tol = 1e-06)$mu)
  med.L = c(med.E, median(Y.L[,t], na.rm=TRUE))
}
mean.L = colMeans(Y.L, na.rm=TRUE)

plot(huber.L, type="l");lines(huber.unscale.L, col="green")
lines(med.L, col="red"); lines(mean.L, col="blue")

########
n.E = nrow(Y.E); n.L=nrow(Y.L); n=n.E+n.L
theta.E = huber.E; theta.L = huber.L

theta.mean=(n.E*theta.E + n.L*theta.L)/n
#################################
# Robust iwt
################################
#load("EHSKMTK.comb.res.040523.RData")



##################################
# Robust effect size
##################################

############
# estimat covariance matrix using bootstrp samples
boot.rep=3000

huber.E.mat=c();
huber.L.mat=c();
#boot.freq = freq[-c]
for (i in 1051:boot.rep){
  
  # resample on raw data
  boot.E=Y.E[sample(1:nrow(Y.E),nrow(Y.E),replace=TRUE), ]
  boot.L=Y.L[sample(1:nrow(Y.L),nrow(Y.L),replace=TRUE), ]
  
  id.E=which( apply(boot.E, 2, mad, na.rm=TRUE)>0)
  id.L=which(apply(boot.L, 2, mad, na.rm=TRUE)>0)

  ##################
  # huber estimation - scaled
  huber.E.boot=c()
  for (te in id.E){
    #huber.E.boot=c(huber.E.boot, huber(boot.E[,t], k = huber.k, tol = 1e-06)$mu)
    huber.E.boot[te]=huber(boot.E[,te], k = huber.k, tol = 1e-06)$mu
  }
  
  huber.L.boot=c()
  for (tl in id.L){
    #huber.L.boot=c(huber.L.boot, huber(boot.L[,t], k = huber.k, tol = 1e-06)$mu)
    huber.L.boot[tl]=huber(boot.L[,tl], k = huber.k, tol = 1e-06)$mu
  }
  
  huber.E.mat=rbind(huber.E.mat,huber.E.boot)
  huber.L.mat=rbind(huber.L.mat,huber.L.boot)
  
  if(i %%300 == 0) print(i)
}

##############################################
# SSR
SSR=(n.L*(theta.L-theta.mean)^2 + n.E*(theta.E-theta.mean)^2)/n

#############################################
# SSE
boot.mean.E=colMeans(huber.E.mat, na.rm=TRUE)
boot.mean.L=colMeans(huber.L.mat, na.rm=TRUE)

mean.mat.E=t(matrix(boot.mean.E,ncol=nrow(huber.E.mat), nrow=length(freq)))
mean.mat.L=t(matrix(boot.mean.L,ncol=nrow(huber.L.mat), nrow=length(freq)))

# cov matrix
tmp.E=rbind((huber.E.mat-mean.mat.E),(huber.L.mat-mean.mat.L))
tmp.E[is.na(tmp.E)] <- mean(tmp.E, na.rm = TRUE)
cov.E=t(tmp.E)%*%tmp.E/nrow(tmp.E)

# pointwise effect size
fSNR=(SSR)/(diag(cov.E))

#plot(freq, sqrt(fSNR), type="l")

size.mat=matrix(nrow=length(freq), ncol=length(freq))
#size.mat[p,]<- sqrt(diff/sigma.2)
size.mat[1,]<- sqrt(fSNR)

for(pp in 1:(length(freq)-1)){
  for(k in 1:length(freq)){
    tmp.id = c((k-pp):(k+pp))
    inc.id = tmp.id[tmp.id >= 1 & tmp.id <= length(freq)]
    size.mat[ (pp+1),k] =sqrt(mean(fSNR[inc.id]))
  }
  if(pp%%30==0) print(pp)
}

#save(fSNR, size.mat, file="EHSKMTK.robustES.RData")
#################
# plot
load("EHSKMTK.comb.res.040523.RData")
load("EHSKMTK.robustES.RData")

#heatmap(size.mat, scale="none", Rowv = NA, Colv = NA, labRow="",labCol=freq )
nlevel=10
#col_heatmap = rainbow(nlevel, start = 0.15, end = 0.67)[nlevel:1]
col_heatmap = rev(hcl.colors(20, "terrain"))

#####
dev.new()

layout(rbind(cbind(1:2, 0),3:4), widths = c(9, 1)) #, heights = c(4, 4, 4))
#par(mar = c(4, 3, 4, 3))
par(mar = c(5, 5, 2.5, 2))
# raw data
plot(freq, Y.E[1,], ylim=c(min(c(Y.E,Y.L),na.rm=TRUE), max(c(Y.E, Y.L), na.rm=TRUE)),xaxs="i", type="n", 
     cex=1.5, cex.axis=1, cex.lab=1.2, xlab="Frequency (MHz)", ylab="BSC(dB)")
for(i in 1:nrow(Y.E)){
  lines(freq, Y.E[i,], col="gray24", lty=2)
}
for(j in 1:nrow(Y.L)){
  lines(freq, Y.L[j,], col="lightblue3",lty=2)
}
#lines(freq,mean.E, lwd=3)
#lines(freq, mean.L, lwd=3, col="blue")
lines(freq,theta.E, lwd=3)
lines(freq,theta.L, lwd=3, col="royalblue1")
rect(min(res.f)-0.8,-50,max(res.f)-0.8,-5,col = rgb(0.5,0.5,0.5,1/4), border = F)
legend("bottomright",col=c("gray24","black","lightblue3","royalblue1"), lty=c(2,1,2,1), lwd=c(1,2,1,2), 
       legend = c("EHS tumor","M-estimator of EHS","LMTK tumor","M-estimator of LMTK" ),cex=1.2,bg="transparent")

# p-val
plot(freq[-c(1,2)], min.p, type="l", xaxs="i",xlab="Frequency (MHz)", ylab="p-value", 
     cex=1.5, cex.axis=1, cex.main=1.5, cex.lab=1.2, main="Adjusted p-values")
lines(freq[-c(1,2)], res.comb$corrected.pval, lty=2, col="purple", lwd=2)
abline(h=0.0253, col="red",lty=3) # for two min
rect(min(res.f),0,max(res.f),3,col = rgb(0.5,0.5,0.5,1/4), border = F)
legend("topleft",col=c("black","purple"), lwd=2,lty=c(1,2), 
       legend = c("min p-values", "adjusted p-values" ),cex=1.5,bg="transparent")


# effect size map
#mar.left=4
#par(mar = c(4, mar.left, 4, 1.5), mgp = c(2.5, 1, 0), xpd = FALSE, las = 1)
image(freq, freq, t(size.mat), col=col_heatmap,
      ylab = expression(paste("Interval length ", "(",Delta,")")), xlab="Frequency (MHz)", main="Effect size heatmap", cex.main=1.5, cex.lab=1.2)
box()
par(mar = c(4, 0.5, 2, 3.5))
image(1, seq.int(nlevel) - (nlevel + 1)/2, t(as.matrix(seq.int(nlevel))), 
      axes = FALSE, xlab = "", ylab = "", col = col_heatmap)
axis(side = 4, at = seq(0, nlevel, length.out = 6) - 
       nlevel/2, labels = format(round(seq(0, max(size.mat), length.out = 6)), 2))
box()



###########################################


