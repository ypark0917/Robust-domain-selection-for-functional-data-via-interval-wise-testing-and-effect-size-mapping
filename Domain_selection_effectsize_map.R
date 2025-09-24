library(KernSmooth)
library(reshape2)
library(mvtnorm)
library(MASS)

# illustration using one simulated set of partially observed t(3) process under sigma.1 = 3
source("Rselection_source_functions.R")
load("sim_example_t3_sigma3.RData")
freq <- seq(0, 1, length = 97) 

# pre-smoothing
h.cv0 <- 1 # find the optimal
dat.ker0 <- matrix(nrow=nrow(dat),ncol=ncol(dat))
for (i in 1:nrow(dat)) {
  tmp <- locpoly(x = freq[which(delta[i,]==1)], y = par.dat[i,which(delta[i,]==1)], 
                 kernel = 'epanech', bandwidth = h.cv0, degree = 1,  
                 range.x = range(freq), 
                 gridsize = length(freq))$y
  dat.ker0[i, which(delta[i,]==1)] <- tmp[which(delta[i,]==1)]
}

## l=1 (pre-smoothing)
h.cv1 <- 1.5 # find the optimal
dat.ker1 <- matrix(nrow=nrow(dat),ncol=ncol(dat))
for (i in 1:nrow(dat)) {
  tmp <- locpoly(x = freq[which(delta[i,]==1)], y = par.dat[i,which(delta[i,]==1)], 
                 kernel = 'epanech', bandwidth = h.cv1, degree = 1, drv = 1,
                 range.x = range(freq), 
                 gridsize = length(freq))$y
  dat.ker1[i, which(delta[i,]==1)] <- tmp[which(delta[i,]==1)]
}

n=50

######################################################################
######################################################################
#
#   1. Robust domain selection
#
######################################################################
######################################################################

# order 0 
dat.1 <- dat.ker0[1:n,]
dat.2<- dat.ker0[(n+1):(2*n),]

  dat.1.xy <- data.frame(cbind(rep(freq, n),as.vector(t(dat.1))))
  dat.1.xy <- dat.1.xy[!is.na(dat.1.xy[,2]), ]
  dat.1.xy <- dat.1.xy[order(dat.1.xy[,1]), ]
  
  dat.2.xy <- data.frame(cbind(rep(freq, n),as.vector(t(dat.2))))
  dat.2.xy <- dat.2.xy[!is.na(dat.2.xy[,2]), ]
  dat.2.xy <- dat.2.xy[order(dat.2.xy[,1]), ]
  
  fit.huber.1<- robust_smooth_spline(dat.1.xy, df = 10,max_iter = 50, delta = 1.345, spar = spar.set, tol = 1e-4)
  fit.huber.2 <- robust_smooth_spline(dat.2.xy, df = 10,max_iter = 50, delta = 1.345, spar = spar.set,  tol = 1e-4)

  huber.1.sp <- fit.huber.1$fitted.values
  huber.2.sp <- fit.huber.2$fitted.values
  
  approx.1=approxfun(dat.1.xy[,1], huber.1.sp)
  approx.2=approxfun(dat.2.xy[,1], huber.2.sp)
  
  huber.1=approx.1(freq)
  huber.2 = approx.2(freq)
  
  huber.1.0 = huber.1
  huber.2.0 = huber.2
  huber.mean <- (huber.1.0 + huber.2.0)/2
    
  diff=(huber.1.0 - huber.2.0)^2
  
  diff.mean.0=c(); 
  for(j in 1:(length(diff)-1)){
    diff.mean.0=c(diff.mean.0, mean(diff[j:(j+1)]))  # interval 
  }


B=1000
boot.diff.mat.0=c()
for(b in 1:B){
  boot.dat = dat[sample(1:nrow(dat), size=nrow(dat)) , ]
  
  boot.dat.1 <- boot.dat[1:nrow(dat.1),];
  boot.dat.2 <- boot.dat[(nrow(dat.1)+1):nrow(dat),];
  
  ##############################################
    ## penalized spline with robust loss function
    ################################################
    dat.1.xy.boot <- data.frame(cbind(rep(freq, n),as.vector(t(boot.dat.1))))
    dat.1.xy.boot <- dat.1.xy.boot[!is.na(dat.1.xy.boot[,2]), ]
    dat.1.xy.boot <- dat.1.xy.boot[order(dat.1.xy.boot[,1]), ]
    
    dat.2.xy.boot<- data.frame(cbind(rep(freq, n),as.vector(t(boot.dat.2))))
    dat.2.xy.boot <- dat.2.xy.boot[!is.na(dat.2.xy.boot[,2]), ]
    dat.2.xy.boot <- dat.2.xy.boot[order(dat.2.xy.boot[,1]), ]
    
    fit.huber.1.boot<- robust_smooth_spline(dat.1.xy.boot, df = 10, max_iter = 50, delta = 1.345, spar = spar.set, tol = 1e-4)
    fit.huber.2.boot <- robust_smooth_spline(dat.2.xy.boot, df = 10,max_iter = 50, delta = 1.345, spar = spar.set,  tol = 1e-4)
    
    huber.1.sp.boot <- fit.huber.1.boot$fitted.values
    huber.2.sp.boot <- fit.huber.2.boot$fitted.values
    
    approx.1.boot=approxfun(dat.1.xy.boot[,1], huber.1.sp.boot)
    approx.2.boot=approxfun(dat.2.xy.boot[,1], huber.2.sp.boot)
    
    boot.huber.1=approx.1.boot(freq)
    boot.huber.2 = approx.2.boot(freq)
    
    boot.huber.mean <- ( boot.huber.1 +  boot.huber.2)/2
    boot.diff=(boot.huber.1 - boot.huber.2)^2

    boot.diff.mean=c();
    for(jj in 1:(length(boot.diff)-1)){
      boot.diff.mean=c(boot.diff.mean, mean(boot.diff[jj:(jj+1)]))
    }
    
    boot.diff.mat.0 = rbind(boot.diff.mat.0, boot.diff.mean)

# unadjusted p-val function for order 0
pt.pval.0=c()
for(j in 1:length(diff.mean.0)){
  pt.pval.0=c(pt.pval.0, mean(boot.diff.mat.0[,j] >= diff.mean.0[j] ))
}

####################################################################
# order 1
dat.1 <- dat.ker1[1:n,]
dat.2<- dat.ker1[(n+1):(2*n),]

 dat.1.xy <- data.frame(cbind(rep(freq, n),as.vector(t(dat.1))))
  dat.1.xy <- dat.1.xy[!is.na(dat.1.xy[,2]), ]
  dat.1.xy <- dat.1.xy[order(dat.1.xy[,1]), ]
  
  dat.2.xy <- data.frame(cbind(rep(freq, n),as.vector(t(dat.2))))
  dat.2.xy <- dat.2.xy[!is.na(dat.2.xy[,2]), ]
  dat.2.xy <- dat.2.xy[order(dat.2.xy[,1]), ]
  
  fit.huber.1<- robust_smooth_spline(dat.1.xy, df = 10,max_iter = 50, delta = 1.345, spar = spar.set, tol = 1e-4)
  fit.huber.2 <- robust_smooth_spline(dat.2.xy, df = 10,max_iter = 50, delta = 1.345, spar = spar.set,  tol = 1e-4)
  
  huber.1.sp <- fit.huber.1$fitted.values
  huber.2.sp <- fit.huber.2$fitted.values
  
  approx.1=approxfun(dat.1.xy[,1], huber.1.sp)
  approx.2=approxfun(dat.2.xy[,1], huber.2.sp)
  
  huber.1=approx.1(freq)
  huber.2 = approx.2(freq)
  
  huber.1.1 = huber.1
  huber.2.1 = huber.2
    
  huber.mean <- (huber.1.1 + huber.2.1)/2
  
  diff=(huber.1.1-huber.2.1)^2

  
  diff.mean.1=c(); diff.mean.1.rs=c();  diff.mean.1.sd=c()
  for(j in 1:(length(diff)-1)){
    diff.mean.1=c(diff.mean.1, mean(diff[j:(j+1)]))  # interval 
  }


B=1000
boot.diff.mat.1=c()
for(b in 1:B){
  boot.dat = dat[sample(1:nrow(dat), size=nrow(dat)) , ]
  
  boot.dat.1 <- boot.dat[1:nrow(dat.1),];
  boot.dat.2 <- boot.dat[(nrow(dat.1)+1):nrow(dat),];
  
  ##############################################
    ## penalized spline with robust loss function
    ################################################
    dat.1.xy.boot <- data.frame(cbind(rep(freq, n),as.vector(t(boot.dat.1))))
    dat.1.xy.boot <- dat.1.xy.boot[!is.na(dat.1.xy.boot[,2]), ]
    dat.1.xy.boot <- dat.1.xy.boot[order(dat.1.xy.boot[,1]), ]
    
    dat.2.xy.boot<- data.frame(cbind(rep(freq, n),as.vector(t(boot.dat.2))))
    dat.2.xy.boot <- dat.2.xy.boot[!is.na(dat.2.xy.boot[,2]), ]
    dat.2.xy.boot <- dat.2.xy.boot[order(dat.2.xy.boot[,1]), ]
    
    fit.huber.1.boot<- robust_smooth_spline(dat.1.xy.boot, df = 10, max_iter = 50, delta = 1.345, spar = spar.set, tol = 1e-4)
    fit.huber.2.boot <- robust_smooth_spline(dat.2.xy.boot, df = 10,max_iter = 50, delta = 1.345, spar = spar.set,  tol = 1e-4)
    
    huber.1.sp.boot <- fit.huber.1.boot$fitted.values
    huber.2.sp.boot <- fit.huber.2.boot$fitted.values
    
    approx.1.boot=approxfun(dat.1.xy.boot[,1], huber.1.sp.boot)
    approx.2.boot=approxfun(dat.2.xy.boot[,1], huber.2.sp.boot)
    
    boot.huber.1=approx.1.boot(freq)
    boot.huber.2 = approx.2.boot(freq)
        
    boot.huber.mean <- ( boot.huber.1 +  boot.huber.2)/2
    boot.diff=(boot.huber.1 - boot.huber.2)^2

    
    boot.diff.mean=c(); 
    for(jj in 1:(length(boot.diff)-1)){
      boot.diff.mean=c(boot.diff.mean, mean(boot.diff[jj:(jj+1)]))
    }
    
    boot.diff.mat.1=rbind(boot.diff.mat.1, boot.diff.mean)

}

# unadjusted p-val function for order 1
pt.pval.1=c()
for(j in 1:length(diff.mean.1)){
  pt.pval.1=c(pt.pval.1, mean(boot.diff.mat.1[,j] >= diff.mean.1[j] ))
}

###########################################
# combined
###########################################

min.p <- apply(cbind(pt.pval.0[-1], pt.pval.1), 1, min)

res.comb= pval.cor.res.21(p=length(pt.pval.1), pval=pt.pval.0[-1], pval.deriv=pt.pval.1, 
                          T0=diff.mean.0[-1], T_coeff = boot.diff.mat.0[,-1], 
                          T0.deriv=diff.mean.1, T_coeff.deriv = boot.diff.mat.1, B= 1000)

res.ind = which(res.comb$corrected.pval<0.025)

# print out domain with separable group behavior
# true separable interval: [0.34, 1]
(res.f=freq[res.ind]) 


######################################################################
######################################################################
#
#   2. Effect size heatmap
#
######################################################################
######################################################################

load("sim_example_t3_sigma3.RData")
dat.1 <- dat[1:50,]
dat.2<- dat[51:100, ]

huber.k=1.3 # robust tuning parameter for M-estimation
huber.1<-c()
for (t in 1:(ncol(dat.1))){
  tmp.s=mad(dat.1[,t],na.rm=TRUE)
  huber.1=c(huber.1, huber(dat.1[,t], k = huber.k, tol = 1e-06)$mu)
}

huber.2<-c()
for (t in 1:(ncol(dat.2))){
  tmp.s=mad(dat.2[,t],na.rm=TRUE)
  huber.2=c(huber.2, huber(dat.2[,t], k = huber.k, tol = 1e-06)$mu)
}

########
n.1 = nrow(dat.1); n.2=nrow(dat.2); n=n.1+n.2
theta.1 = huber.1; theta.2 = huber.2

theta.mean=(n.1*theta.1 + n.2*theta.2)/n

# estimate covariance matrix of functional M-estimate using bootstrp samples
B=1000
huber.1.mat <- huber.2.mat <- c();

for (i in 1:B){
  
  # bootstrap samples
  boot.1=dat.1[sample(1:nrow(dat.1),nrow(dat.1),replace=TRUE), ]
  boot.2=dat.2[sample(1:nrow(dat.2),nrow(dat.2),replace=TRUE), ]
  
  id.1=which( apply(boot.1, 2, mad, na.rm=TRUE)>0)
  id.2=which(apply(boot.2, 2, mad, na.rm=TRUE)>0)
  
  huber.1.boot=c()
  for (te in id.1){
    huber.1.boot[te]=huber(boot.1[,te], k = huber.k, tol = 1e-06)$mu
  }
  
  huber.2.boot=c()
  for (tl in id.2){
    huber.2.boot[tl]=huber(boot.2[,tl], k = huber.k, tol = 1e-06)$mu
  }
  
  huber.1.mat=rbind(huber.1.mat,huber.1.boot)
  huber.2.mat=rbind(huber.2.mat,huber.2.boot)
  
  #if(i %%300 == 0) print(i)
}

##############################################
# SSR
SSR=(n.2*(theta.2-theta.mean)^2 + n.1*(theta.1-theta.mean)^2)/n

#############################################
# SSE
boot.mean.1=colMeans(huber.1.mat, na.rm=TRUE)
boot.mean.2=colMeans(huber.2.mat, na.rm=TRUE)

mean.mat.1=t(matrix(boot.mean.1, ncol=nrow(huber.1.mat), nrow=length(freq)))
mean.mat.2=t(matrix(boot.mean.2, ncol=nrow(huber.2.mat), nrow=length(freq)))

# cov matrix of functional M-estimator
tmp.1=rbind((huber.1.mat-mean.mat.1),(huber.2.mat-mean.mat.2))
tmp.1[is.na(tmp.1)] <- mean(tmp.1, na.rm = TRUE)
cov.1=t(tmp.1)%*%tmp.1/nrow(tmp.1)

# pointwise effect size
fSNR=(SSR)/(diag(cov.1))

# robust effect size matrix - each row corresponds to effect size from each scale
size.mat=matrix(nrow=length(freq), ncol=length(freq))
size.mat[1,]<- sqrt(fSNR)

for(pp in 1:(length(freq)-1)){
  for(k in 1:length(freq)){
    tmp.id = c((k-pp):(k+pp))
    inc.id = tmp.id[tmp.id >= 1 & tmp.id <= length(freq)]
    size.mat[ (pp+1),k] =sqrt(mean(fSNR[inc.id]))
  }
  #if(pp%%30==0) print(pp)
}


#################################
# generate effect size heatmap
#################################

nlevel=10
col_heatmap = rev(hcl.colors(20, "terrain"))

dev.new()
layout(matrix(c(1,2),nrow=1), widths = c(9, 1))
par(mar = c(5, 5, 2.5, 2))
image(freq, freq, t(size.mat), col=col_heatmap,
      ylab = expression(paste("Interval length ", "(",Delta,")")), xlab="Frequency (MHz)", main="Effect size heatmap", cex.main=1.8, cex.2ab=1.5)
box()
par(mar = c(4, 0.1, 4, 3))
image(1, seq.int(nlevel) - (nlevel + 1)/2, t(as.matrix(seq.int(nlevel))), 
      axes = FALSE, xlab = "", ylab = "", col = col_heatmap)
axis(side = 4, at = seq(0, nlevel, length.out = 6) - 
       nlevel/2, labels = format(round(seq(0, max(size.mat), length.out = 6)), 2))
box()
