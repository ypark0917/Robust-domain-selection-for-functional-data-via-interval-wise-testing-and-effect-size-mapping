
# setwd(...)

library(KernSmooth)
library(reshape2)
library(mvtnorm)
library(MASS)

######################################
# mu_1(t) and mu_2(t) in Section 3 
load("group_mean_functions.RData")
freq <- seq(0, 1, length = 97) 

# illustration of mu_1(t) and mu_2(t) displaying distinct behaviors over [0.34, 1]
plot(freq, gr1.mean, type="l")
lines(freq, gr2.mean, col="blue")
abline(v=0.34, lty=2, col="red")

######################################
# setting
type=1 # type=1: Gaussian process/ type=2: t(3) process/ type=3: curve outlier/ type=4: local outlier
n=50   # group size
sigma.e=2  # noise level for Gaussian process or t(3) process
partial=2  # 1=full sampling; 2=partial sampling
con.prob = 0.05  # contamination rate for type 3 and type 4
con.size = 3    # outlying parameter to generate type 3 and type 4. For type 3, we set the con.size = 3; for type 4, we set the con.size = 4. The larger the severe outlyingness
rep=100  # number of repetitions


# cov function 
f.range <- 20
cor.mat<-matrix(nrow=length(freq), ncol=length(freq))
for(i in 1: length(freq)){
  for(j in 1:length(freq)){
    cor.mat[i,j] <- exp(-abs(i-j)/f.range)
  }
}

# generate simulation data
## par.dat : raw trajectories where missing are replaced by NA
## dat.ker : kernel smoothing of par.dat over observed domain

for(ii in 1:rep){

  set.seed(12*ii) 
  
  ###########################################
  # type of outliers
  ###########################################
 
  if(type==1){  # Gaussian error process
    error.gr1 <- (rmvnorm(n, sigma=sigma.e^2*cor.mat))
    error.gr2 <- (rmvnorm(n, sigma=sigma.e^2*cor.mat))
    
    gr1 <- matrix(gr1.mean, ncol=ncol(error.gr1), nrow=nrow(error.gr1), byrow = T) + error.gr1
    gr2 <- matrix(gr2.mean, ncol=ncol(error.gr2), nrow=nrow(error.gr2), byrow = T) + error.gr2
    
  }else if(type==2){  # t(3) error process
    error.gr1 <- (rmvt(n, sigma=sigma.e^2*cor.mat, df=3))
    error.gr2 <- (rmvt(n, sigma=sigma.e^2*cor.mat, df=3))
    
    gr1 <- matrix(gr1.mean, ncol=ncol(error.gr1), nrow=nrow(error.gr1), byrow = T) + error.gr1
    gr2 <- matrix(gr2.mean, ncol=ncol(error.gr2), nrow=nrow(error.gr2), byrow = T) + error.gr2
    
  }else if(type==3){  # curve outlier
    error.gr1 <- (rmvnorm(n, sigma=sigma.e^2*cor.mat))
    error.gr2 <- (rmvnorm(n, sigma=sigma.e^2*cor.mat))
    
    set.seed(12*ii + 3)
    tmp1 <- matrix(rep(rbinom(n, 1, prob= con.prob), each=length(freq)), nrow=n, byrow = TRUE)
    tmp2 <- matrix(rep(rbinom(n, 1, prob= con.prob), each=length(freq)), nrow=n, byrow = TRUE)
    
    error.con.gr1 <- tmp1 * matrix(rep(con.size*sigma.e*rt(n, df=3), each=length(freq)), nrow=n, byrow = TRUE)
    error.con.gr2 <- tmp2 * matrix(rep(con.size*sigma.e*rt(n, df=3), each=length(freq)), nrow=n, byrow = TRUE)
    
    gr1 <- matrix(gr1.mean, ncol=ncol(error.gr1), nrow=nrow(error.gr1), byrow = T) + error.gr1 + error.con.gr1
    gr2 <- matrix(gr2.mean, ncol=ncol(error.gr2), nrow=nrow(error.gr2), byrow = T) + error.gr2 + error.con.gr2
    
  }else if(type==4){  # local outlier
    error.gr1 <- (rmvnorm(n, sigma=sigma.e^2*cor.mat))
    error.gr2 <- (rmvnorm(n, sigma=sigma.e^2*cor.mat))
    
    set.seed(12*ii + 1)
    tmp1 <- matrix(rbinom(n*length(freq), 1, prob= con.prob), nrow=n)
    tmp2 <- matrix(rbinom(n*length(freq), 1, prob= con.prob), nrow=n)
    
    error.con.gr1 <- tmp1 * matrix(con.size*sigma.e*rt(n*length(freq), df=3), nrow=n)
    error.con.gr2 <- tmp2 * matrix(con.size*sigma.e*rt(n*length(freq), df=3), nrow=n)
    
    gr1 <- matrix(gr1.mean, ncol=ncol(error.gr1), nrow=nrow(error.gr1), byrow = T) + error.gr1 + error.con.gr1
    gr2 <- matrix(gr2.mean, ncol=ncol(error.gr2), nrow=nrow(error.gr2), byrow = T) + error.gr2 + error.con.gr2
    
  }
  
  dat0 <- rbind(gr1, gr2)
  
  ################################
  # illustration of original data
  
  plot(freq, gr1.mean, type="n", ylim=c(-60, -10))
  for(i in 1:n){
    lines(freq, dat0[i,]); lines(freq, dat0[50+i,], col="blue")
  }
  
  ############################################################
  # sampling framework (fully observed of partially observed)
  ############################################################
  delta=matrix(1, nr=nrow(dat0), nc=ncol(dat0))
  
  if(partial==1){  # full sampling structure
    par.dat0=dat0*delta  
    par.dat = par.dat0
    par.dat[par.dat==0]=NA       
    
    # smoothing
    h.cv <- 0.01  # kernel smoothing with the given bandwidth over the entire domain
    dat.ker <- matrix(nrow=nrow(dat0),ncol=ncol(dat0))
    for (i in 1:nrow(dat0)) {
      tmp <- locpoly(x = freq[which(delta[i,]==1)], y = par.dat[i,which(delta[i,]==1)], 
                     kernel = 'epanech', bandwidth = h.cv, degree = 1, 
                     range.x = range(freq), 
                     gridsize = length(freq))$y
      dat.ker[i, ] <- tmp
    }
     
  }else if(partial==2){ # partial sampling structure
    m.prob=0.4; dd=1.2; f=0.3
    if(m.prob > 0){  
      
      for(jj in 1:nrow(delta)){
        if(rbinom(n=1, size=1, prob=m.prob)==1){
          
          set.seed(45+jj)
          u1=(runif(1,1,length(freq))); u2=runif(1,1,length(freq))
          Ci=round(dd*u1); Ei =round(f*u2)
          if((Ci-Ei)<1){
            l.b=1
          }else if((Ci-Ei)> length(freq)){
            l.b=length(freq)
          } else {
            l.b=Ci-Ei 
          }

          if((Ci+Ei)>= length(freq)){
            u.b=length(freq)
          }else if((Ci+Ei)<1){
            u.b=1
          } else {u.b=Ci+Ei }
          
          if(l.b==u.b & l.b==length(freq)){
            delta[jj,] = 1
          }else if(l.b==u.b & u.b==1){
            delta[jj,] = 1
          }else delta[jj, c(l.b:u.b)] = 0
          # }  
          
        }
      }
    } 
    par.dat0=dat0*delta          # par.dat0: missing is replaced by 0
    par.dat = par.dat0
    par.dat[par.dat==0]=NA       # par.dat: missing is prelaced by NA
    
    # smoothing
    h.cv <- 0.01  # kernel smoothing with the given bandwidth over individual-specific domain
    dat.ker <- matrix(nrow=nrow(dat0),ncol=ncol(dat0))
    for (i in 1:nrow(dat0)) {
      tmp <- locpoly(x = freq[which(delta[i,]==1)], y = par.dat[i,which(delta[i,]==1)], 
                     kernel = 'epanech', bandwidth = h.cv, degree = 1, 
                     range.x = range(freq), 
                     gridsize = length(freq))$y
      dat.ker[i, which(delta[i,]==1)] <- tmp[which(delta[i,]==1)]
    }
    
  }
}

