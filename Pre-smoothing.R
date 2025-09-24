
load("sim_example_t3_sigma3.RData")
freq <- seq(0, 1, length = 97) 

######################################
#### fully observed data
######################################
## l=0 (pre-smoothing)
h.cv0 <- 1 # find the optimal
dat.ker0 <- matrix(nrow=nrow(dat),ncol=ncol(dat))
for (i in 1:nrow(dat)) {
  tmp <- locpoly(x = freq[which(delta[i,]==1)], y = par.dat[i,which(delta[i,]==1)], 
                 kernel = 'epanech', bandwidth = h.cv0, degree = 0, 
                 range.x = range(freq), 
                 gridsize = length(freq))$y
  dat.ker0[i, ] <- tmp
}

## l=1 (pre-smoothing)
h.cv1 <- 1.5 # find the optimal
dat.ker1 <- matrix(nrow=nrow(dat),ncol=ncol(dat))
for (i in 1:nrow(dat)) {
  tmp <- locpoly(x = freq[which(delta[i,]==1)], y = par.dat[i,which(delta[i,]==1)], 
                 kernel = 'epanech', bandwidth = h.cv1, degree = 1, drv = 1, 
                 range.x = range(freq), 
                 gridsize = length(freq))$y
  dat.ker1[i, ] <- tmp
}

######################################
#### partially observed data
######################################

## l=0 (pre-smoothing)
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

