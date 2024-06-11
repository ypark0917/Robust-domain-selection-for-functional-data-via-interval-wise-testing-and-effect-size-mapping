
###########################################################
##  A part of source code borrowed from
## "Interval-wise testing for functional data." by Pini and Vantini (2017). Journal of Nonparametric Statistics 29(2), 407–424.

pval.correct <- function(pval.matrix){
  matrice_pval_2_2x <- cbind(pval.matrix,pval.matrix)
  p <- dim(pval.matrix)[2]
  matrice_pval_2_2x <- matrice_pval_2_2x[,(2*p):1]
  corrected.pval <- numeric(p)
  for(var in 1:p){
    pval_var <- matrice_pval_2_2x[p,var]
    inizio <- var
    fine <- var #inizio fisso, fine aumenta salendo nelle righe
    for(riga in (p-1):1){
      fine <- fine + 1
      pval_cono <- matrice_pval_2_2x[riga,inizio:fine]
      pval_var <- max(pval_var,pval_cono)
    }
    corrected.pval[var] <- pval_var
  }
  corrected.pval <- corrected.pval[p:1]
  return(corrected.pval)
}


##########################################################
## source function to calculate \tilde{p}_m(t) under L=1
#########################################################
pval.cor.res.L1 <- function(p, pval, pval.deriv, T0, T_coeff, T0.deriv, T_coeff.deriv, B=5000){
  
  # p: the number of grids
  # pval: calculated unadjusted p-val of test statistic at each grid from the raw data (l=0)
  # pval.deriv: calculated unadjusted p-val of test statistic at each grid from the 1st derivated data (l=1)
  # T0: test-statistic at each grid from the raw data (l=0)
  # T0.deriv: test-statistic at each grid from the 1st derivated data (l=1)
  # T_coeff: matrix(ncol=p,nrow=B) -  calculated test statistic based on bootstrap sampled (l=0)
  # T_coeff.deriv: matrix(ncol=p,nrow=B) -  calculated test statistic based on bootstrap sampled (l=1)
  
  #asymmetric combination matrix:
  matrice_pval_asymm <- matrix(nrow=p,ncol=p)
  matrice_pval_asymm[p,] <- apply(cbind(pval[1:p], pval.deriv[1:p]),1,min)
  T0_2x <- c(T0,T0)
  T_coeff_2x <- as.matrix(cbind(T_coeff,T_coeff))

  T0_2x.deriv <- c(T0.deriv,T0.deriv)
  T_coeff_2x.deriv <- as.matrix(cbind(T_coeff.deriv,T_coeff.deriv))
  
  ####################################################################

  for(i in (p-1):1){
    for(j in 1:p){
      inf <- j
      sup <- (p-i)+j
      T0_temp <- sum(T0_2x[inf:sup])
      T_temp <- rowSums(T_coeff_2x[,inf:sup])
      
      T0_temp.deriv <- sum(T0_2x.deriv[inf:sup])
      T_temp.deriv <- rowSums(T_coeff_2x.deriv[,inf:sup])
      
      pval_temp <- min(sum(T_temp>=T0_temp)/B, sum(T_temp.deriv >=T0_temp.deriv)/B)
      matrice_pval_asymm[i,j] <- pval_temp
      #print(c(inf,sup))
    }
    print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
  }
  
  ###################################
  # versione in Rcpp
  #matrice_pval_asymm <- combineMatrix(L_2x, pval_2x, matrice_pval_asymm, p,B)
  
  #symmetric combination matrix
  matrice_pval_symm <- matrix(nrow=p,ncol=4*p)
  for(i in 0:(p-1)){
    for(j in 1:(2*p)){
      matrice_pval_symm[p-i,j+i+p] <- matrice_pval_asymm[p-i,(j+1)%/%2]
      if(j+i>2*p-i){
        matrice_pval_symm[p-i,j+i-p] <- matrice_pval_asymm[p-i,(j+1)%/%2]
      }
    }
  }
  corrected.pval <- pval.correct(matrice_pval_asymm)
  
  ITPresult <- list(pval=pval,
                    pval.matrix=matrice_pval_asymm,
                    corrected.pval=corrected.pval,
                    heatmap.matrix=matrice_pval_symm)
  #class(ITPresult) = 'ITPpointwise'
  return(ITPresult)
}

