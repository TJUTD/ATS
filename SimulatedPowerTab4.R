library(MASS)
library(magic)

###########################################################
# utility fuction
# trace of matrix
tr <- function (m) {
  if (!is.matrix(m) | (dim(m)[1] != dim(m)[2])) 
    stop("m must be a square matrix")
  return(sum(diag(m), na.rm = TRUE))
}

###########################################################
# Table 4. Simulated Power of the Rank Statistic Q_N(C) and F_N(M)
nppwrcmp <- function(a, b, n, delta, tau, reps = 5000) {

  # heteroscedasticity
  mu <- kronecker(matrix(delta,1,a), matrix(1,1,b)) + kronecker(matrix(1,1,a),matrix(tau,1,b))
  d = a * b
  NT <- sum(n)
  rid <- c(0,cumsum(n))
  de <- c(0,0.1,0.2,0.25,0.3,0.4,0.5)
  nd <- length(de)
  
  # projection matrix P = I - Jn/n
  Ja <- rep(1, a) %*% t(rep(1, a))
  Jb <- rep(1, b) %*% t(rep(1, b))
  Pa <- diag(rep(1, a)) - Ja/a
  Pb <- diag(rep(1, b)) - Jb/b
  MaW <- kronecker(Pa, t(rep(1/b, b)))
  MbW <- kronecker(t(rep(1/a, a)), Pb)
  Ma <- kronecker(Pa, Jb/b)
  Mb <- kronecker(Ja/a, Pb)
  Mab <- kronecker(Pa, Pb)  
  Da <- diag(diag(Ma))
  Db <- diag(diag(Mb))
  Dab <- diag(diag(Mab))
  Lambda <- diag(1/(n-1))

  # data and moments generation 
  dat <- matrix(0,NT,reps)
  dev <- matrix(0,NT,nd)
  for (i in 1:d) {
    dat[(rid[i]+1):rid[i+1],] <- matrix(rnorm(n[i]*reps, mu[i], 1), n[i], reps)
    if ((i%%a) == 1) {
      dev[(rid[i]+1):rid[i+1],] <- kronecker(matrix(de,1,nd),matrix(1,n[i],1))
    }
  }
  repsnd <- reps*nd
  rdat <- matrix(0,NT,repsnd)
  for (j in 1:nd) {
    for (i in 1:reps) {
      rdat[,(i+(j-1)*reps)] <- rank(dat[,i]+dev[,j])
    }
  }
  xbar <- matrix(0, d, repsnd)
  SN <- matrix(0, d, repsnd)
  for (i in 1:d) {
    M <- matrix(1/n[i],n[i],n[i])
    MX <- M %*% rdat[(rid[i]+1):rid[i+1],] 
    xbar[i,] <- (MX[1,] - 0.5) / NT
    SN[i,] <- M[1,] %*% ((rdat[(rid[i]+1):rid[i+1],]  - MX) * (rdat[(rid[i]+1):rid[i+1],]  - MX)) / NT / (n[i]-1)
  }

  
  Qa <- matrix(0, reps, nd)
  Qb <- matrix(0, reps, nd)
  Qab <- matrix(0, reps, nd)
  ATSa <- matrix(0, reps, nd)
  ATSb <- matrix(0, reps, nd)
  ATSab <- matrix(0, reps, nd)
  fa <- matrix(0, reps, nd)
  f0a <- matrix(0, reps, nd)
  fb <- matrix(0, reps, nd)
  f0b <- matrix(0, reps, nd)
  fab <- matrix(0, reps, nd)
  f0ab <- matrix(0, reps, nd)
  
  for (i in 1:reps) {
    for (j in 1:nd) {
      x <- xbar[,(i+(j-1)*reps)]
      S <- diag(SN[,(i+(j-1)*reps)])
    
      # Wald statistic
      Qa[i,j] <- NT * t(x) %*% t(MaW) %*% ginv(MaW %*% S %*% t(MaW)) %*% MaW %*% x
      Qb[i,j] <- NT * t(x) %*% t(MbW) %*% ginv(MbW %*% S %*% t(MbW)) %*% MbW %*% x
      Qab[i,j] <- NT * t(x) %*% t(Mab) %*% ginv(Mab %*% S %*% t(Mab)) %*% Mab %*% x
    
      # FN(M) = N / tr(DM*SN) * X'*M*X ~ F(f,f0), f = tr(DM*SN)^2 / tr(M*SN*M*SN) f0 = tr(DM*SN)^2 / tr(DM^2*SN^2*Lambda)    
      ATSa[i,j] <- NT / tr(Da %*% S) * t(x) %*% Ma %*% x 
      ATSb[i,j] <- NT / tr(Db %*% S) * t(x) %*% Mb %*% x 
      ATSab[i,j] <- NT / tr(Dab %*% S) * t(x) %*% Mab %*% x 
      fa[i,j] <- tr(Da %*% S)^2 / tr(Ma %*% S %*% Ma %*% S) 
      fb[i,j] <- tr(Db %*% S)^2 / tr(Mb %*% S %*% Mb %*% S) 
      fab[i,j] <- tr(Dab %*% S)^2 / tr(Mab %*% S %*% Mab %*% S) 
      f0a[i,j] <- tr(Da %*% S)^2 / tr(Da %*% Da %*% S %*% S %*% Lambda)    
      f0b[i,j] <- tr(Db %*% S)^2 / tr(Db %*% Db %*% S %*% S %*% Lambda)    
      f0ab[i,j] <- tr(Dab %*% S)^2 / tr(Dab %*% Dab %*% S %*% S %*% Lambda)    
    }
  }
  
  # chi(f)
  qQa <- qchisq(.95, df = a - 1)
  qQb <- qchisq(.95, df = b - 1)
  qQab <- qchisq(.95, df = (a - 1) * (b-1))
  pwrQ <- matrix(0,nd,3)
  for (j in 1:nd) {
      pwrQ[j,1] <- mean(Qa[,j] >= qQa)
      pwrQ[j,2] <- mean(Qb[,j] >= qQb)
      pwrQ[j,3] <- mean(Qab[,j] >= qQab)
  }
  
  # F(f, f0)
  qFa <- matrix(0, reps, nd)
  qFb <- matrix(0, reps, nd)
  qFab <- matrix(0, reps, nd)
  for (i in 1:reps) {
    for (j in 1:nd) {
    qFa[i,j] <- qf(.95, df1 = fa[i,j], df2 = f0a[i,j])
    qFb[i,j] <- qf(.95, df1 = fb[i,j], df2 = f0b[i,j])
    qFab[i,j] <- qf(.95, df1 = fab[i,j], df2 = f0ab[i,j])
    }
  }
  pwrF <- matrix(0,nd,3)
  for (j in 1:nd) {
    pwrF[j,1] <- mean(ATSa[,j] >= qFa[,j])
    pwrF[j,2] <- mean(ATSb[,j] >= qFb[,j])
    pwrF[j,3] <- mean(ATSab[,j] >= qFab[,j])
  }
  
  comparisontable <- cbind(de, pwrQ[,1], pwrF[,1], pwrQ[,2], pwrF[,2], pwrQ[,3], pwrF[,3])
  colnames(comparisontable) <- c("delta", "Qa", "ATSa", "Qb", "ATSb", "Qab", "ATSab")
  return(comparisontable)
}

# two-way ANOVA y_ijk = alpha_i + beta_j + (alphabeta)_ij + e_ijk, alpha_i = 0, beta_j = 0, (alphabeta)_ij = 0
reps <- 5000
a <- 3; b <- 2
n <- rep(50, 6)
set.seed(201711)
delta <- c(0,0,0); tau <- c(0,0)
table4.1 <- nppwrcmp(a,b,n,delta,tau,reps)
delta <- c(0,0,0); tau <- c(0,1)
table4.2 <- nppwrcmp(a,b,n,delta,tau,reps)

table4.1
table4.2





