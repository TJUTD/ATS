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
# Table 3. Simulated Acceptance Probabilities for the Rank Statistic Q_N(C) and F_N(M)
npBoxType <- function(a, b, n, delta, tau, reps = 5000) {

  
  d = a * b
  NT <- sum(n)

  # data and moments generation 
  # heteroscedasticity
  mu <- kronecker(matrix(delta,1,a), matrix(1,1,b)) + kronecker(matrix(1,1,a),matrix(tau,1,b))
  rid <- c(0,cumsum(n))
  dat <- matrix(0,NT,reps)
  for (i in 1:d) {
    dat[(rid[i]+1):rid[i+1],] <- matrix(rnorm(n[i]*reps, mu[i], 1), n[i], reps)
  }

  # transfer to rank
  for (i in 1:reps) {
    dat[,i] <- rank(dat[,i])
  }
  xbar <- matrix(0, d, reps)
  SN <- matrix(0, d, reps)
  for (i in 1:d) {
    M <- matrix(1/n[i],n[i],n[i])
    MX <- M %*% dat[(rid[i]+1):rid[i+1],] 
    xbar[i,] <- (MX[1,] - 0.5) / NT
    SN[i,] <- M[1,] %*% ((dat[(rid[i]+1):rid[i+1],]  - MX) * (dat[(rid[i]+1):rid[i+1],]  - MX)) / NT / (n[i]-1)
  }

  
  
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
  
  results <- matrix(0, reps, 12)
  for (i in 1:reps) {
    x <- xbar[,i]
    S <- diag(SN[,i])
    
    # Wald statistic
    Qa <- NT * t(x) %*% t(MaW) %*% ginv(MaW %*% S %*% t(MaW)) %*% MaW %*% x
    Qb <- NT * t(x) %*% t(MbW) %*% ginv(MbW %*% S %*% t(MbW)) %*% MbW %*% x
    Qab <- NT * t(x) %*% t(Mab) %*% ginv(Mab %*% S %*% t(Mab)) %*% Mab %*% x
    
    # FN(M) = N / tr(DM*SN) * X'*M*X ~ F(f,f0), f = tr(DM*SN)^2 / tr(M*SN*M*SN) f0 = tr(DM*SN)^2 / tr(DM^2*SN^2*Lambda)    
    ATSa <- NT / tr(Da %*% S) * t(x) %*% Ma %*% x 
    ATSb <- NT / tr(Db %*% S) * t(x) %*% Mb %*% x 
    ATSab <- NT / tr(Dab %*% S) * t(x) %*% Mab %*% x 
    fa <- tr(Da %*% S)^2 / tr(Ma %*% S %*% Ma %*% S) 
    fb <- tr(Db %*% S)^2 / tr(Mb %*% S %*% Mb %*% S) 
    fab <- tr(Dab %*% S)^2 / tr(Mab %*% S %*% Mab %*% S) 
    f0a <- tr(Da %*% S)^2 / tr(Da %*% Da %*% S %*% S %*% Lambda)    
    f0b <- tr(Db %*% S)^2 / tr(Db %*% Db %*% S %*% S %*% Lambda)    
    f0ab <- tr(Dab %*% S)^2 / tr(Dab %*% Dab %*% S %*% S %*% Lambda)    
    
    results[i,] <- c(Qa, Qb, Qab, ATSa, ATSb, ATSab, fa, fb, fab, f0a, f0b, f0ab)
  }
  
  colnames(results) <- c("Qa", "Qb", "Qab", "ATSa", "ATSb", "ATSab", "fa", "fb", "fab", "f0a", "f0b", "f0ab")
  
  # chi(f)
  qQa <- qchisq(c(.9,.95,.99), df = a - 1)
  qQb <- qchisq(c(.9,.95,.99), df = b - 1)
  qQab <- qchisq(c(.9,.95,.99), df = (a - 1) * (b-1))
  pQ <- matrix(0,3,3)
  # colnames(pQ) <- c("Qa","Qb", "Qab")
  for (i in 1:3) {
    pQ[i,1] <- sum(results[, "Qa"] <= qQa[i])/reps
    pQ[i,2] <- sum(results[, "Qb"] <= qQb[i])/reps
    pQ[i,3] <- sum(results[, "Qab"] <= qQab[i])/reps
  }
  
  # F(f, f0)
  qFa <- matrix(0, reps, 3)
  qFb <- matrix(0, reps, 3)
  qFab <- matrix(0, reps, 3)
  for (i in 1:reps) {
    qFa[i,] <- qf(c(.9,.95,.99), df1 = results[i, "fa"], df2 = results[i, "f0a"])
    qFb[i,] <- qf(c(.9,.95,.99), df1 = results[i, "fb"], df2 = results[i, "f0b"])
    qFab[i,] <- qf(c(.9,.95,.99), df1 = results[i, "fab"], df2 = results[i, "f0ab"])
  }
  pATS <- matrix(0,3,3)
  # colnames(pATS) <- c("ATSa", "ATSb", "ATSab")
  for (i in 1:3) {
    pATS[i,1] <- mean(results[, "ATSa"] <= qFa[,i])
    pATS[i,2] <- mean(results[, "ATSb"] <= qFb[,i])
    pATS[i,3] <- mean(results[, "ATSab"] <= qFab[,i])
  }
  
  comparisontable <- cbind(pQ[,1], pATS[,1], pQ[,2], pATS[,2], pQ[,3], pATS[,3])
  colnames(comparisontable) <- c("Qa", "ATSa", "Qb", "ATSb", "Qab", "ATSab")
  return(comparisontable)
}

# two-way ANOVA y_ijk = alpha_i + beta_j + (alphabeta)_ij + e_ijk, alpha_i = 0, beta_j = 0, (alphabeta)_ij = 0
reps <- 5000
a <- 3; b <- 6
n <- rep(7,18)
set.seed(201711)
delta <- c(0,0,0); tau <- c(0,0,0,0,0,0)
table3.1 <- npBoxType(a,b,n,delta,tau,reps)
delta <- c(0,0,0); tau <- c(0,1,1,1,1,1)
table3.2 <- npBoxType(a,b,n,delta,tau,reps)
table3.1
table3.2






