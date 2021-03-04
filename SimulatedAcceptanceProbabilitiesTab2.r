library(MASS)
library(magic)

###########################################################
# trace from "psych"
tr <- function (m) {
  if (!is.matrix(m) | (dim(m)[1] != dim(m)[2])) 
    stop("m must be a square matrix")
  return(sum(diag(m), na.rm = TRUE))
}

###########################################################
# Table 2. Simulated Acceptance Probabilities for the Wald Statistic Q_N(C) and the ANOVA-Type Statistics F_N(M)
BTA <- function(a,b,n,reps = 5000) {

  # heteroscedasticity
  sigma <- (1 + kronecker(1:a, 1:b)/2)^2
  sd <- sqrt(sigma)
  d = a * b
  NT <- sum(n)
  
  SS <- NT * diag(sigma / n)

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
  xbar <- matrix(0, d, reps)
  SN <- matrix(0, d, reps)
  #  s0 <- matrix(0, d, reps)
  for (i in 1:d) {
      dat <- matrix(rnorm(n[i]*reps, 0, sd[i]), n[i], reps)
      M <- matrix(1/n[i],n[i],n[i])
      MX <- M %*% dat
      xbar[i,] <- MX[1,] 
      SN[i,] <- NT / (n[i]-1) * M[1,] %*% ((dat - MX) * (dat - MX))
      #      s0[i,] <- matrix(1/(n[i]-1),1,n[i]) %*% ((dat - MX) * (dat - MX)) 
  }

  results <- matrix(0, reps, 18)
  for (i in 1:reps) {
    x <- xbar[,i]
    S <- diag(SN[,i])
    
    # Wald statistic
    results[i,1] <- NT * t(x) %*% t(MaW) %*% ginv(MaW %*% S %*% t(MaW)) %*% MaW %*% x
    results[i,2] <- NT * t(x) %*% t(MbW) %*% ginv(MbW %*% S %*% t(MbW)) %*% MbW %*% x
    results[i,3] <- NT * t(x) %*% t(Mab) %*% ginv(Mab %*% S %*% t(Mab)) %*% Mab %*% x
  
    # FN(M) = N / tr(DM*SN) * X'*M*X ~ F(f,f0), f = tr(DM*SN)^2 / tr(M*SN*M*SN) f0 = tr(DM*SN)^2 / tr(DM^2*SN^2*Lambda)    
    results[i,4] <- NT / tr(Da %*% SS) * t(x) %*% Ma %*% x 
    results[i,5] <- NT / tr(Db %*% SS) * t(x) %*% Mb %*% x 
    results[i,6] <- NT / tr(Dab %*% SS) * t(x) %*% Mab %*% x 
    results[i,7] <- NT / tr(Da %*% S) * t(x) %*% Ma %*% x 
    results[i,8] <- NT / tr(Db %*% S) * t(x) %*% Mb %*% x 
    results[i,9] <- NT / tr(Dab %*% S) * t(x) %*% Mab %*% x 
    results[i,10] <- tr(Da %*% SS)^2 / tr(Ma %*% SS %*% Ma %*% SS) 
    results[i,11] <- tr(Db %*% SS)^2 / tr(Mb %*% SS %*% Mb %*% SS) 
    results[i,12] <- tr(Dab %*% SS)^2 / tr(Mab %*% SS %*% Mab %*% SS) 
    results[i,13] <- tr(Da %*% S)^2 / tr(Ma %*% S %*% Ma %*% S) 
    results[i,14] <- tr(Db %*% S)^2 / tr(Mb %*% S %*% Mb %*% S) 
    results[i,15] <- tr(Dab %*% S)^2 / tr(Mab %*% S %*% Mab %*% S) 
    results[i,16] <- tr(Da %*% S)^2 / tr(Da %*% Da %*% S %*% S %*% Lambda)    
    results[i,17] <- tr(Db %*% S)^2 / tr(Db %*% Db %*% S %*% S %*% Lambda)    
    results[i,18] <- tr(Dab %*% S)^2 / tr(Dab %*% Dab %*% S %*% S %*% Lambda)    
    
  }

  colnames(results) <- c("Qa", "Qb", "Qab", "ATSaS", "ATSbS", "ATSabS", "ATSa", "ATSb", "ATSab", "faS","fbS", "fabS", "fa", "fb", "fab", "f0a", "f0b", "f0ab")

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

  # chi(f) / f
  qATSaS <- matrix(0, reps, 3)
  qATSbS <- matrix(0, reps, 3)
  qATSabS <- matrix(0, reps, 3)
  for (i in 1:reps) {
    qATSaS[i,] <- qchisq(c(.9,.95,.99), df = results[i, "faS"]) / results[i, "faS"]
    qATSbS[i,] <- qchisq(c(.9,.95,.99), df = results[i, "fbS"]) / results[i, "fbS"]
    qATSabS[i,] <- qchisq(c(.9,.95,.99), df = results[i, "fabS"]) / results[i, "fabS"]
  }
  pATSS <- matrix(0, 3 ,3)
  # colnames(pATSS) <- c("ATSaS", "ATSbS", "ATSabS")
  for (i in 1:3) {
    pATSS[i,1] <- mean(results[, "ATSaS"] <= qATSaS[,i])
    pATSS[i,2] <- mean(results[, "ATSbS"] <= qATSbS[,i])
    pATSS[i,3] <- mean(results[, "ATSabS"] <= qATSabS[,i])
  }
  
  qATSaSn <- matrix(0, reps, 3)
  qATSbSn <- matrix(0, reps, 3)
  qATSabSn <- matrix(0, reps, 3)
  for (i in 1:reps) {
    qATSaSn[i,] <- qchisq(c(.9,.95,.99), df = results[i, "fa"]) / results[i, "fa"]
    qATSbSn[i,] <- qchisq(c(.9,.95,.99), df = results[i, "fb"]) / results[i, "fb"]
    qATSabSn[i,] <- qchisq(c(.9,.95,.99), df = results[i, "fab"]) / results[i, "fab"]
  }
  pATSSn <- matrix(0, 3 ,3)
  # colnames(pATSS) <- c("ATSaS", "ATSbS", "ATSabS")
  for (i in 1:3) {
    pATSSn[i,1] <- mean(results[, "ATSa"] <= qATSaSn[,i])
    pATSSn[i,2] <- mean(results[, "ATSb"] <= qATSbSn[,i])
    pATSSn[i,3] <- mean(results[, "ATSab"] <= qATSabSn[,i])
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
  
  comparisontable <- cbind(pQ[,1], pATSS[,1], pATSSn[,1], pATS[,1], pQ[,2], pATSS[,2], pATSSn[,2], pATS[,2], pQ[,3], pATSS[,3], pATSSn[,3], pATS[,3])
  colnames(comparisontable) <- c("Qa", "ATSaS", "ATSaSn", "ATSa", "Qb", "ATSbS", "ATSbSn", "ATSb", "Qab", "ATSSab", "ATSabSn", "ATSab")
  return(comparisontable)
}

# two-way ANOVA y_ijk = alpha_i + beta_j + (alphabeta)_ij + e_ijk, alpha_i = 0, beta_j = 0, (alphabeta)_ij = 0
reps <- 5000
a = 2; b = 5
set.seed(201710)
# case 1
n1 <- rep(7,10)
table2.1 <- BTA(a,b,n1,reps)

# case 2 
n2 <- rep(c(11, 10, 9, 8, 7), 2)
table2.2 <- BTA(a,b,n2,reps)

table2.1
table2.2





