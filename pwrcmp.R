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
BoxType <- function(a, b, n, reps = 5000, numpt = 31) {
  
  # heteroscedasticity
  sigma <- (1 + kronecker(1:a, 1:b)/2)^2
  sd <- sqrt(sigma)
  d = a * b
  NT <- sum(n)
  
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
  for (i in 1:d) {
    dat <- matrix(rnorm(n[i]*reps, 0, sd[i]), n[i], reps)
    M <- matrix(1/n[i],n[i],n[i])
    MX <- M %*% dat
    xbar[i,] <- MX[1,] 
    SN[i,] <- NT / (n[i]-1) * M[1,] %*% ((dat - MX) * (dat - MX))
  }

 
  
  #delta <- seq(-20,20,length = numpt)
  delta <- seq(0,20,length = numpt)
  Qa <- matrix(0, reps, numpt)
  Qb <- matrix(0, reps, numpt)
  Qab <- matrix(0, reps, numpt)
  ATSa <- matrix(0, reps, numpt)
  ATSb <- matrix(0, reps, numpt)
  ATSab <- matrix(0, reps, numpt)
  f <- matrix(0, reps, 3)
  f0 <- matrix(0, reps, 3)
  
  for (i in 1:reps) {
    
    S <- diag(SN[,i])
    f[i,1] <- tr(Da %*% S)^2 / tr(Ma %*% S %*% Ma %*% S) 
    f[i,2] <- tr(Db %*% S)^2 / tr(Mb %*% S %*% Mb %*% S) 
    f[i,3] <- tr(Dab %*% S)^2 / tr(Mab %*% S %*% Mab %*% S) 
    f0[i,1] <- tr(Da %*% S)^2 / tr(Da %*% Da %*% S %*% S %*% Lambda)    
    f0[i,2] <- tr(Db %*% S)^2 / tr(Db %*% Db %*% S %*% S %*% Lambda)    
    f0[i,3] <- tr(Dab %*% S)^2 / tr(Dab %*% Dab %*% S %*% S %*% Lambda)  
  
    for (j in 1:numpt) {
      x <- xbar[,i]
      x[d] <- x[d] + delta[j]
      # Wald statistic
      Qa[i,j] <- NT * t(x) %*% t(MaW) %*% ginv(MaW %*% S %*% t(MaW)) %*% MaW %*% x
      Qb[i,j] <- NT * t(x) %*% t(MbW) %*% ginv(MbW %*% S %*% t(MbW)) %*% MbW %*% x
      Qab[i,j] <- NT * t(x) %*% t(Mab) %*% ginv(Mab %*% S %*% t(Mab)) %*% Mab %*% x
  
      # FN(M) = N / tr(DM*SN) * X'*M*X ~ F(f,f0), f = tr(DM*SN)^2 / tr(M*SN*M*SN) f0 = tr(DM*SN)^2 / tr(DM^2*SN^2*Lambda)    
      ATSa[i,j] <- NT / tr(Da %*% S) * t(x) %*% Ma %*% x 
      ATSb[i,j] <- NT / tr(Db %*% S) * t(x) %*% Mb %*% x 
      ATSab[i,j] <- NT / tr(Dab %*% S) * t(x) %*% Mab %*% x 
    
    }
  }

  # chi(f)
  qQa <- qchisq(c(.9,.95,.99), df = a - 1)
  qQb <- qchisq(c(.9,.95,.99), df = b - 1)
  qQab <- qchisq(c(.9,.95,.99), df = (a - 1) * (b-1))
  pQ <- matrix(0,3,3)
  pwrQa <- matrix(0,numpt,3)
  pwrQb <- matrix(0,numpt,3)
  pwrQab <- matrix(0,numpt,3)
  # colnames(pQ) <- c("Qa","Qb", "Qab")
  for (i in 1:3) {
    pQ[i,1] <- mean(Qa[,1] <= qQa[i])
    pQ[i,2] <- mean(Qb[,1] <= qQb[i])
    pQ[i,3] <- mean(Qab[,1] <= qQab[i])
    for (j in 1:numpt) {
      pwrQa[j,i] <- mean(Qa[,j] >= qQa[i])
      pwrQb[j,i] <- mean(Qb[,j] >= qQb[i])
      pwrQab[j,i] <- mean(Qab[,j] >= qQab[i])
    }
  }

  # chi(f) / f
  qATSaN <- matrix(0, reps, 3)
  qATSbN <- matrix(0, reps, 3)
  qATSabN <- matrix(0, reps, 3)
  for (i in 1:reps) {
    qATSaN[i,] <- qchisq(c(.9,.95,.99), df = f[i,1]) / f[i,1]
    qATSbN[i,] <- qchisq(c(.9,.95,.99), df = f[i,2]) / f[i,2]
    qATSabN[i,] <- qchisq(c(.9,.95,.99), df = f[i,3]) / f[i,3]
  }
  pATSN <- matrix(0, 3 ,3)
  pwrATSaN <- matrix(0,numpt,3)
  pwrATSbN <- matrix(0,numpt,3)
  pwrATSabN <- matrix(0,numpt,3)
  # colnames(pATSS) <- c("ATSaS", "ATSbS", "ATSabS")
  for (i in 1:3) {
    pATSN[i,1] <- mean(ATSa[,1] <= qATSaN[,i])
    pATSN[i,2] <- mean(ATSb[,1] <= qATSbN[,i])
    pATSN[i,3] <- mean(ATSab[,1] <= qATSabN[,i])
    for (j in 1:numpt) {
      pwrATSaN[j,i] <- mean(ATSa[,j] >= qATSaN[,i])
      pwrATSbN[j,i] <- mean(ATSb[,j] >= qATSbN[,i])
      pwrATSabN[j,i] <- mean(ATSab[,j] >= qATSabN[,i])
    }
  }


  # F(f, f0)
  qFa <- matrix(0, reps, 3)
  qFb <- matrix(0, reps, 3)
  qFab <- matrix(0, reps, 3)
  for (i in 1:reps) {
    qFa[i,] <- qf(c(.9,.95,.99), df1 = f[i,1], df2 = f0[i,1])
    qFb[i,] <- qf(c(.9,.95,.99), df1 = f[i,2], df2 = f0[i,2])
    qFab[i,] <- qf(c(.9,.95,.99), df1 = f[i,3], df2 = f0[i,3])
  }
  pATS <- matrix(0,3,3)
  pwrATSa <- matrix(0,numpt,3)
  pwrATSb <- matrix(0,numpt,3)
  pwrATSab <- matrix(0,numpt,3)
  # colnames(pATS) <- c("ATSa", "ATSb", "ATSab")
  for (i in 1:3) {
    pATS[i,1] <- mean(ATSa[,1] <= qFa[,i])
    pATS[i,2] <- mean(ATSb[,1] <= qFb[,i])
    pATS[i,3] <- mean(ATSab[,1] <= qFab[,i])
    for (j in 1:numpt) {
      pwrATSa[j,i] <- mean(ATSa[,j] >= qFa[,i])
      pwrATSb[j,i] <- mean(ATSb[,j] >= qFb[,i])
      pwrATSab[j,i] <- mean(ATSab[,j] >= qFab[,i])
    }
  }
  
  delta <- c(-rev(delta[-1]),delta)
  powQa <- c(rev(pwrQa[-1,2]),pwrQa[,2])
  powQb <- c(rev(pwrQb[-1,2]),pwrQb[,2])
  powQab <- c(rev(pwrQab[-1,2]),pwrQab[,2])
  powATSaN <- c(rev(pwrATSaN[-1,2]),pwrATSaN[,2])
  powATSbN <- c(rev(pwrATSbN[-1,2]),pwrATSbN[,2])
  powATSabN <- c(rev(pwrATSabN[-1,2]),pwrATSabN[,2])
  powATSa <- c(rev(pwrATSa[-1,2]),pwrATSa[,2])
  powATSb <- c(rev(pwrATSb[-1,2]),pwrATSb[,2])
  powATSab <- c(rev(pwrATSab[-1,2]),pwrATSab[,2])
  
  dev.new(width=12, height=4)
  
  par(mfrow=c(2,3))
  #matplot(delta, cbind(pwrQa[,2],pwrATSaN[,2],pwrATSa[,2]), type = c("l"), xlab = expression(delta), ylab = "Power", main="Power Comparison on A", pch = 1:3, col = 1:3)
  matplot(delta, cbind(powQa,powATSaN,powATSa), type = c("l"), xlab = expression(delta), ylab = "Power", main="Power Comparison on A", pch = 1:3, col = 1:3)
  legend("top", c("Q","Fn","F"), pch = 1:3, col = 1:3)
  #matplot(delta, cbind(pwrQb[,2],pwrATSbN[,2],pwrATSb[,2]), type = c("l"), xlab = expression(delta), ylab = "Power",main="Power Comparison on B", pch = 1:3, col = 1:3)
  matplot(delta, cbind(powQb,powATSbN,powATSb), type = c("l"), xlab = expression(delta), ylab = "Power", main="Power Comparison on B", pch = 1:3, col = 1:3)
  legend("top", c("Q","Fn","F"), pch = 1:3, col = 1:3)
  #matplot(delta, cbind(pwrQab[,2],pwrATSabN[,2],pwrATSab[,2]), type = c("l"), xlab = expression(delta), ylab = "Power",main="Power Comparison on AB", pch = 1:3, col = 1:3)
  matplot(delta, cbind(powQab,powATSabN,powATSab), type = c("l"), xlab = expression(delta), ylab = "Power", main="Power Comparison on AB", pch = 1:3, col = 1:3)
  legend("top", c("Q","Fn","F"), pch = 1:3, col = 1:3)
  #matplot(delta[(numpt/2):(numpt/2+20)], cbind(pwrQa[,2],pwrATSaN[,2],pwrATSa[,2])[(numpt/2):(numpt/2+20),], type = c("b"), xlab = expression(delta), ylab = "Power", main="Power Comparison on A", pch = 1:3, col = 1:3)
  matplot(delta[(numpt-1):(numpt+20)], cbind(powQa,powATSaN,powATSa)[(numpt-1):(numpt+20),], type = c("b"), xlab = expression(delta), ylab = "Power", main="Power Comparison on A", pch = 1:3, col = 1:3)
  legend("topleft", c("Q","Fn","F"), pch = 1:3, col = 1:3)
  #matplot(delta[(numpt/2):(numpt/2+20)], cbind(pwrQb[,2],pwrATSbN[,2],pwrATSb[,2])[(numpt/2):(numpt/2+20),], type = c("b"), xlab = expression(delta), ylab = "Power",main="Power Comparison on B", pch = 1:3, col = 1:3)
  matplot(delta[(numpt-1):(numpt+20)], cbind(powQb,powATSbN,powATSb)[(numpt-1):(numpt+20),], type = c("b"), xlab = expression(delta), ylab = "Power",main="Power Comparison on B", pch = 1:3, col = 1:3)
  legend("topleft", c("Q","Fn","F"), pch = 1:3, col = 1:3)
  #matplot(delta[(numpt/2):(numpt/2+20)], cbind(pwrQab[,2],pwrATSabN[,2],pwrATSab[,2])[(numpt/2):(numpt/2+20),], type = c("b"), xlab = expression(delta), ylab = "Power",main="Power Comparison on AB", pch = 1:3, col = 1:3)
  matplot(delta[(numpt-1):(numpt+20)], cbind(powQab,powATSabN,powATSab)[(numpt-1):(numpt+20),], type = c("b"), xlab = expression(delta), ylab = "Power",main="Power Comparison on AB", pch = 1:3, col = 1:3)
  legend("topleft", c("Q","Fn","F"), pch = 1:3, col = 1:3)
  
  comparisontable <- cbind(pQ[,1], pATSN[,1], pATS[,1], pQ[,2], pATSN[,2], pATS[,2], pQ[,3], pATSN[,3], pATS[,3])
  colnames(comparisontable) <- c("Qa", "ATSaN", "ATSa", "Qb", "ATSbN", "ATSb", "Qab", "ATSabN", "ATSab")
  return(comparisontable)
}

# two-way ANOVA y_ijk = alpha_i + beta_j + (alphabeta)_ij + e_ijk, alpha_i = 0, beta_j = 0, (alphabeta)_ij = 0
reps <- 5000
a = 2; b = 5
set.seed(201710)
# case 1
n1 <- rep(7,10)
table2.1 <- BoxType(a,b,n1,reps)


# case 2 
n2 <- rep(c(11, 10, 9, 8, 7), 2)
table2.2 <- BoxType(a,b,n2,reps)

table2.1
table2.2





