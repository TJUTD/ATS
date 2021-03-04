library(MASS)
library(magic)

###########################################################
# nonparametric two way Anova Type Statistics example
###########################################################
npATS2Way <- function(formula,data) {

  rawdat <- model.frame(formula, data)
  for(i in ncol(rawdat):2) {
	  rawdat[,i] <-as.factor(rawdat[,i]) }
  nf <- ncol(rawdat)-1
  fname <- colnames(rawdat[,2:ncol(rawdat)])
  
  for(i in ncol(rawdat):2) {
	  rawdat <- rawdat[order(rawdat[,i]),] }
  
  
  n <- as.vector(table(rawdat[,ncol(rawdat):2]))
  rid <- c(0,cumsum(n))
  
  levela <- levels(rawdat[,2])
  levelb <- levels(rawdat[,3])
  a <- nlevels(rawdat[,2])
  b <- nlevels(rawdat[,3])
  d <- a * b
  NT <- sum(n)
  psf <- rep(1:d, n)


  # transfer to rank
  dat <- rank(rawdat[,1])
 
  pd <- matrix(0, d, 1)
  Sd <- matrix(0,d,1)
  VN <- matrix(0, d, 1)
  for (i in 1:d) {
    M <- matrix(1/n[i],n[i],n[i])
    MX <- M %*% dat[(rid[i]+1):rid[i+1]] 
    pd[i] <- (MX[1] - 0.5) / NT
    VN[i] <- M[1,] %*% ((dat[psf==i]  - MX) * (dat[psf==i]  - MX))  / NT / (n[i]-1)
    Sd[i] <- n[i] * M[1,] %*% ((dat[psf==i]  - MX) * (dat[psf==i]  - MX))  / (NT^2) / (n[i]-1)
  }

  # utility fuction
  # trace of matrix
  tr <- function (m) {
    if (!is.matrix(m) | (dim(m)[1] != dim(m)[2])) 
      stop("m must be a square matrix")
    return(sum(diag(m), na.rm = TRUE))
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

  # Wald test
  Wald <- function(p,M,Var,df){
    Q =t(M %*% p) %*% ginv(M %*% Var %*% t(M)) %*% M %*% p
    pValue = 1-pchisq(Q,df)
    return(c(Q, df, pValue))
  }

  # ANOVA-type test
  ATS <- function(p,M,D,Var){
    F <- 1 / tr(D %*% Var) * t(p) %*% M %*% p
    f <- tr(D %*% Var)^2 / tr(M %*% Var %*% M %*% Var)
    f0 <- tr(D %*% Var)^2 / tr(D %*% D %*% Var %*% Var %*% diag(1/(n-1)))
    pValue = 1-pf(F,f,f0)
    return(c(F, f, f0, pValue))
  }  

  # result
  S <- matrix(0,d,d)
  diag(S) <- VN 
  # Wald statistic
  WaldA <- Wald(pd, MaW, S/NT, a-1)
  WaldB <- Wald(pd, MbW, S/NT, b-1)
  WaldAB <- Wald(pd, Mab, S/NT, (a - 1) * (b-1))

  # ATS
  ATSA <- ATS(pd, Ma, Da, S/NT)
  ATSB <- ATS(pd, Mb, Db, S/NT)
  ATSAB <- ATS(pd, Mab, Dab, S/NT)  
 
  wts <- rbind(WaldA, WaldB, WaldAB)
  colnames(wts) <- c("Statistic", "df", "p-Value")
  ats <- rbind(ATSA, ATSB, ATSAB)
  colnames(ats) <- c("Statistic", "f", "f0", "p-Value")
  descriptive <- data.frame(cbind(rep(levela, rep(b,a)), rep(levelb, a), n, pd, Sd))
  colnames(descriptive) <- c(fname, "size", "p", "Sn")

  result <- list(Descriptive = descriptive, Wald = wts, ANOVAType = ats)
  return(result)
}


file <- "~/kidney.txt"
kdata <- read.table(file, header = T)
# weight ~ Dose * sex
npATS2Way(weight ~ sex * Dose , kdata)
npATS2Way(weight ~ Dose * sex , kdata)