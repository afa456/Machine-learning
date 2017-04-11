rm(list=ls())
myData = read.csv("myData(manipulate)-1.csv")
# Remove row names
myData = myData[,2:dim(myData)[2]]
###################### ITERATIONS ######################
myX = as.matrix(myData[,2:dim(myData)[2]])
myY = myData[,1]
dn = dim(myX)[1]
dp = dim(myX)[2]
meany =mean(myY)
Aerror=c()
Noinfo=c()
# Standardize the predictors
myX = scale(myX)
myY = scale(myY)
nbsteps = 150
k=1
XA = NULL
activeSet = rep(F, dp)
mu = rep(0,dn)
beta = rep(0,dp)
BETA = matrix(0,dp,nbsteps)
C_remove = F
logic_end = T
while(k<=nbsteps && logic_end) {
  error = myY-mu
  if(!C_remove) {
    errCor = cor(myX, error)
    errCor[activeSet] = 0
    index = which.max(abs(errCor))
    activeSet[index] = T
  } else {
    index = which(rmVar==names(as.data.frame(myX)))
    activeSet[index] = F
  }
  sign = sign(cor(myX, error))
  signId = diag(dp)
  for(i in 1:dp) {
    signId[i,i] =sign[i]
  }  
  XA = (myX %*% signId)[,activeSet]
  XA = as.matrix((myX %*% signId)[,activeSet])
  oneN = rep(1,sum(activeSet))
  cHat = t(myX) %*% (error)
  cCap =max(abs(t(XA) %*% error))
  a = t(myX) %*% XA %*% solve(t(XA) %*% XA) %*% oneN
  gamma = c(((cCap-cHat)/(1-a))[!activeSet], 
             ((cCap+cHat)/(1+a))[!activeSet])
  gamma = min(gamma[gamma>0.00001])
  #lasso modification
  bTilde = solve(t(XA) %*% XA) %*% oneN
  cStar = -(BETA[activeSet,k-1] / (sign[activeSet]*bTilde))
  cj=Inf
  if(length(cStar[cStar>0.00001])==0) {
    cj=Inf
    indexStar = 0
    rmVar = NULL
  } else {
    cj = min(cStar[cStar>0.00001])
    indexStar = which(cStar==cj)
    rmVar = names(as.data.frame(myX))[activeSet][indexStar]
  }
  C_remove = cj < gamma
  if(C_remove) { gamma = cj }
  
  uTilde = gamma * (XA %*% solve(t(XA) %*% XA) %*% oneN)
  mu = mu + uTilde
  b = gamma * (solve(t(XA) %*% XA) %*% oneN)
  beta[activeSet] = beta[activeSet] + b*sign[activeSet]
  BETA[,k] = beta
  error = myY-mu
  cHat = t(myX) %*% error
  k = k+1
  logic_end = dim(XA)[2] < dim(XA)[1]-1
  if(C_remove) {logic_end = T}
  Aerror=c(Aerror,sum(abs((myY-meany-myX %*%beta)/29)))
  summ = 0
  for(i in 1:29){
    response = matrix(myY[i],29,1)
    summ = summ + sum(abs(response-meany-mu))
  }
  Noinfo=c(Noinfo,summ/29^2)
}
k = k-1
BETA = BETA[,1:k]
t = colSums(abs(BETA))
plot(t, rep(0,k), type="l", ylim=c(-1.2, 1.2), xlab="t value", ylab="coefficient estimate")
title("Lasso Path No.1", cex.main = 1, font.main= 1)
for(i in 1:dp) { lines(t, BETA[i,], lwd=2, col=sample(colours(), 1)) }
plot(t,Aerror)
plot(t,Noinfo)