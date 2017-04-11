rm(list=ls())
myData = read.csv("myData-3.csv")
#myData = read.csv("myData(manipulate)-1.csv")
# Remove row names
myData = myData[,2:dim(myData)[2]]

###################### ITERATIONS ######################

myX = as.matrix(myData[,2:dim(myData)[2]])
myY = myData[,1]
dn = dim(myX)[1]
dp = dim(myX)[2]

# Standardize the predictors
myX = scale(myX, drop(rep(1,dn) %*% myX)/dn, FALSE)
myX = scale(myX, FALSE, sqrt(drop(rep(1,dn) %*% (myX^2))))
myY = drop(myY - mean(myY))

nbsteps = 28
k=1
XA = NULL
activeSet = rep(F, dp)
mu = rep(0,dn)
beta = rep(0,dp)
BETA = matrix(0,dp,nbsteps)

while(k<=nbsteps) {
  error = myY-mu
  errCor = cor(myX, error)
  errCor[activeSet] = 0
  index = which.max(abs(errCor))
  activeSet[index] = T

  sign = sign(cor(myX, error))
  signId = diag(dp)
  for(i in 1:dp) {
    signId[i,i] =sign[i]
  }  
  XA = (myX %*% signId)[,activeSet]
  XA = as.matrix((myX %*% signId)[,activeSet])
  
  oneN = rep(1,sum(activeSet))
  
  cHat = t(myX) %*% (error)
  cCap =max(abs(cHat))
  a = t(myX) %*% XA %*% solve(t(XA) %*% XA) %*% oneN
  
  gamma = c(((cCap-cHat)/(1-a))[!activeSet], 
             ((cCap+cHat)/(1+a))[!activeSet])
  gamma = min(gamma[gamma>0.0001])
  
  uTilde = gamma * (XA %*% solve(t(XA) %*% XA) %*% oneN)
  mu = mu + uTilde
  b = gamma * (solve(t(XA) %*% XA) %*% oneN)
  beta[activeSet] = beta[activeSet] + b*sign[activeSet]
  BETA[,k] = beta
  k = k+1
}
beta = solve(t(XA) %*% XA) %*% t(XA) %*% myY