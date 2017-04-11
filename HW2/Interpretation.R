################Interpretation################
rm(list=ls())
library(Hmisc)
source("Lars Al.R")
Results=vector("list",3)

myData = read.csv("myData-1.csv")
myData = myData[,2:dim(myData)[2]]
myData2 = read.csv("myData-2.csv")
myData2 = myData2[,2:dim(myData2)[2]]
myData3 = read.csv("myData-3.csv")
myData3 = myData3[,2:dim(myData3)[2]]
myY=as.vector(myData[,1])

Data=myData
Data0=apply(Data[,2:{dim(Data)[2]}],2,function(xx){
return((xx-mean(xx,na.rm=TRUE))/sd(xx,na.rm=TRUE))
})
Data[,2:dim(Data)[2]]=Data0
Data=as.matrix(Data)
Results[[1]]=myLARS(Data,"lasso")

Data2=myData2
Data0=apply(Data2[,2:{dim(Data2)[2]}],2,function(xx){
return((xx-mean(xx,na.rm=TRUE))/sd(xx,na.rm=TRUE))
})
Data2[,2:dim(Data2)[2]]=Data0
Data2=as.matrix(Data2)
Results[[2]]=myLARS(Data2,"lasso")

Data3=myData3
Data0=apply(Data3[,2:{dim(Data3)[2]}],2,function(xx){
return((xx-mean(xx,na.rm=TRUE))/sd(xx,na.rm=TRUE))
})
Data3[,2:dim(Data3)[2]]=Data0
Data3=as.matrix(Data3)
Results[[3]]=myLARS(Data3,"lasso")

Ts=Results[[1]]$Ts
Estimates=matrix(0,150,1)

T=2
for(j in 1:3){
Coef2=matrix(0,150,1)
Betas=Results[[j]]$Betas
Bet=as.matrix(Betas)
Coef=t(sapply(Bet, unlist))
for (c in 1:150){
Coef2[c]=approxExtrap(as.matrix(Ts),as.matrix(Coef[,c]),T)$y}
Estimates=Estimates+Coef2
}
Estimates=Estimates/3
Error=myY-Data[,2:dim(Data)[2]]%*%Estimates
plot(Error)
# selected variables
Selected=attributes(Data[,2:{dim(Data)[2]}])[[2]][[2]][which(abs(Estimates)>1e-5)]
# Largest coefficient
LargeCoef=which(abs(Estimates)==max(abs(Estimates)))
LCoef=attributes(Data[,2:{dim(Data)[2]}])[[2]][[2]][LargeCoef]
Data2=Data[,2:{dim(Data)[2]}]
Data2[,LargeCoef]=Data[,LargeCoef+1]+1
change=(Data2-Data[,2:dim(Data)[2]])%*%Estimates