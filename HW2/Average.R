rm(list=ls())
library(Hmisc)
###############Average 3 Dataset###############
Error1 = read.csv("Error632-1p.csv")
Error2 = read.csv("Error632-2p.csv")
Error3 = read.csv("Error632-3p.csv")
#myData = read.csv("myData(manipulate)-1.csv")
# Remove row names
Error1 = as.matrix(Error1[,2:dim(Error1)[2]])
Error2 = as.matrix(Error2[,2:dim(Error2)[2]])
Error3 = as.matrix(Error3[,2:dim(Error3)[2]])

MatErr=vector("list",3)
MatErr[[1]]=Error1
MatErr[[2]]=Error2
MatErr[[3]]=Error3
Ts=MatErr[[1]][,1]
Overall=matrix(0,length(Ts),2)
Overall[,1]=Ts
for(t in 1:length(Ts)){
err=MatErr[[1]][t,2]
for(j in 2:3){
err=err+approxExtrap(as.matrix(MatErr[[j]][,1]),as.matrix(MatErr[[j]][,2]),Ts[t])$y
}
Overall[t,2]=err/3
}
par(mfrow=c(1,1))
plot(Overall, xlab="t value", ylab="predictive error")
title("average 0.632+ bootstrap estimate", cex.main = 1, font.main= 1)