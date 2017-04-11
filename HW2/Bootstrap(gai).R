rm(list=ls())
library(Hmisc)
#myData = read.csv("myData-3.csv")
myData = read.csv("myData(manipulate)-1.csv")
# Remove row names
myData = myData[,2:dim(myData)[2]]
source("Lars Al.R")


###############Set up the bootStrap###############
B=25
myY=as.vector(myData[,1])
Mean=mean(myY)
dn=length(myData[,1])

Data=myData
Data0=apply(Data[,2:{dim(Data)[2]}],2,function(xx){
return((xx-mean(xx,na.rm=TRUE))/sd(xx,na.rm=TRUE))
})
Data[,2:dim(Data)[2]]=Data0
Data=as.matrix(Data)
boots=matrix(0,B,dn)
for(i in 1:B){
boots[i,]=sample(dn,replace=TRUE)
}
Bootstrps=vector("list",B)
for(k in 1:B){

Data2=Data[boots[k,],]
Data0=apply(Data2[,2:{dim(Data2)[2]}],2,function(xx){
return((xx-mean(xx,na.rm=TRUE))/sd(xx,na.rm=TRUE))
})
Data2[,2:dim(Data2)[2]]=Data0
Data2=as.matrix(Data2)
Bootstrps[[k]]=myLARS(Data2,"lasso",500)
}

###############Leave-One-Out Bootstrap###############
L=c()
for (b in 1:B){
	L=c(L,length(Bootstrps[[b]]$Ts))
}
index=which(L==max(L))[1]
Ts=Bootstrps[[index]]$Ts
Err=matrix(0,length(Ts),2)
for (t in 1:length(Ts)){
err=0
N=0
for(i in 1:dn){
Ci=apply(boots,1,function(xx){return(!any(xx==i))})
if(sum(Ci)>0){
bootstrp=Bootstrps[Ci]
N=N+1
er=0
for (bi in bootstrp){
Tsx=bi$Ts
Betay=as.matrix(bi$Betas)
Coef1=t(sapply(Betay, unlist))
Coef2=matrix(0,150,1)
for (c in 1:150){
Coef2[c]=approxExtrap(as.matrix(bi$Ts),as.matrix(Coef1[,c]),Ts[t])$y
#Coef2[c]=approx(as.matrix(bi$Ts),as.matrix(Coef[,c]),Ts[t],rule=2)$y
}
er=er+abs(myY[i]-Mean-Data[i,2:dim(Data)[2]]%*%Coef2)
}}
er=er/sum(Ci)
err=err+er
}
err=err/N
Err[t,1]=Ts[t]
Err[t,2]=err
}

###############0.632+ bootstrap###############
Results=myLARS(Data,"lasso")

AbsErr=Results$AppErrs
NoInErr=Results$NoInfErrs
Tss1=Err[,1]
Tss2=Results$Ts
Err632plus=matrix(0,length(Tss1),2)
Err632=matrix(0,length(Tss1),2)
for (t in 1:length(Tss1)){
errhat=Err[t,2]
abs_err=approxExtrap(as.matrix(Tss2),as.matrix(AbsErr),Tss1[t])$y
noinf_er=approxExtrap(as.matrix(Tss2),as.matrix(NoInErr),Tss1[t])$y

R=(errhat-abs_err)/(noinf_er-abs_err)
omega=0.632/(1-0.368*R)
Err632plus[t,1]=Tss1[t]
Err632[t,1]=Tss1[t]
Err632plus[t,2]=(1-omega)*abs_err+omega*errhat
Err632[t,2]=0.632*abs_err+0.368*errhat
}

plot(Err632plus, xlab="t value", ylab="predictive error")
title("0.632+ bootstrap estimate", cex.main = 1, font.main= 1)
plot(Err632, xlab="t value", ylab="predictive error")
title("0.632+ bootstrap estimate", cex.main = 1, font.main= 1)

write.csv(Err632plus,file="Error632p-3.csv")
write.csv(Err632,file="Error632-3.csv")
