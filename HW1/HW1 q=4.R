# ISyE 6740 ยก Take Home Exam #1
library(mvtnorm)

# Read handwritten digits data
myData=read.csv("semeion.csv",header=FALSE)
myX=data.matrix(myData[,1:256])
myLabel=apply(myData[,257:266],1,function(xx){
return(which(xx=="1")-1)
})

# Number of rows
NR=dim(myX)[1]
# Number of columns
NC=dim(myX)[2]
# Number of clusters
Nclu=10
# Number of principal component
q=4

#Cluster data by using kmeans method, K=10
myCluster=kmeans(myX,10,iter.max=20,nstart=10)

#Initialization: assignments of data
gamma=matrix(0,nrow=NR,ncol=Nclu)
for(i in 1:NR) {
	Clc=myCluster$cluster[i]
	gamma[i, Clc]=1
}

#For likelihood
likeli=rep(0,20)


#Iterations
for(it in 1:20){
    
    N=matrix(0,1,10)
    for(i in 1:10) {
	    N[i]=sum(gamma[,i])	
    }

    Mu=matrix(0,nrow=Nclu,ncol=NC)
    pi=matrix(0,10,1) 

#Initialization of covariance matrices
    sigma=array(dim=c(256,256,10))
    px=matrix(0,NR,Nclu)

#------------------------M-Step------------------------

    for(k in 1:Nclu){
        mu_k=rep(0,256)
        for(n in 1:NR){
            mu_k=mu_k+gamma[n,k]*myX[n, ] 
        }
        Mu[k,]=mu_k/N[k]
    }

    pi=colSums(gamma)/NR
    
    for (k in 1:Nclu){
        Covar_k=matrix(0,256,256)

        for(n in 1:NR){
            Xi_Bar=myX[n, ]-Mu[k, ]
            Covk_temp=(Xi_Bar %*% t(Xi_Bar))*gamma[n,k]
            Covar_k=Covar_k+Covk_temp
        }
        Covar_k=Covar_k/N[k]
        myeigen=eigen(Covar_k,symmetric=TRUE)
        Vq=myeigen$vectors[,1:q]
        
        sigma_2=sum(myeigen$values[q+1:NC],na.rm=TRUE)/(NC-q)
        diag_q=diag(q)

        for(nq in 1:q){
            diag_q[nq,nq]=sqrt(myeigen$values[nq]-sigma_2)	
        }
        
        Wq=Vq %*% diag_q
        sigma[ , ,k]=Wq %*% t(Wq)+(sigma_2*diag(NC))
    }
#------------------------E-Step----------------------
    for (k in 1:Nclu){
	    px[ ,k]=pi[k]*dmvnorm(myX,Mu[k, ],sigma[ , ,k],log=FALSE) 
    }
    for (i in 1:NR){
        for(k in 1:Nclu){
        gamma[i,k]=px[i,k]/sum(px[i, ])
        }
    }   
#------------------------Likelihood------------------------
likeli[it]=sum(log(rowSums(px)))
print(c(it,likeli[it]))
}

#------------------------Computer AIC----------------------
AIC= -2*likeli[20]+2*(NC*q+1-q*(q-1)/2)
likeli
AIC
# Plot likelihood VS iter
dev.new(width=6,height=4)
par(mai=c(0.5,0.45,0.35,0.05),cex=0.8)
plot(1:20,likeli, pch=19,axes=TRUE)
title("q=6 Likelihood")

# Visulazation pictures
dev.new(width=6,height=10)
par(mai=c(0,0,0,0),cex=0.8,mfrow=c(10,6))
for(k in 1:Nclu){
    image(t(matrix(Mu[k, ],byrow=TRUE,16,16)[16:1, ]),col=gray(0:256/256),axes=FALSE)
    box( )
    for(j in 1:5){
        temp=rmvnorm(1,mean=Mu[k, ],sigma[ , ,k])
        image(t(matrix(temp,byrow=TRUE,16,16)[16:1, ]),col=gray(0:256/256),axes=FALSE)
    box ( )
    }
}

# Accuracy Assessment







