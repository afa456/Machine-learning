############ Data Manipulation ##############
rm(list=ls())
myData=read.csv("AngleClosure.csv",header=TRUE,na.strings=c("NA","."))

myData1 = myData[,-c(1,15,16)]
myData1[,-21] = data.matrix(myData1[,-21])
myData1[,21] = factor(myData1[,21])
myData2 = myData1[,1:11]
myData3 = na.omit(myData2)
y = myData1[,21]
newData = cbind(y, myData3)

############ Stacking  ##############
library(e1071)
library(nnet)
library(randomForest)
library(ada)
library(pROC)
library(quadprog)

nIter=25
X = matrix(0,,5)
Y = matrix(0,,1)

for(iter in 1:nIter){
	index = sample(nrow(newData))[1:round(nrow(newData)/10)]
  	trainData = newData[-index,]
  	testData = newData[index,]
	fusionX = matrix(0,dim(testData)[1],)	

    #SVM Model
    svm_model = svm(y~., data=trainData, gamma=0.001, cost=3.162, probability=TRUE)
    svm_myPred = predict(svm_model, testData, probability=TRUE)
    svm_myProb = attr(svm_myPred, "probabilities")[,1]
	fusionX = cbind(fusionX,svm_myProb)    
   
    #NNET Model
    nnet_model = nnet(y~., data=trainData, size=3, decay=0.1)
    nnet_myPred = as.numeric(predict(nnet_model, testData, probability = TRUE))
	fusionX = cbind(fusionX,nnet_myPred) 

    #RF Model
    rf_model  = randomForest(y~., data=trainData, ntree=150)
    rf_myPred = predict(rf_model, testData, type="prob")[,2]
	fusionX = cbind(fusionX,rf_myPred) 

    #ADA Model
    ada_model = ada(y~.,  data=trainData, nu=0.1)
    ada_myPred = predict(ada_model, testData, type="prob")[,2]
	fusionX = cbind(fusionX,ada_myPred) 

    #LG Model
	lm = glm(y~.,data=trainData,family=binomial(link=probit))
	step = step(lm, direction = "both", steps = 2)
	newtrainData = step $ model
	lg_model = glm(y~.,data=newtrainData,family=binomial(link=probit))
    lg_myPred = predict(lg_model, testData, family=binomial(link=probit), type = "response")
	fusionX = cbind(fusionX,lg_myPred)

    temp = as.numeric(testData$y=='YES')
	Y = rbind(Y,as.matrix(temp))
	X = rbind(X, fusionX[,2:6])
}
Miu = X
Response = Y
Miu = Miu[2:dim(Miu)[1],]
Response = Response[2:dim(Response)[1],]

D = t(Miu) %*% Miu
d = t(Response) %*% Miu
A = cbind(rep(1,5), diag(5))
b = c(1, 0, 0, 0, 0, 0)

solution = solve.QP(D, d, A, b,meq =1)
weightsConstrained = solution$solution
# weightsConstrained
# 1.457315e-01 2.662107e-01 4.312247e-19 1.891759e-01 3.988819e-01
weightsUnConstrained = solution$unconstrained.solution
#weightsUnConstrained
#  0.15072539  0.27214284 -0.06550852  0.21353243  0.41952495