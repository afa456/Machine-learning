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

############ Develop Prediction Models ##############
############ #1 SVM Model ##############
library(e1071)
library(pROC)
library(lattice)
nIter = 25
gamma = 10^seq(-6,0,0.5)
cost = 10^seq(-6,2,0.5)
svm.auc = array(0, dim=c(nIter,length(gamma),length(cost)))

for(iter in 1:nIter){
  index = sample(nrow(newData))[1:round(nrow(newData)/10)]
  trainData = newData[-index,]
  testData = newData[index,]
  for(g in gamma){ 
  	for(c in cost){
    	model = svm(y~., data=trainData, gamma=g, cost=c, probability=TRUE)
    	myPred = predict(model, testData, probability=TRUE)
    	myProb = attr(myPred, "probabilities")[,1]
    	myroc = roc(testData$y, myProb)
    	svm.auc[iter,which(gamma==g),which(cost==c)] = myroc$auc
  }}
  print(iter)
}
svmtest = apply(svm.auc,c(2,3),mean)
levelplot(svmtest)
# AUC max 0.95903
# Cost = 3.162, Gamma = 0.001

############ #2 Neural Network Model ##############
library(nnet)
library(pROC)
library(lattice)
nIter = 25
size = seq(3,24,3)
decay = 10^seq(-6,1,0.5)
nnet.auc = array(0, dim=c(nIter,length(size),length(decay)))

for(iter in 1:nIter){
  index = sample(nrow(newData))[1:round(nrow(newData)/10)]
  trainData = newData[-index,]
  testData = newData[index,]
  for(s in size){ 
  	for(d in decay){
    	model = nnet(y~., data=trainData, size=s, decay=d)
    	myPred = predict(model, testData, probability = TRUE)
    	myroc = roc(testData$y, myPred)
    	nnet.auc[iter,which(size==s),which(decay==d)] = myroc$auc
  }}
  print(iter)
}
nnettest = apply(nnet.auc,c(2,3),mean)
levelplot(nnettest)
# AUC max 0.95957
# Size = 3, Decay = 0.1

############ #3 Random Forest Model ##############
library(randomForest)
library(pROC)
library(ggplot2)
nIter = 25
ntree = seq(10,300,10)
rf.auc = matrix(NA, nIter, length(ntree))

for(iter in 1:nIter){
  index = sample(nrow(newData))[1:round(nrow(newData)/10)]
  trainData = newData[-index,]
  testData = newData[index,]
  for(n in ntree){
    model  = randomForest(y~., data=trainData, ntree=n)
    myPred = predict(model, testData, type="prob")[,2]
    myroc = roc(testData$y, myPred)
    rf.auc[iter,which(ntree==n)] = as.numeric(myroc$auc)
  }
  print(iter)
}
rfplot = data.frame(ntree, colMeans(rf.auc))
names(rfplot) = c("ntree", "average_auc")
ggplot(rfplot, aes(x=ntree, y=average_auc, colour=factor(ntree))) + 
geom_point() + geom_smooth()
# AUC max 0.950198
# ntree = 150

############ #4 Boosted Model ##############
library(ada)
library(pROC)
library(ggplot2)
nIter = 25
nu = 10^seq(-6,0,0.25)
ada.auc = matrix(NA, nIter, length(nu))

for(iter in 1:nIter){
  index = sample(nrow(newData))[1:round(nrow(newData)/10)]
  trainData = newData[-index,]
  testData = newData[index,]
  for(n in nu){
    model = ada(y~.,  data=trainData, nu=n)
    myPred = predict(model, testData, type="prob")[,2]
    myroc = roc(testData$y, myPred)
    ada.auc[iter,which(nu==n)] = as.numeric(myroc$auc)
  }
  print(iter)
}
adaplot = data.frame(nu, colMeans(ada.auc))
names(adaplot) = c("nu", "average_auc")
ggplot(adaplot, aes(x=log(nu), y=average_auc, colour=factor(nu))) + 
geom_point() + geom_smooth()
# AUC max 0.94445
# nu = 0.1

############ #5 Logistic Regression Model ##############
library(pROC)
library(ggplot2)
nIter = 25
mysteps = seq(1,20,1)
glm.auc = matrix(NA, nIter, length(mysteps))

for(iter in 1:nIter){
  index = sample(nrow(newData))[1:round(nrow(newData)/10)]
  trainData = newData[-index,]
  testData = newData[index,]
  for(s in mysteps){
	test = glm(y~.,data=trainData,family=binomial(link=probit))
	step = step(test, direction = "both", steps = s)
	newtrainData = step $ model
	model = glm(y~.,data=newtrainData,family=binomial(link=probit))
    myPred = predict(model, testData, family=binomial(link=probit), type = "response")
    myroc = roc(testData$y, myPred)
    glm.auc[iter,which(mysteps==s)] = as.numeric(myroc$auc)
  }
  print(iter)
}
lgplot = data.frame(mysteps, colMeans(glm.auc))
names(lgplot) = c("steps", "average_auc")
ggplot(lgplot, aes(x=steps, y=average_auc, colour=factor(steps))) + 
geom_point() + geom_smooth()
# AUC max 0.95808
# mystep = 2



