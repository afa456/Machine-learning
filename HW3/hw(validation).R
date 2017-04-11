############ Cross Validation ##############
rm(list=ls())
myData=read.csv("AngleClosure.csv",header=TRUE,na.strings=c("NA","."))

myData1 = myData[,-c(1,15,16)]
myData1[,-21] = data.matrix(myData1[,-21])
myData1[,21] = factor(myData1[,21])
myData2 = myData1[,1:11]
myData3 = na.omit(myData2)
y = myData1[,21]
newData = cbind(y, myData3)

attnames = names(newData)[-1]
caseData = read.csv("AngleClosure_ValidationCases.csv")
controlData = read.csv("AngleClosure_ValidationControls.csv")

myCase.right = caseData[,c(19,21,23,24,25,26,27,30,31,32,36)]
myCase.left = caseData[,c(7,9,11,12,13,14,15,30,31,32,36)]

myControl.right = controlData[,c(18,20,22,23,24,25,26,29,30,31,35)]
myControl.left = controlData[,c(6,8,10,11,12,13,14,29,30,31,35)]

# for cases data
# delete rows which have any missing values. 
# preferentially take right eye data
logic_r = which(complete.cases(myCase.right))
logic_l = which(complete.cases(myCase.left))
logic_l = logic_l[which(!logic_l %in% logic_r)]

temp_r = myCase.right[logic_r,]
temp_l = myCase.left[logic_l,]
names(temp_r) = names(temp_l) = attnames
myCase = rbind(temp_r, temp_l)
y = rep("YES", nrow(myCase))
myCase = cbind(y, myCase)

# for controls data
logic_r = which(complete.cases(myControl.right))
# No missing value rows in control for right eye
#logic_l = which(complete.cases(myControl.left))
#logic_l = logic_l[which(!logic_l %in% logic_r)]
temp_r = myControl.right[logic_r,]
#temp_l = myControl.left[logic_l,]
names(temp_r) = attnames
#names(temp_l) = attnames
myControl = temp_r
y = rep("NO", nrow(myControl))
myControl = cbind(y, myControl)

valid_data = rbind(myCase, myControl)
row.names(valid_data) = NULL
############ Validation ##############
library(e1071)
library(nnet)
library(randomForest)
library(ada)
library(pROC)

# SVM Model
svm_model = svm(y~., data=newData, gamma=0.001, cost=3.162, probability=TRUE)
svm_myPred = predict(svm_model, valid_data, probability=TRUE)
svm_myProb = attr(svm_myPred, "probabilities")[,1]
svm_roc = roc(as.numeric(valid_data[,1]),as.numeric(svm_myProb), plot=T)
svm_auc = svm_roc$auc
# Area under the curve: 0.9519

#NNET Model
nnet_model = nnet(y~., data=newData, size=3, decay=0.1)
nnet_myPred = as.numeric(predict(nnet_model, valid_data, probability = TRUE))
nnet_roc = roc(valid_data$y, nnet_myPred, plot=T)
nnet_auc = nnet_roc$auc
#Area under the curve: 0.9649

#RF Model
rf_model  = randomForest(y~., data=newData, ntree=150)
rf_myPred = predict(rf_model, valid_data, type="prob")[,2]
rf_roc = roc(valid_data$y, rf_myPred, plot=T)
rf_auc = rf_roc$auc
#Area under the curve: 0.952

#ADA Model
ada_model = ada(y~.,  data=newData, nu=0.1)
ada_myPred = predict(ada_model, valid_data, type="prob")[,2]
ada_roc = roc(valid_data$y, ada_myPred, plot=T)
ada_auc = ada_roc$auc
#Area under the curve: 0.961

#LG Model
lm = glm(y~.,data=newData,family=binomial(link=probit))
step = step(lm, direction = "both", steps = 2)
newtrainData = step $ model
lg_model = glm(y~.,data=newtrainData,family=binomial(link=probit))
lg_myPred = predict(lg_model, valid_data, family=binomial(link=probit), type = "response")
lg_roc = roc(valid_data$y, lg_myPred, plot=T)
lg_auc = lg_roc$auc
#Area under the curve: 0.9536

P = cbind (svm_myProb, nnet_myPred, rf_myPred, ada_myPred, lg_myPred)
# stacked constrained
stackedC_myPred = P %*% weightsConstrained
Con_roc = roc(as.numeric(valid_data[,1]),as.numeric(stackedC_myPred), plot=T)
Con_auc = Con_roc$auc
# Area under the curve: 0.9597
# stacked unconstrained
stackedU_myPred = P %*% weightsUnConstrained
Uncon_roc = roc(as.numeric(valid_data[,1]),as.numeric(stackedU_myPred), grid=TRUE, plot=T) 
Uncon_auc = Uncon_roc$auc
# Area under the curve: 0.959