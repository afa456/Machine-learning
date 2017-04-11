rm(list=ls())
# Read in labelname and response data
myXData=read.csv("texture.csv",header=TRUE)
myYData=read.csv("textureResponse.csv",header=TRUE)
# Restrict attention to particular variables

myLogical=sapply(attributes(myXData)$names,function(varname){
	if(varname %in% c("PatID","TimePoint","PFS")){
		return(TRUE)
	}else{
		return(((length(grep(pattern="v1_",varname))>0))&
			  (!((length(grep(pattern="Xa",varname))>0)|
			(length(grep(pattern="DeltaX",varname))>0))))
	}
})

myXData=myXData[,myLogical]
# Rearrange data

myXtemp=data.matrix(myXData[,6:dim(myXData)[2]])

patientIDsTemp=sort(unique(c(myXData$PatID,myYData$Pat.No)))

timePointsTemp=sort(unique(as.character(myXData$TimePoint)))

PFS=matrix(NA,length(patientIDsTemp),1)

predictors=matrix(NA,length(patientIDsTemp),(dim(myXtemp)[2])*length(timePointsTemp))

counter=1

for(pt in patientIDsTemp){
	index=which(myYData$Pat.No==pt)
	if(length(index)==1){
		PFS[counter]=log(myYData$PFS[index])
}
for(jj in 1:length(timePointsTemp)){
	index=which((myXData$TimePoint==timePointsTemp[jj])&
				(myXData$PatID==pt))
	if(length(index)==1){
		predictors[counter,((jj-1)*dim(myXtemp)[2]+1):(jj*dim(myXtemp)[2])]=
		myXtemp[index,]
	}
}
counter=counter+1
}
predictorsStar=apply(predictors,2,function(xx){
	if(sum(!is.na(xx))>=2){
		return((xx-mean(xx,na.rm=TRUE))/sd(xx,na.rm=TRUE))
	}else{
		return(xx)
	}
})
myData=data.frame(PFS,predictorsStar)
myData$PFS=myData$PFS-mean(myData$PFS,na.rm=TRUE)
# Label columns
myNames=c("PFS",
	paste(attributes(myXData)$names[6:dim(myXData)[2]],"_baseline1",sep=""),
	paste(attributes(myXData)$names[6:dim(myXData)[2]],"_baseline2",sep=""),
	paste(attributes(myXData)$names[6:dim(myXData)[2]],"_followUp",sep=""))
attributes(myData)$names=myNames


############Data Manipulation##############

# Remove completely missing columns.
missCol=logical(length = dim(myData)[2])
for(ii in 1:dim(myData)[2]){
	if(sum(abs(myData[,ii])>matrix(0.0001,dim(myData)[1],1), na.rm = TRUE)>0){
	missCol[ii] = T
	}else{
	missCol[ii] = F
	}
}
myData1=myData[,missCol]

# Remove rows with completely missing predictors.
missRow = logical(length = dim(myData1)[1])
for(ii in 1:dim(myData1)[1]){
	if(sum(abs(myData1[ii,2:dim(myData1)[2]])>matrix(0.0001,1,dim(myData1)[2]), na.rm = TRUE)>0){
	missRow[ii] = T
	}else{
	missRow[ii] = F
	}
}
myData2=myData1[missRow,]

# Remove columns with constant predictors.
myData3=myData2[,apply(myData2,2,function(xx){
	if(var(xx,na.rm=TRUE)<0.0001){
	return(F)
	}else{
	return(T)}
	}
)]

# Remove (second or later occurrance of) columns which are collinear (correlation 1) with another.
correlation = cor(myData3[,-1], use ="pairwise.complete.obs")
temp = matrix(data = FALSE, nrow = dim(myData3)[2]-1, ncol = dim(myData3)[2]-1)
for(ii in 1:(dim(myData3)[2]-1)){
	for(jj in ii:(dim(myData3)[2]-1)){
		if (correlation[ii,jj]>0.9999 & ii!=jj){
			temp[ii,jj] = T
		}
	}
}
corrtemp = colSums(temp)
corrtemp = (corrtemp==0)
corrtemp = append(corrtemp,1,after = 0)
myData4 = myData3[,as.logical(corrtemp)]

############ Multiple Imputation ##############
newmyData=myData4
NA_logical =logical(length = dim(newmyData)[2])
for(i in 1:dim(newmyData)[2]){
  NA_logical[i]=any(which(is.na(newmyData[,i])))
}  
  tempData = newmyData
  for (i in 1:dim(newmyData)[2]){
     if(NA_logical[i] == T){
   	   row_index = which((is.na(tempData[,i])) == T)
       row_No_NA = which((is.na(tempData[,i])) == F)
       tempData[row_index,i] = sample(tempData[row_No_NA,i],size=length(row_index),replace=T)
	   }
  }

  for (i in 1:dim(newmyData)[2]) {
    if(NA_logical[i] == T){
	    row_NA = which((is.na(newmyData[,i])) == T)
      lineareg = tempData[-row_NA,i]
      mini=lm(lineareg~1,data = tempData[-row_NA,-i])
	    maxi=lm(lineareg~.,data = tempData[-row_NA,-i])
	    selectvar = step(mini,scope = list(upper=maxi,lower = mini),direction="forward",steps = 10,k = 2)
	    coefficient = coef(summary(selectvar))
      currentList = names(selectvar$coefficients)
	    regdata = c()
	    for(j in 1:length(currentList)){
	      if(currentList[j] == "(Intercept)"){
	         regdata = cbind(regdata,tempData[,i])
	      }else{
	         regdata = cbind(regdata,tempData[,which(attributes(tempData)$name == currentList[j])])
	      }
	    }
      RegData = regdata %*% as.matrix(coefficient[,1])
      noisyplus = RegData[row_NA] + rnorm(length(row_NA),0,var(coefficient[,2])) 
	    for(j in 1:length(row_NA)){
	       indexfinal = which.min(abs(RegData[-row_NA]-noisyplus[j]))
         newmyData[row_NA[j],i] = newmyData[-row_NA[j],i][indexfinal]
	    }
  	}
  }	
############ After Imputation ##############
Data=newmyData
FUlabel=grep("followUp",attributes(Data)$names)
BL1label=grep("baseline1",attributes(Data)$names)
BL2label=grep("baseline2",attributes(Data)$names)

avg_dif=Data[1]
for (i in 1:75){
	labelname=attributes(Data)$names[FUlabel[i]]
	labelname1=substr(labelname, 1,(nchar(labelname)-8))
	indexbs1=BL1label[which(grepl(labelname1,attributes(Data)$names[BL1label]))]
	indexbs2=BL2label[which(grepl(labelname1,attributes(Data)$names[BL2label]))]
	a=b=0
	if(length(indexbs1)>0){
		a=Data[,indexbs1]
	}
	if(length(indexbs2)>0){
		b=Data[,indexbs2]
	}
	
	avg=paste(labelname1,"_Average",sep="")
	diff=paste(labelname1,"_Diff",sep="")
	if(length(a)>1 && length(b)>1){
		avg_dif[avg]=(a+b)/2
		avg_dif[diff]=Data[,FUlabel[i]]-(a+b)/2
	}
	if(length(a)>1 && length(b)==1){
		avg_dif[avg]=a
		avg_dif[diff]=Data[,FUlabel[i]]-a
	}
	if(length(a)==1 && length(b)>1){
	avg_dif[avg]=b
	avg_dif[diff]=Data[,FUlabel[i]]-b
	}
}
write.csv(avg_dif,"myData(manipulate)-1.csv")

