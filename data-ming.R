setwd("~/Desktop/Data Mining/project")
library(kernlab)
library(e1071)
library(ROCR)
library(performanceEstimation)
install.packages("caret")
library(caret)
library(help="caret")
install.packages("caret", dependencies = c("Depends", "Suggests"))
install.packages("doMPI")
install.packages("Rmpi")
install.packages("mlbench")
library(mlbench)
install.packages("AppliedPredictiveModeling")
library(AppliedPredictiveModeling)
#######
install.packages("rsubgroup")
library(rsubgroup)
data(credit.data)
#load("~/Desktop/Data Mining/project/TREE.RData")
setwd("~/Desktop/Data Mining/project")
GER=read.csv("german_credit.csv")
attach(GER)
factor_vars <- c("Creditability","Account.Balance","Payment.Status.of.Previous.Credit","Purpose","Value.Savings.Stocks","Length.of.current.employment",
                 "Sex...Marital.Status","Guarantors","Most.valuable.available.asset","Concurrent.Credits","Type.of.apartment",
                 "Occupation","Telephone","Foreign.Worker")
GER[factor_vars] <- lapply(GER[factor_vars], function(x) as.factor(x))
sapply(GER, class)
head(GER)

#response variable distribution
y <- GER$Creditability
cbind(freq=table(y), percentage=prop.table(table(y))*100)
##no NAs
install.packages("Amelia")
library(Amelia)
missmap(GER, col=c("black", "grey"), legend=FALSE)
##Scatter plot Matrix By Class
pairs(Creditability~., data=GER[,1:5], col=GER$Creditability)
pairs(Creditability~., data=GER[,c(1,6:10)], col=GER$Creditability)
pairs(Creditability~., data=GER[,c(1,11:15)], col=GER$Creditability)
pairs(Creditability~., data=GER[,c(1,16:21)], col=GER$Creditability)




cost=c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000)

library(caret)
install.packages('klaR')
library(klaR)
set.seed(1234)
trainIndex=sample(1:1000,800,replace=FALSE)
#trainIndex <- createDataPartition(GER$Creditability, p=0.80, list=FALSE)
trainset <- GER[ trainIndex,]
testset <- GER[-trainIndex,]
#set.seed(123456)
#V = 10



library(caret)
library(mlbench)

set.seed(7)
#trainControl <- trainControl(method="repeatedcv", number=10,repeats = 1,search="grid")
#grid <- expand.grid(.C=10**seq(-4, 4, by=1))                    
#fit.svmRadial <- train(Creditability~., data=trainset, method="svmLinear",
                       #metric="Accuracy",tuneLength=10, trControl=trainControl,tuneGrid=grid)


###!!!We find that when C is large, the accuracy is low, so we avoid c>=1.
trainControl <- trainControl(method="repeatedcv", number=10,repeats = 1,search="grid")
grid <- expand.grid(.C=10**seq(-4, 0, by=1))                    
fit.svmRadial <- train(Creditability~., data=trainset, method="svmLinear",
                       preProc=c("center", "scale", "BoxCox"),metric="Accuracy",tuneLength=10, trControl=trainControl,tuneGrid=grid)

# summarize fit
print(fit.svmRadial)
plot(fit.svmRadial)
names(getModelInfo())

sapply(trainset,class)

min(trainset$No.of.dependents)
index = split(sample(dim(trainset)[1]), rep(1:V, length = dim(trainset)[1]))
library(performanceEstimation) # Loading our infra-structure
library(e1071)
res <- performanceEstimation(
  PredTask(Creditability ~ .,data=trainset),
  workflowVariants(learner="svm",learner.pars=list(cost=0.0001)),
  EstimationTask(metrics="err",method=CV()))


plot(res)



score=function(kernel,kpar)
{
  errors=array(NA,dim=list(5,2,9))
  for (i in 1:9){
    for (ii in 1:5){
      model.linear1 = ksvm(V15~.,  data = trainset[-index[[ii]], ],kernel =kernel,kpar=kpar,C= cost[i])
      predict.train = predict(model.linear1,newdata =trainset[-index[[ii]],1:24 ])
      predict.test = predict(model.linear1,newdata =trainset[index[[ii]],1:24 ])
      errors[ii,1,i]=mean(predict.train!=trainset[-index[[ii]],25])
      errors[ii,2,i]=mean(predict.test!=trainset[index[[ii]],25])
    }
  }
  meanerror=matrix(NA,nrow=9,ncol=2)
  for (k in 1:9){
    meanerror[k,]=  apply(errors[,,k],2,mean)
  }
  meanerror
}

linearscore=score(kernel="vanilladot",kpar=list())
polydf2=score(kernel="polydot",kpar=list(degree=2))
polydf3=score(kernel="polydot",kpar=list(degree=3))
rbf0.1=score(kernel = "radial",kpar = list(sigma =0.1))
rbf1=score(kernel = "radial",kpar = list(sigma =1))
rbf10=score(kernel = "radial",kpar = list(sigma =10))
