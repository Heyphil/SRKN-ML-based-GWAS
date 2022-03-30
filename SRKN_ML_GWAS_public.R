library(readxl)
library(readr)
library(dplyr)
library(pls)
library(plsVarSel)
library(ggplot2)
library (e1071)
library(randomForest)
library(neuralnet)
library(gbm)
library(caret)
library(mltools)

#----- prepare the dataset---------------------
# set the working dictionary
setwd("D:/")
# read the SNP dataset
data.dummy <- read.delim("genotype_numerical.txt")
# read the dependent dataset 
data.score = read_csv('SRKN_Galls1721.csv') %>% data.frame()

# remove those only have one class 
data.factor = as.matrix(data.dummy[,-1]) %>% as.factor()
View(data.factor)
nonsenseInd = NULL
for (i in 2:6001) {
  temp = data.dummy[,i]
  if (sd(temp) == 0) {
    nonsenseInd = c(nonsenseInd,i)
  }
}
data.plsr = data.dummy[,-nonsenseInd]

# put the SNP and dependent data togather

for (i in 1:nrow(data.plsr)) {
  name = data.plsr$taxa[i]
  ind = data.score$Check[data.score$Name == name]
  if (length(ind) != 0){
    data.plsr$score[i] = ind}
}

#---------- PLSR for variable importance ---------
data.plsr = data.plsr[,-1]
set.seed(1)
pls.fit = plsr(score~.,data = data.plsr,scale = FALSE,validation = "CV")
#summary(pls.fit)

# calculate the variable importance with VIP
comp <- which.min(pls.fit$validation$PRESS) # optimal number of components of PLS model
vip <- VIP(pls.fit, comp) %>% data.frame()# variable importance in projection
vipname = row.names(vip)
vip = vip[grepl('Gm',vipname),]
vipname = vipname[grepl('Gm',vipname)]
vip = cbind.data.frame('SNP' = vipname,'vip' = vip)
plot(vip)
write.csv(vip, file = "PLSRVIP.csv")

#---------- Remove variables with high correlation ---------
# correlation between variables
var.correlation = round(cor(data.plsr[,-which(names(data.plsr) %in% c('taxa','score'))]),2)

# remove high VIP with correlations
vip <- read_csv("PLSRVIP.csv") %>% data.frame()
vip = vip[order(-vip$vip),]
newvip = vip
vipname = newvip$SNP

for (i in 1:nrow(newvip)){
  
  name = vipname[i]
  if (is.na(name)){
    break
  }
  if (grepl('Gm',name)){
    col.ind = grep(name, colnames(var.correlation))
    col.corr = var.correlation[,col.ind]
    
    col.corr = sort(col.corr,decreasing = TRUE) %>% as.data.frame()
    
    snpList = row.names(col.corr)
    
    corr.ind = (abs(col.corr) > 0.5)
    
    if (length(which(corr.ind)) > 1) {
      col.corr = col.corr[corr.ind,] %>% as.data.frame()
      snpList = snpList[corr.ind]
      removeInd = which(newvip$SNP %in% snpList[-1])
      if (sum(removeInd) != 0){
        newvip = newvip[-removeInd,]
        vipname = newvip$SNP
      }
    }
  }
}

write.csv(newvip, file = "VIP_removeCorrelated.csv")


###-------------- visulize the VIP scores
vip <- read_csv("PLSRVIP.csv") %>% data.frame()
vip = vip[order(vip$SNP),]
vip$num = seq(1,nrow(vip))

VIP_removeCorrelated <- read_csv("VIP_removeCorrelated.csv")
View(VIP_removeCorrelated)

VIP_selected = as.data.frame(VIP_removeCorrelated[VIP_removeCorrelated$vip>2,])
num = matrix(0,nrow = nrow(VIP_selected),ncol = 1)
for (i in 1:nrow(VIP_selected)){
  snp = VIP_selected$SNP[i]
  vip.ind = grep(snp, vip$SNP)
  VIP_selected$num[i] = vip.ind
}

basic = ggplot() +
  geom_point(data = vip, aes(x = num, y = vip), color = "black") + # must include argument label "data"
  geom_point(data = VIP_selected, aes(x = num, y = vip), color = "red")

basic + theme(
  plot.background = element_rect(fill = "white"), 
  panel.background = element_rect(fill = "white", colour="white")
)

plot(vip$num,vip$vip)
plot(VIP_selected$num,VIP_selected$newvip,col = 'red')


#----------- Stepforward variable selection----------------------------
# put the SNP and dependent data togather
for (i in 1:nrow(data.dummy)) {
  name = data.dummy$taxa[i]
  ind = data.score$Check[data.score$Name == name]
  if (length(ind) != 0){
    data.dummy$score[i] = ind }
}

###------ these are user functions --------------------###
select.models = function(train.set,test.set,fmla,fct.name = 'RF'){
  
  if (fct.name == 'SVM'){
    # using SVM
    set.seed(1)
    tune.out=tune(svm,fmla,data=train.set, kernel ="radial",scale = FALSE, 
                  ranges =list(cost=c(0.01,0.1,1,10,100,1000),
                               gamma=c(0.0001,0.001,0.01,0.5,1)))
    svm.fit = tune.out$best.model
    # prediction accuracy
    ypred= predict(svm.fit,test.set)
    #remove(test.set,train.set,tune.out,svm.fit,ypred,cv.table)
  } else if (fct.name == 'RF'){
    rf.fit =randomForest(fmla,data = train.set)
    ypred = predict(rf.fit,test.set,"response")
    #} else if (fct.name == 'Boosting'){
    #boost.fit =gbm(fmla,data =train.set,
    #               distribution="gaussian",
    #               n.trees = 100,
    #               shrinkage = 0.01)
    #ypred = predict(boost.fit,test.set,n.trees = 100,"response")
    
  } else if (fct.name == 'ANN'){
    #fmla.elements = all.vars(fmla)
    set.seed(1)
    nn.fit = neuralnet(fmla,data = train.set, hidden=5,
                       act.fct = "logistic", linear.output = FALSE,
                       stepmax = 1e+010)
    ypred = predict(nn.fit,test.set,rep = 1,all.units = FALSE)
    data.summary = summary(test.set$class3)
    if (length(data.summary) == 3){
      ypred=matrix(ifelse(max.col(ypred[ ,1:3]) == 1, 1,
                          ifelse(max.col(ypred[ ,1:3]) == 2, 2, 3)))
    } else if (length(data.summary) == 2){
      ypred=matrix(ifelse(max.col(ypred[ ,1:2]) == 1, 1, 2))
    }
  }
  return(ypred)
}

classification.5foldcv = function(data.dummy.select,fmla,s,fct.name){
  
  cv.single = matrix(0,nrow = 5)
  for (j in 1:5) {
    test.set <- data.dummy.select[seq_len(length(s))[(s == j)],] #test data
    train.set <- data.dummy.select[seq_len(length(s))[(s != j)],] #training data
    ypred = select.models(train.set,test.set,fmla,fct.name)
    cv.table = table(predict=ypred, truth=test.set$class3)
    if (nrow(cv.table) == ncol(cv.table)){
      cv.single[j] = sum(diag(cv.table))/length(ypred)
    } else{
      diff = (ypred - as.numeric(test.set$class3)) == 0 
      cv.single[j] = length(diff[diff == 'TRUE'])/length(ypred)
    }
  }
  return(cv.single)
} 

classification.accuracy = function(data.dummy.select,fmla,s,fct.name){
  #set.seed(1)
  cv.single = matrix(0,nrow = 5)
  for (j in 1:5) {
    test.set <- data.dummy.select[seq_len(length(s))[(s == j)],] #test data
    train.set <- data.dummy.select[seq_len(length(s))[(s != j)],] #training data
    ypred = select.models(train.set,test.set,fmla,fct.name)
    #cv.table = table(predict=ypred, truth=test.set$class3)
  }
  return(cbind.data.frame('prediction' = ypred,'truth' = test.set$class3))
} 

creatCVlist = function(max.var,variable.pool){
  names = NULL
  if (is.null(max.var)){  # the first round without preference
    names = variable.pool
  } else {
    for (var1 in max.var){
      for (var2 in variable.pool) {
        if (!grepl(var2,var1)){  # to determine if the variable has been used in previous formula
          name = paste(var1,var2,sep = '+')
          names  = rbind(names,name)
        } 
      }
    }
  }
  
  return(names)
}

roundTrain = function(variable.pool,max.var,data.dummy.select,s,fct.name){
  
  names = creatCVlist(max.var,variable.pool)
  cv.round = matrix(0,nrow = length(names))
  rownames(cv.round) = names
  for (i in 1:length(names)){
    fmla = as.formula(paste("class3 ~ ", paste(names[i])))
    cv.single = classification.5foldcv(data.dummy.select,fmla,s,fct.name)
    cv.round[i] = mean(cv.single)
    remove(cv.single,fmla)
  }
  return(cv.round)
}

roundTrain_gain = function(variable.pool,max.var,data.dummy.select,s,fct.name){
  names = creatCVlist(max.var,variable.pool)
  cv.round = NULL
  acc.gain.tolerance = 0
  acc.gain = 1
  max.acc = 0
  num.SNP = 0
  while (acc.gain > acc.gain.tolerance){
    num.SNP = num.SNP+1    
    fmla = as.formula(paste("class3 ~ ", paste(names[num.SNP])))
    cv.single = classification.5foldcv(data.dummy.select,fmla,s,fct.name)
    cv.round = rbind(cv.round,mean(cv.single))
    acc.gain = max(cv.round) - max.acc
    max.acc = max(cv.round)    
    remove(cv.single,fmla)
  }
  cv.round = as.matrix(cv.round)
  rownames(cv.round) = names[1:num.SNP]
  return(cv.round)
}
###----------------------------------------------------###

# Initiate a few parameters
set.seed(1)
s <- sample(rep(1:5,ceiling(nrow(data.dummy)/5)), nrow(data.dummy)) 
acc.gain.tolerance = 0
acc.gain = 1
max.acc = 0
num.SNP = 1
master.list = vector(mode = "list", length = 0)
max.var = NULL  # can put a preference here, for example, max.var = 'Gm10_1232205_A_G'

# Specific a classification function here. It could be 'SVM', 'ANN', or 'RF'. 
# Note: when choose RF, it's better to specify a preference by, e.g., max.var = 'Gm10_1232205_A_G'. just to help RG model converge. Because there is a bug when RF took single SNPs as predictors. 
cla.function = 'SVM'
# max.var = 'Gm10_1232205_A_G'

# start selection from the variables with high VIP
Variables = read_csv('VIP_removeCorrelated.csv') %>% as.data.frame()
variableNames =Variables$SNP[Variables$vip > 2]
data.dummy.select = data.dummy[variableNames]

data.dummy.select$class3 = 2
data.dummy.select$class3[data.dummy$score>4] = 3
data.dummy.select$class3[data.dummy$score<2] = 1
data.dummy.select$class3 = as.factor(data.dummy.select$class3)
summary(data.dummy.select$class3)

# stop criteria: acc.gain > acc.gain.tolerance
resultTable = NULL
while (num.SNP<=29){
  cv.round = roundTrain(variableNames,max.var,data.dummy.select,s,cla.function)
  master.list = append(master.list,list(cv.round))
  names(master.list)[num.SNP] = paste('Round',num.SNP,sep = "_")
  cv.round = master.list[[paste('Round',num.SNP,sep = "_")]]
  max.ind = which.max(cv.round)
  max.var = row.names(cv.round)[max.ind]
  num.SNP = num.SNP+1
  acc.gain = max(cv.round) - max.acc
  max.acc = max(cv.round)
  remove(cv.round)
  # return the classification matrices
  fmla = as.formula(paste("class3 ~ ", max.var))
  cv.pred = classification.accuracy(data.dummy.select,fmla,s,cla.function)
  cv.table = confusionMatrix(cv.pred$prediction, cv.pred$truth)
  byclass = cv.table[["byClass"]]
  byclass2 = byclass[,c(11,5,2,6)]    
  
  temp = c(cv.table[["overall"]][["Accuracy"]],
           mcc(cv.pred$prediction, cv.pred$truth),
           as.vector(t(byclass2)))
  resultTable = rbind(resultTable,temp)
  remove(fmla,cv.pred, cv.table, byclass, byclass2, temp)
}

write.csv(resultTable, file = "Metrics_SVM_29_MCC.csv")
save(master.list,file = 'classificationACC_ANN_29.RData')
