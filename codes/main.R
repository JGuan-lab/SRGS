library("spls")
library("reshape2")
library("parallel")
library("pROC")

source("calcROC.R")
source("SRGS.R")
source("parSRGS.R")

######   Run SRGS with the nodropout data   ######
######   Prepare data   ######
data<-read.table('/data/Simulate/Without\ dropout/Size50/data/sample500/Ecoli1.txt',sep = ',')
net<-read.table('/data/Simulate/Without\ dropout/Size50/gold\ standards/InSilicoSize50-Ecoli1_goldstandard.tsv')

######   Run SRGS with the nodropout data   ######
predictNetwork <- SRGS(data, FALSE, iter = 10, k = 1, stepsize = 0.01)
######   run with parallel   ######
predictNetwork <- parSRGS(data, FALSE, iter = 10, k = 1, stepsize = 0.01, num.cores = 5)
######   Calculate ROC   ######
AUCROC<-calcROC(predictNetwork,net)


######   Run SRGS with the dropout data   ######
######   Prepare data   ######
data<-read.table('/data/Simulate/With\ dropout/Size50/data/50/Ecoli1.txt',sep = ',')
net<-read.table('/data/Simulate/With\ dropout/Size50/gold\ standards/InSilicoSize50-Ecoli1_goldstandard.tsv')

######   Run SRGS with the dropout data   ######
set.seed(30)
predictNetwork <- SRGS(data, TRUE, p = 0.6, iter = 10, k = 1, stepsize = 0.01)
######   run with parallel   ######
predictNetwork <- parSRGS(data, TRUE, p = 0.6,iter = 10, k = 1, stepsize = 0.01, num.cores = 6)
######   Calculate ROC   ######
AUCROC<-calcROC(predictNetwork,net)







