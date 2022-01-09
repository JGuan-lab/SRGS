# SRGS </br> 
## 1. Introduction  

  
## 2. Installation
Depends: 

       R (>= 3.6.3)   

Requirements: 

      library("spls")
      library("reshape2")
      library("parallel")
      library("pROC")
      
## 3. Quick start

Run `main.R`. The parameters can be changed as below.

### 3.1 Prepare data
The SGRS() function takes as input argument a gene expression dataframe or matrix data.
Each row of that matrix must correspond to a sample and each column must correspond to a gene.
The gene names must be specified in colnames(exprMatr).

      data<-read.table('/data/Simulate/Without\ dropout/Size50/data/sample500/Ecoli1.txt',sep = ',')
      net<-read.table('/data/Simulate/Without\ dropout/Size50/gold\ standards/InSilicoSize50-Ecoli1_goldstandard.tsv')

      > data[1:5,1:5]
            G1        G2        G3        G4        G5
      1 0.6232604 0.7038660 0.5866520 0.7534090 0.6247685
      2 0.6607892 0.4835473 0.2811439 0.5049232 0.7115181
      3 0.5951789 0.6791481 0.2123983 0.6271043 0.7488240
      4 0.6397211 0.5582082 0.1787238 0.6244249 0.5817057
      5 0.6271739 0.4982657 0.3603228 0.5544687 0.7246676
      
      > net[1:5,]
        V1 V2 V3
      1 G1 G2  1
      2 G3 G2  1
      3 G4 G2  1
      4 G4 G5  1
      5 G4 G6  1
      
### 3.2.1 Run SRGS with the nodropout data
The resulting linkList matrix contains the ranking of links. Each row corresponds to a regulatory link. The first column shows the regulator, the second column shows the target gene, and the last column indicates the weight of the link.
    
    library(SRGS)
    
    data<-read.table('/data/Simulate/Without\ dropout/Size50/data/sample500/Ecoli1.txt',sep = ',')
    net<-read.table('/data/Simulate/Without\ dropout/Size50/gold\ standards/InSilicoSize50-Ecoli1_goldstandard.tsv')
    predictNetwork <- SRGS(data, FALSE, iter = 10, k = 1, stepsize = 0.01, num.cores=1)
    ###run with parallel###
    predictNetwork <- parSRGS(data, FALSE, iter = 10, k = 1, stepsize = 0.01, num.cores = 5)
    > predictNetwork[1:5,]
      regulatoryGene targetGene     weight
    1             G2         G1 0.04123711
    2             G3         G1 0.76288660
    3             G4         G1 0.76288660
    4             G5         G1 0.21649485
    5             G6         G1 0.43298969
    
### 3.2.2 Calculate ROC 

    AUCROC<-calcROC(n,net)
    
    > AUCROC
    Area under the curve: 0.8051
    
### 3.3.1 Run SRGS with the dropout data
The resulting linkList matrix contains the ranking of links. Each row corresponds to a regulatory link. The first column shows the regulator, the second column shows the target gene, and the last column indicates the weight of the link.
    
    library(SRGS)
    
    data<-read.table('/data/Simulate/With\ dropout/Size50/data/50/Ecoli1.txt',sep = ',')
    net<-read.table('/data/Simulate/With\ dropout/Size50/gold\ standards/InSilicoSize50-Ecoli1_goldstandard.tsv')
    set.seed(30)
    predictNetwork <- SRGS(data, TRUE, 0.6, iter = 10, k = 1, stepsize = 0.01, num.cores=1)
    ###run with parallel###
    predictNetwork <- parSRGS(data, TRUE, iter = 10, k = 1, stepsize = 0.01, num.cores = 5)
    > predictNetwork[1:5,]
      regulatoryGene targetGene     weight
    1             G2         G1 0.04123711
    2             G3         G1 0.76288660
    3             G4         G1 0.76288660
    4             G5         G1 0.21649485
    5             G6         G1 0.43298969
    
### 3.2.2 Calculate ROC 
  
    AUCROC<-calcROC(n,net)
    
    > AUCROC
    Area under the curve: 0.6639
