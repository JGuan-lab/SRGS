# SRGS </br> 

建模：针对基因调控网络构建中基因耦合关系复杂等复杂动态系统建模难题，提出了复杂耦合关联系统的建模方法，建立了融合数据和知识的复杂动态网络系统模型，提升了基因调控网络动态系统模型的准确性和鲁棒性。

## 1. Introduction

SRGS, SPLS (sparse partial least squares)-based recursive gene selection, can be used for gene regulatory network inference from bulk or single-cell expression data. Based on SPLS, SRGS can achieve the purpose of regression and feature selection simultaneously. SRGS recursively selects and scores the genes which may have regulations on the considered target gene. To consider the characteristic of single-cell data, we randomly scramble samples, set some expression values to zeroes, and generate multiple copies of data through multiple iterations, making SRGS more robust.

SRGS corresponds to the following paper:
Guan, J., Wang, Y., Wang, Y. et al. SRGS: sparse partial least squares-based recursive gene selection for gene regulatory network inference. BMC Genomics 23, 782 (2022). https://doi.org/10.1186/s12864-022-09020-7
  
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

The SRGS() function takes a gene expression dataframe or matrix data as input.
Each row of the matrix must correspond to a sample and each column must correspond to a gene.
The gene names must be specified.

      data<-read.table('/data/Simulate/Without dropout/Size50/data/sample500/Ecoli1.txt',sep = ',')
      net<-read.table('/data/Simulate/Without dropout/Size50/gold standards/InSilicoSize50-Ecoli1_goldstandard.tsv')

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
   
    library(SRGS)
    
    data<-read.table('/data/Simulate/Without dropout/Size50/data/sample500/Ecoli1.txt',sep = ',')
    net<-read.table('/data/Simulate/Without dropout/Size50/gold standards/InSilicoSize50-Ecoli1_goldstandard.tsv')
    
    predictNetwork <- SRGS(data, FALSE, iter = 10, k = 1, stepsize = 0.01, num.cores=1) # when SRGS is tested on a dense expression data, the parameter 'dropout' could be set as FALSE.
    ### If one would like to run in parallel, SRGS can provide parallel computing by using the function of parSRGS instead ###
    predictNetwork <- parSRGS(data, FALSE, iter = 10, k = 1, stepsize = 0.01, num.cores = 5)
    
    > predictNetwork[1:5,]
      regulatoryGene targetGene     weight
    1             G2         G1 0.04123711
    2             G3         G1 0.76288660
    3             G4         G1 0.76288660
    4             G5         G1 0.21649485
    5             G6         G1 0.43298969
The result contains the links. Each row corresponds to a regulatory link. The first column shows the regulator, the second column shows the target gene, and the last column indicates the weight of the link.

### 3.2.2 Calculate AUC 

    AUCROC<-calcROC(n,net)
    
    > AUCROC
    Area under the curve: 0.8051
    
### 3.3.1 Run SRGS with the dropout data
  
    library(SRGS)
    
    data<-read.table('/data/Simulate/With dropout/Size50/data/50/Ecoli1.txt',sep = ',')
    net<-read.table('/data/Simulate/With dropout/Size50/gold standards/InSilicoSize50-Ecoli1_goldstandard.tsv')
    
    set.seed(30)
    predictNetwork <- SRGS(data, TRUE, 0.6, iter = 10, k = 1, stepsize = 0.01, num.cores=1) # when SRGS is tested on a sparse expression data, the parameter 'dropout' should be set as TRUE.
    
    ### If one would like to run in parallel, SRGS can provide parallel computing by using the function of parSRGS instead ###
    predictNetwork <- parSRGS(data, TRUE, iter = 10, k = 1, stepsize = 0.01, num.cores = 5)
    
    > predictNetwork[1:5,]
      regulatoryGene targetGene     weight
    1             G2         G1 0.04123711
    2             G3         G1 0.76288660
    3             G4         G1 0.76288660
    4             G5         G1 0.21649485
    5             G6         G1 0.43298969
The result contains the links. Each row corresponds to a regulatory link. The first column shows the regulator, the second column shows the target gene, and the last column indicates the weight of the link.   

### 3.3.2 Calculate AUC 
  
    AUCROC<-calcROC(n,net)
    
    > AUCROC
    Area under the curve: 0.6639
