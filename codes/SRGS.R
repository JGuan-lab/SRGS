SRGS<-function(data, dropout, p = 0.8, iter = 10, k = 1, stepsize = 0.01){
  
  scoreMatrix<-matrix(nrow=0,ncol=ncol(data))
  colnames(scoreMatrix)<-colnames(data)
  scoreMatrix<-data.frame(scoreMatrix)
  
  for ( i in 1:ncol(data))
  {
    
    dataX<-data[,-i]
    dataY<-data[,i]
    Fscore<-data.frame(name=colnames(data),score=0)
    
    res<-lapply(seq(0.01,0.99,stepsize),function(x,dataX,dataY,k,dropout,p,stepsize,iter){
      try({
        if(dropout==FALSE)
        {
          f<- spls(dataX, dataY, K=k, eta=x)
          be1<-data.frame(name=rownames(f$betahat),be = f$betahat)
          be1[be1$be!=0,]$be<-1
          scoreList<-be1$name[which(be1$be!=0)]
          return(scoreList)
        }else{
          for(it in 1:iter)
          {
            ran<-sample(1:nrow(dataX),nrow(dataX))
            dataXOld<-dataX[ran,]
            dataYNew<-dataY[ran]
            a<-matrix(rbinom(nrow(dataXOld)*ncol(dataXOld),1,p),nrow = nrow(dataXOld),ncol = ncol(dataXOld))
            dataXNew<-a*dataXOld
            colnames(dataXNew)<-colnames(dataXOld)
            f<- spls(dataXNew, dataYNew, K=k, eta=x)
            be1<-data.frame(name=rownames(f$betahat),be = f$betahat)
            be1[be1$be!=0,]$be<-1
            scoreList<-be1$name[which(be1$be!=0)]
            if(length(scoreList)>0 && length(scoreList)<ncol(dataXNew)){
              Fscore$score[which(Fscore$name %in% be1$name)]=Fscore$score[which(Fscore$name %in% be1$name)]+be1$be
            }
          }
          return(Fscore)
        }
      },silent = TRUE)
    },dataX=dataX,dataY=dataY,k=k,dropout=dropout,p=p,stepsize=stepsize,iter=iter)
    
    
    
    if(dropout == FALSE){
      for(q in 1:length(res)){
        if(length(res[[q]])>0 && length(res[[q]])<ncol(dataX)){
          Fscore$score[which(Fscore$name %in% res[[q]])]=Fscore$score[which(Fscore$name %in% res[[q]])]+1
        }
      }
    }else{
      for(q in 1:length(res)){
        Fscore$score[which(Fscore$name %in% res[[q]]$name)]=Fscore$score[which(Fscore$name %in% res[[q]]$name)]+res[[q]]$score
      }
    }
    
    
    maxScore<-max(Fscore$score)
    minScore<-min(Fscore[-i,]$score)
    maxMin<-max(Fscore$score)-min(Fscore[-i,]$score)
    Fscore$score<-(Fscore$score-minScore)/maxMin
    Fscore[i,]$score<-0
    
    row<-matrix(Fscore$score,nrow=1,ncol=nrow(Fscore))
    colnames(row)<-Fscore$name
    row<-data.frame(row)
    
    scoreMatrix<-rbind(scoreMatrix,row)
    
  }
  
  scoreMatrix<-data.frame(t(scoreMatrix))
  colnames(scoreMatrix)<-rownames(scoreMatrix)
  
  aa<-melt(as.matrix(scoreMatrix))
  n<-aa[-which(aa$Var1 == aa$Var2),]
  colnames(n) <- c('regulatoryGene', 'targetGene', 'weight')
  return(n)
}
