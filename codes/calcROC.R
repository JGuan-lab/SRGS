calcROC<-function(predNet,trueNet){
  
  m<-data.frame()
  for (i in 1:nrow(trueNet)) {
    
    d<-predNet[intersect(which(predNet[,1]==trueNet[i,1]),which(predNet[,2]==trueNet[i,2])),]
    m <- rbind(m,d)
    
  }
  
  roc1<- roc(trueNet[,3], m[,3],levels = c(0,1))
  return(roc1$auc)
}
