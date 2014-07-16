
distance <-function ( x, xk, p=2 ){
  
  nk<-dim(xk)[1]
  
  m<-dim(xk)[2]
  
  dist=numeric(nk) 
  
  dist<-rep(0,nk) 
  
  for( i in 1: nk) {
    
    ssqd <- 0
    
    for( j in 1:m){
      
      ssqd<-ssqd+(x[j]-xk[i,j])^p}
    
    dist[i]<-sqrt(ssqd)
  }
  
  return(dist)
  
}

