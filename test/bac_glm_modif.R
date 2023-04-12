BAC_GLM = function(effects,train, maxiter = 100, eps = 0.001, coeff.init=NULL, family='poisson'){
  
  coeffs = list()
  new_coeffs = list()
  for (i in 1:length(effects)){
    if(is.null(coeff.init))
    {
      coeffs[[i]] = rep(0,ncol(effects[[i]]))
    }
    else
    {
      coeffs[[i]] = coeff.init[[i]]
    }

    new_coeffs[[i]] = rep(0,ncol(effects[[i]]))
  }
  
  new_coeffs_stock <- as.array(new_coeffs)
  
  j=1
  while(j<=maxiter){
    print(paste("Iteration:",j))
    norm_coeffs = list()
    res = list()
    for (k in sample(1:length(effects)) ){
      
      offset = rep(0,nrow(effects[[1]]))
      for (l in seq(1,length(effects))[-k]){
        offset = offset + as.matrix(effects[[l]])%*%coeffs[[l]]
      }
      
      g <- glm(event~effects[[k]]-1, data=train, family=family, 
          offset=offset)

      new_coeffs[[k]] = as.vector(g$coefficients)
      print(new_coeffs)
      res[[k]] = sqrt(sum((new_coeffs[[k]] - coeffs[[k]])^2)) 
      norm_coeffs[[k]] = sqrt(sum((coeffs[[k]])^2)) 
      
      
    }
    
    coeffs = new_coeffs
    
    print(unlist(res)>unlist(norm_coeffs)*eps)
    if(sum(unlist(res)>unlist(norm_coeffs)*eps) == 0){
      break
    }
    lambda = rep(0,nrow(effects[[1]]))
    for (h in 1:length(effects)){
      lambda = lambda  + (effects[[h]]%*%coeffs[[h]])
    }
    lambda = as.vector(exp(lambda))
    print(paste("#############################", mean(lambda)))
    
    j = j+1
    if(j==maxiter){print("Warning: reached maximum number of iterations")}
    
  }
  
  
  lambda = rep(0,nrow(effects[[1]]))
  for (h in 1:length(effects)){
    lambda = lambda  + (effects[[h]]%*%coeffs[[h]])
  }
  lambda = as.vector(exp(lambda))
  
  return(list("fullcoef"=coeffs,"lambdafit"=lambda,"iterations"=j-1, "new_coeffs_stock"=new_coeffs_stock))
  
}





