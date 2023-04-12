##### Lambda Sim #####
lambda_sim = function(t, step=1){
  lambda = 100*exp(dnorm( (t*step) %% 1440,9*60,1*60) + dnorm((t*step) %% 1440,16*60,2*60) + 5*dnorm((t*step) %% (1440*5),1440*0.5,12*60) - 5*dnorm((t*step) %% (1440*5),1440*4.5,12*60) - 1.5)
  lambda = (lambda-min(lambda))*max(lambda)/(max(lambda)-min(lambda))
  return(lambda)
}

thin_sim = function(lambda, lambda_plus, T_end, seed=24){
  # Simulate the homogenous process
  set.seed(seed)
  m = round(3*T_end*lambda_plus) # Taking a number of events much larger than the expected value in order to ensure we overshoot the boundary
  u = runif(m,0,1)
  t = -1/lambda_plus*log(u) # inter arrival times
  s = cumsum(t) # arrival times
  s = s[s<=T_end] # trimming events that happen after then end of our window
  nstar = length(s)
  
  # Selecting which event we are keeping with the appropriate probability
  w = runif(nstar,0,1)
  Ind = (w <= lambda[ceiling(s)]/lambda_plus)
  
  # Returning the arrival times of the non homogenous poisson process
  return(round(s[Ind]))
  
}

valid_nhpp = function(df_all,lambda,step="15 mins"){
  simu = thin_sim(lambda, 1,nrow(df_all))
  df_all = df_all %>% mutate(time = seq(1,nrow(df_all)))
  final = right_join(tibble(time=simu,pred_event=1),df_all,by="time") %>%
    mutate(pred_event = ifelse(is.na(pred_event),0,pred_event),
           bucket = ceiling_date(date,step)) %>%
    group_by(bucket) %>%
    summarise(event = sum(event),pred_event = sum(pred_event))
  return(tibble(rmse=rmse(final$event,final$pred_event),
                mae=mae(final$event,final$pred_event),
                res_mean = mean(final$pred_event - final$event),
                res_var = var(final$pred_event - final$event),
  ))
}

mae = function(y_true,y_pred){
  return(mean(abs(y_true-y_pred)))
}

rmse = function(y_true,y_pred){
  return(sqrt(mean((y_true-y_pred)**2)))
}

##### Algo 0 - Fit it all at once #####

AAO = function(covariates,posE,pen=0){
  
  minuslogl = function(beta)
  {
    mllikpois = -sum(as.matrix(covariates[posE,])%*%beta) + sum(exp(as.matrix(covariates)%*%beta )) + pen*sum(abs(beta))
    return(as.double(mllikpois))
  }
  
  out = nlminb(start=rep(0,ncol(covariates)), objective=minuslogl)
  # out = nlm(f=minuslogl, p=rep(0,ncol(covariates)))
  # out = optim(fn=minuslogl, par=rep(0,ncol(covariates)), method="BFGS")
  
  return(list("fullcoef"=out$par,"lambdafit"=as.vector(exp(covariates %*% out$par))))
  # return(list("fullcoef"=out$estimate,"lambdafit"=as.vector(exp(covariates %*% out$estimate))))
}


##### Algo 1 - Backfitting #####

BAC = function(effects,posE, maxiter = 100, eps = 0.001){

  coeffs = list()
  new_coeffs = list()
  for (i in 1:length(effects)){
    coeffs[[i]] = rep(0,ncol(effects[[i]]))
    new_coeffs[[i]] = rep(0,ncol(effects[[i]]))
  }
  
  j=1
  while(j<=maxiter){
    print(paste("Iteration:",j))
    norm_coeffs = list()
    res = list()
    for (k in sample(1:length(effects)) ){
      
      
      minuslogl = function(beta)
      {
        first_term = as.matrix(effects[[k]][posE,])%*%beta
        second_term = as.matrix(effects[[k]])%*%beta 
        for (l in seq(1,length(effects))[-k]){
          first_term = c(first_term,as.matrix(effects[[l]][posE,])%*%coeffs[[l]])
          second_term = c(second_term,as.matrix(effects[[l]])%*%coeffs[[l]] )
        }
        mllikpois = -sum(first_term) + sum(exp(second_term))
        return(as.double(mllikpois))
      }
      
      new_coeffs[[k]] = nlminb(start=coeffs[[k]], objective=minuslogl)$par
      res[[k]] = sqrt(sum((new_coeffs[[k]] - coeffs[[k]])^2)) 
      norm_coeffs[[k]] = sqrt(sum((coeffs[[k]])^2)) 
    }
    
    coeffs = new_coeffs
    
    print(unlist(res)>unlist(norm_coeffs)*eps)
    if(sum(unlist(res)>unlist(norm_coeffs)*eps) == 0){
      break
    }
    
    j = j+1
    if(j==maxiter){print("Warning: reached maximum number of iterations")}
  }
  
  
  lambda = rep(0,nrow(effects[[1]]))
  for (h in 1:length(effects)){
    lambda = lambda  + (effects[[h]]%*%coeffs[[h]])
  }
  lambda = as.vector(exp(lambda))
  
  return(list("fullcoef"=coeffs,"lambdafit"=lambda))
  
}

BAC_GLM = function(effects,train, maxiter = 100, eps = 0.001){
  
  coeffs = list()
  new_coeffs = list()
  for (i in 1:length(effects)){
    coeffs[[i]] = rep(0,ncol(effects[[i]]))
    new_coeffs[[i]] = rep(0,ncol(effects[[i]]))
  }
  
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
      
      new_coeffs[[k]] = as.vector(glm(event~effects[[k]]-1, data=train, family='poisson', 
                            offset=offset)$coefficients)
      print(new_coeffs)
      res[[k]] = sqrt(sum((new_coeffs[[k]] - coeffs[[k]])^2)) 
      norm_coeffs[[k]] = sqrt(sum((coeffs[[k]])^2)) 
    }
    
    coeffs = new_coeffs
    
    print(unlist(res)>unlist(norm_coeffs)*eps)
    if(sum(unlist(res)>unlist(norm_coeffs)*eps) == 0){
      break
    }
    
    j = j+1
    if(j==maxiter){print("Warning: reached maximum number of iterations")}
  }
  
  
  lambda = rep(0,nrow(effects[[1]]))
  for (h in 1:length(effects)){
    lambda = lambda  + (effects[[h]]%*%coeffs[[h]])
  }
  lambda = as.vector(exp(lambda))
  
  return(list("fullcoef"=coeffs,"lambdafit"=lambda))
    
}

BAC_GLMnet = function(effects,train, maxiter = 100, eps = 0.001){
  
  coeffs = list()
  new_coeffs = list()
  cvfits = list()
  for (i in 1:length(effects)){
    coeffs[[i]] = rep(0,ncol(effects[[i]]))
    new_coeffs[[i]] = rep(0,ncol(effects[[i]]))
    cvfits[[i]] = 0 
  }
  
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
      
      cvfit = cv.glmnet(x=effects[[k]], y=train$event, nfolds = 10, family = "poisson", 
                        offset=offset,
                        intercept = F,
                        type.measure="deviance", lambda=exp(seq(log(1000), log(0.001), length.out = 20)))
      
      cvfits[[k]] = cvfit
      
      # glm_fit = glmnet(x=effects[[k]], y=train$event, family='poisson',
      #                  intercept = F,
      #                  offset=offset,lambda=cvfit$lambda.min)
      # new_coeffs[[k]] = as.vector(glm_fit$beta)
      
      new_coeffs[[k]] = as.vector(coef(cvfit, s="lambda.min")[-1])
      print(new_coeffs)
      res[[k]] = sqrt(sum((new_coeffs[[k]] - coeffs[[k]])^2)) 
      norm_coeffs[[k]] = sqrt(sum((coeffs[[k]])^2)) 
    }
    
    coeffs = new_coeffs
    
    print(unlist(res)>unlist(norm_coeffs)*eps)
    if(sum(unlist(res)>unlist(norm_coeffs)*eps) == 0){
      break
    }
    
    j = j+1
    if(j==maxiter){print("Warning: reached maximum number of iterations")}
  }
  
  
  lambda = rep(0,nrow(effects[[1]]))
  for (h in 1:length(effects)){
    lambda = lambda  + (effects[[h]]%*%coeffs[[h]])
  }
  lambda = as.vector(exp(lambda))
  
  return(list("fullcoef"=coeffs,"lambdafit"=lambda,"cvfits"=cvfits))
  
}

##### Algo 2 - Backfitting with order of updates based on correlations #####


##### Algo 3 - Backfitting with LARS #####


##### Algo 4 - Boosting #####

# LCP = function(train){
#   
# }
# 
# 
# GBA = function(){
#   
# }

##### END