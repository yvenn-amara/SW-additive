##### Plotting Wavelet Basis #####
plot_wav = function(x,numLevels, filterNumber = 1, ylim=c(-5,5), coef_is_null=NULL){
  #vector of zeros and ones of the same size as 2^numlevels
  require("lattice")
  n = (ceiling(max(x))-floor(min(x)))*10
  x_seq=seq(min(x),max(x), length=n)
  Z = WavD(x_seq, numLevels=numLevels, filterNumber = filterNumber)
  B = cbind(rep(1,n),Z)
  
  if(is.null(coef_is_null)==FALSE){
    B = t(t(B)* coef_is_null)
  }
  
  df = data.frame(x = x_seq,B)
  # print(names(df))
  
  form = as.formula(paste(paste(names(df)[- 1],  collapse = ' + '),'x',  sep = '~'))
  return(xyplot(form,data = df,
                type = 'l', grid=TRUE,outer =TRUE,
                ylab="DaubLeAsymm",scales=list(tick.number=10)))
  
}

##### Lambda Sim #####
l_smooth = function(data, step=1){
  lambda = exp(-sqrt(abs(data$temp))-(data$tod/max(data$tod)-1)^2+data$traffic^3/max(data$traffic^3) + abs(data$day-4)/5 +4)
  # lambda = (lambda-min(lambda))*max(lambda)/(max(lambda)-min(lambda))
  return(lambda)
}

l_peaky = function(data, step=1){
  lambda = exp(-data$temp/max(data$temp)-sin(data$tod-mean(data$tod))+4*data$traffic/max(data$traffic) + abs(data$day-4)/5+1)
  # lambda = (lambda-min(lambda))*max(lambda)/(max(lambda)-min(lambda))
  return(lambda)
}

l_both = function(data, step=1){
  lambda = exp(-data$temp/max(data$temp)-sin(data$tod-mean(data$tod))+data$traffic/max(data$traffic) + abs(data$day-4)/5 + 2)
  # lambda = (lambda-min(lambda))*max(lambda)/(max(lambda)-min(lambda))
  return(lambda)
}




lambda_sim = function(t, step=1){
  lambda = exp(dnorm( (t*step) %% 1440,9*60,1*60) + dnorm((t*step) %% 1440,16*60,2*60) + 5*dnorm((t*step) %% (1440*5),1440*0.5,12*60) - 5*dnorm((t*step) %% (1440*5),1440*4.5,12*60) +1.5)
  lambda = (lambda-min(lambda))*max(lambda)/(max(lambda)-min(lambda))
  return(lambda)
}

lambda_sim2 = function(t,step=1){
  return(exp(sin(2*pi*t/(1440/step)) + abs(((floor(t/(1440/step))%%5)-2))/5 ))
}

lambda_sim3 = function(t,step=1,peaks=1){
  return(exp(2*(t%%(1440/step/peaks))/(1440/step) + abs(((floor(t/(1440/step))%%5)-2))/5 ))
  # return(exp( ((t%%(1440/step/peaks))/(1440/step))^3 -(((t+1440/2)%%(1440/step/2))/(1440/step))^2 + abs(((floor(t/(1440/step))%%5)-2))/5 ))
  # return(exp( -(((t+1440/2)%%(1440/step/peaks))/(1440/step))^2 + abs(((floor(t/(1440/step))%%5)-2))/5 ))
}

lambda_mix = function(t,step=1,peaks=1){
  return(exp((log(lambda_sim2(t,step))+log(lambda_sim3(t,step,peaks)))/3))
}

lambda_peaky = function(t,step=1,peaks=1){
  return(exp(sin((1440/step)/(2*pi*(t%%(1440/step) +1) )) + (t%%(1440/step/peaks))/(1440/step) + abs(((floor(t/(1440/step))%%5)-2))/5 ))
}

thin_sim = function(lambda, lambda_plus, T_end, seed=NULL){
  # Simulate the homogeneous process
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

AAO = function(effects,train){
  covariates = do.call(cbind,effects)
  out = glm(event~covariates-1,family="poisson",data=train)
  
  return(list("fullcoef"=out$coefficients,"lambdafit"=as.vector(exp(covariates %*% out$coefficients))))
}

AA0_pen = function(effects,train){
  covariates = do.call(cbind,effects)
  
  
  cvfit = cv.glmnet(x=covariates, y=train$event, nfolds = 10, family = "poisson", 
                    intercept = F,
                    type.measure="deviance")#, lambda=exp(seq(log(100), log(0.001), length.out = 50)))
  
  coeffs = as.vector(coef(cvfit, s="lambda.min"))[-1]
  
  return(list("fullcoef"=coeffs,"lambdafit"=as.vector(exp(covariates %*% coeffs)),"cvfit"=cvfit))
  
}

##### Algo 1 - Fit One by One #####

OBO = function(effects,train,initial_offset=0){
  offset= rep(0,nrow(effects[[1]])) + initial_offset
  coeffs=list()
  metrics = NULL
  
  for (i in 1:length(effects)){
    mod = glm(event~effects[[i]]-1,data=train,offset = offset)
    coeffs = append(coeffs,list(as.vector(mod$coefficients)))
    print(coeffs)
    offset = offset + as.matrix(effects[[i]])%*%coeffs[[i]]
    metrics = rbind(metrics,tibble(rmse=rmse(train$lambda,offset),
                                   mae=mae(train$lambda,offset)))
    
  }
  
  
  
  return(list("fullcoef"=coeffs,
              "lambdafit"=as.vector(exp(offset)),
              "metrics"=metrics))
}

OBO_pen = function(effects,train,initial_offset=0){
  offset= rep(0,nrow(effects[[1]])) + initial_offset
  coeffs=list()
  metrics = NULL
  
  for (i in 1:length(effects)){
    cvfit = cv.glmnet(x=effects[[i]], y=train$event, offset=offset,
                      nfolds = 10, family = "poisson", 
                      intercept = F,
                      type.measure="deviance")
    coeffs = append(coeffs,list(coef(cvfit, s="lambda.min")[-1]) )
    offset = offset + as.matrix(effects[[i]])%*%coeffs[[i]]
    
    metrics = rbind(metrics,tibble(rmse=rmse(train$lambda,offset),
                                   mae=mae(train$lambda,offset)))
  }
  
  return(list("fullcoef"=coeffs,
              "lambdafit"=as.vector(exp(offset)),
              "metrics"=metrics))
}

##### Algo 2 - Backfitting #####

rig_BAC = function(effects,train, maxiter = 100, eps = 0.001, init_offset = 0, ignore_errors = FALSE){
  coeffs = list()
  new_coeffs = list()
  fk = list()
  
  for (i in 1:length(effects)){
    coeffs[[i]] = rep(0,ncol(effects[[i]]))
    new_coeffs[[i]] = rep(0,ncol(effects[[i]]))
    fk[[i]] = rep(0, nrow(train))
  }
  
  check = rep(TRUE, length(effects))
  metrics = NULL
  
  intercept = mean(log(train$event+0.1))
  
  j=1
  while(j<=maxiter){
    print(paste("Iteration:",j))
    norm_coeffs = list()
    res = list()
    # for (k in sample((1:length(effects))[check]) )
    for (k in sample((1:length(effects))) )
    {
      
      offset = rep(0,nrow(effects[[1]]))+intercept+init_offset
      for (l in seq(1,length(effects))[-k]){
        offset = offset + fk[[l]]
      }
      
      if (ignore_errors){
        tryCatch({
          new_coeffs[[k]] = as.vector(glm(event~effects[[k]]-1, data=train, family='poisson', 
                                          offset=offset)$coefficients)
        }, error=function(e){
          cat("ERROR :",conditionMessage(e), "\n")
        })
      } else{
        new_coeffs[[k]] = as.vector(glm(event~effects[[k]]-1, data=train, family='poisson', 
                                        offset=offset)$coefficients)
      }

      
      res[[k]] = sqrt(sum((new_coeffs[[k]] - coeffs[[k]])^2)) 
      norm_coeffs[[k]] = sqrt(sum((coeffs[[k]])^2))
      
      fk[[k]] = as.matrix(effects[[k]])%*%new_coeffs[[k]]
      fk[[k]] = fk[[k]] - mean(fk[[k]])
      
    }
    print(new_coeffs)
    coeffs = new_coeffs
    
    
    lambda = rep(0,nrow(effects[[1]]))
    for (h in 1:length(effects)){
      lambda = lambda  + (effects[[h]]%*%coeffs[[h]])
    }
    lambda = as.vector(exp(lambda))
    print(paste("#############################", mean(lambda)))
    
    # print(res)
    print(paste("norm_coeff:",norm_coeffs))
    res[sapply(res, is.null)] = 0
    norm_coeffs[sapply(norm_coeffs, is.null)] = 0
    
    check = unlist(res)>unlist(norm_coeffs)*eps
    # print(check)
    # if(sum(check) == 0 | mean(lambda) < lambda_threshold)
    if(sum(check) == 0){
      break
    }
    
    metrics = rbind(metrics,tibble(rmse=rmse(train$lambda,lambda),
                                   mae=mae(train$lambda,lambda)))
    
    j = j+1
    if(j==maxiter){print("Warning: reached maximum number of iterations")}
  } 
  
  lambda = intercept
  for (h in 1:length(effects)){
    lambda = lambda  + fk[[h]]
  }
  lambda = as.vector(exp(lambda))
  return(list("fullcoef"=coeffs,
              "lambdafit"=lambda,
              "iterations"=j-1,
              "metrics"=metrics
  ))
}

rig_BAC_pen = function(effects,train, maxiter = 100, eps = 0.001, init_offset = 0, ignore_errors = FALSE){
  coeffs = list()
  new_coeffs = list()
  cvfits = list()
  fk = list()
  for (i in 1:length(effects)){
    coeffs[[i]] = rep(0,ncol(effects[[i]]))
    new_coeffs[[i]] = rep(0,ncol(effects[[i]]))
    cvfits[[i]] = 0 
    fk[[i]] = rep(0, nrow(train))
  }
  
  intercept = mean(log(train$event+0.1))
  metrics = NULL
  j=1
  while(j<=maxiter){
    
    print(paste("Iteration:",j))
    norm_coeffs = list()
    res = list()
    
    
    for (k in sample(1:length(effects)) ){
      
      offset = rep(0,nrow(effects[[1]])) + intercept + init_offset
      for (l in seq(1,length(effects))[-k]){
        offset = offset + fk[[l]]
      }
      
      
      if (ignore_errors){
      tryCatch({
        cvfit = cv.glmnet(x=effects[[k]], y=train$event, nfolds = 5, family = "poisson",
                          offset=offset,
                          nlambda=10,
                          intercept = F,
                          type.measure="deviance")# lambda=exp(seq(log(1000), log(0.001), length.out = 20)))
        # cvfit = glmnet(x=effects[[k]], y=train$event, family = "poisson",
        #                   offset=offset,
        #                   nlambda=10,
        #                   intercept = F,
        #                   type.measure="deviance")
      }, error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
      })}
      else {
        cvfit = cv.glmnet(x=effects[[k]], y=train$event, nfolds = 5, family = "poisson",
                          offset=offset,
                          nlambda=10,
                          intercept = F,
                          type.measure="deviance")
        }
      
      cvfits[[k]] = cvfit
      
      # glm_fit = glmnet(x=effects[[k]], y=train$event, family='poisson',
      #                  intercept = F,
      #                  offset=offset,lambda=cvfit$lambda.min)
      # new_coeffs[[k]] = as.vector(glm_fit$beta)
      
      # print(coef(cvfit, s="lambda.min"))
      new_coeffs[[k]] = as.vector(coef(cvfit, s="lambda.min")[-1])
      res[[k]] = sqrt(sum((new_coeffs[[k]] - coeffs[[k]])^2)) 
      norm_coeffs[[k]] = sqrt(sum((coeffs[[k]])^2)) 
      
      fk[[k]] = as.matrix(effects[[k]])%*%new_coeffs[[k]]
      fk[[k]] = fk[[k]] - mean(fk[[k]])
    }
    print(new_coeffs)
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
    metrics = rbind(metrics,tibble(rmse=rmse(train$lambda,lambda),
                                   mae=mae(train$lambda,lambda)))
    
    j = j+1
    if(j==maxiter){print("Warning: reached maximum number of iterations")}
  }
  
  
  lambda = intercept
  for (h in 1:length(effects)){
    lambda = lambda  + fk[[h]]
  }
  lambda = as.vector(exp(lambda))
  
  return(list("fullcoef"=coeffs,"lambdafit"=lambda,"cvfits"=cvfits,"iterations"=j-1))
  
}

BAC_GLM = function(effects,train, maxiter = 100, eps = 0.001, init_offset = 0, lambda_threshold=1e-8){
  
  coeffs = list()
  new_coeffs = list()
  for (i in 1:length(effects)){
    coeffs[[i]] = rep(0,ncol(effects[[i]]))
    new_coeffs[[i]] = rep(0,ncol(effects[[i]]))
  }
  
  check = rep(TRUE, length(effects))
  metrics = NULL
  
  j=1
  while(j<=maxiter){
    print(paste("Iteration:",j))
    norm_coeffs = list()
    res = list()
    # for (k in sample((1:length(effects))[check]) )
      for (k in sample((1:length(effects))) )
        {
      
      offset = rep(0,nrow(effects[[1]]))+init_offset
      for (l in seq(1,length(effects))[-k]){
        offset = offset + as.matrix(effects[[l]])%*%coeffs[[l]]
      }
      
      
      tryCatch({
        new_coeffs[[k]] = as.vector(glm(event~effects[[k]]-1, data=train, family='poisson', 
                                        offset=offset)$coefficients)
      }, error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
        })

      res[[k]] = sqrt(sum((new_coeffs[[k]] - coeffs[[k]])^2)) 
      norm_coeffs[[k]] = sqrt(sum((coeffs[[k]])^2)) 
    }
    print(new_coeffs)
    coeffs = new_coeffs
    
    
    lambda = rep(0,nrow(effects[[1]]))
    for (h in 1:length(effects)){
      lambda = lambda  + (effects[[h]]%*%coeffs[[h]])
    }
    lambda = as.vector(exp(lambda))
    print(paste("#############################", mean(lambda)))
    
    # print(res)
    print(paste("norm_coeff:",norm_coeffs))
    res[sapply(res, is.null)] = 0
    norm_coeffs[sapply(norm_coeffs, is.null)] = 0

    check = unlist(res)>unlist(norm_coeffs)*eps
    # print(check)
    # if(sum(check) == 0 | mean(lambda) < lambda_threshold)
    if(sum(check) == 0){
      break
    }
    
    metrics = rbind(metrics,tibble(rmse=rmse(train$lambda,lambda),
                                   mae=mae(train$lambda,lambda)))
    
    j = j+1
    if(j==maxiter){print("Warning: reached maximum number of iterations")}
  }
  
  
  lambda = rep(0,nrow(effects[[1]]))
  for (h in 1:length(effects)){
    lambda = lambda  + (effects[[h]]%*%coeffs[[h]])
  }
  lambda = as.vector(exp(lambda))
  return(list("fullcoef"=coeffs,
              "lambdafit"=lambda,
              "iterations"=j-1,
              "metrics"=metrics
              ))
    
}

BAC_GLMnet = function(effects,train, maxiter = 100, eps = 0.001, init_offset = 0){
  
  coeffs = list()
  new_coeffs = list()
  cvfits = list()
  for (i in 1:length(effects)){
    coeffs[[i]] = rep(0,ncol(effects[[i]]))
    new_coeffs[[i]] = rep(0,ncol(effects[[i]]))
    cvfits[[i]] = 0 
  }
  
  metrics = NULL
  j=1
  while(j<=maxiter){
    
    print(paste("Iteration:",j))
    norm_coeffs = list()
    res = list()

    
    for (k in sample(1:length(effects)) ){
      
      offset = rep(0,nrow(effects[[1]])) + init_offset
      for (l in seq(1,length(effects))[-k]){
        offset = offset + as.matrix(effects[[l]])%*%coeffs[[l]]
      }
      
      tryCatch({
        cvfit = cv.glmnet(x=effects[[k]], y=train$event, nfolds = 5, family = "poisson",
                          offset=offset,
                          nlambda=10,
                          intercept = F,
                          type.measure="deviance")# lambda=exp(seq(log(1000), log(0.001), length.out = 20)))
        # cvfit = glmnet(x=effects[[k]], y=train$event, family = "poisson",
        #                   offset=offset,
        #                   nlambda=10,
        #                   intercept = F,
        #                   type.measure="deviance")
      }, error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
      })

      cvfits[[k]] = cvfit
      
      # glm_fit = glmnet(x=effects[[k]], y=train$event, family='poisson',
      #                  intercept = F,
      #                  offset=offset,lambda=cvfit$lambda.min)
      # new_coeffs[[k]] = as.vector(glm_fit$beta)
      
      # print(coef(cvfit, s="lambda.min"))
      new_coeffs[[k]] = as.vector(coef(cvfit, s="lambda.min")[-1])
      res[[k]] = sqrt(sum((new_coeffs[[k]] - coeffs[[k]])^2)) 
      norm_coeffs[[k]] = sqrt(sum((coeffs[[k]])^2)) 
    }
    print(new_coeffs)
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
    metrics = rbind(metrics,tibble(rmse=rmse(train$lambda,lambda),
                                   mae=mae(train$lambda,lambda)))
    
    j = j+1
    if(j==maxiter){print("Warning: reached maximum number of iterations")}
  }
  
  
  lambda = rep(0,nrow(effects[[1]]))
  for (h in 1:length(effects)){
    lambda = lambda  + (effects[[h]]%*%coeffs[[h]])
  }
  lambda = as.vector(exp(lambda))
  
  return(list("fullcoef"=coeffs,"lambdafit"=lambda,"cvfits"=cvfits,"iterations"=j-1))
  
}

GAM_then_wav = function(formula_gam,effects_wav,train, maxiter = 100, eps = 0.001){
  
  mod_gam = gam(formula_gam, data=train,family=poisson)
  
  return(OBO_pen(effects_wav,train,initial_offset = as.vector(predict(mod_gam,type="response")) ))
  
}

BAC_GAM_then_wav = function(formula_gam,effects_wav,train, maxiter = 100, eps = 0.001){
  
   mod_gam = gam(formula_gam, data=train,family=poisson)
  
   return(BAC_GLM(effects_wav,train,init_offset = as.vector(predict(mod_gam,type="response")) ))
  
}

BAC_GLMnet_glob = function(effects_glm,effects_gam,train, maxiter = 100, eps = 0.001){
  
  if(length(effects_gam)>0){
    colnames_gam = c(paste0("effect_",seq(1,ncol(effects_gam))),"event")
    gam_data = tibble(effects_gam) %>%
      mutate(event=train$event)
    colnames(gam_data) = colnames_gam
    # print(gam_data)
  }  
  
  effects=c(effects_glm,effects_gam)
  boundary=length(effects_glm)
  coeffs = list()
  new_coeffs = list()
  fits = list("cvfits"=list(),"gams"=list())
  for (i in 1:length(effects)){
    if(i<=boundary){
      fits[[1]][[i]] = 0
      new_coeffs[[i]] = rep(0,ncol(effects[[i]]))
      coeffs[[i]] = rep(0,ncol(effects[[i]]))
    } else if (i>boundary){
      coeffs[[i]] = 0
      new_coeffs[[i]] = 0
      fits[[2]][[i-boundary]] = 0
    }
  }
  
  j=1
  while(j<=maxiter){
    
    print(paste("Iteration:",j))
    norm_coeffs = list()
    res = list()
    
    
    for (k in sample(1:length(effects)) ){
      
      offset = rep(0,nrow(effects[[1]]))
      for (l in seq(1,length(effects))[-k]){
        if(l <= boundary){
          offset = offset + as.matrix(effects[[l]])%*%coeffs[[l]]
        } else if (l > boundary){
          if(typeof(fits[[2]][[l-boundary]])=="list"){
            offset = offset + as.matrix(predict(fits[[2]][[l-boundary]],type="terms"))
          }
        }
      }
      
      if( k <= boundary){
        cvfit = cv.glmnet(x=effects[[k]], y=train$event, nfolds = 10, family = "poisson", 
                          offset=offset,
                          intercept = F,
                          type.measure="deviance")# lambda=exp(seq(log(1000), log(0.001), length.out = 20)))
        
        fits[[1]][[k]] = cvfit
        new_coeffs[[k]] = as.vector(coef(cvfit, s="lambda.min")[-1])
        
      } else if( k > boundary){
        
        fmla = as.formula(paste0("event ~ s(",colnames_gam[[k-boundary]] ,") - 1"))
        mod = gam(fmla, 
                  data=gam_data, family = "poisson",
                  offset=offset)
        
        fits[[2]][[k-boundary]] = mod
        new_coeffs[[k]] = as.vector(coef(mod))
      }
      res[[k]] = sqrt(sum((new_coeffs[[k]] - coeffs[[k]])^2)) 
      norm_coeffs[[k]] = sqrt(sum((coeffs[[k]])^2)) 
    }
    print(new_coeffs)
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
    if(h<=boundary){
      lambda = lambda  + (effects[[h]]%*%coeffs[[h]])
    } else if(h>boundary){
      lambda = lambda  + predict(fits[[2]][[h-boundary]],type="terms")
    }
    
  }
  lambda = as.vector(exp(lambda))
  
  return(list("fullcoef"=coeffs,"lambdafit"=lambda,"cvfits"=fits,"iterations"=j-1))
  
}

rig_BAC_final = function(fmla,train, maxiter = 100, eps = 0.001, init_offset = 0, ignore_errors = FALSE){
  
  full_effects = extract_effects(fmla,data,res_max = 8)
  effects = compact(compact(c(full_effects$linear,full_effects$wavelets)))
  sterms =full_effects$splines
  
  coeffs = list()
  new_coeffs = list()
  fk = list()
  
  #### Adding splines terms
  boundary=length(effects)
  effects = c(effects,sterms)
  
  for (i in 1:length(effects)){
    if (i<=boundary){
      coeffs[[i]] = rep(0,ncol(effects[[i]]))
      new_coeffs[[i]] = rep(0,ncol(effects[[i]]))
      fk[[i]] = rep(0, nrow(train))
    }
    else if (i>boundary){
      coeffs[[i]] = 0
      new_coeffs[[i]] = 0
      fk[[i]] = rep(0, nrow(train))
    }
  }
  
  check = rep(TRUE, length(effects))
  metrics = NULL
  
  intercept = mean(log(train$event+0.1))
  
  j=1
  while(j<=maxiter){
    print(paste("Iteration:",j))
    norm_coeffs = list()
    res = list()
    # for (k in sample((1:length(effects))[check]) )
    for (k in sample((1:length(effects))) )
    {
      
      offset = rep(0,nrow(train))+intercept+init_offset
      for (l in seq(1,length(effects))[-k]){
        offset = offset + fk[[l]]
      }
      
      if(k<=boundary){
                if (ignore_errors){
                  tryCatch({
                    new_coeffs[[k]] = as.vector(glm(event~effects[[k]]-1, data=train, family='poisson', 
                                                    offset=offset)$coefficients)
                  }, error=function(e){
                    cat("ERROR :",conditionMessage(e), "\n")
                  })
                } else{
                  new_coeffs[[k]] = as.vector(glm(event~effects[[k]]-1, data=train, family='poisson', 
                                                  offset=offset)$coefficients)}
                fk[[k]] = as.matrix(effects[[k]])%*%new_coeffs[[k]]
        
      } else if(k>boundary){
        fmla = as.formula(paste0("event ~ ", effects[[k]] ,"- 1"))
        mod = gam(fmla, 
                  data=train, family = "poisson",
                  offset=offset)
        new_coeffs[[k]] = as.vector(coef(mod))
        
        fk[[k]] = as.matrix(predict(mod))
      }
      
      res[[k]] = sqrt(sum((new_coeffs[[k]] - coeffs[[k]])^2)) 
      norm_coeffs[[k]] = sqrt(sum((coeffs[[k]])^2))
      fk[[k]] = fk[[k]] - mean(fk[[k]])

    }
    print(new_coeffs)
    coeffs = new_coeffs
    
    lambda = intercept
    for (h in 1:length(effects)){
      lambda = lambda  + fk[[h]]
    }
    lambda = as.vector(exp(lambda))
    print(paste("#############################", mean(lambda)))
    
    # print(res)
    print(paste("norm_coeff:",norm_coeffs))
    res[sapply(res, is.null)] = 0
    norm_coeffs[sapply(norm_coeffs, is.null)] = 0
    
    check = unlist(res)>unlist(norm_coeffs)*eps
    print(check)
    # if(sum(check) == 0 | mean(lambda) < lambda_threshold)
    if(sum(check) == 0){
      break
    }
    
    metrics = rbind(metrics,tibble(rmse=rmse(train$lambda,lambda),
                                   mae=mae(train$lambda,lambda)))
    
    j = j+1
    if(j==maxiter){print("Warning: reached maximum number of iterations")}
  } 
  
  lambda = intercept
  for (h in 1:length(effects)){
    lambda = lambda  + fk[[h]]
  }
  lambda = as.vector(exp(lambda))
  return(list("fullcoef"=coeffs,
              "lambdafit"=lambda,
              "iterations"=j-1,
              "metrics"=metrics
  ))
}



##### Algo 2 - Backfitting 2 #####

BAC_GLM2 = function(effects,train, maxiter = 100, eps = 0.001, init_offset = 0, lambda_threshold=1e-8){
  
  coeffs = list()
  new_coeffs = list()
  for (i in 1:length(effects)){
    coeffs[[i]] = rep(0,ncol(effects[[i]]))
    new_coeffs[[i]] = rep(0,ncol(effects[[i]]))
  }
  
  check = rep(TRUE, length(effects))
  metrics = NULL
  
  j=1
  while(j<=maxiter){
    print(paste("Iteration:",j))
    norm_coeffs = list()
    res = list()
    
    offset = rep(0,nrow(effects[[1]]))+init_offset
    for (l in seq(1,length(effects))){
      offset = offset + as.matrix(effects[[l]])%*%coeffs[[l]]
    }
  
    # for (k in sample((1:length(effects))[check]) )
    for (k in sample((1:length(effects))) )
    {
      
      tryCatch({
        new_coeffs[[k]] = as.vector(glm(event~effects[[k]]-1, data=train, family='poisson', 
                                        offset=offset)$coefficients)
      }, error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
      })
      
      res[[k]] = sqrt(sum((new_coeffs[[k]] - coeffs[[k]])^2)) 
      norm_coeffs[[k]] = sqrt(sum((coeffs[[k]])^2)) 
    }
    print(new_coeffs)
    coeffs = new_coeffs
    
    
    lambda = rep(0,nrow(effects[[1]]))
    for (h in 1:length(effects)){
      lambda = lambda  + (effects[[h]]%*%coeffs[[h]])
    }
    lambda = as.vector(exp(lambda))
    print(paste("#############################", mean(lambda)))
    
    # print(res)
    # print(norm_coeffs)
    res[sapply(res, is.null)] = 0
    norm_coeffs[sapply(norm_coeffs, is.null)] = 0
    
    check = unlist(res)>unlist(norm_coeffs)*eps
    # print(check)
    # if(sum(check) == 0 | mean(lambda) < lambda_threshold)
    if(sum(check) == 0){
      break
    }
    
    metrics = rbind(metrics,tibble(rmse=rmse(train$lambda,lambda),
                                   mae=mae(train$lambda,lambda)))
    
    j = j+1
    if(j==maxiter){print("Warning: reached maximum number of iterations")}
  }
  
  
  lambda = rep(0,nrow(effects[[1]]))
  for (h in 1:length(effects)){
    lambda = lambda  + (effects[[h]]%*%coeffs[[h]])
  }
  lambda = as.vector(exp(lambda))
  return(list("fullcoef"=coeffs,
              "lambdafit"=lambda,
              "iterations"=j-1,
              "metrics"=metrics
  ))
  
}

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