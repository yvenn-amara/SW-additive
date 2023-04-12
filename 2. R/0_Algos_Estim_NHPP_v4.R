library(fastDummies)
library(gamwave)
library(glmnet)
library(mgcv)

extract_effects = function(fmla,data,forced_res_list=NULL,res_max=8,wav_intercept=FALSE,get_full_rank=TRUE){
  tf = terms.formula(fmla,specials=c("s","w","factor","wp"))
  if(attr(tf,"response")>0){
    lterms = rownames(attr(tf,'factors'))[-c(1,unlist(attr(tf, "specials")))]
  }else{
    lterms = rownames(attr(tf,'factors'))[-unlist(attr(tf, "specials"))]
  }
  sterms = rownames(attr(tf,'factors'))[attr(tf, "specials")$s]
  wterms = c(rownames(attr(tf,'factors'))[attr(tf, "specials")$w],rownames(attr(tf,'factors'))[attr(tf, "specials")$wp])
  fterms = rownames(attr(tf,'factors'))[attr(tf, "specials")$factor]
  ft = sub("\\).*", "", sub(".*\\(", "", fterms)) 
  
  if(length(lterms)>0){
    effect_lin = as.matrix(data[lterms])
  } else {
    effect_lin = NULL
  }
  
  effect_fac = NULL
  if (length(ft)>0){
    for (z in 1:length(ft)){
      effect_fac = cbind(effect_fac,as.matrix(dummy_cols(data[ft[[z]]],remove_first_dummy = TRUE,remove_selected_columns = TRUE)))
    }
  }
  
  effect_lin_full = cbind(effect_lin,effect_fac)
  
  
  wt = sub("\\).*", "", sub(".*\\(", "", wterms)) 
  res_list=NULL
  if (length(wt>0)){
    effects_wav = list()
    res_list = list()
    for (i in 1:length(wt)){
      var = sub("\\).*", "", sub(".*\\(", "", wterms[i])) 
      series = data[[var]]

      if (is.null(forced_res_list)){
        res = res_max
      } else {
        res = forced_res_list[[i]] 
        get_full_rank = FALSE
      }
      
      if(wav_intercept){
        wav_effect = cbind(1,WavD(series, numLevels=res, filterNumber = 1))
        colnames(wav_effect) = paste0("wav_",var,"_",as.character(seq(1,2^res)))
      }else{
        wav_effect = WavD(series, numLevels=res, filterNumber = 1)
        colnames(wav_effect) = paste0("wav_",var,"_",as.character(seq(1,2^res-1)))
      }
      
      effects_wav[[i]] = wav_effect
      X=as.data.frame(wav_effect)
      
      if(get_full_rank){
        while(ncol(X) != qr(X)$rank){
          print(res)
          res = res - 1
          
          if(wav_intercept){
            wav_effect = cbind(1,WavD(series, numLevels=res, filterNumber = 1))
            colnames(wav_effect) = paste0("wav_",var,"_",as.character(seq(1,2^res)))
          }else{
            wav_effect = WavD(series, numLevels=res, filterNumber = 1)
            colnames(wav_effect) = paste0("wav_",var,"_",as.character(seq(1,2^res-1)))
          }
          
          effects_wav[[i]] = wav_effect
          X=as.data.frame(wav_effect)
        }
      }

      
      res_list[[i]] = res
      
    }
    X = as.data.frame(effects_wav)
    if(ncol(X) == qr(X)$rank){print("Full Rank")}else{print("Not Full Rank")}
  } else{effects_wav = NULL}
  
  w = length(rownames(attr(tf,'factors'))[attr(tf, "specials")$w])
  wp = length(rownames(attr(tf,'factors'))[attr(tf, "specials")$wp])
  
  if(w>0){
    wav =effects_wav[1:w]
  } else {
    wav = NULL
  }
  
  if(wp>0){
    wav_pen =effects_wav[w+1:wp]
  } else {
    wav_pen = NULL
  }
  
  return(list("linear"=list(effect_lin_full),
              "splines"=sterms,
              "wavelets"=wav,
              "wavelets_pen"=wav_pen,
              "res_list"=res_list))
}

rig_BAC_final = function(fmla,train, forced_res = NULL, res_max = 8, maxiter = 100, eps = 0.001, init_offset = 0, ignore_errors = FALSE,wav_intercept=FALSE,get_full_rank=TRUE,simu=TRUE){
  start_time  = Sys.time()
  full_effects = extract_effects(fmla,train,forced_res=forced_res,res_max = res_max,wav_intercept=wav_intercept,get_full_rank=get_full_rank)
  # effects = compact(compact(c(full_effects$linear,full_effects$wavelets)))
  effects_lin = compact(full_effects$linear)
  sterms =full_effects$splines
  effects_wav = compact(full_effects$wavelets)
  effects_pen = full_effects$wavelets_pen
  
  coeffs = list()
  new_coeffs = list()
  fk = list()
  res_fk = list()
  old_fk = list()
  
  norm_coeffs = list()
  res = list()
  models = list()
  
  #### Adding splines terms
  boundary=length(effects_lin)
  boundary2=length(sterms)
  boundary3=length(effects_wav)
  boundary4=length(effects_pen)
  effects = c(effects_lin,sterms,effects_wav,effects_pen)
  
  for (i in 1:length(effects)){
    if (i>boundary & i <= boundary+boundary2){
      coeffs[[i]] = 0
      new_coeffs[[i]] = 0
      norm_coeffs[[i]] = 0
      res[[i]] = 0
    } else {
      coeffs[[i]] = rep(0,ncol(effects[[i]]))
      new_coeffs[[i]] = rep(0,ncol(effects[[i]]))
      norm_coeffs[[i]] = rep(0,ncol(effects[[i]]))
      res[[i]] = rep(0,ncol(effects[[i]]))
    }
    fk[[i]] = rep(0, nrow(train))
    res_fk[[i]] = rep(0, nrow(train))
    old_fk[[i]] = rep(0, nrow(train))
    models[[i]] = 0
  }
  
  check = rep(TRUE, length(effects))
  metrics = NULL
  
  intercept = mean(log(train$event+0.1))
  
  j=1
  while(j<=maxiter){
    print(paste("======= Iteration",j," :",j))
    
    still_to_fit = (1:length(effects))[check]
    
    if(length(still_to_fit)==1){
      index = still_to_fit
    } else{
      index = sample(still_to_fit)
    }
    
    # for (k in sample((1:length(effects))[check]) )
    for (k in index)
    # for (k in ((1:length(effects))[check]) )
    # for (k in sample((1:length(effects))) )
    {
      offset = rep(0,nrow(train))+intercept+init_offset
      for (l in seq(1,length(effects))[-k]){
        offset = offset + fk[[l]]
      }
      
      if(k<=boundary | ((k > boundary + boundary2) & (k <= boundary + boundary2 + boundary3)) ){
        # print(paste("GLM:",k))
        if (ignore_errors){
          tryCatch({
            models[[k]] = glm(event~effects[[k]]-1, data=train, family='poisson', 
                              offset=offset)
            new_coeffs[[k]] = as.vector(models[[k]]$coefficients)
          }, error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
          })
        } else{
          models[[k]] = glm(event~effects[[k]]-1, data=train, family='poisson', 
                            offset=offset)
          new_coeffs[[k]] = as.vector(models[[k]]$coefficients)}
        
        fk[[k]] = as.matrix(effects[[k]])%*%new_coeffs[[k]]
        
      } else if (k>boundary+boundary2+boundary3){
        # print(paste("WAV_pen:",k))
        if (ignore_errors){
          tryCatch({
            cvfit = cv.glmnet(x=effects[[k]], y=train$event, nfolds = 5, family = "poisson",
                              offset=offset,
                              nlambda=100,
                              alpha=1,
                              lambda.min.ratio=0.001,
                              intercept = F,
                              type.measure="deviance")
            new_coeffs[[k]] = as.vector(coef(cvfit, s="lambda.min")[-1])
          }, error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
          })
        } else{
          cvfit = cv.glmnet(x=effects[[k]], y=train$event, nfolds = 5, family = "poisson",
                            offset=offset,
                            nlambda=100,
                            alpha=1,
                            lambda.min.ratio=0.001,
                            intercept = F,
                            type.measure="deviance")
          new_coeffs[[k]] = as.vector(coef(cvfit, s="lambda.min")[-1])}
        models[[k]] = cvfit
        fk[[k]] = as.matrix(effects[[k]])%*%new_coeffs[[k]]
        
        } else if(k>boundary & k <= boundary + boundary2){
          # print(paste("GAM:",k))
        fmla_gam = as.formula(paste0("event ~ ", effects[[k]] ,"- 1"))
        mod = gam(fmla_gam, 
                  data=train, family = "poisson",
                  offset=offset)
        new_coeffs[[k]] = as.vector(coef(mod))
        
        fk[[k]] = as.matrix(predict(mod))
        models[[k]] = mod
      }
      
      # res[[k]] = sqrt(sum((new_coeffs[[k]] - coeffs[[k]])^2)) 
      res[[k]] = sqrt(mean((new_coeffs[[k]] - coeffs[[k]])^2)) 
      
      # norm_coeffs[[k]] = sqrt(sum((coeffs[[k]])^2))
      norm_coeffs[[k]] = sqrt(mean((coeffs[[k]])^2))
      
      fk[[k]] = fk[[k]] - mean(fk[[k]])
      
      res_fk[[k]] = sqrt(mean((old_fk[[k]] - fk[[k]])^2))
      old_fk[[k]] = fk[[k]] 
    }
    # print(new_coeffs)
    coeffs = new_coeffs
    
    lambda = intercept
    for (h in 1:length(effects)){
      lambda = lambda  + fk[[h]]
    }
    lambda = as.vector(exp(lambda))
    # print(paste("#############################", mean(lambda)))
    
    # print(res)
    # print(paste("norm_coeff:",norm_coeffs))
    
    # res[sapply(res, is.null)] = 0
    # norm_coeffs[sapply(norm_coeffs, is.null)] = 0
    
    # check = unlist(res)>unlist(norm_coeffs)*eps
    # print(paste("res:",unlist(res) ))
    # print(paste("res2:",eps*unlist(norm_coeffs)))
    # check = unlist(res)>eps*unlist(norm_coeffs)
    # check = unlist(res)>eps
    
    check = unlist(res_fk) > eps
    # if(length(effects_lin)>0){check[1]=FALSE}
    print(check)
    # if(sum(check) == 0 | mean(lambda) < lambda_threshold)
    if(sum(check) == 0){
      break
    }
    
    # metrics = rbind(metrics,tibble(rmse=rmse(train$lambda,lambda),
                                   # mae=mae(train$lambda,lambda)))
    
    j = j+1
    if(j==maxiter){print("Warning: reached maximum number of iterations")}
  } 
  
  lin_pred = intercept
  spl_pred = 0
  wav_pred = 0
  lambdafit = intercept
  for (g in 1:length(effects)){
    if(g<=boundary | ((g > boundary + boundary2) & (g <= boundary + boundary2 + boundary3)) ){
      lin_pred = lin_pred + fk[[g]]
    } else if (g>boundary+boundary2+boundary3){
      spl_pred = spl_pred + fk[[g]]
    } else if (g>boundary & g <= boundary + boundary2){
      wav_pred = wav_pred + fk[[g]]
    }
    lambdafit = lambdafit  + fk[[g]]
  }
  lambdafit = as.vector(exp(lambdafit))
  
  
  
  ### Degrees of freedom
  dof = 1
  for (z in 1:length(effects)){
    if (z>boundary & z <= boundary+boundary2){
     dof = dof + sum(influence(models[[z]]))
    } else if (z<=boundary | ((z > boundary + boundary2) & (z <= boundary + boundary2 + boundary3))){
      dof = dof + sum(influence(models[[z]])$hat)
    } else {
      dof = dof + sum(coeffs[[z]]!=0)
    }
  }
  
  if(simu){
    final = train %>% dplyr::select(date,ymd_date,lambda) %>% dplyr::rename("y_true"="lambda") %>% mutate(y_pred=lambdafit )
  }else{
    final = train %>% dplyr::select(date,ymd_date,event) %>% dplyr::rename("y_true"="event") %>% mutate(y_pred=round(lambdafit ))
    
  }
  
  return(list("runtime"=difftime(Sys.time(),start_time,units = "secs"),
              "simu"=simu,
              "res_list"=full_effects$res_list,
              "formula"=fmla,
              "intercept"=intercept,
              "wav_intercept"=wav_intercept,
              "get_full_rank"=get_full_rank,
              "res_max"=res_max,
              "fullcoef"=coeffs,
              "lambdafit"=lambdafit,
              "iterations"=j-1,
              # "metrics"=metrics,
              "models"=models,
              "num_param"=length(unlist(coeffs)),
              "num_param_non_zeros"=sum(unlist(coeffs)!=0),
              "dof"= dof,
              "final"=final,
              "lin_pred"=lin_pred,
              "spl_pred"=spl_pred,
              "wav_pred"=wav_pred
              ))
}

OBO_pen = function(fmla,train, forced_res = NULL, res_max = 8, maxiter = 100, eps = 0.001, init_offset = 0, ignore_errors = FALSE,wav_intercept=FALSE,get_full_rank=TRUE,simu=TRUE){
  start_time  = Sys.time()
  full_effects = extract_effects(fmla,train,forced_res=forced_res,res_max = res_max,wav_intercept=wav_intercept,get_full_rank=get_full_rank)
  # effects = compact(compact(c(full_effects$linear,full_effects$wavelets)))
  effects_lin = compact(full_effects$linear)
  sterms =full_effects$splines
  effects_wav = compact(full_effects$wavelets)
  effects_pen = full_effects$wavelets_pen
  
  coeffs = list()
  new_coeffs = list()
  fk = list()
  res_fk = list()
  old_fk = list()
  
  norm_coeffs = list()
  res = list()
  models = list()
  
  #### Adding splines terms
  boundary=length(effects_lin)
  boundary2=length(sterms)
  boundary3=length(effects_wav)
  boundary4=length(effects_pen)
  effects = c(effects_lin,sterms,effects_wav,effects_pen)
  
  for (i in 1:length(effects)){
    if (i>boundary & i <= boundary+boundary2){
      coeffs[[i]] = 0
      new_coeffs[[i]] = 0
      norm_coeffs[[i]] = 0
      res[[i]] = 0
      
    } else {
      coeffs[[i]] = rep(0,ncol(effects[[i]]))
      new_coeffs[[i]] = rep(0,ncol(effects[[i]]))
      norm_coeffs[[i]] = rep(0,ncol(effects[[i]]))
      res[[i]] = rep(0,ncol(effects[[i]]))

    }
    
    fk[[i]] = rep(0, nrow(train))
    res_fk[[i]] = rep(0, nrow(train))
    old_fk[[i]] = rep(0, nrow(train))
    models[[i]] = 0
  }
  
  check = rep(TRUE, length(effects))
  metrics = NULL
  
  intercept = mean(log(train$event+0.1))
  
  offset = rep(0,nrow(train))+intercept+init_offset
  
  for (k in 1:length(effects)){
    if(k<=boundary | ((k > boundary + boundary2) & (k <= boundary + boundary2 + boundary3)) ){
      if (ignore_errors){
        tryCatch({
          models[[k]] = glm(event~effects[[k]]-1, data=train, family='poisson', 
                            offset=offset)
        }, error=function(e){
          cat("ERROR :",conditionMessage(e), "\n")
        })
      } else{
        models[[k]] = glm(event~effects[[k]]-1, data=train, family='poisson', 
                          offset=offset)}
      # print("GLM")
      new_coeffs[[k]] = as.vector(models[[k]]$coefficients)
      fk[[k]] = as.matrix(effects[[k]])%*%new_coeffs[[k]]
      
    } else if (k>boundary+boundary2+boundary3){
      if (ignore_errors){
        tryCatch({
          cvfit = cv.glmnet(x=effects[[k]], y=train$event, nfolds = 5, family = "poisson",
                            offset=offset,
                            nlambda=100,
                            alpha=1,
                            lambda.min.ratio=0.001,
                            intercept = F,
                            type.measure="deviance")
        }, error=function(e){
          cat("ERROR :",conditionMessage(e), "\n")
        })
      } else{
        # print("GLMNET")
        cvfit = cv.glmnet(x=effects[[k]], y=train$event, nfolds = 5, family = "poisson",
                          offset=offset,
                          nlambda=100,
                          alpha=1,
                          lambda.min.ratio=0.001,
                          intercept = F,
                          type.measure="deviance")
      }
      new_coeffs[[k]] = as.vector(coef(cvfit, s="lambda.min")[-1])
      fk[[k]] = as.matrix(effects[[k]])%*%new_coeffs[[k]]
      models[[k]] = cvfit
    } else if(k>boundary & k <= boundary + boundary2){
      # print("GAM")
      fmla_gam = as.formula(paste0("event ~ ", effects[[k]] ,"- 1"))
      mod = gam(fmla_gam, 
                data=train, family = "poisson",
                offset=offset)
      new_coeffs[[k]] = as.vector(coef(mod))
      
      fk[[k]] = as.matrix(predict(mod))
      models[[k]] = mod
    }
    
    res[[k]] = sqrt(sum((new_coeffs[[k]] - coeffs[[k]])^2)) 
    norm_coeffs[[k]] = sqrt(sum((coeffs[[k]])^2))
    fk[[k]] = fk[[k]] - mean(fk[[k]])
    
    offset = offset + fk[[k]]
    
    # metrics = rbind(metrics,tibble(rmse=rmse(train$lambda,offset),
    #                                mae=mae(train$lambda,offset)))
  }
  coeffs = new_coeffs
  
  ####
  lin_pred = intercept
  spl_pred = 0
  wav_pred = 0
  lambdafit = intercept
  for (g in 1:length(effects)){
    if(g<=boundary | ((g > boundary + boundary2) & (g <= boundary + boundary2 + boundary3)) ){
      lin_pred = lin_pred + fk[[g]]
    } else if (g>boundary+boundary2+boundary3){
      spl_pred = spl_pred + fk[[g]]
    } else if (g>boundary & g <= boundary + boundary2){
      wav_pred = wav_pred + fk[[g]]
    }
    lambdafit = lambdafit  + fk[[g]]
  }
  lambdafit = as.vector(exp(lambdafit))
  
  ### Degrees of freedom
  dof = 1
  
  for (z in 1:length(effects)){
    if (z>boundary & z <= boundary+boundary2){
      dof = dof + sum(influence(models[[z]]))
    } else if (z<=boundary | ((z > boundary + boundary2) & (z <= boundary + boundary2 + boundary3))){
      dof = dof + sum(influence(models[[z]])$hat)
    } else {
      dof = dof + sum(coeffs[[z]]!=0)
    }
  }
  
  if(simu){
    final = train %>% dplyr::select(date,ymd_date,lambda) %>% dplyr::rename("y_true"="lambda") %>% mutate(y_pred=as.vector(exp(offset)))
  }else{
    final = train %>% dplyr::select(date,ymd_date,event) %>% dplyr::rename("y_true"="event") %>% mutate(y_pred=round(as.vector(exp(offset)) ))
    
  }
  
  return(list("runtime"=difftime(Sys.time(),start_time,units = "secs"),
              "simu"=simu,
              "res_list"=full_effects$res_list,
              "formula"=fmla,
              "intercept"=intercept,
              "wav_intercept"=wav_intercept,
              "get_full_rank"=get_full_rank,
              "res_max"=res_max,
              "fullcoef"=coeffs,
              "lambdafit"=lambdafit,
              "iterations"=1,
              # "metrics"=metrics,
              "models"=models,
              "num_param"=length(unlist(coeffs)),
              "num_param_non_zeros"=sum(unlist(coeffs)!=0),
              "dof"= dof,
              "final"=final,
              "lin_pred"=lin_pred,
              "spl_pred"=spl_pred,
              "wav_pred"=wav_pred
  ))
  
  
  # return(list("runtime"=difftime(Sys.time(),start_time,units = "secs"),
  #             "fullcoef"=coeffs,
  #             "lambdafit"=as.vector(exp(offset)),
  #             #"metrics"=metrics,
  #             "models"=models,
  #             "num_param"=length(unlist(coeffs)),
  #             "num_param_non_zeros"=sum(unlist(coeffs)!=0),
  #             "dof"= dof,
  #             "final"= final))

}

predict.sw = function(out, newdata){
  full_effects = extract_effects(out$formula,newdata,forced_res_list=out$res_list,res_max = out$res_max,wav_intercept=out$wav_intercept,get_full_rank=FALSE)
  
  effects_lin = compact(full_effects$linear)
  sterms =full_effects$splines
  effects_wav = compact(full_effects$wavelets)
  effects_pen = full_effects$wavelets_pen
  
  # fk = list()
  
  #### Adding splines terms
  boundary=length(effects_lin)
  boundary2=length(sterms)
  boundary3=length(effects_wav)
  boundary4=length(effects_pen)
  effects = c(effects_lin,sterms,effects_wav,effects_pen)
  # print(effects)
  lambda = out$intercept
  for (k in 1:length(effects)){
    print(k)
    # print(out$models[[k]])
    if(k<=boundary | ((k > boundary + boundary2) & (k <= boundary + boundary2 + boundary3)) ){
      # print('glm')
      lambda = lambda + as.matrix(effects[[k]])%*%out$fullcoef[[k]]
      
    } else if (k>boundary+boundary2+boundary3){
      # print('glmnet')
      # print(effects[[k]])
      lambda = lambda + as.matrix(effects[[k]])%*%out$fullcoef[[k]]
      
    } else if(k>boundary & k <= boundary + boundary2){
      # print('gam')
      lambda = lambda + as.matrix(predict(out$models[[k]],newdata))
    }
    
  }
  lambda_pred = as.vector(exp(lambda))
  
  if(out$simu){
    final = newdata %>% dplyr::select(date,ymd_date,lambda) %>% dplyr::rename("y_true"="lambda") %>% mutate(y_pred=lambda_pred )
  }else{
    final = newdata %>% dplyr::select(date,ymd_date,event) %>% dplyr::rename("y_true"="event") %>% mutate(y_pred=round(lambda_pred ))
  }
  return(final)
}
