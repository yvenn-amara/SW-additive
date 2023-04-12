##### Plotting Wavelet Basis #####
plot_wav = function(x,numLevels, filterNumber = 1, ylim=c(-5,5), coef_is_null=NULL){
  #vector of zeros and ones of the same size as 2^numlevels
  require("lattice")
  n = (ceiling(max(x))-floor(min(x)))*10
  x_seq=seq(min(x),max(x), length=n)
  Z = WavD(x_seq, numLevels=numLevels, filterNumber = filterNumber)
  # B = cbind(rep(1,n),Z)
  B = Z
  
  if(is.null(coef_is_null)==FALSE){
    B = t(t(B) * coef_is_null)
  }
  
  df = data.frame(x = x_seq,B)
  # print(names(df))
  
  form = as.formula(paste(paste(names(df)[- 1],  collapse = ' + '),'x',  sep = '~'))
  return(xyplot(form,data = df,
                type = 'l', grid=TRUE,outer =TRUE,
                ylab="DaubLeAsymm",scales=list(tick.number=10)))
  
}

plot_wav2 = function(x,numLevels, filterNumber = 1, ylim=c(-5,5), coef_is_null=NULL){
  #vector of zeros and ones of the same size as 2^numlevels
  require("lattice")
  n = (ceiling(max(x))-floor(min(x)))*10
  x_seq=seq(min(x),max(x), length=n)
  Z = WavD(x_seq, numLevels=numLevels, filterNumber = filterNumber)
  # B = cbind(rep(1,n),Z)
  B = Z
  
  if(is.null(coef_is_null)==FALSE){
    B = t(t(B) * coef_is_null)
  }
  
  df = tibble(x = x_seq,wav=apply(B,1,sum))
  # print(names(df))
  
  return(ggplot(data=df, aes(x=x))+geom_line(aes(y=wav))+theme_Publication()+xlab("Time of day [hours]")+ylab('Response'))
}

##### Basic functions
blocks = function(x){
  x = (x-min(x))/(max(x)-min(x))
  t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 
         0.78, 0.81)
  h1 <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
  res = 0
  for (i in seq(1, length(h1))) {
    res = res + (h1[i] * (1 + sign(x - t[i])))/2
  }
  return(res +sign(min(res))*min(res))
}

heavisine = function(x){
  x = (x-min(x))/(max(x)-min(x))
  res = 4 * sin(4 * pi * x) - sign(x - 0.3) - sign(0.72 - x)
  return(res + sign(min(res))*min(res))
}

bumps = function(x){
  x = (x-min(x))/(max(x)-min(x))
  x=as.matrix(x)
  t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 
         0.78, 0.81)
  h2 <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
  w <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 
         0.005, 0.008, 0.005)
  res = 0
  for (i in seq(1, length(h2))) {
    res = res + h2[i] * pmax(0, (1 - abs((x - t[i])/w[i])))^4
  }
  return(res)
}

ya_doppler = function(x){
  x = (x-min(x))/(max(x)-min(x))
  eps <- 0.05
  res = sqrt(x * (1 - x)) * sin((2 * pi * (1 - eps))/(x + eps))
  return(res+sign(min(res))*min(res))
}

create_effects = function(num_effects = 10,data){
  for(i in 1:num_effects){
    data[paste0('effects_',i)] = rnorm(n=nrow(data),mean=0,sd=1)
    # data[paste0('effects_',i)] = runif(n=nrow(data),min=0,max=1)
  }
  return(data)
}

##### Lambda Sim #####
lambda_donoho = function(data,type,num_effects){
  ### 4 types, blocks, heavisine, bumps and ya_doppler
  res = 0
  if(type == "blocks"){
    for (i in 1:num_effects){
      res = res + blocks(data[paste0('effects_',i)])
    }
  } else if(type == "heavisine"){
    for (i in 1:num_effects){
      res = res + heavisine(data[paste0('effects_',i)])
    }
  } else if(type == "bumps"){
    for (i in 1:num_effects){
      res = res + bumps(data[paste0('effects_',i)])
    }
  } else if(type == "ya_doppler"){
    for (i in 1:num_effects){
      res = res + ya_doppler(data[paste0('effects_',i)])
    }
  }
  res = (res+ abs(data$day-4)/5)
  res = 2*res/max(res) 
  # print(res)
  return(as.matrix(exp(res)))
}

# plot(lambda_donoho(data,"heavisine"),type='l')

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

# mae = function(y_true,y_pred){
#   return(mean(abs(y_true-y_pred)))
# }
# 
# rmse = function(y_true,y_pred){
#   return(sqrt(mean((y_true-y_pred)**2)))
# }

mae = function(data){
  return(mean(abs(data$y_true-data$y_pred)))
}

rmse = function(data){
  return(sqrt(mean((data$y_true-data$y_pred)**2)))
}

p_rmse = function(data){
  df = data %>% 
    dplyr::select(ymd_date,y_true) %>%
    group_by(ymd_date) %>%
    summarise(y_true = max(y_true))
  
  dff = data %>% 
    dplyr::select(ymd_date,y_true,y_pred) %>%
    left_join(df,.,by=c("ymd_date"="ymd_date","y_true"="y_true")) %>%
    distinct(ymd_date,.keep_all = TRUE)
  return(rmse(dff))
}

deviance = function(data){
  return(2*sum(data$y_true*ifelse(data$y_true==0,0,log(data$y_true/data$y_pred)) - (data$y_true - data$y_pred) ))
}

### Plotting theme
theme_Publication = function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            #panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            #legend.key = element_rect(colour = NA),
            # legend.position = c(.95, .95),
            # legend.justification = c("right", "top"),
            # legend.box.just = "right",
            legend.position = "right",
            legend.direction = "vertical",
            #legend.key.size= unit(0.2, "cm"),
            #legend.key.size= unit(0.75, "cm"),
            legend.key.size=unit(3,"line"),
            # legend.margin = margin(6, 6, 6, 6),
            legend.title = element_text(face="italic"),
            panel.grid.major = element_line(colour = "black", size = 0.2,linetype="dashed"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

### END