library(tidyverse)
library(lubridate)
library(riem)
library(zoo)
source("2. R/0_Algos_Estim_NHPP_v4.R")
source("2. R/utility.R")

step = 60 # in mins from 1 to 60
step_string = paste(as.character(step),"mins")

df = read.csv("1. Preprocessed Data/domestics.csv")
# df$Start = as.POSIXct(df$Start,tz="US/Pacific")
df$Start = as.POSIXct(df$Start,tz="Europe/London")


#### Top 10 Largest cities in the UK in 1,000 (April 2022)
london = 8900 #EGLC
birmingham = 1150 #EGBB
glasgow = 612 #EGPF
liverpool = 579 #EGGP
bristol = 572 #EGGD
manchester = 554 #EGCC
# sheffield = 544 #Couldn't find
leeds = 503 #EGNM
edinburgh = 508 #EGPH
# leicester = 471 #Couldn't find

total_pop = london + birmingham + glasgow + liverpool + bristol + manchester + leeds + edinburgh

stations = list("EGLC","EGBB","EGPF","EGGP","EGGD","EGCC","EGNM","EGPH")
weights = list(london/total_pop,birmingham/total_pop,glasgow/total_pop,liverpool/total_pop,bristol/total_pop,manchester/total_pop,leeds/total_pop,edinburgh/total_pop)


generate_toy = function(start,end,step){
  year_start = floor_date(start, "year")
  year_end = ceiling_date(end, "year")
  dates = tibble(date=seq(year_start,year_end,step))
  toy = NULL
  for (i in unique(year(dates$date))){
    temp = dates %>%
      filter(year(date) == i) %>%
      nrow()
    toy = c(toy,seq(0,1,length=temp))
  }
  return(dates %>% mutate(toy=toy))
}

prep = function(df,stations,trend=FALSE,weights=NULL,tzone){
  
  #### Arrivals
  dates_step = tibble(Start = seq(floor_date(min(df$Start),unit = "day"),
                                   ceiling_date(max(df$Start),unit = "day")-60*step,
                                   step_string))
  arrivals = df %>% 
    mutate(Start = floor_date(Start,step_string)) %>%
    group_by(Start) %>%
    summarise(event = n()) %>%
    left_join(dates_step,.) %>% 
    mutate(event = ifelse(is.na(event),0,event)) %>% 
    arrange(Start)
  
  #### Weather Data
  
  if(length(weights) != length(stations) & length(stations)>1){
    print("Error: Multiple stations given without specifying weights")
    break
  }
  
  if(is.null(weights) & length(stations)==1){
    weights = list(1)
  }
  
  # dates_1min = tibble(Start = seq(min(measurements$valid),max(measurements$valid),"1 min"))
  
  dates_1min = tibble(Start = seq(min(arrivals$Start),max(arrivals$Start),"1 min"))
  
  measurements = riem_measures(station = stations[[1]], 
                               date_start = format(min(with_tz(arrivals$Start,tzone="UTC") ), format = "%Y-%m-%d") , 
                               date_end = format(max(with_tz(arrivals$Start,tzone="UTC") ) +24*3600, format = "%Y-%m-%d") )
  
  weather = measurements %>% 
    dplyr::select(valid,tmpf) %>%
    left_join(dates_1min,.,by=c("Start"="valid")) %>%
    mutate(tmpf = na.spline(tmpf),
           Start = with_tz(Start, tzone = tzone))
  weather$tmpf = weather$tmpf*weights[[1]]
  
  if (length(stations) > 1){
    for (i in 2:length(stations)){
      measurements = riem_measures(station = stations[[i]], 
                                   date_start = format(min(with_tz(arrivals$Start,tzone="UTC") ), format = "%Y-%m-%d") , 
                                   date_end = format(max(with_tz(arrivals$Start,tzone="UTC") )+24*3600, format = "%Y-%m-%d") )
      
      #### Remove outliers
      measurements = measurements %>% filter(tmpf < 150) %>% filter(tmpf > -150)
      
      temp_weather = measurements %>% 
        dplyr::select(valid,tmpf) %>%
        left_join(dates_1min,.,by=c("Start"="valid")) %>%
        mutate(tmpf = na.spline(tmpf),
               Start = with_tz(Start, tzone = tzone))
      weather$tmpf = weather$tmpf + temp_weather$tmpf*weights[[i]]
      
    }
  }
  
  #### Getting toy
  toy = generate_toy(min(df$Start),max(df$Start),step*60)
  
  #### Merging datasets and adding calendar variables
  
  final = left_join(arrivals,weather) %>%
    mutate(tod = floor((hour(Start)/24 + minute(Start)/(60*24)) * (60*24/step)),
           dow = wday(Start),
           temp = (tmpf-32) * 5/9,
           ymd_date = format(Start,format="%Y-%m-%d")) %>%
    rename("date"="Start") %>%
    left_join(.,toy)
    
  
  # #### Remove trend if specified
  if(trend){
    final$event = final$event%/%gam(event ~ s(toy), data=final)$fitted.values
  }
  
  return(final)
  
}

save_metrics = function(OBOs,OBOw,OBOsw,BACs,BACw,BACsw,test){
  #### Saving results
  metrics_train = tibble(model=c("OBOs","OBOw","OBOsw","BACs","BACw","BACsw"),
                   mae=c(mae(OBOs$final),
                         mae(OBOw$final),
                         mae(OBOsw$final),
                         mae(BACs$final),
                         mae(BACw$final),
                         mae(BACsw$final)),
                   rmse=c(rmse(OBOs$final),
                          rmse(OBOw$final),
                          rmse(OBOsw$final),
                          rmse(BACs$final),
                          rmse(BACw$final),
                          rmse(BACsw$final)),
                   p_rmse=c(p_rmse(OBOs$final),
                            p_rmse(OBOw$final),
                            p_rmse(OBOsw$final),
                            p_rmse(BACs$final),
                            p_rmse(BACw$final),
                            p_rmse(BACsw$final)),
                   deviance=c(deviance(OBOs$final),
                              deviance(OBOw$final),
                              deviance(OBOsw$final),
                              deviance(BACs$final),
                              deviance(BACw$final),
                              deviance(BACsw$final)),
                   dof=c(OBOs$dof,
                         OBOw$dof,
                         OBOsw$dof,
                         BACs$dof,
                         BACw$dof,
                         BACsw$dof),
                   runtime=c(OBOs$runtime,
                             OBOw$runtime,
                             OBOsw$runtime,
                             BACs$runtime,
                             BACw$runtime,
                             BACsw$runtime),
                   param=c(OBOs$num_param,
                           OBOw$num_param,
                           OBOsw$num_param,
                           BACs$num_param,
                           BACw$num_param,
                           BACsw$num_param),
                   param_non_zeros=c(OBOs$num_param_non_zeros,
                                     OBOw$num_param_non_zeros,
                                     OBOsw$num_param_non_zeros,
                                     BACs$num_param_non_zeros,
                                     BACw$num_param_non_zeros,
                                     BACsw$num_param_non_zeros),
                   iter = c(OBOs$iterations,
                            OBOw$iterations,
                            OBOsw$iterations,
                            BACs$iterations,
                            BACw$iterations,
                            BACsw$iterations)
                   
  ) %>%
    mutate(runtime = as.numeric(runtime))
  
  #### Test
  OBOw_test = predict.sw(OBOw,test)
  OBOs_test = predict.sw(OBOs,test)
  OBOsw_test = predict.sw(OBOsw,test)
  BACw_test = predict.sw(BACw,test)
  BACs_test = predict.sw(BACs,test)
  BACsw_test = predict.sw(BACsw,test)
  
  metrics_test = tibble(model=c("OBOs","OBOw","OBOsw","BACs","BACw","BACsw"),
                         mae=c(mae(OBOs_test),
                               mae(OBOw_test),
                               mae(OBOsw_test),
                               mae(BACs_test),
                               mae(BACw_test),
                               mae(BACsw_test)),
                         rmse=c(rmse(OBOs_test),
                                rmse(OBOw_test),
                                rmse(OBOsw_test),
                                rmse(BACs_test),
                                rmse(BACw_test),
                                rmse(BACsw_test)),
                         p_rmse=c(p_rmse(OBOs_test),
                                  p_rmse(OBOw_test),
                                  p_rmse(OBOsw_test),
                                  p_rmse(BACs_test),
                                  p_rmse(BACw_test),
                                  p_rmse(BACsw_test)),
                         deviance=c(deviance(OBOs_test),
                                    deviance(OBOw_test),
                                    deviance(OBOsw_test),
                                    deviance(BACs_test),
                                    deviance(BACw_test),
                                    deviance(BACsw_test))
                         
  )
  
  return(list(metrics_train,metrics_test))
}

# palo_alto = prep(df,station="PAO",tzone="US/Pacific")
domestics = prep(df,trend=TRUE,stations=stations,weights=weights,tzone="Europe/London")

# filtered = palo_alto %>% 
#   filter(hour(date)>=8 & hour(date)<22) # %>%
  # filter(date <= "2016-07-31" & date >= "2016-01-01")

### dom filter
filtered = domestics %>% 
  filter(dow %in% c(2,3,4,5,6)) # %>%
  # filter(hour(date)>=8 & hour(date)<22)

# plot(filtered$temp,filtered$event)
# plot(filtered$dow,filtered$event)
# plot(filtered$tod,filtered$event)
# plot(filtered$toy,filtered$event)

train = filtered %>%
  filter(date < "2017-09-04")
test = filtered %>%
  filter(date <= "2017-09-11" & date >= "2017-09-04")

####
# fmla_w = as.formula("event ~ factor(dow) + wp(tod) + wp(toy) + wp(temp)")
# fmla_s = as.formula("event ~ factor(dow) + s(tod,k=10) + s(toy,k=10) + s(temp,k=10)")
# fmla_sw = as.formula("event ~ factor(dow) + s(tod,k=10) + wp(tod) + s(toy,k=10) + wp(toy) + s(temp,k=10) + wp(temp)")

fmla_w = as.formula("event ~ factor(dow) + wp(tod) + wp(temp)")
fmla_s = as.formula("event ~ factor(dow) + s(tod,k=10) + s(temp,k=10)")
fmla_sw = as.formula("event ~ factor(dow) + s(tod,k=10) + wp(tod) + s(temp,k=10) + wp(temp)")


train_dates = list('2017-09-04','2017-09-11','2017-09-18','2017-09-25',
                   '2017-10-02','2017-10-09','2017-10-16','2017-10-23',
                   '2017-10-30','2017-11-06','2017-11-13','2017-11-20',
                   '2017-11-27','2017-12-04','2017-12-11','2017-12-18',
                   '2017-12-25')
test_dates = list('2017-09-11','2017-09-18','2017-09-25',
                  '2017-10-02','2017-10-09','2017-10-16','2017-10-23',
                  '2017-10-30','2017-11-06','2017-11-13','2017-11-20',
                  '2017-11-27','2017-12-04','2017-12-11','2017-12-18',
                  '2017-12-25','2018-01-01')

# train_dates = list('2017-02-04')
# test_dates = list('2017-02-11')

training_routine = function(train_dates,test_dates){
  metrics_train = NULL
  metrics_test = NULL
  for(i in 1:length(train_dates)){
    train = filtered %>%
      filter(date < train_dates[[i]])
    test = filtered %>%
      filter(date < test_dates[[i]] & date >= train_dates[[i]])
    
    BACw = rig_BAC_final(fmla_w, train,simu=FALSE,get_full_rank = TRUE, res_max=6)
    BACs = rig_BAC_final(fmla_s, train,simu=FALSE,get_full_rank = TRUE, res_max=6)
    BACsw = rig_BAC_final(fmla_sw, train,simu=FALSE,get_full_rank = TRUE, res_max=6)
    
    OBOw = OBO_pen(fmla_w, train,simu=FALSE, eps=1e-3, maxiter=100,get_full_rank=TRUE, res_max=6)
    OBOs = OBO_pen(fmla_s, train,simu=FALSE, eps=1e-3, maxiter=100,get_full_rank=TRUE, res_max=6)
    OBOsw = OBO_pen(fmla_sw, train,simu=FALSE, eps=1e-3, maxiter=100,get_full_rank=TRUE, res_max=6)
    
    final = save_metrics(OBOs,OBOw,OBOsw,BACs,BACw,BACsw,test)
    metrics_train = bind_rows(metrics_train,final[[1]])
    metrics_test = bind_rows(metrics_test,final[[2]])
    
    ### plot pred
    pred_table = bind_rows(predict.sw(BACw,test) %>% mutate(model = "BACw") %>% rename("pred"="y_pred") %>% dplyr::select(date,ymd_date,pred,model),
                           predict.sw(BACs,test) %>% mutate(model = "BACs") %>% rename("pred"="y_pred") %>% dplyr::select(date,ymd_date,pred,model),
                           predict.sw(BACsw,test) %>% mutate(model = "BACsw") %>% rename("pred"="y_pred") %>% dplyr::select(date,ymd_date,pred,model),
                           predict.sw(OBOw,test) %>% mutate(model = "OBOw") %>% rename("pred"="y_pred") %>% dplyr::select(date,ymd_date,pred,model),
                           predict.sw(OBOs,test) %>% mutate(model = "OBOs") %>% rename("pred"="y_pred") %>% dplyr::select(date,ymd_date,pred,model),
                           predict.sw(OBOsw,test) %>% mutate(model = "OBOsw") %>% rename("pred"="y_pred") %>% dplyr::select(date,ymd_date,pred,model),
                           predict.sw(OBOs,test) %>% mutate(model = "observed") %>% rename("pred"="y_true") %>% dplyr::select(date,ymd_date,pred,model))
    
    pdf(file=paste0("5. Outputs/general/","pred_plot",i,".pdf"),width=10,height=7)
    print(ggplot(data = pred_table, aes(x=date,y=pred)) +
      geom_line(aes(color = model),size=1)+
      ylab('Number of Arrivals') +
      xlab('Date [hours]') +
      theme_Publication())
    dev.off()
    
  }
  return(list(metrics_train,metrics_test))
}

met = training_routine(train_dates,test_dates)

write.csv(met[[1]], "5. Outputs/general/metrics_train.csv",row.names=FALSE)
write.csv(met[[2]], "5. Outputs/general/metrics_test.csv",row.names=FALSE)

# w = predict.sw(out_w,train)
# s = predict.sw(out_s,train)
# sw = predict.sw(out_sw,train)
# w = predict.sw(out_w,test)
# s = predict.sw(out_s,test)
# sw = predict.sw(out_sw,test)

# plot(out_sw$models[[1]])
# plot(out_sw$models[[2]])
# plot(out_sw$models[[3]])
# plot(out_sw$models[[4]])
# plot(out_sw$models[[4]]$glmnet.fit)
# plot(out_sw$models[[5]]$glmnet.fit)
# plot_wav(x = train$tod, numLevels = log2(length(out_sw$fullcoef[[4]])+1),filterNumber = 1, coef_is_null = out_sw$fullcoef[[4]])
# plot_wav(x = train$toy, numLevels = log2(length(out_sw$fullcoef[[5]])+1),filterNumber = 1, coef_is_null = out_sw$fullcoef[[5]])

# aa = gam(fmla,data=filtered,family=poisson)
# sqrt( mean((aa$y-aa$fitted.values)^2))

# plot(x=w$date,y=w$y_true,type='l')
# lines(w$date,w$y_pred,col='green')
# 
# plot(x=s$date,y=s$y_true,type='l')
# lines(s$date,s$y_pred,col='green')
# 
# plot(x=sw$date,y=sw$y_true,type='l')
# lines(sw$date,sw$y_pred,col='green')

#### One-off PLOTS FOR PAPER
# w = 15
# a = 1 + (24*5)*w 
# b = a + 24*5-1
# 
# fpe = tibble(date = out_sw$final$date[a:b],
#              y_true = out_sw$final$y_true[a:b],
#              lin_part = out_sw$lin_pred[a:b],
#              spl_part = out_sw$spl_pred[a:b],
#              wav_part = out_sw$wav_pred[a:b])

# Sys.setlocale("LC_TIME", "English")
# eff = ggplot(data = fpe, aes(x=date))+
#   geom_line(aes(y=y_true,colour="observed"),size=1) +
#   geom_line(aes(y=exp(lin_part), colour="lin"),size=1) +
#   geom_line(aes(y=exp(lin_part+spl_part), colour="lin+spl"),size=1) +
#   geom_line(aes(y=exp(lin_part+spl_part+wav_part), colour="lin+spl+wav"),size=1)+
#   scale_colour_manual("", 
#                       breaks = c("observed","lin", "lin+spl", "lin+spl+wav"),
#                       values = c("black","grey", "red", "steelblue")) +
#   xlab("Date [hours]") +
#   scale_y_continuous("Number of arrivals") +
#   theme_Publication()
# 
# pdf(file=paste0("5. Outputs/general/","random_week_fit.pdf"),width=10,height=7)
# eff
# dev.off()
# 
# pdf(file=paste0("5. Outputs/general/","wav_before_fit.pdf"),width=10,height=7)
# plot_wav(x = train$tod, numLevels = log2(length(out_sw$fullcoef[[5]])+1),filterNumber = 1)
# dev.off()
# 
# pdf(file=paste0("5. Outputs/general/","wav_after_fit.pdf"),width=10,height=7)
# plot_wav(x = train$tod, numLevels = log2(length(out_sw$fullcoef[[5]])+1),filterNumber = 1, coef_is_null = out_sw$fullcoef[[5]])
# dev.off()
# 
# pdf(file=paste0("5. Outputs/general/","wavfinal_after_fit.pdf"),width=10,height=7)
# plot_wav2(x = train$tod, numLevels = log2(length(out_sw$fullcoef[[5]])+1),filterNumber = 1, coef_is_null = out_sw$fullcoef[[5]])
# dev.off()
#### END PLOTS FOR PAPER


# plot_wav(x = train$toy*365, numLevels = log2(length(out_sw$fullcoef[[6]])+1),filterNumber = 1)
# plot_wav(x = train$toy*365, numLevels = log2(length(out_sw$fullcoef[[6]])+1),filterNumber = 1, coef_is_null = out_sw$fullcoef[[6]])


#### END