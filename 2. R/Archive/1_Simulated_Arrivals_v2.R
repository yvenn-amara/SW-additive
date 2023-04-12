##### Simulating Arrivals based on a known intensity function #####

##### 1. Importing functions libraries
library(tidyverse)
library(lubridate)
library(wavethresh)
library(gamwave)
library(glmnet)
# args=commandArgs(TRUE)


# points_simu = args[1]
# step = args[2]
# lambda_num = args[3]
step=15
points_simu=(7*1440/step)
lambda_num = 1

tod_res = 5
# posweek_res = 2


dataset = "simulated_arrivals"
start_date = "2022-02-21"
source("2. R/0_Algos_Estim_NHPP_v2.R")

#####

start_date = as.POSIXct(start_date, tz='UTC') 

data = tibble(date = seq(start_date,by = paste(as.character(step),"mins"), length.out=points_simu)) %>%
  filter(wday(date) %in% c(2,3,4,5,6)) %>%
  mutate(t=seq(0,nrow(.)-1),
         tod = floor((hour(date)/24 + minute(date)/(60*24)) * (2^tod_res)),
         day = wday(date),
         # posweek = floor((((t*step) %% 7200)/7200) * (2^posweek_res)),
         indicator_mon = ifelse(wday(date) == 2,1,0),
         indicator_tue = ifelse(wday(date) == 3,1,0),
         indicator_wed = ifelse(wday(date) == 4,1,0),
         indicator_thu = ifelse(wday(date) == 5,1,0),
         indicator_fri = ifelse(wday(date) == 6,1,0)
         )

lambdas = list(lambda_sim(data$t,step=step),
               lambda_sim2(data$t,step=step),
               lambda_sim3(data$t,step=step, peaks=2),
               lambda_sim3(data$t,step=step, peaks=5),
               lambda_sim3(data$t,step=step, peaks=10),
               lambda_peaky(data$t,step=step),
               lambda_mix(data$t,step=step))
data$lambda = lambdas[[lambda_num]]
plot(data$date,data$lambda,type='l')

basis_tod = WavD(sort(unique(data$tod)), numLevels=tod_res, filterNumber = 1)
colnames(basis_tod) = paste0("wav_tod_",as.character(seq(1,2^tod_res-1)))
# matplot(basis_tod, type='l')

#####

# simu = tibble(date = min(data$date) + step*minutes(thin_sim(data$lambda, ceiling(max(data['lambda'])), nrow(data))),
#               event = 1) %>%
#   mutate(date=ceiling_date(date, paste(as.character(step),"mins"))) %>%
#   group_by(date) %>%
#   summarise_all(sum)

simu = tibble(t_simu = thin_sim(data$lambda, ceiling(max(data['lambda'])), nrow(data)),
              event = 1) %>%
  inner_join(data%>%dplyr::select(date,t),.,by=c("t"="t_simu")) %>%
  dplyr::select(date,event) %>%
  group_by(date) %>%
  summarise_all(sum)

data = left_join(data,simu,by="date") %>%
  mutate(event = ifelse(is.na(event),0,event))

plot(data$date,data$lambda,type='l',ylim=c(min(data$lambda,data$event),max(data$lambda,data$event)))
points(data$date,data$event)

train = left_join(data,
                  tibble(as.data.frame(basis_tod)) %>% mutate(tod = sort(unique(data$tod))),
                  by='tod'
                  ) %>% 
  # left_join(.,                  
  #           tibble(as.data.frame(basis_posweek)) %>% mutate(posweek = sort(unique(data$posweek))),
  #           by='posweek') %>%
  mutate(pos = as.numeric(difftime(date,min(.$date),units="mins")),
         intercept = 1)
  
# valid_nhpp(data,data$lambda)

#####

effects = list()
# effects[[1]] = as.matrix(train %>% dplyr::select(intercept))
effects[[1]] = as.matrix(train %>% dplyr::select(starts_with("wav_tod")))
effects[[2]] = as.matrix(train %>% dplyr::select(starts_with("indicator")))
# effects[[3]] = as.matrix(train %>% dplyr::select(indicator_tue))
# effects[[4]] = as.matrix(train %>% dplyr::select(indicator_wed))
# effects[[5]] = as.matrix(train %>% dplyr::select(indicator_thu))
# effects[[6]] = as.matrix(train %>% dplyr::select(indicator_fri))

#### 0. ALL AT ONCE

start_time  = Sys.time()
out_AAO = AAO(effects,train)
AAO_time  = difftime(Sys.time(),start_time,units = "secs")

AAO_tibble = data %>%
  dplyr::select(date) %>%
  mutate(lambda = out_AAO$lambdafit,
         model = "AAO")


start_time  = Sys.time()
out_AAO_pen= AA0_pen(effects,train)
AAO_pen_time  = difftime(Sys.time(),start_time,units = "secs")

AAO_pen_tibble = data %>%
  dplyr::select(date) %>%
  mutate(lambda = out_AAO_pen$lambdafit,
         model = "AAO_pen")

# plot(data$lambda, type='l')
# lines(out_AAO$lambdafit,type='l', col='red')
# lines(out_AAO_pen$lambdafit,type='l', col='green')
# 
# mae(data$lambda,out_AAO$lambdafit)
# mae(data$lambda,out_AAO_pen$lambdafit)

# valid_nhpp(data,out_AAO$lambdafit, step='15 min')
# valid_nhpp(data,out_AAO_pen$lambdafit, step='15 min')

#### 1. BACKFITTING

start_time  = Sys.time()
out_BAC_GLM = BAC_GLM(effects,train)
BAC_GLM_time  = difftime(Sys.time(),start_time,units = "secs")

BAC_GLM_tibble = data %>%
  dplyr::select(date) %>%
  mutate(lambda = out_BAC_GLM$lambdafit,
         model = "BAC")

start_time  = Sys.time()
out_BAC_GLMnet = BAC_GLMnet(effects,train)
BAC_GLMnet_time  = difftime(Sys.time(),start_time,units = "secs")

BAC_GLM_tibble = data %>%
  dplyr::select(date) %>%
  mutate(lambda = out_BAC_GLM$lambdafit,
         model = "BAC")

# plot(data$lambda, type='l')
# lines(out_BAC_GLM$lambdafit,type='l', col='green')
# lines(out_BAC_GLMnet$lambdafit,type='l', col='blue')

# mae(data$lambda,out_BAC_GLM$lambdafit)
# mae(data$lambda,out_BAC_GLMnet$lambdafit)

plot(out_BAC_GLMnet$cvfits[[1]])
plot(out_BAC_GLMnet$cvfits[[2]])

plot(out_BAC_GLMnet$cvfits[[1]]$glmnet.fit)
# abline(v = out_BAC_GLMnet$cvfits[[1]]$lambda.min)

plot(out_BAC_GLMnet$cvfits[[2]]$glmnet.fit)
# abline(v = out_BAC_GLMnet$cvfits[[2]]$lambda)

out_BAC_GLMnet$fullcoef

#####

start_time  = Sys.time()
out_GAM = gam(event ~ s(tod) + factor(day), data=train,family=poisson)
GAM_time  = difftime(Sys.time(),start_time,units = "secs")

plot(data$lambda, type='l')
lines(predict(out_GAM,type="response"),type='l', col='green')
mae(data$lambda,predict(out_GAM,type="response"))
out_GAM$coefficients

#####

start_time  = Sys.time()
out_glob = BAC_GLMnet_glob(effects_glm = effects,effects_gam = train %>% dplyr::select(tod),train)
glob_time  = difftime(Sys.time(),start_time,units = "secs")

plot(data$lambda, type='l')
lines(out_glob$lambdafit,type='l', col='green')
mae(data$lambda,out_glob$lambdafit)

out_glob$fullcoef


#### Saving results

metrics = tibble(model=c("AAO","AAO_pen","BAC","BAC_pen","GAM","BAC_GW"),
       mae=c(mae(data$lambda,out_AAO$lambdafit),
             mae(data$lambda,out_AAO_pen$lambdafit),
             mae(data$lambda,out_BAC_GLM$lambdafit),
             mae(data$lambda,out_BAC_GLMnet$lambdafit),
             mae(data$lambda,predict(out_GAM,type="response")),
             mae(data$lambda,out_glob$lambdafit)),
       rmse=c(rmse(data$lambda,out_AAO$lambdafit),
              rmse(data$lambda,out_AAO_pen$lambdafit),
              rmse(data$lambda,out_BAC_GLM$lambdafit),
              rmse(data$lambda,out_BAC_GLMnet$lambdafit),
              rmse(data$lambda,predict(out_GAM,type="response")),
              rmse(data$lambda,out_glob$lambdafit)),
       runtime=c(AAO_time,
                 AAO_pen_time,
                 BAC_GLM_time,
                 BAC_GLMnet_time,
                 GAM_time,
                 glob_time)
       
       ) %>%
  mutate(runtime = as.numeric(runtime))


forplot = data %>%
  dplyr::select(date,lambda) %>%
  mutate(model="observed") %>%
  
  


ggplot(data = )

##### END
