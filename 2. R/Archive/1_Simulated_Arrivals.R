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
points_simu=480*5
step=15

tod_res = 5
# posweek_res = 2


dataset = "simulated_arrivals"
start_date = "2022-02-21"
source("2. R/0_Algos_Estim_NHPP.R")

#####

start_date = as.POSIXct(start_date, tz='UTC') 

data = tibble(date = seq(start_date,by = paste(as.character(step),"mins"), length.out=points_simu)) %>%
  filter(wday(date) %in% c(2,3,4,5,6)) %>%
  mutate(t=seq(0,nrow(.)-1),
         tod = floor((hour(date)/24 + minute(date)/(60*24)) * (2^tod_res)),
         # posweek = floor((((t*step) %% 7200)/7200) * (2^posweek_res)),
         indicator_mon = ifelse(wday(date) == 2,1,0),
         indicator_tue = ifelse(wday(date) == 3,1,0),
         indicator_wed = ifelse(wday(date) == 4,1,0),
         indicator_thu = ifelse(wday(date) == 5,1,0),
         indicator_fri = ifelse(wday(date) == 6,1,0)
         )

data$lambda = lambda_sim(data$t,step=step)
plot(data$date,data$lambda,type='l')
# decomposition = wd(data$lambda)

basis_tod = WavD(sort(unique(data$tod)), numLevels=tod_res)
colnames(basis_tod) = paste0("wav_tod_",as.character(seq(1,2^tod_res-1)))
matplot(basis_tod, type='l')


# basis_posweek = WavD(sort(unique(data$posweek)), numLevels=posweek_res)
# colnames(basis_posweek) = paste0("wav_posweek_",as.character(seq(1,2^posweek_res-1)))
# matplot(basis_posweek, type='l')

#####

simu = tibble(date = min(data$date) + step*minutes(thin_sim(data$lambda, ceiling(max(data['lambda'])), nrow(data))),
              event = 1) %>%
  mutate(date=ceiling_date(date, paste(as.character(step),"mins"))) %>%
  group_by(date) %>%
  summarise_all(sum)

data = left_join(data,simu,by="date") %>%
  mutate(event = ifelse(is.na(event),0,event))

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

#### 0. ALL AT ONCE
# covar = as.matrix(train %>% dplyr::select(indicator,intercept))
# pos = (train %>% filter(event >0))$t+1
# 
# out_AAO = AAO(covar,pos)
# out_AAO_pen= AAO(covar,pos,pen=10)
# 
# plot(data$lambda, type='l')
# lines(out_AAO$lambdafit,type='l', col='red')
# lines(out_AAO_pen$lambdafit,type='l', col='green')
# 
# mae(data$lambda,out_AAO$lambdafit)
# mae(data$lambda,out_AAO_pen$lambdafit)

# valid_nhpp(data,out_AAO$lambdafit, step='15 min')
# valid_nhpp(data,out_AAO_pen$lambdafit, step='15 min')

#### 1. BACKFITTING

# X = as.matrix(train %>% dplyr::select(starts_with("indicator")))
# cvfit = cv.glmnet(x=X, y=train$event, nfolds = 10, family = "poisson", type.measure="deviance", lambda=seq(0.0001,10, length=100))
# cvfit$lambda.min
# plot(cvfit)
# 
# glmnetfit = glmnet(x=X, y=train$event, lambda=cvfit$lambda.min)
# glmnetfit$beta

# glm(event ~ indicator_mon + indicator_tue + indicator_wed + indicator_thu + indicator_fri-1, 
#     data=train)

effects = list()

# effects[[1]] = as.matrix(train %>% dplyr::select(intercept))
effects[[1]] = as.matrix(train %>% dplyr::select(starts_with("wav_tod")))
effects[[2]] = as.matrix(train %>% dplyr::select(starts_with("indicator")))
# effects[[3]] = as.matrix(train %>% dplyr::select(indicator_tue))
# effects[[4]] = as.matrix(train %>% dplyr::select(indicator_wed))
# effects[[5]] = as.matrix(train %>% dplyr::select(indicator_thu))
# effects[[6]] = as.matrix(train %>% dplyr::select(indicator_fri))

out_BAC_GLM = BAC_GLM(effects,train)
plot(data$lambda, type='l')
lines(out_BAC_GLM$lambdafit,type='l', col='green')

mae(data$lambda,out_BAC_GLM$lambdafit)


out_BAC_GLMnet = BAC_GLMnet(effects,train)
plot(data$lambda, type='l')
lines(out_BAC_GLMnet$lambdafit,type='l', col='blue')
mae(data$lambda,out_BAC_GLMnet$lambdafit)

plot(out_BAC_GLMnet$cvfits[[1]])
plot(out_BAC_GLMnet$cvfits[[2]])

plot(out_BAC_GLMnet$cvfits[[1]]$glmnet.fit)
# abline(v = out_BAC_GLMnet$cvfits[[1]]$lambda.min)

plot(out_BAC_GLMnet$cvfits[[2]]$glmnet.fit)
# abline(v = out_BAC_GLMnet$cvfits[[2]]$lambda)

out_BAC_GLMnet$fullcoef

#####
