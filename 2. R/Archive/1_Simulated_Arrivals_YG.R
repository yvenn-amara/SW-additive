##### Simulating Arrivals based on a known intensity function #####

##### 1. Importing functions libraries
library(tidyverse)
library(lubridate)
library(wavethresh)
library(gamwave)


#install.packages("caret")
# ncvreg’, ‘grpSLOPE’, ‘grpreg’ 
#install.packages("~/Documents/PackagesR/stackeR_0.5.tar.gz", repos = NULL, type = "source")
#install.packages("~/Documents/PackagesR/gamwave_1.0.tar.gz", repos = NULL, type = "source")
#install.packages("wavethresh")
#library(gamwave)
# args=commandArgs(TRUE)


# points_simu = args[1]
# step = args[2]
points_simu=1440
step=1

dataset = "simulated_arrivals"
start_date = "2022-02-21"
source("2. R/0_Algos_Estim_NHPP.R")
#benchmark-master/2. R/2_NHPP_wavelets_per_dow.R
#####

start_date = as.POSIXct(start_date, tz='UTC') 

data = tibble(date = seq(start_date,by = paste(as.character(step),"mins"), length.out=points_simu)) %>%
  filter(wday(date) %in% c(2,3,4,5,6)) %>%
  mutate(t=seq(0,nrow(.)-1),
         tod = floor((hour(date)/24 + minute(date)/(60*24)) * 32),
         posweek = floor(((t*step %% 7200)/7200) * 32),
         indicator1 = ifelse(t %in% seq(450,650),1,0),
         indicator2 = ifelse(t %in% seq(850,1050),1,0)
         # tod32 = ((hour(date) + minute(date)/60)*100) %/% 75,
         # posweek32 = floor((t*step %% 7200)/225)
         # posweek32 = floor(rep(seq(0,7199),nrow(data)/7200)/225
         )

data$lambda = lambda_sim(data$t,step=step)
plot(data$date,data$lambda,type='l')
# decomposition = wd(data$lambda)

basis_tod = WavD(sort(unique(data$tod)), numLevels=5)
colnames(basis_tod) = paste0("wav_tod_",as.character(seq(1,31)))
matplot(basis_tod, type='l')

basis_posweek = WavD(sort(unique(data$posweek)), numLevels=5)
colnames(basis_posweek) = paste0("wav_posweek_",as.character(seq(1,31)))
matplot(basis_posweek, type='l')

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
  left_join(.,                  
            tibble(as.data.frame(basis_posweek)) %>% mutate(posweek = sort(unique(data$posweek))),
            by='posweek') %>%
  mutate(pos = as.numeric(difftime(date,min(.$date),units="mins")),
         intercept = 1)
  
# valid_nhpp(data,data$lambda)

#####

#### 0. ALL AT ONCE
# covar = as.matrix(train %>% dplyr::select(starts_with("wav"),intercept))
covar = as.matrix(train %>% dplyr::select(intercept, indicator1, indicator2))
# pos = (train %>% filter(event >0))$pos/step
# pos = diff((train %>% filter(event >0))$t)

# train_poi = train %>% filter(event > 0)
# pos = as.numeric(difftime(train_poi$date,min(data$date),units="mins"))/step
pos = (train %>% filter(event >0))$t+1

out_AAO = AAO(covar,pos)
out_AAO_pen= AAO(covar,pos,pen=10)

plot(data$lambda, type='l')
lines(out_AAO$lambdafit,type='l', col='red')
#lines(out_AAO_pen$lambdafit,type='l', col='green')

mae(data$lambda,out_AAO$lambdafit)
mae(data$lambda,out_AAO_pen$lambdafit)

test <-  covar%*%out_AAO$fullcoef
lines(exp(test), type='l')
range(exp(test)-out_AAO$lambdafit)

out_AAO_1 = AAO(covar[,1],pos)

mean(data$lambda)

exp(out_AAO$fullcoef[1])


# valid_nhpp(data,out_AAO$lambdafit, step='15 min')
# valid_nhpp(data,out_AAO_pen$lambdafit, step='15 min')

#### 1. BACKFITTING

effects = list()
#effects[[1]] = covar

effects[[1]] = as.matrix(train %>% dplyr::select(intercept))
effects[[2]] = as.matrix(train %>% dplyr::select(intercept, indicator1))
effects[[3]] = as.matrix(train %>% dplyr::select(intercept, indicator2))

effects[[1]] = as.matrix(train %>% dplyr::select(intercept,indicator1))
effects[[2]] = as.matrix(train %>% dplyr::select(indicator2))


effects[[1]] = as.matrix(train %>% dplyr::select(intercept))
effects[[2]] = as.matrix(train %>% dplyr::select(indicator1,indicator2))

effects = list()
effects[[1]] = as.matrix(train %>% dplyr::select(intercept))
effects[[2]] = as.matrix(train %>% dplyr::select(intercept, indicator1))
effects[[3]] = as.matrix(train %>% dplyr::select(intercept, indicator2))

effects = list()
effects[[1]] = as.matrix(train %>% dplyr::select(intercept, indicator1))
effects[[2]] = as.matrix(train %>% dplyr::select(intercept, indicator2))


# effects[[1]] = as.matrix(train %>% dplyr::select(intercept))
# effects[[2]] = as.matrix(train %>% dplyr::select(intercept,indicator2))



#effects[[2]] = as.matrix(train %>% dplyr::select(indicator))
# effects[[2]] = as.matrix(train %>% dplyr::select(tod) %>% mutate(tod=tod*tod))
#effects[[1]] = as.matrix(train %>% dplyr::select(intercept))


# effects[[1]] = as.matrix(train %>% dplyr::select(starts_with("wav_tod")))
# effects[[2]] = as.matrix(train %>% dplyr::select(starts_with("wav_posweek")))
# effects[[3]] = as.matrix(train %>% dplyr::select(intercept))

out_BAC = BAC(effects,pos)
plot(data$lambda, type='l')
lines(out_BAC$lambdafit,type='l', col='purple')

test <-  covar%*%unlist(out_BAC$fullcoef)


range(exp(test)-out_BAC$lambdafit)

out_AAO$fullcoef
out_BAC$fullcoef%>%unlist



cor(covar)

mae(data$lambda,out_BAC$lambdafit)
#####

head(data)
?glm


reg <- lm(lambda~indicator1+indicator2, data=data)
plot(data$lambda, type='l')
lines(reg$fitted.values, col='red')


plot(data$event, type='b')


?filter
K <- 60
filtre <- array(1/K, dim=K)
smooth_events <- stats::filter(data$event, filter=filtre, method='convolution', circular = T)
plot(data$lambda, type='l')
lines(smooth_events, col='red')

reg <- lm(smooth_events~indicator1+indicator2, data=data)
plot(data$lambda, type='l')
lines(reg$fitted.values, col='red')


plot(data$lambda, type='l')


poisson_reg <- glm(event~indicator1+indicator2, data=data, family='poisson')

gam_reg <- mgcv::gam(event~s(tod, k=5), data=data, family='poisson')

poisson_reg2 <- glm(event~as.factor(tod), data=data, family='poisson')



plot(data$lambda, type='l')
lines(out_AAO$lambdafit,type='l', col='red')
lines(poisson_reg$fitted.values, col='green')
lines(reg$fitted.values, col='purple')
lines(gam_reg$fitted.values, col='blue')


mae(data$lambda, reg$fitted.values)
mae(data$lambda, poisson_reg$fitted.values)
mae(data$lambda, out_AAO$lambdafit)
mae(data$lambda, gam_reg$fitted.values)

?glm

poisson_reg1<- glm(event~1, data=data, family='poisson')

out_AAO$fullcoef
poisson_reg$coefficients



#### #### #### #### #### #### #### #### #### Bac avec glm


####  sans intercept
poisson_reg1 <- glm(event~1, data=data, family='poisson')
bac_glm1 <- glm(event~indicator1-1, data=data, family='poisson', offset=poisson_reg1$coefficients[1]%>%rep(nrow(data)))
bac_glm2 <- glm(event~indicator2-1, data=data, family='poisson', offset=
                  poisson_reg1$coefficients[1]%>%rep(nrow(data)) + data$indicator1*bac_glm1$coefficients[1])

plot(data$lambda, type='l')
lines(poisson_reg1$fitted.values)
lines(poisson_reg$fitted.values, col='green')
lines(bac_glm2$fitted.values, col='purple')


# ####  avec offset
# poisson_reg1 <- glm(event~1, data=data, family='poisson')
# bac_glm1 <- glm(event~indicator1, data=data, family='poisson', offset=poisson_reg1$coefficients[1]%>%rep(nrow(data)))
# bac_glm2 <- glm(event~indicator2, data=data, family='poisson', offset=
#                   poisson_reg1$coefficients[1]%>%rep(nrow(data)) + data$indicator1*bac_glm1$coefficients[2])
# 
# plot(data$lambda, type='l')
# lines(poisson_reg$fitted.values, col='green')
# lines(bac_glm2$fitted.values, col='purple')
# 



#### #### #### #### #### #### #### #### #### #### 
####  en centrant
#### #### #### #### #### #### #### #### #### #### 

indicator1_center <- data$indicator1-mean(data$indicator1)
indicator2_center <- data$indicator2-mean(data$indicator2)

poisson_reg1 <- glm(event~1, data=data, family='poisson')
bac_glm1 <- glm(event~indicator1_center-1, data=data, family='poisson', offset=poisson_reg1$coefficients[1]%>%rep(nrow(data)))
bac_glm2 <- glm(event~indicator2_center-1, data=data, family='poisson', offset=
                  poisson_reg1$coefficients[1]%>%rep(nrow(data)) + indicator1_center*bac_glm1$coefficients[1])

plot(data$lambda, type='l')
lines(poisson_reg$fitted.values, col='green')
lines(bac_glm2$fitted.values, col='purple')


















