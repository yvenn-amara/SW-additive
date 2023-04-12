##### Simulating Arrivals based on a known intensity function #####

##### 1. Importing functions libraries
library(tidyverse)
library(lubridate)
library(wavethresh)
library(gamwave)
library(glmnet)
args=commandArgs(TRUE)

step = as.numeric(args[1])
num_weeks = as.numeric(args[2])
lambda_num = as.numeric(args[3])
iter = as.numeric(args[4])
points_simu = (7*1440/step)*num_weeks
# step=15
# points_simu=(7*1440/step)*2
# lambda_num = 1
# iter = 1

tod_res = 5
# posweek_res = 2


dataset = "simulated_arrivals"
start_date = "2022-02-21"
source("2. R/0_Algos_Estim_NHPP_v2.R")
path = paste0("5. Outputs/",num_weeks,"_",lambda_num,"/")
dir.create(path)
#####

start_date = as.POSIXct(start_date, tz='UTC') 

data = tibble(date = seq(start_date,by = paste(as.character(step),"mins"), length.out=points_simu)) %>%
  filter(wday(date) %in% c(2,3,4,5,6)) %>%
  mutate(t=seq(0,nrow(.)-1),
         t_day = t %% (1440/step),
         tod = floor((hour(date)/24 + minute(date)/(60*24)) * (2^tod_res)),
         temp = -7.5*sin((2*pi/(1440/step))*t+pi/6)+10.5,
         # traffic = -80*sin((2*pi/(1440/step))*(t+1000))+100,
         t_traffic = (t_day-min(t_day))/(max(t_day)-min(t_day))*6-3,
         traffic=100*exp(-abs(abs(t_traffic)-1)),
         day = wday(date),
         indicator_mon = ifelse(wday(date) == 2,1,0),
         indicator_tue = ifelse(wday(date) == 3,1,0),
         indicator_wed = ifelse(wday(date) == 4,1,0),
         indicator_thu = ifelse(wday(date) == 5,1,0),
         indicator_fri = ifelse(wday(date) == 6,1,0)
         )

lambdas = list(l_smooth(data,step=1),
               l_peaky(data,step=1),
               l_both(data,step=1),
               lambda_sim(data$t,step=step),
               lambda_sim2(data$t,step=step),
               lambda_sim3(data$t,step=step, peaks=2),
               lambda_sim3(data$t,step=step, peaks=5),
               lambda_sim3(data$t,step=step, peaks=10),
               lambda_peaky(data$t,step=step),
               lambda_mix(data$t,step=step))
data$lambda = lambdas[[lambda_num]]
plot(data$date,data$lambda,type='l')
plot(data$date,data$temp,type='l')
plot(data$date,data$traffic,type='l')

plot(data$date,lambdas[[1]],type='l')
plot(data$date,lambdas[[2]],type='l')
plot(data$date,lambdas[[3]],type='l')

pdf(file=paste0(path,"intensity",num_weeks,"_",lambda_num,".pdf"),width=9,height=5)
ggplot(data = data, aes(x=date, y = lambda)) + 
  geom_line()+
  ylab('Lambda') +
  xlab("Date") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("")
dev.off()

basis_tod = WavD(sort(unique(data$tod)), numLevels=tod_res, filterNumber = 1)
colnames(basis_tod) = paste0("wav_tod_",as.character(seq(1,2^tod_res-1)))
# matplot(basis_tod, type='l')

#####

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
  mutate(pos = as.numeric(difftime(date,min(.$date),units="mins")),
         intercept = 1)
  
# valid_nhpp(data,data$lambda)

#####

effects = list()
effects[[1]] = as.matrix(train %>% dplyr::select(starts_with("wav_tod")))
effects[[2]] = as.matrix(train %>% dplyr::select(starts_with("indicator")))

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

BAC_GLMnet_tibble = data %>%
  dplyr::select(date) %>%
  mutate(lambda = out_BAC_GLMnet$lambdafit,
         model = "BAC_pen")

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

GAM_tibble = data %>%
  dplyr::select(date) %>%
  mutate(lambda = predict(out_GAM,type="response"),
         model = "GAM")

# plot(data$lambda, type='l')
# lines(predict(out_GAM,type="response"),type='l', col='green')
# mae(data$lambda,predict(out_GAM,type="response"))
# out_GAM$coefficients

#####

start_time  = Sys.time()
out_glob = BAC_GLMnet_glob(effects_glm = effects,effects_gam = train %>% dplyr::select(tod),train)
glob_time  = difftime(Sys.time(),start_time,units = "secs")

GLOB_tibble = data %>%
  dplyr::select(date) %>%
  mutate(lambda = out_glob$lambdafit,
         model = "BAC_GW")

# plot(data$lambda, type='l')
# lines(out_glob$lambdafit,type='l', col='green')
# mae(data$lambda,out_glob$lambdafit)
# 
# out_glob$fullcoef


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
                 glob_time),
         param=c(sum(unlist(out_AAO$fullcoef)!=0),
                 sum(unlist(out_AAO_pen$fullcoef)!=0),
                 sum(unlist(out_BAC_GLM$fullcoef)!=0),
                 sum(unlist(out_BAC_GLMnet$fullcoef)!=0),
                 length(unlist(coef(out_GAM))!=0 ),
                 sum(unlist(out_glob$fullcoef)!=0)
                 ),
       iter = c(NA,
                NA,
                out_BAC_GLM$iterations,
                out_BAC_GLMnet$iterations,
                NA,
                out_glob$iterations)
       
       ) %>%
  mutate(runtime = as.numeric(runtime))

write.csv(metrics,paste0("5. Outputs/",num_weeks,"_",lambda_num,"/metrics_step",step,"_nweeks",num_weeks,"_nlambda",lambda_num,"_iter",iter,".csv"), row.names = F)

forplot = data %>%
  dplyr::select(date,lambda) %>%
  mutate(model="observed") %>%
  bind_rows(.,AAO_tibble,AAO_pen_tibble,BAC_GLM_tibble,BAC_GLMnet_tibble,GAM_tibble,GLOB_tibble) %>%
  filter(date <= min(date)+1440*5*60)
  
if (iter ==1){
  pdf(file=paste0("5. Outputs/",num_weeks,"_",lambda_num,"/lineplot_step",step,"_nweeks",num_weeks,"_nlambda",lambda_num,"_iter",iter,".pdf"),width=9,height=5)
  ggplot(data = forplot, aes(x=date, y = lambda)) + 
    geom_line(aes(color = model)) +
    theme_bw()
  dev.off()
}


##### END
