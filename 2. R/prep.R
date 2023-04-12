source("2. R/utility.R")

#### Prep data

start_date = as.POSIXct(start_date, tz='UTC')
set.seed(20201)
data = tibble(date = seq(start_date,by = paste(as.character(step),"mins"), length.out=points_simu)) %>%
  filter(wday(date) %in% c(2,3,4,5,6)) %>%
  mutate(t=seq(0,nrow(.)-1),
         t_day = t %% (1440/step),
         tod = floor((hour(date)/24 + minute(date)/(60*24)) * (2^tod_res)),
         # tod = 0.5*rnorm(nrow(.),9,1) + 0.5*rnorm(nrow(.),17,3),
         # temp = -7.5*sin((2*pi/(1440/step))*t+pi/6)+10.5,
         temp = rnorm(nrow(.),15,3),
         ## traffic = -80*sin((2*pi/(1440/step))*(t+1000))+100,
         # t_traffic = (t_day-min(t_day))/(max(t_day)-min(t_day))*6-3,
         # traffic=100*exp(-abs(abs(t_traffic)-1)),
         traffic=rnorm(nrow(.),100,20),
         day = wday(date),
         ymd_date = as.Date(as.character(date))
  ) %>% create_effects(data=.,num_effects = num_effects)

lambdas = list("blocks"=lambda_donoho(data,"blocks",num_effects = num_effects),
               "doppler"=lambda_donoho(data,"ya_doppler",num_effects = num_effects),
               "heavisine"=lambda_donoho(data,"heavisine",num_effects = num_effects),
               "bumps"=lambda_donoho(data,"bumps",num_effects = num_effects),
               l_smooth(data,step=1),
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
# plot(data$date,data$lambda,type='l')
# plot(data$date,data$temp,type='l')
# plot(data$date,data$traffic,type='l')
#
plot(data$date,data$lambda,type='l')
# plot(data$date,lambdas[[1]],type='l')
# plot(data$date,lambdas[[2]],type='l')
# plot(data$date,lambdas[[3]],type='l')
# plot(data$date,lambdas[[4]],type='l')

# plot(data$temp,data$lambda)
# plot(data$traffic, data$lambda)
# plot(data$tod, data$lambda)
# plot(data$lambda,type='l')

wav_tod = WavD(data$tod, numLevels=tod_res, filterNumber = 1)
colnames(wav_tod) = paste0("wav_tod_",as.character(seq(1,2^tod_res-1)))
plot_wav(data$tod,numLevels = tod_res)

#####

simu = tibble(t_simu = thin_sim(data$lambda, ceiling(max(data['lambda'])), nrow(data), seed=NULL),
              event = 1) %>%
  left_join(data%>%dplyr::select(date,t),.,by=c("t"="t_simu")) %>%
  dplyr::select(date,event) %>%
  group_by(date) %>%
  summarise_all(sum) %>%
  mutate(event=ifelse(is.na(event),0,event))

plot(simu$date,simu$event)


train = left_join(data,simu,by="date") %>%
  mutate(event = ifelse(is.na(event),0,event)) %>%
  # cbind(.,wav_tod,wav_temp,wav_traffic) %>%
  mutate(pos = as.numeric(difftime(date,min(.$date),units="mins")),
         intercept = 1)

# plot(train$date,train$event,type='l')

