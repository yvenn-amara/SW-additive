##### 1. Importing functions libraries #####
library(tidyverse)
library(lubridate)
library(NHPoisson)
library(mgcv)
library(plotly)

dataset = "simulated_data"
start_date = "2021-03-01"
end_date = "2021-05-02"
# end_date = "2021-08-01"
# end_date = "2022-03-01"

window=60
tz = "UTC"
source("2. R/0_Utility.R")

##### 2. Creating the intensity function of the process #####

# Covariates data 2 months
start_date = as.POSIXct(start_date, tz='UTC') 
end_date = as.POSIXct(end_date, tz='UTC')

data = tibble(date = seq(start_date,end_date,'1 mins')) %>%
  mutate(dow=wday(date),
         tod=hour(date),
         tod32 = ((hour(date) + minute(date)/60)*100) %/% 75)

# Effects
end = 10000
plot(data$date[1:end],lambda_tod32(data$tod32, data$dow)[1:end],type='l')
data['lambda'] = lambda_tod32(data$tod32, data$dow)
##### 3. Simulating the process #####
# simu = simNHP.fun(lambda,fixed.seed=204)
# simu_arrival = tibble(date=min(data$date)+simu$posNH*60)

# simu = thin_sim(lambda_tod32, 1, min(df_all$Start),nrow(df_all))
# simu_arrival = tibble(date=min(data$date)+simu*60)
# 
# write.csv(simu_arrival,paste("2. R/1. Outputs/simulated_data/sim_data_arrivals.csv",sep=''), row.names = F)

##### 4.Estimating uncompressible theoritical error #####

time.span = nrow(df_all) #All time

# df_all = df_all %>%
#   mutate(real_f_1 = dnorm(tod,9,1),
#          real_f_2 = dnorm(tod,12.5,1.5),
#          real_f_3 = dnorm(tod,19,2),
#          real_g_1 = dnorm(Day,2,1),
#          real_g_2 = dnorm(Day,4,1),
#          real_g_3 = dnorm(Day,6,1)
#          )

df_all = df_all %>%
  mutate(real_f_1 = dnorm(tod32,12,1),#dnorm(tod,9,1),
         #real_f_2 = dnorm(tod32,16,1.5),#dnorm(tod,12.5,1.5),
         real_f_3 = dnorm(tod32,25,2)#dnorm(tod,19,2)
         ) %>%
  mutate(real_week = ifelse(wday(Start) %in% c(2,3,4,5,6),1,0)#,
         #real_weekend = ifelse(wday(Start) %in% c(1,7),1,0)
         ) #%>%
  #bind_cols(.,fastDummies::dummy_cols(df_all %>% dplyr::select(Day),remove_first_dummy = T) %>% select(-Day) %>% rename_with( ~ paste0("real_", .x)))

#
covar = cbind(#df_all$Day[1:time.span],
  (df_all %>% dplyr::select(starts_with("real_")
                            #starts_with("wav_Day")
  ) %>% as.matrix())[1:time.span,])
colnames(covar) = NULL

# Check if covar is ill-conditioned
test = Matrix::rankMatrix(covar)[[1]] != ncol(covar)
test
if (test){
  covar = WeightIt::make_full_rank(covar,with.intercept = T)
}
# covar_corr = round(cor(covar),2)
Matrix::rankMatrix(covar)

start = as.list(rep(0,ncol(covar) + 1))
names(start) = paste("b",seq(0,ncol(covar)),sep="")
out = fitPP.fun(covariates = covar, posE = df_poisson$events[df_poisson$events < time.span],
                start=start, modSim=T, tind=T, minfun="nlminb"#, fixed = list("b0"=-2,"b1"=0.5)#,method='L-BFGS'
)

# saveRDS(out, paste("2. R/1. Outputs/",dataset,"/nhp/models/real",sep=''))

##### 5. Comparing artificial intensity with fitted intensities #####
# mod = readRDS(paste("2. R/1. Outputs/simulated_data/nhp/models/NHPP_wavtod_Day_fit_DaubEx1_1",sep=''))
mod = readRDS(paste("2. R/1. Outputs/simulated_data/nhp/models/my_NHPP_wavtod_Day_fit_DaubEx1_1",sep=''))
# mod = readRDS(paste("2. R/1. Outputs/simulated_data/nhp/models/NHPP_wavtod_Day_fit_DaubEx5_1",sep=''))
# mod = readRDS(paste("2. R/1. Outputs/simulated_data/nhp/models/NHPP_wavtod_Day_fit_DaubEx10_1",sep=''))
# mod = readRDS(paste("2. R/1. Outputs/simulated_data/nhp/models/NHPP_splinestod_Day_fit_cc_1",sep=''))
# mod = readRDS(paste("2. R/1. Outputs/",dataset,"/nhp/models/real",sep=''))

# data_filtered = data %>% filter(date %in% unique(df_all$Start))

# mae(data_filtered$lambda,mod@lambdafit)
# mape(data_filtered$lambda,mod@lambdafit)
# rmse(data_filtered$lambda,mod@lambdafit)

# lambdafit = mod@lambdafit[1:nrow(data)]#+0.0125

lambdafit = mod$lambdafit#[1:nrow(data)]
lambda = data$lambda[1:length(lambdafit)]


# mod@fullcoef = optimal_coef(data,mod)
# lambdafit = as.vector(exp(mod@covariates %*% mod@fullcoef[2:length(mod@fullcoef)] + mod@fullcoef[1]))[1:nrow(data)]

p = tibble(date = data$date[1:length(lambdafit)],actual = lambda,fitted = lambdafit) %>%
  ggplot(aes(x=date)) +
  geom_line(aes(y=actual,col='Actual')) +
  geom_line(aes(y=fitted,col='Fitted')) +
  ylab('Intensity') +
  xlab('Date') +
  ggtitle(paste("Real -",
                "MAE:",round(mae(lambda,lambdafit),3),
                "MAPE:",round(mape(lambda,lambdafit),2),
                "RMSE:",round(rmse(lambda,lambdafit),3))) +
  theme_bw()

ggplotly(p)
### END