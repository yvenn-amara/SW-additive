##### Simulating Arrivals based on a known intensity function #####

##### 1. Importing functions libraries
library(tidyverse)
library(lubridate)
library(wavethresh)
library(gamwave)
library(glmnet)
library(ggpubr)
library(gridExtra)
library(mgcv)
require(fastDummies)
args=commandArgs(TRUE)

# step = as.numeric(args[1])
# num_weeks = as.numeric(args[2])
# lambda_num = as.numeric(args[3])
  # iter = as.numeric(args[4])
step=15
# num_weeks=1
num_weeks = as.numeric(args[1])
# lambda_num = "blocks"
lambda_num = args[2]
# iter = 1
iter = as.numeric(args[3])
tod_res = 6
points_simu = (7*1440/step)*num_weeks
# num_effects = 1
num_effects = as.numeric(args[4])

dataset = "simulated_arrivals"
start_date = "2022-02-21"
source("2. R/0_Algos_Estim_NHPP_v4.R")
path = paste0("5. Outputs/",num_weeks,"_",lambda_num,"_",num_effects,"/")
dir.create(path)
source("2. R/prep.R")
##### plots

pdf(file=paste0(path,"intensity",num_weeks,"_",lambda_num,"_",num_effects,".pdf"),width=9,height=5)
ggplot(data = data[1:480,], aes(x=date, y = lambda)) +
  geom_line()+
  ylab('Lambda') +
  xlab("Date") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("")
dev.off()

pdf(file=paste0(path,"simulation",num_weeks,"_",lambda_num,"_",num_effects,".pdf"),width=9,height=5)
ggplot(data = simu[1:480,], aes(x=date, y = event)) +
  geom_line()+
  ylab('Arrivals') +
  xlab("Date") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("")
dev.off()

######

# fmla = as.formula("event ~ factor(day) + s(tod, k=30) + s(traffic, k=30) + s(temp, k=30) + w(tod) + w(traffic) + w(temp)")
# fmla1 = as.formula("event ~ factor(day) + s(tod, k=10) + s(traffic, k=10) + s(temp, k=10)")
# fmla = as.formula("event ~ factor(day) + s(tod, k=10) + s(traffic, k=10) + s(temp, k=10) + wp(tod) + wp(traffic) + wp(temp)")
# fmla = as.formula("event ~ factor(day) + wp(tod) + wp(traffic) + wp(temp)")
# fmla1 = as.formula("event ~ factor(day) + s(tod, k=30) + s(traffic, k=30) + s(temp, k=30)")
# fmla = as.formula("event ~ factor(day) + s(tod, k=30) + s(traffic, k=30) + s(temp, k=30)")

# fmla1 = as.formula("event ~ factor(day) + s(tod, k=30)")
# fmla = as.formula("event ~ factor(day) + s(tod, k=30) + wp(tod)")

# fmla = as.formula("event ~ factor(day) + wp(traffic) + wp(temp)")
# fmla = as.formula("event ~ factor(day) + s(traffic, k=10) + s(temp, k=10) + wp(traffic) + wp(temp)")
# fmla1 = as.formula("event ~ factor(day) + s(traffic, k=10) + s(temp, k=10)")


fmla_s = "event ~ factor(day)"
for(i in 1:num_effects){fmla_s = paste0(fmla_s,'+s(effects_',i,',k=50)')}                  
fmla_s = as.formula(fmla_s)

fmla_w = "event ~ factor(day)"
for(i in 1:num_effects){fmla_w = paste0(fmla_w,'+wp(effects_',i,')')}                  
fmla_w = as.formula(fmla_w)

fmla_sw = "event ~ factor(day)"
for(i in 1:num_effects){fmla_sw = paste0(fmla_sw,'+s(effects_',i,',k=30)','+wp(effects_',i,')')}                  
fmla_sw = as.formula(fmla_sw)

BACs = rig_BAC_final(fmla = fmla_s, train = train, eps=1e-3, maxiter=100,get_full_rank=TRUE)
BACw = rig_BAC_final(fmla = fmla_w, train = train, eps=1e-3, maxiter=100,get_full_rank=TRUE)
BACsw = rig_BAC_final(fmla = fmla_sw, train = train, eps=1e-3, maxiter=100,get_full_rank=TRUE)


# gam.check(out$models[[3]])
# plot(out$models[[5]])
# out$fullcoef

pdf(file=paste0(path,"plotwav_last_effect_beforefit",num_weeks,"_",lambda_num,"_",num_effects,"_",iter,".pdf"),width=9,height=5)
plot_wav(train[paste0('effects_',num_effects)],
         numLevels = log2(length(BACsw$fullcoef[[length(BACsw$fullcoef)]])+1),
         filterNumber = 1)
dev.off()

pdf(file=paste0(path,"plotwav_last_effect_afterfit",num_weeks,"_",lambda_num,"_",num_effects,"_",iter,".pdf"),width=9,height=5)
plot_wav(train[paste0('effects_',num_effects)],
         numLevels = log2(length(BACsw$fullcoef[[length(BACsw$fullcoef)]])+1),
         filterNumber = 1,
         coef_is_null = BACsw$fullcoef[[length(BACsw$fullcoef)]])
dev.off()

# plot(train$lambda, type='l')
# plot(train$event, type='l')
# lines(out$lambdafit, col='blue')
# mae(train$lambda,out$lambdafit)
# rmse(train$lambda,out$lambdafit)
# 
# mae(train$lambda,out1$lambdafit)
# rmse(train$lambda,out1$lambdafit)
# 
# mae(train$lambda,out2$lambdafit)
# rmse(train$lambda,out2$lambdafit)

#### OBO

OBOs = OBO_pen(fmla = fmla_s, train = train, eps=1e-3, maxiter=100,get_full_rank=TRUE)
OBOw = OBO_pen(fmla = fmla_w, train = train, eps=1e-3, maxiter=100,get_full_rank=TRUE)
OBOsw = OBO_pen(fmla = fmla_sw, train = train, eps=1e-3, maxiter=100,get_full_rank=TRUE)

# out = OBO_pen(fmla = fmla, train = train, eps=1e-8, maxiter=1000)
# plot(train$lambda, type='l')
# lines(out$lambdafit, col='blue')
mae(OBOsw$final)
rmse(OBOsw$final)
p_rmse(OBOsw$final)

#### Saving results

metrics = tibble(model=c("OBOs","OBOw","OBOsw","BACs","BACw","BACsw"),
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
                 iter = c(NA,
                          NA,
                          NA,
                          BACs$iterations,
                          BACw$iterations,
                          BACsw$iterations)
                 
) %>%
  mutate(runtime = as.numeric(runtime))

write.csv(metrics,paste0(path,"metrics",num_weeks,"_",lambda_num,"_",num_effects,"_",iter,".csv"), row.names = F)

# forplot = data %>%
#   dplyr::select(date,lambda) %>%
#   mutate(model="observed") %>%
#   bind_rows(.,AAO_tibble,AAO_pen_tibble,BAC_GLM_tibble,BAC_GLMnet_tibble,GAM_tibble,GLOB_tibble) %>%
#   filter(date <= min(date)+1440*5*60)
# 
# if (iter ==1){
#   pdf(file=paste0("5. Outputs/",num_weeks,"_",lambda_num,"/lineplot_step",step,"_nweeks",num_weeks,"_nlambda",lambda_num,"_iter",iter,".pdf"),width=9,height=5)
#   ggplot(data = forplot, aes(x=date, y = lambda)) + 
#     geom_line(aes(color = model)) +
#     theme_bw()
#   dev.off()
# }


##### END