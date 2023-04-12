library(tidyverse)
library(ggpubr)

####
# args=commandArgs(TRUE)
# num_weeks = as.numeric(args[1])
# lambda_num = args[2]
# num_effects = as.numeric(args[3])
# num_weeks = 1
# lambda_num = "blocks"
# num_effects = 2

dataset="palo_alto"
path = paste0("5. Outputs/general/",dataset)
train = read.csv(paste0(path,"/metrics_train.csv"))
test = read.csv(paste0(path,"/metrics_test.csv"))

#### Performance

aa = ggplot(data = train, aes(x=model, y = mae)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("MAE")

ab = ggplot(data = train, aes(x=model, y = rmse)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("RMSE")

ac = ggplot(data = train, aes(x=model, y = p_rmse)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Peak RMSE")

ad = ggplot(data = train, aes(x=model, y = deviance)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Deviance")

ba = ggplot(data = train, aes(x=model, y = param)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Number of parameters")

bb = ggplot(data = train, aes(x=model, y = param_non_zeros)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Number of non-zero parameters")

bc = ggplot(data = train, aes(x=model, y = runtime)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Runtime")

bd = ggplot(data = train, aes(x=model, y = dof)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Model d.o.f")

be = ggplot(data = train %>% filter(model %in% c("BACs","BACw","BACsw")), aes(x=model, y = iter)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Number of iterations")

#### Test
ca = ggplot(data = test, aes(x=model, y = mae)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("MAE")

cb = ggplot(data = test, aes(x=model, y = rmse)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("RMSE")

cc = ggplot(data = test, aes(x=model, y = p_rmse)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Peak RMSE")

cd = ggplot(data = test, aes(x=model, y = deviance)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Deviance")

#### Saving figures
pdf(file=paste0(path,"/boxplots_perf_train.pdf"),width=9,height=7)
ggarrange(aa, ab, ac, ad,
          ncol = 2, nrow = 2)
dev.off()

pdf(file=paste0(path,"/boxplots_parc_train.pdf"),width=9,height=7)
ggarrange(ba, bb, bc, bd, be,
          ncol = 2, nrow = 3)
dev.off()

pdf(file=paste0(path,"/boxplots_perf_test.pdf"),width=9,height=7)
ggarrange(ca, cb, cc, cd,
          ncol = 2, nrow = 2)
dev.off()
