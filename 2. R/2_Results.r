library(tidyverse)
library(ggpubr)

####
args=commandArgs(TRUE)
num_weeks = as.numeric(args[1])
lambda_num = args[2]
num_effects = as.numeric(args[3])
# num_weeks = 1
# lambda_num = "blocks"
# num_effects = 2

path = paste0("5. Outputs/",num_weeks,"_",lambda_num,"_",num_effects,"/")
temp = list.files(path=path,pattern="*.csv")
df = bind_rows(lapply(paste0(path,temp), read.csv))

#### Performance

aa = ggplot(data = df, aes(x=model, y = mae)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("MAE")

ab = ggplot(data = df, aes(x=model, y = rmse)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("RMSE")

ac = ggplot(data = df, aes(x=model, y = p_rmse)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Peak RMSE")

ad = ggplot(data = df, aes(x=model, y = deviance)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Deviance")

ba = ggplot(data = df, aes(x=model, y = param)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Number of parameters")

bb = ggplot(data = df, aes(x=model, y = param_non_zeros)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Number of non-zero parameters")

bc = ggplot(data = df, aes(x=model, y = runtime)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Runtime")

bd = ggplot(data = df, aes(x=model, y = dof)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Model d.o.f")

be = ggplot(data = df %>% filter(model %in% c("BACs","BACw","BACsw")), aes(x=model, y = iter)) + 
  geom_boxplot(aes(group = model))+
  ylab('') +
  xlab("Model") +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Number of iterations")

pdf(file=paste0(path,"boxplots_perf",num_weeks,"_",lambda_num,"_",num_effects,".pdf"),width=9,height=7)
ggarrange(aa, ab, ac, ad,
          ncol = 2, nrow = 2)
dev.off()

pdf(file=paste0(path,"boxplots_parc",num_weeks,"_",lambda_num,"_",num_effects,".pdf"),width=9,height=7)
ggarrange(ba, bb, bc, bd, be,
          ncol = 2, nrow = 3)
dev.off()