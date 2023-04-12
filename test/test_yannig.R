summary(train)

effects = list()
effects[[1]] = as.matrix(train %>% dplyr::select(starts_with("indicator")))
effects[[2]] = as.matrix(train %>% dplyr::select(starts_with("wav_tod")))
effects[[3]] = as.matrix(train %>% dplyr::select(starts_with("wav_temp")))
effects[[4]] = as.matrix(train %>% dplyr::select(starts_with("wav_traffic")))


X=as.data.frame(effects)
ncol(X)
qr(X)$rank
ncol(X) == qr(X)$rank

#### 1. BACKFITTING

plot(data$lambda, type='l')


out_BAC_GLM = BAC_GLM(effects[[2]],train)



BAC_GLM_time  = difftime(Sys.time(),start_time,units = "secs")


#########@
p <- read.csv("https://stats.idre.ucla.edu/stat/data/poisson_sim.csv")
p <- within(p, {
  prog <- factor(prog, levels=1:3, labels=c("General", "Academic", 
                                            "Vocational"))
  id <- factor(id)
})
summary(p)
p$event <- p$num_awards
p$math = (p$math - mean(p$math))/sd(p$math)

m1 <- glm(num_awards ~ prog + math, family="poisson", data=p)
m1

###1er test
effects <- list()
prog <- model.matrix(~prog, p)
effects[[1]] = as.matrix(prog)
effects[[2]] = as.matrix(p %>% dplyr::select(math))
out_BAC_GLM = BAC_GLM(effects,p)


m1$coefficients
out_BAC_GLM$fullcoef  ###résultats proches


###2e test
###on a les NA à l'tération 13
effects <- list()
prog <- model.matrix(~prog, p)
effects[[1]] = as.matrix(prog[,1])
effects[[2]] = as.matrix(prog[,2])
effects[[3]] = as.matrix(prog[,3])
effects[[4]] = as.matrix(p %>% dplyr::select(math))
out_BAC_GLM = BAC_GLM(effects,p, maxiter = 100)

##en initialisant correctement

coeff.init[[1]] <- m1$coefficients[1]
coeff.init[[2]] <- m1$coefficients[2]
coeff.init[[3]] <- m1$coefficients[3]
coeff.init[[4]] <- m1$coefficients[4]
out_BAC_GLM = BAC_GLM(effects,p, maxiter = 13, coeff.init=coeff.init)
out_BAC_GLM$fullcoef

m1$coefficients
out_BAC_GLM$new_coeffs_stock


####ça ne converge pas:
coeff.init <- list()
for(i in c(1:4))
{
  delta <- abs(m1$coefficients[i])*0.1
  coeff.init[[i]] <- runif(1, min=m1$coefficients[i]-delta, max=m1$coefficients[i]+delta)
}

out_BAC_GLM   <- BAC_GLM(effects,p, coeff.init=coeff.init)
m1$coefficients


####en changeant le coeff 1 (intercept) de 1% seulement  ça marche (mais pas tout le temps)

delta <- abs(m1$coefficients[i])*0.01
coeff.init[[1]] <- runif(1, min=m1$coefficients[1]-delta, max=m1$coefficients[1]+delta)
coeff.init[[2]] <- m1$coefficients[2]
coeff.init[[3]] <- m1$coefficients[3]
coeff.init[[4]] <- m1$coefficients[4]

out_BAC_GLM   <- BAC_GLM(effects,p, coeff.init=coeff.init)

out_BAC_GLM$fullcoef
m1$coefficients

####en changeant le coeff 2 (intercept) de 5% seulement  ça marche (mais pas tout le temps)
coeff.init[[1]] <- m1$coefficients[1]
delta <- abs(m1$coefficients[i])*0.05
coeff.init[[2]] <- runif(1, min=m1$coefficients[2]-delta, max=m1$coefficients[2]+delta)
coeff.init[[3]] <- m1$coefficients[3]
coeff.init[[4]] <- m1$coefficients[4]

out_BAC_GLM   <- BAC_GLM(effects,p, coeff.init=coeff.init, family="quasipoisson")

out_BAC_GLM$fullcoef
m1$coefficients


###exemple de ce type d'erreur:
#https://www.statology.org/glm-fit-fitted-probabilities-numerically-0-or-1-occurred/#:~:text=This%20warning%20occurs%20when%20you,message%20and%20not%20an%20error.




