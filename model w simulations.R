
#clear everything from R
if(!is.sull(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

#load libraries
library(rjags)
library(tidyverse)
load.module("mix") #this is for the normal mixture models, for sensitivity

set.seed(12345)


# #set working directory to unzipped folder with data
# setwd("~/Desktop/Typhoid Underreporting/Nepal") 

# nep_results <- read.csv(file = "Nepal Results DS 8.22.19.csv")

#model
jcode.s <-"
model{
p_BCpos[1] ~ dnormmix(mu1, tau1, pi1)
#p_BCpos[1] ~ dnorm(mu1, tau1)
p_BCpos[2] ~ dnormmix(mu2, tau2, pi2)
p_BCpos[3] ~ dnormmix(mu3, tau3, pi3)

RR_BCpos <-  exp(logRR_BCpos)
logRR_BCpos ~ dnorm(log(1.87),1/0.3731117^2)

for (a in 1:3){
#probability of receiving a blood culture (without adjustment)
p_BC[a] ~ dbeta(alpha_bc[a],beta_bc[a])

#obs rate --> true rate (not per 100,000 pyo yet), adjusted for sensitivity, pr(BC)
n.BCpos[a] ~ dpois(lambda.obs[a]) # positive BC results ~ poisson
lambda.obs[a] <- lambda.true[a]*p_BCpos[a]*(1-((1-p_BC[a])*(1/RR_BCpos)))
# lambda.obs[a] <- lambda.true[a]*p_BCpos[a]*p_BC[a]
log(lambda.true[a]) <- beta0[a] + log(persontime[a])  #adjust for persontime
true.cases[a] ~ dpois(lambda.true[a])
beta0[a] ~ dnorm(0, 1/100000000)  #weakly informative prior for the intercept
true.rate[a] <- exp(beta0[a])  #true rate (per persontime, but not /100,000 pyo yet)
}


# distribution of relative risk (odds ratio) of risk factor for having typhoid
#unshared vs. shared latrine
RR_TF <-  exp(logRR_TF)
logRR_TF ~ dnorm(log(5.71),1/.4719435)

for (a in 1:3){
#probabilities used for adjusting for healthcare seeking (see tree diagram)
p_HCfev[a] ~dbeta(alpha_HCfev[a], beta_HCfev[a])
p1[a] ~ dbeta(alpha_p1[a], beta_p1[a])
p0[a] ~ dbeta(alpha_p0[a], beta_p0[a])
S1[a] ~ dbeta(alpha_S1[a], beta_S1[a])
S0[a] ~ dbeta(alpha_S0[a], beta_S0[a])
f1[a] ~ dbeta(alpha_f1[a], beta_f1[a])
f0[a] ~ dbeta(alpha_f0[a], beta_f0[a])

RR_F[a] <- f1[a]/f0[a] #relative risk for fever (risk factor)

#probability of typhoid without risk factor OR incidence rate of typhoid without risk factor
lambda_0[a] <- true.rate[a]*((1/(p0[a]+RR_F[a]*p1[a]))/p_HCfev[a])*(((S1[a]*RR_F[a]*p1[a])+(S0[a]*p0[a]))/((S1[a]*RR_TF*p1[a])+(S0[a]*p0[a])))

#incidence of typhoid adjusting for sensitivity, pr(BC), pr(healthcare seeking)
final.inc[a] <- (RR_TF*lambda_0[a]*p1[a]+lambda_0[a]*p0[a])*100000

#probability of healthcare seeking
p_HC[a] <- (true.rate[a]*100000)/final.inc[a]

}

}"


#make data  
true_count <- c(900, 900, 900)
#persontime 
persontime <- rep(100000,3) 

pr_HCfev <- c(.2, .5, .8)
pr_p1    <- c(.2, .5, .8)
pr_p0    <- c(.2, .5, .8)
pr_S1    <- c(.2, .5, .8)
pr_S0    <- c(.2, .5, .8)
pr_f1    <- c(.2, .5, .8)
pr_f0    <- c(.2, .5, .8)
pr_BC    <- c(.2, .5, .8)
# pr_BCpos <- c(.2, .5, .8)

n=735
# n=2000
# n=1000
alpha_p1 <- c(rbinom(1, n, pr_p1[1]),
              rbinom(1, n, pr_p1[2]),
              rbinom(1, n, pr_p1[3]))
beta_p1 <- c(n-alpha_p1[1],
             n-alpha_p1[2],
             n-alpha_p1[3])
alpha_p0 <- beta_p1
beta_p0 <- alpha_p1

alpha_f1 <- c(rbinom(1, alpha_p1[1], pr_f1[1]),
              rbinom(1, alpha_p1[2], pr_f1[2]),
              rbinom(1, alpha_p1[3], pr_f1[3]))
beta_f1 <- c(alpha_p1[1]-alpha_f1[1],
             alpha_p1[2]-alpha_f1[2],
             alpha_p1[3]-alpha_f1[3])
alpha_f0 <- c(rbinom(1, alpha_p0[1], pr_f0[1]),
              rbinom(1, alpha_p0[2], pr_f0[2]),
              rbinom(1, alpha_p0[3], pr_f0[3]))
beta_f0 <- c(alpha_p0[1]-alpha_f0[1],
             alpha_p0[2]-alpha_f0[2],
             alpha_p0[3]-alpha_f0[3])

alpha_S1 <- c(rbinom(1, alpha_p1[1], pr_S1[1]),
              rbinom(1, alpha_p1[2], pr_S1[2]),
              rbinom(1, alpha_p1[3], pr_S1[3]))
beta_S1 <- c(alpha_p1[1]-alpha_S1[1],
             alpha_p1[2]-alpha_S1[2],
             alpha_p1[3]-alpha_S1[3])
alpha_S0 <- c(rbinom(1, alpha_p0[1], pr_S0[1]),
              rbinom(1, alpha_p0[2], pr_S0[2]),
              rbinom(1, alpha_p0[3], pr_S0[3]))
beta_S0 <- c(alpha_p0[1]-alpha_S0[1],
             alpha_p0[2]-alpha_S0[2],
             alpha_p0[3]-alpha_S0[3])
alpha_HCfev <- c(rbinom(1, alpha_f1[1]+alpha_f0[1], pr_HCfev[1]),
              rbinom(1, alpha_f1[2]+alpha_f0[2], pr_HCfev[2]),
              rbinom(1, alpha_f1[3]+alpha_f0[3], pr_HCfev[3]))
beta_HCfev <- c(alpha_f1[1]+alpha_f0[1]-alpha_HCfev[1],
                alpha_f1[2]+alpha_f0[2]-alpha_HCfev[2],
                alpha_f1[3]+alpha_f0[3]-alpha_HCfev[3])


true_inc <- (true_count/persontime)
true_incp100000 <- true_inc*100000

# typh_norisk <- c(rbinom(1, true_count[1], pr_p0[1]),
#                  rbinom(1, true_count[2], pr_p0[2]),
#                  rbinom(1, true_count[3], pr_p0[3]),
#                  rbinom(1, true_count[4], pr_p0[4]),
#                  rbinom(1, true_count[5], pr_p0[5]))
# 
# typh_inc_norisk <- typh_norisk/persontime
# true_lambda0 <- typh_inc_norisk

# logRR_TF <- rnorm(1,log(5.71),.4719435)
# RR_TF <-  exp(logRR_TF)
RR_TF <-  5.71
# obs_p1 <- c(rbeta(1,alpha_p1[1], beta_p1[1]),
#             rbeta(1,alpha_p1[2], beta_p1[2]),
#             rbeta(1,alpha_p1[3], beta_p1[3]),
#             rbeta(1,alpha_p1[4], beta_p1[4]),
#             rbeta(1,alpha_p1[5], beta_p1[5]))
# obs_p0 <- c(rbeta(1,alpha_p0[1], beta_p0[1]),
#             rbeta(1,alpha_p0[2], beta_p0[2]),
#             rbeta(1,alpha_p0[3], beta_p0[3]),
#             rbeta(1,alpha_p0[4], beta_p0[4]),
#             rbeta(1,alpha_p0[5], beta_p0[5]))

RR_F <- pr_f1/pr_f0

true_lambda0 <- true_inc/(RR_TF*pr_p1 + pr_p0)

true.rate <- true_lambda0/(((1/(pr_p0+RR_F*pr_p1))/pr_HCfev)*((pr_S1*RR_F*pr_p1+pr_S0*pr_p0)/(pr_S1*RR_TF*pr_p1+pr_S0*pr_p0)))


pr_HC <- round(true.rate/true_inc,2)

lambda_true <- true.rate
lambda_true*100000
lambda_true*persontime

# logRR_BCpos <- rnorm(3,log(1.87),0.3731117)
# RR_BCpos <-  exp(logRR_BCpos)
RR_BCpos <- 1.87

alpha_bc2 <- c(rbinom(1, round((lambda_true*persontime),0)[1],  ifelse(pr_BC[1]*RR_BCpos>1,1,pr_BC[1]*RR_BCpos)),
              rbinom(1, round((lambda_true*persontime),0)[2],  ifelse(pr_BC[2]*RR_BCpos>1,1,pr_BC[2]*RR_BCpos)),
               rbinom(1, round((lambda_true*persontime),0)[3],  ifelse(pr_BC[3]*RR_BCpos>1,1,pr_BC[3]*RR_BCpos)))
# alpha_bc <- c(rbinom(1, round((lambda_true*persontime),0)[1],  ifelse(pr_BC[1]*2>1,1,pr_BC[1]*2)),
#               rbinom(1, round((lambda_true*persontime),0)[2],  ifelse(pr_BC[2]*2>1,1,pr_BC[2]*2)),
#               rbinom(1, round((lambda_true*persontime),0)[3],  ifelse(pr_BC[3]*2>1,1,pr_BC[3]*2)))
alpha_bc <- c(rbinom(1, round((lambda_true*persontime),0)[1],  pr_BC[1]),
              rbinom(1, round((lambda_true*persontime),0)[2], pr_BC[2]),
              rbinom(1, round((lambda_true*persontime),0)[3],  pr_BC[3]))

alpha_bc
alpha_bc2

beta_bc <- c(round((lambda_true*persontime),0)[1]-alpha_bc[1],
             round((lambda_true*persontime),0)[2]-alpha_bc[2],
             round((lambda_true*persontime),0)[3]-alpha_bc[3])
beta_bc2 <- c(round((lambda_true*persontime),0)[1]-alpha_bc2[1],
             round((lambda_true*persontime),0)[2]-alpha_bc2[2],
             round((lambda_true*persontime),0)[3]-alpha_bc2[3]+.5)




bc_vol1 <- rnorm(n = alpha_bc2[1], mean = 4, sd = .1)
bc_vol2 <- rnorm(n = alpha_bc2[2], mean = 4, sd = .1)
bc_vol3 <- rnorm(n = alpha_bc2[3], mean = 4, sd = .1)

p_abx <- c(.1, .2, .3)

abx_1 <- rbinom(n= alpha_bc2[1], size=1, prob = p_abx[1])
abx_2 <- rbinom(n= alpha_bc2[2], size=1, prob = p_abx[2])
abx_3 <- rbinom(n= alpha_bc2[3], size=1, prob = p_abx[3])

sens1 <- (exp(-0.732088+0.030566*bc_vol1))*(1-.34*abx_1)
sens2 <- (exp(-0.732088+0.030566*bc_vol2))*(1-.34*abx_2)
sens3 <- (exp(-0.732088+0.030566*bc_vol3))*(1-.34*abx_3)

bc_yes <- alpha_bc2
bc_no <- beta_bc2



mean(sens1[which(abx_1==1)])
sd(sens1[which(abx_1==1)])
length(sens1[which(abx_1==1)])/length(sens1)
mean(sens1[which(abx_1==0)])
sd(sens1[which(abx_1==0)])
length(sens1[which(abx_1==0)])/length(sens1)
sens.tbl1 <- as.data.frame(matrix(c(mean(sens1[which(abx_1==0)]), mean(sens1[which(abx_1==1)]),
                                          1/(sd(sens1[which(abx_1==0)]))^2, 1/(sd(sens1[which(abx_1==1)]))^2,
                                          length(sens1[which(abx_1==0)])/length(sens1),length(sens1[which(abx_1==1)])/length(sens1)),nrow=3,ncol=2,byrow = T))
sens.tbl1

mean(sens2[which(abx_2==1)])
sd(sens2[which(abx_2==1)])
length(sens2[which(abx_2==1)])/length(sens2)
mean(sens2[which(abx_2==0)])
sd(sens2[which(abx_2==0)])
length(sens2[which(abx_2==0)])/length(sens2)
sens.tbl2 <- as.data.frame(matrix(c(mean(sens2[which(abx_2==0)]), mean(sens2[which(abx_2==1)]),
                                    1/(sd(sens2[which(abx_2==0)]))^2,1/(sd(sens2[which(abx_2==1)]))^2,
                                    length(sens2[which(abx_2==0)])/length(sens2),length(sens2[which(abx_2==1)])/length(sens2)),nrow=3,ncol=2,byrow = T))

mean(sens3[which(abx_3==1)])
sd(sens3[which(abx_3==1)])
length(sens3[which(abx_3==1)])/length(sens3)
mean(sens3[which(abx_3==0)])
sd(sens3[which(abx_3==0)])
length(sens3[which(abx_3==0)])/length(sens3)
sens.tbl3 <- as.data.frame(matrix(c(mean(sens3[which(abx_3==0)]), mean(sens3[which(abx_3==1)]),
                                    1/(sd(sens3[which(abx_3==0)]))^2,1/(sd(sens3[which(abx_3==1)]))^2,
                                    length(sens3[which(abx_3==0)])/length(sens3),length(sens3[which(abx_3==1)])/length(sens3)),nrow=3,ncol=2,byrow = T))

mu1  <- sens.tbl1[1,]
tau1 <- sens.tbl1[2,]
pi1  <- sens.tbl1[3,]

# mu1  <- sens.tbl1[1,1]
# tau1 <- sens.tbl1[2,1]

mu2  <- sens.tbl2[1,]
tau2 <- sens.tbl2[2,]
pi2  <- sens.tbl2[3,]

mu3  <- sens.tbl3[1,]
tau3 <- sens.tbl3[2,]
pi3  <- sens.tbl3[3,]

BCpos1 <- rep(NA, length(bc_vol1))
for (i in 1:length(BCpos1)){
  BCpos1[i] <- as.numeric(rbernoulli(1,sens1[i]))
}

BCpos2 <- rep(NA, length(bc_vol2))
for (i in 1:length(BCpos2)){
  BCpos2[i] <- as.numeric(rbernoulli(1,sens2[i]))
}

BCpos3 <- rep(NA, length(bc_vol3))
for (i in 1:length(BCpos3)){
  BCpos3[i] <- as.numeric(rbernoulli(1,sens3[i]))
}

n.BCpos <- c(sum(BCpos1),sum(BCpos2),sum(BCpos3))

pr_BCpos <- c(median(sens1),median(sens2),median(sens3))

# pr_BC2 <- c((pr_BC[1]+(1-pr_BC[1])*(1/RR_BCpos[1])),
#            (pr_BC[2]+(1-pr_BC[2])*(1/RR_BCpos[2])),
#            (pr_BC[3]+(1-pr_BC[3])*(1/RR_BCpos[3])))

# lambda_obs <- lambda_true*pr_BCpos*pr_BC2
# lambda_obs*100000
# lambda_obs*persontime

#check crude incidence
(n.BCpos/persontime)*100000

#list of all input variables
jdat.s <- list(n.BCpos=n.BCpos, persontime=persontime, mu1=mu1, mu2=mu2, mu3=mu3, 
               tau1=tau1, tau2=tau2, tau3=tau3, pi2=pi2, pi3=pi3, pi1=pi1,
               alpha_bc=alpha_bc2, beta_bc=beta_bc2, alpha_HCfev=alpha_HCfev, beta_HCfev=beta_HCfev,
               alpha_p1=alpha_p1, alpha_p0=alpha_p0, beta_p1=beta_p1, beta_p0=beta_p0,
               alpha_S0=alpha_S0, alpha_S1=alpha_S1, beta_S0=beta_S0, beta_S1=beta_S1,
               alpha_f0=alpha_f0, alpha_f1=alpha_f1, beta_f0=beta_f0, beta_f1=beta_f1)

#initialize model
jmod.s <- jags.model(textConnection(jcode.s), data=jdat.s, n.chains=3)

#set burn-in period
update(jmod.s,10000)


#fit full model, and decide which parameters to track
#thin 3, run 100000 iterations
jpost.s <- coda.samples(jmod.s, thin=3, c('p_BC','p_BCpos','p_HC','p_HCfev', 
                                          'p0','p1', 'f0', 'f1',
                                          'S0','S1','RR_F','RR_TF', 
                                          'lambda_0', 'final.inc', 'true.rate'
                                          ), n.iter=100000) 

#check for convergence (quickly); do more thoroughly
#might have to run a few times to get full convergence (mixture models are difficult)
(gelman.diag(jpost.s, multivariate = FALSE))

#extract posterior summaries
sum.post.s <- summary(jpost.s)

#view posterior summery of quantiles
View(round(sum.post.s$quantiles,0))

#overall adjustment ratios
#round(sum.post.s$quantiles[27:31,3]/((n.BCpos/persontime)*100000),1)

#final adjusted rates
round(sum.post.s$quantiles[17:19,c(3,1,5)],0)
# final.inc[1] 493  283  1092

# final.inc[2] 900  748  1089
# final.inc[3] 871  777   973
true_incp100000

View(cbind(true_incp100000, round(sum.post.s$quantiles[17:19,c(3,1,5)],0)))

sim_est <- as.data.frame(cbind(round(true_incp100000), round(sum.post.s$quantiles[17:19,c(3,1,5)],0)))

#posterior probabilities
View(round(sum.post.s$quantiles[c(29:37),c(3,1,5)],2))

cbind(round(c(ifelse(pr_BC*RR_BCpos>1,1,pr_BC*RR_BCpos),pr_BCpos,pr_HC),2),round(sum.post.s$quantiles[c(29:37),c(3,1,5)],2))

View(cbind(round(c(ifelse(pr_BC*RR_BCpos>1,1,pr_BC*RR_BCpos),pr_BCpos,pr_HC),2),round(sum.post.s$quantiles[c(29:37),c(3,1,5)],2)))
# 50% 2.5% 97.5%
# 
# p_BC[2]    0.94 0.96 0.94  0.97
# p_BC[3]    1.00 1.00 1.00  1.00
# 
# p_BCpos[2] 0.54 0.54 0.54  0.55
# p_BCpos[3] 0.54 0.54 0.54  0.55
# 
# p_HC[2]    0.50 0.45 0.39  0.51
# p_HC[3]    0.80 0.78 0.74  0.81

# p_BC[1]    0.37 0.39 0.32  0.47
# 
# p_BCpos[1] 0.54 0.54 0.54  0.55
# 
# p_HC[1]    0.20 0.21 0.15  0.30

sim_probs <- as.data.frame(cbind(round(c(ifelse(pr_BC*RR_BCpos>1,1,pr_BC*RR_BCpos),pr_BCpos,pr_HC),2),round(sum.post.s$quantiles[c(29:37),c(3,1,5)],2)))

# write.csv(sim_probs,"~/Desktop/Typhoid Underreporting/sim_probs_NOmsng735.csv")
# write.csv(sim_est,"~/Desktop/Typhoid Underreporting/sim_est_NOmsng735.csv")

# write.csv(sim_probs,"~/Desktop/Typhoid Underreporting/sim_probs_NOmsng1000.csv")
# write.csv(sim_est,"~/Desktop/Typhoid Underreporting/sim_est_NOmsng1000.csv")

write.csv(sim_probs,"~/Desktop/Typhoid Underreporting/sim_probs_NOmsng2000.csv")
write.csv(sim_est,"~/Desktop/Typhoid Underreporting/sim_est_NOmsng2000.csv")

# write.csv(sim_probs,"~/Desktop/Typhoid Underreporting/sim_probs_msng2000.csv")
# write.csv(sim_est,"~/Desktop/Typhoid Underreporting/sim_est_msng2000.csv")

# write.csv(sim_probs,"~/Desktop/Typhoid Underreporting/sim_probs_msng1000.csv")
# write.csv(sim_est,"~/Desktop/Typhoid Underreporting/sim_est_msng1000.csv")

# write.csv(sim_probs,"~/Desktop/Typhoid Underreporting/sim_probs_msng735.csv")
# write.csv(sim_est,"~/Desktop/Typhoid Underreporting/sim_est_msng735.csv")
