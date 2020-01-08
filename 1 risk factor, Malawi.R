
#clear everything from R
if(!is.bull(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

#load libraries
library(rjags)
load.module("mix") #this is for the normal mixture models, for sensitivity

#set working directory to unzipped folder with data
setwd("~/Desktop/Typhoid Underreporting/Malawi")

set.seed(123456)

#########################
#########MODEL###########
#########################
jcode.m <-"
model{

#BC sensitivity; normal mixture models, created from the data
p_BCpos[1] ~ dnormmix(mu.u5, tau.u5, pi.u5)
p_BCpos[2] ~ dnormmix(mu.5_9, tau.5_9, pi.5_9)
p_BCpos[3] ~ dnormmix(mu.10_14, tau.10_14, pi.10_14)
p_BCpos[4] ~ dnormmix(mu.15_29, tau.15_29, pi.15_29)
p_BCpos[5] ~ dnormmix(mu.30plus, tau.30plus, pi.30plus)
p_BCpos[6] ~ dnormmix(mu.all, tau.all, pi.all)

for (a in 1:6){
#probability of receiving a blood culture; number of positive cultures
p_BC[a] ~ dbeta(alpha.bc[a],beta.bc[a])
n.BCpos[a] ~ dpois(lambda.obs[a])

#obs rate --> true rate (not per 100,000 pyo yet), adjusted for sensitivity, pr(BC)
lambda.obs[a] <- lambda.true[a]*p_BCpos[a]*p_BC[a]
log(lambda.true[a]) <- beta0[a] + log(persontime[a])  #adjust for persontime
true.cases[a] ~ dpois(lambda.true[a])
beta0[a] ~ dnorm(0, 1/100000000)  #weakly informative prior for the intercept
true.rate[a] <- exp(beta0[a]) #true rate (per persontime, but not /100,000 pyo yet)
}

# distribution of relative risk (odds ratio) of risk factor for having typhoid
#soap available vs. not available after defecation
RR_TF <-  exp(logRR_TF)
logRR_TF ~ dnorm(log(.5),1/ 0.343339^2)

for (a in 1:6){
#probabilities used for adjusting for healthcare seeking (see tree diagram)
p_HCfev[a] ~dbeta(alpha.HCfev[a], beta.HCfev[a])
p1[a] ~ dbeta(alpha.p1[a], beta.p1[a])
p0[a] ~ dbeta(alpha.p0[a], beta.p0[a])
S1[a] ~ dbeta(alpha.S1[a], beta.S1[a])
S0[a] ~ dbeta(alpha.S0[a], beta.S0[a])
f1[a] ~ dbeta(alpha.f1[a], beta.f1[a])
f0[a] ~ dbeta(alpha.f0[a], beta.f0[a])

RR_F[a] <- f1[a]/f0[a]

#probability of typhoid without risk factor OR incidence rate of typhoid without risk factor
lambda_0[a] <- true.rate[a]*((1/(p0[a]+RR_F[a]*p1[a]))/p_HCfev[a])*(((S1[a]*RR_F[a]*p1[a])+(S0[a]*p0[a]))/((S1[a]*RR_TF*p1[a])+(S0[a]*p0[a])))

#incidence of typhoid adjusting for sensitivity, pr(BC), pr(healthcare seeking)
final.inc[a] <- (RR_TF*lambda_0[a]*p1[a]+lambda_0[a]*p0[a])*100000

#probability of healthcare seeking
p_HC[a] <- (true.rate[a]*100000)/final.inc[a]
}


}"



#########################
####INPUT PARAMETERS#####
#########################

#load typhoid results and eligibility/enrolled datasets
mal_results <- read.csv("Malawi Results DS 8.22.19.csv")
mal_scrn    <- read.csv("e&e.csv") #malawi screening/ eligible & enrolled

#normal mixture parameters for blood culture sensitivity
mu.u5       <- read.csv("mal.sens.u5.finaltbl.csv")[1,-1]
tau.u5      <- read.csv("mal.sens.u5.finaltbl.csv")[2,-1]
pi.u5       <- read.csv("mal.sens.u5.finaltbl.csv")[3,-1]

mu.5_9      <- read.csv("mal.sens.5_9.finaltbl.csv")[1,-1]
tau.5_9     <- read.csv("mal.sens.5_9.finaltbl.csv")[2,-1]
pi.5_9      <- read.csv("mal.sens.5_9.finaltbl.csv")[3,-1]

mu.10_14    <- read.csv("mal.sens.10_14.finaltbl.csv")[1,-1]
tau.10_14   <- read.csv("mal.sens.10_14.finaltbl.csv")[2,-1]
pi.10_14    <- read.csv("mal.sens.10_14.finaltbl.csv")[3,-1]

mu.15_29    <- read.csv("mal.sens.15_29.finaltbl.csv")[1,-1]
tau.15_29   <- read.csv("mal.sens.15_29.finaltbl.csv")[2,-1]
pi.15_29    <- read.csv("mal.sens.15_29.finaltbl.csv")[3,-1]

mu.30plus   <- read.csv("mal.sens.30plus.finaltbl.csv")[1,-1]
tau.30plus  <- read.csv("mal.sens.30plus.finaltbl.csv")[2,-1]
pi.30plus   <- read.csv("mal.sens.30plus.finaltbl.csv")[3,-1]

mu.all      <- read.csv("mal.sens.all.finaltbl.csv")[1,-1]
tau.all     <- read.csv("mal.sens.all.finaltbl.csv")[2,-1]
pi.all      <- read.csv("mal.sens.all.finaltbl.csv")[3,-1]

#set alpha, beta parameters for beta distributions in model
#pulled from HUS
#Note that HUS only has 3 age groups (<5, 5-14, 15+), so second 2 are used twice
alpha.hc <- c(66, 45, 45, 47, 47, 159)
beta.hc  <- c(41, 40, 40, 44, 44, 124)

alpha.HCfev <- c(65, 65, 65, 59, 59, 189)
beta.HCfev  <- c(44, 17, 17, 32, 32, 93)

alpha.p1 <- c(504, 504, 504, 504, 504, 1512)
beta.p1  <- c(126, 126, 126, 126, 126, 378)

alpha.p0 <- c(126, 126, 126, 126, 126, 378)
beta.p0  <- c(504, 504, 504, 504, 504, 1512)

alpha.S1 <- c(40, 35, 35, 32, 32, 107)
beta.S1  <- c(18, 4,  4,  8,  8,  30)

alpha.S0 <- c(25, 30, 30, 27, 27, 82)
beta.S0  <- c(26, 13, 13, 24, 24, 63)

alpha.f1 <- c(77, 59, 59, 66, 66, 202)
beta.f1  <- c(427, 445, 445, 438, 438, 1310)

alpha.f0 <- c(32, 23, 23, 25, 25, 80)
beta.f0  <- c(94, 103, 103, 101, 101, 298)

#number of blood cultures drawn, by agecat
n.BC <-  c(sum(mal_scrn$enrolled[which(mal_scrn$AgeGp=="0-4")]),
           sum(mal_scrn$enrolled[which(mal_scrn$AgeGp=="5-9yrs")]),
           sum(mal_scrn$enrolled[which(mal_scrn$AgeGp=="10-14yrs")]),
           sum(mal_scrn$enrolled[which(mal_scrn$AgeGp=="15-29")]),
           sum(mal_scrn$enrolled[which(mal_scrn$AgeGp=="30+")]),
           sum(mal_scrn$enrolled[which((mal_scrn$AgeGp=="0-4")|(mal_scrn$AgeGp=="5-9yrs")|(mal_scrn$AgeGp=="10-14yrs")|(mal_scrn$AgeGp=="15-29")|(mal_scrn$AgeGp=="30+"))]))

#number people people who sought healthcare, by agecat
n.HC <- c(sum(mal_scrn$eligible[which(mal_scrn$AgeGp=="0-4")]),
          sum(mal_scrn$eligible[which(mal_scrn$AgeGp=="5-9yrs")]),
          sum(mal_scrn$eligible[which(mal_scrn$AgeGp=="10-14yrs")]),
          sum(mal_scrn$eligible[which(mal_scrn$AgeGp=="15-29")]),
          sum(mal_scrn$eligible[which(mal_scrn$AgeGp=="30+")]),
          sum(mal_scrn$eligible[which((mal_scrn$AgeGp=="0-4")|(mal_scrn$AgeGp=="5-9yrs")|(mal_scrn$AgeGp=="10-14yrs")|(mal_scrn$AgeGp=="15-29")|(mal_scrn$AgeGp=="30+"))]))

#alpha and beta parameters for pr(HC seek) ~Beta
alpha.bc <- n.BC
beta.bc  <- n.HC-n.BC

#persontime
persontime <- c(14417*2,
                12680*2,
                13103*2,
                33179*2,
                26630*2,
                14417*2+12680*2+13103*2+33179*2+26630*2)

#number of blood culture positive results, by agecat
n.BCpos <- c(sum(mal_results$result[which(mal_results$agecat==1)]),
             sum(mal_results$result[which(mal_results$agecat==2)]),
             sum(mal_results$result[which(mal_results$agecat==3)]),
             sum(mal_results$result[which(mal_results$agecat==4)]),
             sum(mal_results$result[which(mal_results$agecat==5)]),
             sum(mal_results$result))

(n.BCpos/persontime)*100000

#########################
#######RUN MODEL#########
#########################

#list of all input variables
jdat.m <- list(n.BCpos=n.BCpos, persontime=persontime, 
               mu.u5=mu.u5,tau.u5=tau.u5,pi.u5=pi.u5,mu.5_9=mu.5_9,tau.5_9=tau.5_9,pi.5_9=pi.5_9,
               mu.10_14=mu.10_14,tau.10_14=tau.10_14,pi.10_14=pi.10_14,mu.15_29=mu.15_29,tau.15_29=tau.15_29,pi.15_29=pi.15_29,
               mu.30plus=mu.30plus,tau.30plus=tau.30plus,pi.30plus=pi.30plus,mu.all=mu.all,tau.all=tau.all,pi.all=pi.all,
               alpha.HCfev=alpha.HCfev, beta.HCfev=beta.HCfev, alpha.bc=alpha.bc, beta.bc=beta.bc, 
               alpha.p1=alpha.p1, alpha.p0=alpha.p0, beta.p1=beta.p1, beta.p0=beta.p0,
               alpha.S0=alpha.S0, alpha.S1=alpha.S1, beta.S0=beta.S0, beta.S1=beta.S1,
               alpha.f0=alpha.f0, alpha.f1=alpha.f1, beta.f0=beta.f0, beta.f1=beta.f1)

#initialize model
jmod.m <- jags.model(textConnection(jcode.m), data=jdat.m, n.chains=3)

#set burn-in period
update(jmod.m,10000)

#fit full model, and decide which parameters to track
#thin 3, run 500000 iterations
jpost.m <- coda.samples(jmod.m, thin=3, c(
  'p_HCfev', 'p0','p1', 'f0', 'f1',
  'S0','S1','RR_F','RR_TF', 
  'lambda_0', 'final.inc','true.rate',
  'p_BC','p_BCpos','p_HC'
), n.iter=500000) 


#########################
######DIAGNOSTICS########
#########################

#check for convergence (quickly); do more thoroughly
#might have to run a few times to get full convergence (mixture models are difficult)
(gelman.diag(jpost.m, multivariate = FALSE))


#########################
#######POSTERIOR#########
#########################

#extract posterior summaries
sum.post.m <- summary(jpost.m)

#view posterior summery of quantiles
View(round(sum.post.m$quantiles,0))

#overall adjustment ratios
round(sum.post.m$quantiles[32:37,3]/((n.BCpos/persontime)*100000),1)

#final adjusted rates
round(sum.post.m$quantiles[32:37,c(3,1,5)],0)

#posterior probabilities
View(round(sum.post.m$quantiles[c(56:73),c(3,1,5)],2))
