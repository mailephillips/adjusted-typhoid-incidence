
#clear everything from R
if(!is.bull(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

#load libraries
library(rjags)
load.module("mix") #this is for the normal mixture models, for sensitivity

set.seed(1)

#set working directory to unzipped folder with data
setwd("~/Desktop/Typhoid Underreporting/Bangladesh")


#########################
#########MODEL###########
#########################
jcode.b <-"
model{

#BC sensitivity; normal mixture models, created from the data
p_BCpos[1] ~ dnormmix(mu.u5, tau.u5, pi.u5)
p_BCpos[2] ~ dnormmix(mu.5_9, tau.5_9, pi.5_9)
p_BCpos[3] ~ dnormmix(mu.10_14, tau.10_14, pi.10_14)
p_BCpos[4] ~ dnormmix(mu.15_29, tau.15_29, pi.15_29)
p_BCpos[5] ~ dnormmix(mu.30plus, tau.30plus, pi.30plus)
p_BCpos[6] ~ dnormmix(mu.all, tau.all, pi.all)

#relative risk for positive BC (whether or not received BC)
RR_BCpos <-  exp(logRR_BCpos)
logRR_BCpos ~ dnorm(log(1.87),1/0.3731117^2)

for (a in 1:6){
#probability of receiving a blood culture (without adjustment)
p_BC[a] ~ dbeta(alpha.bc[a],beta.bc[a])

#obs rate --> true rate (not per 100,000 pyo yet), adjusted for sensitivity, pr(BC)
n.BCpos[a] ~ dpois(lambda.obs[a]) # positive BC results ~ poisson
lambda.obs[a] <- lambda.true[a]*p_BCpos[a]*(p_BC[a]+(1-p_BC[a])*(1/RR_BCpos)) 
log(lambda.true[a]) <- beta0[a] + log(persontime[a])  #adjust for persontime
true.cases[a] ~ dpois(lambda.true[a])
beta0[a] ~ dnorm(0, 1/100000000)  #weakly informative prior for the intercept
true.rate[a] <- exp(beta0[a])  #true rate (per persontime, but not /100,000 pyo yet)
}

# distribution of relative risk (odds ratio) of risk factor for having typhoid
#boiled vs. unboiled drinking water
RR_TF <-  exp(logRR_TF)
logRR_TF ~ dnorm(log(7.6),1/0.6372431^2)

for (a in 1:6){
#probabilities used for adjusting for healthcare seeking (see tree diagram)
p_HCfev[a] ~dbeta(alpha.HCfev[a], beta.HCfev[a])
p1[a] ~ dbeta(alpha.p1[a], beta.p1[a])
p0[a] ~ dbeta(alpha.p0[a], beta.p0[a])
S1[a] ~ dbeta(alpha.S1[a], beta.S1[a])
S0[a] ~ dbeta(alpha.S0[a], beta.S0[a])
f1[a] ~ dbeta(alpha.f1[a], beta.f1[a])
f0[a] ~ dbeta(alpha.f0[a], beta.f0[a])

RR_F[a] <- f1[a]/f0[a] #relative risk for fever (risk factor)

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

#load typhoid results dataset
bang_results <- read.csv(file = "Bangladesh Results DS 8.22.19.csv")

#normal mixture model parameters
mu.u5  <- read.csv("bang.sens.u5.finaltbl.csv")[1,-1]
tau.u5 <- read.csv("bang.sens.u5.finaltbl.csv")[2,-1]
pi.u5  <- read.csv("bang.sens.u5.finaltbl.csv")[3,-1]

mu.5_9  <- read.csv("bang.sens.5_9.finaltbl.csv")[1,-1]
tau.5_9 <- read.csv("bang.sens.5_9.finaltbl.csv")[2,-1]
pi.5_9  <- read.csv("bang.sens.5_9.finaltbl.csv")[3,-1]

mu.10_14  <- read.csv("bang.sens.10_14.finaltbl.csv")[1,-1]
tau.10_14 <- read.csv("bang.sens.10_14.finaltbl.csv")[2,-1]
pi.10_14  <- read.csv("bang.sens.10_14.finaltbl.csv")[3,-1]

mu.15_29  <- read.csv("bang.sens.15_29.finaltbl.csv")[1,-1]
tau.15_29 <- read.csv("bang.sens.15_29.finaltbl.csv")[2,-1]
pi.15_29  <- read.csv("bang.sens.15_29.finaltbl.csv")[3,-1]

mu.30plus  <- read.csv("bang.sens.30plus.finaltbl.csv")[1,-1]
tau.30plus <- read.csv("bang.sens.30plus.finaltbl.csv")[2,-1]
pi.30plus  <- read.csv("bang.sens.30plus.finaltbl.csv")[3,-1]

mu.all  <- read.csv("bang.sens.all.finaltbl.csv")[1,-1]
tau.all <- read.csv("bang.sens.all.finaltbl.csv")[2,-1]
pi.all  <- read.csv("bang.sens.all.finaltbl.csv")[3,-1]

#set alpha, beta parameters for beta distributions in model
#pulled from HUS
#Note that HUS only has 3 age groups (<5, 5-14, 15+), so second 2 are used twice
alpha.bc <- c(4884, 3590, 3590, 4677, 4677, 13151)
beta.bc  <- c(673,  214,  214,  161,  161,  1048)

alpha.HCfev <- c(41, 50,  50,  39,  39,  130)
beta.HCfev  <- c(68, 107, 107, 138, 138, 313)

alpha.p1 <- c(146, 306, 306, 411, 411, 863)
beta.p1  <- c(356, 731, 731, 965, 965, 2052)

alpha.p0 <- c(356,  731, 731, 965, 965, 2052)
beta.p0  <- c(146, 306, 306, 411, 411, 863)

alpha.S1 <- c(10, 16, 16, 9,  9,  35)
beta.S1  <- c(24, 30, 30, 44, 44, 98)

alpha.S0 <- c(31, 34, 34, 30, 30, 95)
beta.S0  <- c(44, 77, 77, 94, 94, 215)

alpha.f1 <- c(34,  46,  46,  53,  53,  133)
beta.f1  <- c(112, 260, 260, 358, 358, 730)

alpha.f0 <- c(75,  111, 111, 124, 124, 310)
beta.f0  <- c(281, 620, 620, 841, 841, 1742)

#persontime
persontime <- c(11262*2,
                10205*2,
                10618*2,
                35773*2,
                43460*2,
                11262*2+10205*2+10618*2+35773*2+43460*2) 

#divide results dataset into age categories
bangu5 <- bang_results[which(bang_results$agecat==1),]
bang5_9 <- bang_results[which(bang_results$agecat==2),]
bang10_14 <- bang_results[which(bang_results$agecat==3),]
bang15_29 <- bang_results[which(bang_results$agecat==4),]
bang30plus <- bang_results[which(bang_results$agecat==5),]
bangall <- bang_results

#number of positive BC results, by agecat
n.BCpos <- c(sum(bang_results$result[which(bang_results$agecat==1)]),
             sum(bang_results$result[which(bang_results$agecat==2)]),
             sum(bang_results$result[which(bang_results$agecat==3)]),
             sum(bang_results$result[which(bang_results$agecat==4)]),
             sum(bang_results$result[which(bang_results$agecat==5)]),
             sum(bang_results$result))

#check crude incidence
(n.BCpos/persontime)*100000


#########################
#######RUN MODEL#########
#########################

#list of all input variables
jdat.b <- list(n.BCpos=n.BCpos, persontime=persontime, 
               mu.u5=mu.u5,tau.u5=tau.u5,pi.u5=pi.u5,mu.5_9=mu.5_9,tau.5_9=tau.5_9,pi.5_9=pi.5_9,
               mu.10_14=mu.10_14,tau.10_14=tau.10_14,pi.10_14=pi.10_14,mu.15_29=mu.15_29,tau.15_29=tau.15_29,pi.15_29=pi.15_29,
               mu.30plus=mu.30plus,tau.30plus=tau.30plus,pi.30plus=pi.30plus,mu.all=mu.all,tau.all=tau.all,pi.all=pi.all,
               alpha.bc=alpha.bc, beta.bc=beta.bc, alpha.HCfev=alpha.HCfev, beta.HCfev=beta.HCfev,
             alpha.p1=alpha.p1, alpha.p0=alpha.p0, beta.p1=beta.p1, beta.p0=beta.p0,
             alpha.S0=alpha.S0, alpha.S1=alpha.S1, beta.S0=beta.S0, beta.S1=beta.S1,
             alpha.f0=alpha.f0, alpha.f1=alpha.f1, beta.f0=beta.f0, beta.f1=beta.f1)

#initialize model
jmod.b <- jags.model(textConnection(jcode.b), data=jdat.b, n.chains=3)

#set burn-in period
update(jmod.b,10000)

#fit full model, and decide which parameters to track
#thin 3, run 100000 iterations
jpost.b <- coda.samples(jmod.b, thin=3, c(
                                      'p_HCfev', 'p0','p1', 'f0', 'f1',
                                      'S0','S1','RR_F','RR_TF', 
                                      'lambda_0', 'final.inc', 'true.rate',
                                      'p_BC','p_BCpos','p_HC'
                                      ), n.iter=100000) 


#########################
######DIAGNOSTICS########
#########################

#check for convergence (quickly); do more thoroughly
#might have to run a few times to get full convergence (mixture models are difficult)
(gelman.diag(jpost.b, multivariate = FALSE))


#########################
#######POSTERIOR#########
#########################

#extract posterior summaries
sum.post.b <- summary(jpost.b)

#view posterior summery of quantiles
View(round(sum.post.b$quantiles,0))

#overall adjustment ratios
round(sum.post.b$quantiles[32:37,3]/((n.BCpos/persontime)*100000),1)

#final adjusted rates
round(sum.post.b$quantiles[32:37,c(3,1,5)],0)

#posterior probabilities
View(round(sum.post.b$quantiles[c(56:73),c(3,1,5)],2))

