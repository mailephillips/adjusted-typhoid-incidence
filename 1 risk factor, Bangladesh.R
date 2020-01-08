

if(!is.bull(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

library(rjags)
library(coda)
library(mcmcplots)
load.module("mix")

set.seed(1)

df.visit <- read.csv(file = "~/Desktop/Typhoid Underreporting/Bangladesh/Bangladesh Results DS 8.22.19.csv")
#df.enrol <- read.csv(file = "~/Desktop/Typhoid Underreporting/Bang_Enroll_2.6.19.csv")

jcode.b <-"
model{

p_BCpos[1] ~ dnormmix(mu.u5, tau.u5, pi.u5)
p_BCpos[2] ~ dnormmix(mu.5_9, tau.5_9, pi.5_9)
p_BCpos[3] ~ dnormmix(mu.10_14, tau.10_14, pi.10_14)
p_BCpos[4] ~ dnormmix(mu.15_29, tau.15_29, pi.15_29)
p_BCpos[5] ~ dnormmix(mu.30plus, tau.30plus, pi.30plus)
p_BCpos[6] ~ dnormmix(mu.all, tau.all, pi.all)

RR_BCpos <-  exp(logRR_BCpos)
logRR_BCpos ~ dnorm(log(1.87),1/0.3731117^2)

for (a in 1:6){
p_BC[a] ~ dbeta(alpha.bc[a],beta.bc[a])

n.BCpos[a] ~ dpois(lambda.obs[a])
lambda.obs[a] <- lambda.true.BCpos.BC[a]*p_BCpos[a]*(p_BC[a]+(1-p_BC[a])*(1/RR_BCpos))
log(lambda.true.BCpos.BC[a]) <- beta0.BCpos.BC[a] + log(persontime[a])
true.cases.BCpos.BC[a] ~ dpois(lambda.true.BCpos.BC[a])
beta0.BCpos.BC[a] ~ dnorm(0, 1/100000000)
true.rate[a] <- exp(beta0.BCpos.BC[a])
}


RR_TF <-  exp(logRR_TF)
logRR_TF ~ dnorm(log(7.6),1/0.6372431^2)

for (a in 1:6){
p_HCfev[a] ~dbeta(alpha.HCfev[a], beta.HCfev[a])
p1[a] ~ dbeta(alpha.p1[a], beta.p1[a])
p0[a] ~ dbeta(alpha.p0[a], beta.p0[a])
S1[a] ~ dbeta(alpha.S1[a], beta.S1[a])
S0[a] ~ dbeta(alpha.S0[a], beta.S0[a])
f1[a] ~ dbeta(alpha.f1[a], beta.f1[a])
f0[a] ~ dbeta(alpha.f0[a], beta.f0[a])

RR_F[a] <- f1[a]/f0[a]

lambda_0[a] <- true.rate[a]*((1/(p0[a]+RR_F[a]*p1[a]))/p_HCfev[a])*(((S1[a]*RR_F[a]*p1[a])+(S0[a]*p0[a]))/((S1[a]*RR_TF*p1[a])+(S0[a]*p0[a])))
final.inc[a] <- (RR_TF*lambda_0[a]*p1[a]+lambda_0[a]*p0[a])*100000
}


}"

mu.u5  <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.u5.finaltbl.csv")[1,-1]
tau.u5 <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.u5.finaltbl.csv")[2,-1]
pi.u5  <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.u5.finaltbl.csv")[3,-1]

mu.5_9  <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.5_9.finaltbl.csv")[1,-1]
tau.5_9 <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.5_9.finaltbl.csv")[2,-1]
pi.5_9  <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.5_9.finaltbl.csv")[3,-1]

mu.10_14  <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.10_14.finaltbl.csv")[1,-1]
tau.10_14 <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.10_14.finaltbl.csv")[2,-1]
pi.10_14  <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.10_14.finaltbl.csv")[3,-1]

mu.15_29  <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.15_29.finaltbl.csv")[1,-1]
tau.15_29 <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.15_29.finaltbl.csv")[2,-1]
pi.15_29  <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.15_29.finaltbl.csv")[3,-1]

mu.30plus  <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.30plus.finaltbl.csv")[1,-1]
tau.30plus <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.30plus.finaltbl.csv")[2,-1]
pi.30plus  <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.30plus.finaltbl.csv")[3,-1]

mu.all  <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.all.finaltbl.csv")[1,-1]
tau.all <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.all.finaltbl.csv")[2,-1]
pi.all  <- read.csv("~/Desktop/Typhoid Underreporting/Bangladesh/bang.sens.all.finaltbl.csv")[3,-1]

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


bang_results <- df.visit

bangu5 <- bang_results[which(bang_results$agecat==1),]
bang5_9 <- bang_results[which(bang_results$agecat==2),]
bang10_14 <- bang_results[which(bang_results$agecat==3),]
bang15_29 <- bang_results[which(bang_results$agecat==4),]
bang30plus <- bang_results[which(bang_results$agecat==5),]
bangall <- bang_results

n.BC <-  c(length(bangu5$result),length(bang5_9$result),length(bang10_14$result),
           length(bang15_29$result),length(bang30plus$result),length(bangall$result))
persontime <- c(11262*2,
                10205*2,
                10618*2,
                35773*2,
                43460*2,
                11262*2+10205*2+10618*2+35773*2+43460*2) #table 4 from paper

BC.obs.u5 <- bangu5$result
BC.obs.5_9 <- bang5_9$result
BC.obs.10_14 <- bang10_14$result
BC.obs.15_29 <- bang15_29$result
BC.obs.30plus <- bang30plus$result
BC.obs.all <- bangall$result

n.BCpos.u5 <- sum(BC.obs.u5)
n.BCpos.5_9 <- sum(BC.obs.5_9)
n.BCpos.10_14 <- sum(BC.obs.10_14)
n.BCpos.15_29 <- sum(BC.obs.15_29)
n.BCpos.30plus <- sum(BC.obs.30plus)
n.BCpos.all <- sum(BC.obs.all)

n.BCpos <- c(n.BCpos.u5,
             n.BCpos.5_9,
             n.BCpos.10_14,
             n.BCpos.15_29,
             n.BCpos.30plus,
             n.BCpos.all)

(n.BCpos/persontime)*100000

jdat.b <- list(n.BCpos=n.BCpos, persontime=persontime, 
               mu.u5=mu.u5,tau.u5=tau.u5,pi.u5=pi.u5,mu.5_9=mu.5_9,tau.5_9=tau.5_9,pi.5_9=pi.5_9,
               mu.10_14=mu.10_14,tau.10_14=tau.10_14,pi.10_14=pi.10_14,mu.15_29=mu.15_29,tau.15_29=tau.15_29,pi.15_29=pi.15_29,
               mu.30plus=mu.30plus,tau.30plus=tau.30plus,pi.30plus=pi.30plus,mu.all=mu.all,tau.all=tau.all,pi.all=pi.all,
               alpha.bc=alpha.bc, beta.bc=beta.bc, alpha.HCfev=alpha.HCfev, beta.HCfev=beta.HCfev,
             alpha.p1=alpha.p1, alpha.p0=alpha.p0, beta.p1=beta.p1, beta.p0=beta.p0,
             alpha.S0=alpha.S0, alpha.S1=alpha.S1, beta.S0=beta.S0, beta.S1=beta.S1,
             alpha.f0=alpha.f0, alpha.f1=alpha.f1, beta.f0=beta.f0, beta.f1=beta.f1)

jmod.b <- jags.model(textConnection(jcode.b), data=jdat.b, n.chains=3)
update(jmod.b,10000)

jpost.b <- coda.samples(jmod.b, thin=3, c(
                                      'p_HCfev', 'p0','p1', 'f0', 'f1',
                                      'S0','S1','RR_F','RR_TF', 
                                      'lambda_0', 'final.inc', 'true.rate'
                                      ), n.iter=100000) 
(gelman.diag(jpost.b, multivariate = FALSE))


sum.post.b2 <- summary(jpost.b)
View(round(sum.post.b2$quantiles,0))
round(sum.post.b2$quantiles[32:37,3]/((n.BCpos/persontime)*100000),1)
round(sum.post.b2$quantiles[32:37,c(3,1,5)],0)
#write.csv(sum.post.b$quantiles,"~/Desktop/Typhoid Underreporting/nep.mod.quants.csv")


# plot.ds <- as.data.frame(round(sum.post.b2$quantiles[32:37,c(3,1,5)],0))
# write.csv(plot.ds, "~/Desktop/Typhoid Underreporting/bang.inc.csv")

