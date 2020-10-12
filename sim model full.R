#clear everything from R
if(!is.bull(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

#load libraries
library(rjags)
load.module("mix") #this is for the normal mixture models, for sensitivity
library(R.matlab)

# set.seed(1) #735
set.seed(01) #1000
# set.seed(12345) #2000

#########################
####INPUT PARAMETERS#####
#########################
#load data
ds <- readMat('strataa_datasim.mat')

names(ds)
# [1] "TFpositive" "bctest"     "fever"      "riskfactor"  "soughtcare"

# define number of samples from HUS to use
# nsamp <- 735
nsamp <- 1000
# nsamp <- 2000

# create indices for random samples (from HUS)
samp.n1 <- sample(1:dim(ds$riskfactor)[1],size = nsamp, replace = F)
samp.n2 <- sample(1:dim(ds$riskfactor)[1],size = nsamp, replace = F)
samp.n3 <- sample(1:dim(ds$riskfactor)[1],size = nsamp, replace = F)
samp.n4 <- sample(1:dim(ds$riskfactor)[1],size = nsamp, replace = F)

# randomly sample from risk factor variable for each simulated scenario
rf1 <- ds$riskfactor[,1][samp.n1]
rf2 <- ds$riskfactor[,2][samp.n2]
rf3 <- ds$riskfactor[,3][samp.n3]
rf4 <- ds$riskfactor[,4][samp.n4]

# create parameters for beta distribution for risk factor probability
alpha.p1 <- c(length(rf1[which(ds$riskfactor[,1][samp.n1]==1)]),
              length(rf2[which(ds$riskfactor[,2][samp.n2]==1)]),
              length(rf3[which(ds$riskfactor[,3][samp.n3]==1)]),
              length(rf4[which(ds$riskfactor[,4][samp.n4]==1)]))
beta.p1  <- c(length(rf1[which(ds$riskfactor[,1][samp.n1]==0)]),
              length(rf2[which(ds$riskfactor[,2][samp.n2]==0)]),
              length(rf3[which(ds$riskfactor[,3][samp.n3]==0)]),
              length(rf4[which(ds$riskfactor[,4][samp.n4]==0)]))
alpha.p0 <- beta.p1
beta.p0  <- alpha.p1


# randomly sample from fever variable for each simulated scenario
fev1.rf1 <- ds$fever[,1][samp.n1][which(ds$riskfactor[,1][samp.n1]==1)]
fev1.rf2 <- ds$fever[,2][samp.n2][which(ds$riskfactor[,2][samp.n2]==1)]
fev1.rf3 <- ds$fever[,3][samp.n3][which(ds$riskfactor[,3][samp.n3]==1)]
fev1.rf4 <- ds$fever[,4][samp.n4][which(ds$riskfactor[,4][samp.n4]==1)]

fev0.rf1 <- ds$fever[,1][samp.n1][which(ds$riskfactor[,1][samp.n1]==0)]
fev0.rf2 <- ds$fever[,2][samp.n2][which(ds$riskfactor[,2][samp.n2]==0)]
fev0.rf3 <- ds$fever[,3][samp.n3][which(ds$riskfactor[,3][samp.n3]==0)]
fev0.rf4 <- ds$fever[,4][samp.n4][which(ds$riskfactor[,4][samp.n4]==0)]

# create parameters for beta distribution for fever probability
alpha.f1 <- c(sum(fev1.rf1),
              sum(fev1.rf2),
              sum(fev1.rf3),
              sum(fev1.rf4))
beta.f1 <- c(length(fev1.rf1)-alpha.f1[1],
             length(fev1.rf2)-alpha.f1[2],
             length(fev1.rf3)-alpha.f1[3],
             length(fev1.rf4)-alpha.f1[4])

alpha.f0 <- c(sum(fev0.rf1),
              sum(fev0.rf2),
              sum(fev0.rf3),
              sum(fev0.rf4))
beta.f0 <- c(length(fev0.rf1)-alpha.f0[1],
             length(fev0.rf2)-alpha.f0[2],
             length(fev0.rf3)-alpha.f0[3],
             length(fev0.rf4)-alpha.f0[4])

# randomly sample from healthcare seeking variable for each simulated scenario
HC1.rf.fev1 <- ds$soughtcare[,1][samp.n1][which((ds$riskfactor[,1][samp.n1]==1)&
                                                  (ds$fever[,1][samp.n1]==1))]
HC1.rf.fev2 <- ds$soughtcare[,2][samp.n2][which((ds$riskfactor[,2][samp.n2]==1)&
                                                  (ds$fever[,2][samp.n2]==1))]
HC1.rf.fev3 <- ds$soughtcare[,3][samp.n3][which((ds$riskfactor[,3][samp.n3]==1)&
                                                  (ds$fever[,3][samp.n3]==1))]
HC1.rf.fev4 <- ds$soughtcare[,4][samp.n4][which((ds$riskfactor[,4][samp.n4]==1)&
                                                  (ds$fever[,4][samp.n4]==1))]
HC0.rf.fev1 <- ds$soughtcare[,1][samp.n1][which((ds$riskfactor[,1][samp.n1]==0)&
                                                  (ds$fever[,1][samp.n1]==1))]
HC0.rf.fev2 <- ds$soughtcare[,2][samp.n2][which((ds$riskfactor[,2][samp.n2]==0)&
                                                  (ds$fever[,2][samp.n2]==1))]
HC0.rf.fev3 <- ds$soughtcare[,3][samp.n3][which((ds$riskfactor[,3][samp.n3]==0)&
                                                  (ds$fever[,3][samp.n3]==1))]
HC0.rf.fev4 <- ds$soughtcare[,4][samp.n4][which((ds$riskfactor[,4][samp.n4]==0)&
                                                  (ds$fever[,4][samp.n4]==1))]

# create parameters for beta distribution for healthcare seeking probability
alpha.S1 <- c(sum(HC1.rf.fev1),
              sum(HC1.rf.fev2),
              sum(HC1.rf.fev3),
              sum(HC1.rf.fev4))
beta.S1  <- c(length(HC1.rf.fev1)-alpha.S1[1],
             length(HC1.rf.fev2)-alpha.S1[2],
             length(HC1.rf.fev3)-alpha.S1[3],
             length(HC1.rf.fev4)-alpha.S1[4])

alpha.S0 <- c(sum(HC0.rf.fev1),
              sum(HC0.rf.fev2),
              sum(HC0.rf.fev3),
              sum(HC0.rf.fev4))
beta.S0  <- c(length(HC0.rf.fev1)-alpha.S0[1],
              length(HC0.rf.fev2)-alpha.S0[2],
              length(HC0.rf.fev3)-alpha.S0[3],
              length(HC0.rf.fev4)-alpha.S0[4])


# randomly sample from healthcare seeking variable among those with fever for each simulated scenario
HC.fev1 <- ds$soughtcare[,1][samp.n1][which(ds$fever[,1][samp.n1]==1)]
HC.fev2 <- ds$soughtcare[,2][samp.n2][which(ds$fever[,2][samp.n2]==1)]
HC.fev3 <- ds$soughtcare[,3][samp.n3][which(ds$fever[,3][samp.n3]==1)]
HC.fev4 <- ds$soughtcare[,4][samp.n4][which(ds$fever[,4][samp.n4]==1)]

# create parameters for beta distribution for healthcare seeking given fever probability
alpha.HCfev <- c(sum(HC.fev1),
                 sum(HC.fev2),
                 sum(HC.fev3),
                 sum(HC.fev4))
beta.HCfev  <- c(length(HC.fev1)-alpha.HCfev[1],
                 length(HC.fev2)-alpha.HCfev[2],
                 length(HC.fev3)-alpha.HCfev[3],
                 length(HC.fev4)-alpha.HCfev[4])

# set persontime
persontime <- rep(100000*2,4) 

# create parameters for beta distribution for receiving blood culture probability
alpha.bc <- apply(ds$bctest, MARGIN = 2, table)[2,]
beta.bc  <- apply(ds$bctest, MARGIN = 2, table)[1,]

# count the "crude" number of blood-culture-positive cases
n.BCpos <- c(sum(ds$TFpositive[,1], na.rm = T),
             sum(ds$TFpositive[,2], na.rm = T),
             sum(ds$TFpositive[,3], na.rm = T),
             sum(ds$TFpositive[,4], na.rm = T))

# create parameters for normal mixed distribution for blood culture sensitivity
#low prior abx use
mu1  <- mu3  <- c(.6,.4)
tau1 <- tau3 <- c(300,300)
pi1  <- pi3  <- c(.8,.2)
#high prior abx use
mu2  <- mu4  <- c(.75,.5)
tau2 <- tau4 <- c(300,300)
pi2  <- pi4  <- c(.2,.8)

#"crude" incidence
(n.BCpos/persontime)*100000
# [1]  64.5  64.5 220.0 194.5

#########################
#########MODEL###########
#########################
jcode <-"
model{

#BC sensitivity; normal mixture models, created from the data
p_BCpos[1] ~ dnormmix(mu1, tau1, pi1)
p_BCpos[2] ~ dnormmix(mu2, tau2, pi2)
p_BCpos[3] ~ dnormmix(mu3, tau3, pi3)
p_BCpos[4] ~ dnormmix(mu4, tau4, pi4)

#relative risk for positive BC (whether or not received BC)
RR_BCpos <-  exp(logRR_BCpos)
logRR_BCpos ~ dnorm(log(1.87),1/0.3731117^2)

for (a in 1:4){
#probability of receiving a blood culture (without adjustment)
p_BC[a] ~ dbeta(alpha.bc[a],beta.bc[a])
p_BC.adj[a] <- (1-(1-p_BC[a])*(1/RR_BCpos)) 

#obs rate --> true rate (not per 100,000 pyo yet), adjusted for sensitivity, pr(BC)
n.BCpos[a] ~ dpois(lambda.obs[a]) # positive BC results ~ poisson
lambda.obs[a] <- lambda.true[a]*p_BCpos[a]*p_BC.adj[a]
log(lambda.true[a]) <- beta0[a] + log(persontime[a])  #adjust for persontime
true.cases[a] ~ dpois(lambda.true[a])
beta0[a] ~ dnorm(0, 1/100000000)  #weakly informative prior for the intercept
true.rate[a] <- exp(beta0[a])  #true rate (per persontime, but not /100,000 pyo yet)
}

# distribution of relative risk (odds ratio) of risk factor for having typhoid
#boiled vs. unboiled drinking water
RR_TF <-  exp(logRR_TF)
logRR_TF ~ dnorm(log(5),10000)

for (a in 1:4){
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

#adjustment factors
adj.fact[a] <- final.inc[a]/((n.BCpos[a]/persontime[a])*100000)
}
}"

#########################
#######RUN MODEL#########
#########################

#list of all input variables
jdat <- list(n.BCpos=n.BCpos, persontime=persontime, 
               mu1=mu1,tau1=tau1,pi1=pi1, mu2=mu2,tau2=tau2,pi2=pi2, mu3=mu3,tau3=tau3,pi3=pi3, mu4=mu4,tau4=tau4,pi4=pi4,
               alpha.bc=alpha.bc, beta.bc=beta.bc, alpha.HCfev=alpha.HCfev, beta.HCfev=beta.HCfev,
               alpha.p1=alpha.p1, alpha.p0=alpha.p0, beta.p1=beta.p1, beta.p0=beta.p0,
               alpha.S0=alpha.S0, alpha.S1=alpha.S1, beta.S0=beta.S0, beta.S1=beta.S1,
               alpha.f0=alpha.f0, alpha.f1=alpha.f1, beta.f0=beta.f0, beta.f1=beta.f1)

#initialize model
jmod <- jags.model(textConnection(jcode), data=jdat, n.chains=3)

#set burn-in period
update(jmod,10000)

#fit full model, and decide which parameters to track
#thin 3, run 100000 iterations
jpost <- coda.samples(jmod, thin=3, c(
  'p_HCfev', 'p0','p1', 'f0', 'f1',
  'S0','S1','RR_F','RR_TF',
  'lambda_0', 'final.inc', 'true.rate',
  'p_BC.adj','p_BCpos','p_HC', 'adj.fact'
), n.iter=100000)

#########################
######DIAGNOSTICS########
#########################

#check for convergence (quickly); do more thoroughly 
#might have to run a few times to get full convergence
(gelman.diag(jpost, multivariate = FALSE))


#########################
#######POSTERIOR#########
#########################

#extract posterior summaries
sum.post <- summary(jpost)

#overall adjustment ratios
# round(sum.post$quantiles[14:17,c(3,1,5)],1)

#final adjusted rates
round(sum.post$quantiles[26:29,c(3,1,5)],0)

#posterior probabilities
round(sum.post$quantiles[c(42:53),c(3,1,5)],2)


# round(sum.post$quantiles[c(6:13,18:25,34:41,54:57),c(3,1,5)],2)

est.inc  <- data.frame(round(sum.post$quantiles[26:29,c(3,1,5)],0))
est.prob <- data.frame(round(sum.post$quantiles[c(42:53),c(3,1,5)],2))
# est.hcprob <- data.frame(round(sum.post$quantiles[c(6:13,18:25,34:41,54:57),c(3,1,5)],2))