#clear everything from R
if(!is.bull(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

#load libraries
library(rjags)
load.module("mix") #this is for the normal mixture models, for sensitivity
library(R.matlab)

# set.seed(0)
set.seed(012)
#########################
####INPUT PARAMETERS#####
#########################
#load data 
ds <- readMat('strataa_datasim.mat')

names(ds)
# [1] "TFpositive" "bctest"     "fever"      "riskfactor"  "soughtcare"

# define number of samples from HUS to use
# nsamp <- 735
# nsamp <- 1000
nsamp <- 2000

# create indices for random samples (from HUS)
samp.n1 <- sample(1:dim(ds$soughtcare)[1],size = nsamp, replace = F)
samp.n2 <- sample(1:dim(ds$soughtcare)[1],size = nsamp, replace = F)
samp.n3 <- sample(1:dim(ds$soughtcare)[1],size = nsamp, replace = F)
samp.n4 <- sample(1:dim(ds$soughtcare)[1],size = nsamp, replace = F)

# randomly sample from healthcare seeking variable among those with fever for each simulated scenario
HC.fev1 <- ds$soughtcare[,1][samp.n1][which(ds$fever[,1][samp.n1]==1)]
HC.fev2 <- ds$soughtcare[,2][samp.n2][which(ds$fever[,2][samp.n2]==1)]
HC.fev3 <- ds$soughtcare[,3][samp.n3][which(ds$fever[,3][samp.n3]==1)]
HC.fev4 <- ds$soughtcare[,4][samp.n4][which(ds$fever[,4][samp.n4]==1)]

# create parameters for beta distribution for healthcare seeking given fever probability
alpha.hc <- c(sum(HC.fev1),
                 sum(HC.fev2),
                 sum(HC.fev3),
                 sum(HC.fev4))
beta.hc  <- c(length(HC.fev1)-alpha.hc[1],
                 length(HC.fev2)-alpha.hc[2],
                 length(HC.fev3)-alpha.hc[3],
                 length(HC.fev4)-alpha.hc[4])

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

# create parameters for (simple) normal distribution for blood culture sensitivity
#mean 0.59, 95%CI: 0.54-0.64
mu  <- c(.59)
tau <- c(1/0.0006507705)

#"crude" incidence
(n.BCpos/persontime)*100000
# [1]  64.5  64.5 220.0 194.5
#########################
#########MODEL###########
#########################
jcode <-"
model{

#BC sensitivity; normal mixture models, created from the data
p_BCpos[1] ~ dnorm(mu, tau)

for (a in 1:4){
#probability of receiving a blood culture (without adjustment)
p_BC[a] ~ dbeta(alpha.bc[a],beta.bc[a])
p_HC[a] ~ dbeta(alpha.hc[a],beta.hc[a])

#obs rate --> true rate (not per 100,000 pyo yet), adjusted for sensitivity, pr(BC)
n.BCpos[a] ~ dpois(lambda.obs[a]) # positive BC results ~ poisson
lambda.obs[a] <- lambda.true[a]*p_BCpos*p_BC[a]*p_HC[a]
log(lambda.true[a]) <- beta0[a] + log(persontime[a])  #adjust for persontime
true.cases[a] ~ dpois(lambda.true[a])
beta0[a] ~ dnorm(0, 1/100000000)  #weakly informative prior for the intercept
true.inc[a] <- exp(beta0[a])*100000  #true rate 

#adjustment factors
adj.fact[a] <- true.inc[a]/((n.BCpos[a]/persontime[a])*100000)

}
}"

#########################
#######RUN MODEL#########
#########################

#list of all input variables
jdat <- list(n.BCpos=n.BCpos, persontime=persontime, 
               mu=mu,tau=tau,
               alpha.bc=alpha.bc, beta.bc=beta.bc, alpha.hc=alpha.hc, beta.hc=beta.hc)

#initialize model
jmod <- jags.model(textConnection(jcode), data=jdat, n.chains=3)

#set burn-in period
update(jmod,10000)

#fit full model, and decide which parameters to track
jpost <- coda.samples(jmod, thin=3, c(
 'true.inc',
  'p_BC','p_BCpos','p_HC', 'adj.fact'
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
round(sum.post$quantiles[14:17,c(3,1,5)],0)

#posterior probabilities
round(sum.post$quantiles[c(5:13),c(3,1,5)],2)


est.inc  <- data.frame(round(sum.post$quantiles[14:17,c(3,1,5)],0))
est.prob <- data.frame(round(sum.post$quantiles[c(5:8,9,9,9,9,10:13),c(3,1,5)],2))