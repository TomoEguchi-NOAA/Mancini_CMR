# fitting von Bertalanffy growth function using rjags

# Tomo Eguchi
# 29 January 2013


rm(list=ls())
library(rjags)
library(coda)

save <- TRUE
## set MCMC parameters
n.adapt <- 10000
n.update <- 50000
n.iter <- 50000
n.thin <- 5
n.chains <- 5

sysInfo <- Sys.info()
if (sysInfo[1] == "Windows") source('~/R/TomosFunctions.R')

#source('~/R/R_Work/TomosFunctions.R')
D <- dirSelector()

runDate <- Sys.Date()

# define model and data file names
modelName <- "Model_FlinfRk_L.txt"
dirName <- paste(D$Dtomo, 'turtles/Australia/data/', sep="")
fileName <- 'BPKAPT_Master5_growth_01-Apr-2014.txt'
dataName <- paste(dirName, fileName, sep = '')

if (save == TRUE) saveFname <- paste(D$Rdir, 
                                     "CMR_Australia/RData/vonBertFlinfRk_", 
                                     runDate, ".RData", sep = "")

dat <- read.table(dataName, 
                  header = FALSE, sep = '\t')

# maximum number of recaptures:
maxN <- 4

## parameters to monitor - when this is changed, make sure to change
## summary statistics indices at the end of this script. 
parameters <- c('CV', 'k', 'A', 
                'Linf',  
                'Shape', 'rate', 'kAlpha', 'kBeta', 
                'deviance')

bugs.data <- list(nIndiv = dim(dat)[1],
                  n = dat[, 1],
                  L = dat[, 2:(maxN+1)],
                  t = dat[, (maxN+2):dim(dat)[2]])

initsFunction <- function(d){
  kAlpha <- runif(1, 0, 100)
  kBeta <- runif(1, 0, 100)
  k <- rbeta(dim(dat)[1], 1, 1)
  CV <- rbeta(1, 1, 1)
  Shape <- runif(1, 0, 100)
  rate <- runif(1, 0, 100)
  A <- list(kAlpha = kAlpha, kBeta = kBeta, k = k, CV = CV,
            Shape = Shape, rate = rate)
  return(A)
}

# Create a list of length nChains with z in initsFunction
#inList <- list(CH, f)
inits <- lapply(c(1:n.chains), initsFunction)

jm <- jags.model(modelName, 
                 data = bugs.data, 
                 inits, 
                 n.chains = length(inits), 
                 n.adapt = n.adapt)

if (save == TRUE) save(list = ls(all = TRUE), file = saveFname)

#parallel.seeds("base::BaseRNG", length(inits))
#update(jm, n.iter = n.update)
update(jm, n.iter = n.iter)
load.module("dic")

# ##############################################
# Using coda.samples - 
zm <- coda.samples(jm,
                   variable.names = parameters, 
                   n.iter = n.iter, 
                   thin = n.thin)

if (save == TRUE) save(list = ls(all = TRUE), file = saveFname)

dicOut <- dic.samples(jm, 
                      n.iter=n.iter, 
                      thin = n.thin,
                      type = "pD")

if (save == TRUE) save(list = ls(all = TRUE), file = saveFname)

# use raftery.diag to figure out the required chain lengths:
r.diag <- raftery.diag(zm)

# test also with other convergence diagnostic tools:
g.diag <- gelman.diag(zm)
h.diag <- heidel.diag(zm)
if (save == TRUE) save(list = ls(all = TRUE), file = saveFname)

#
# zm contains MCMC samples - for each chain. To access
# them, use zm[[1]], etc to the # chains, which is 
# defined by the number of initial values above
# it may be easier to use jags.samples rather than 
# coda.samples - see below
#
# use summary(zm) to look at summary statistics of parameters
# statMat[1] = Mean,SD, statMat[2] = quantiles
statMat <- summary(zm)		
if (save == TRUE) save(list = ls(all = TRUE), file = saveFname)

# plot the results
kAlphaMean <- statMat$statistics[rownames(statMat$statistics) == 'kAlpha', 1]
kBetaMean <- statMat$statistics[rownames(statMat$statistics) == 'kBeta', 1]

kAlphaMedian <- statMat$quantiles[rownames(statMat$quantiles) == 'kAlpha', 3]
kBetaMedian <- statMat$quantiles[rownames(statMat$quantiles) == 'kBeta', 3]

dk <- dbeta(seq(0.01, 1, by = 0.01), kAlphaMedian, kBetaMedian)
plot(seq(0.01, 1, by = 0.01), dk, type = 'l')

ShapeMean <- statMat$statistics[rownames(statMat$statistics) == 'Shape', 1]
rateMean <- statMat$statistics[rownames(statMat$statistics) == 'rate', 1]
dA <- dgamma(seq(0.01, 4, by = 0.01), shape = ShapeMean, rate = rateMean)
plot(seq(0.01, 4, by = 0.01), dA, type = 'l')

Linf <- matrix(data = NA, nrow = 10000, ncol = 4)
for (k in 1:4){
  Linf[,k] <- zm[[k]][,colnames(zm[[k]]) == 'Linf']
}
Linf2 <- as.vector(Linf)
hist(Linf2)

if (save == TRUE) save(list = ls(all = TRUE), file = saveFname)

