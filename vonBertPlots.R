#vonBertPlots


# Tomo Eguchi
# 3 April 2014

#rm(list=ls())
library(rjags)
sysInfo <- Sys.info()
ifelse(sysInfo[1] == "Windows", source('~/R/TomosFunctions.R'),
       source('~/Documents/R/TomosFunctions.R'))

D <- dirSelector()

# first load the datasets
load(paste(D$Rdir, "CMR_Australia/RData/vonBertRlinfRk_2014-04-03.RData", sep=""))

# check for the convergence:
# maxReq sample sizes
max.req.n <- vector(mode = "numeric", length = n.chains)

for (k in 1:n.chains){
  max.req.n[k] <- max(r.diag[[k]]$resmatrix[,2])
}

# scale reduction factor
upper.Rhat <- as.matrix(g.diag$psrf)[,2]
max.Rhat <- max(upper.Rhat)
max.Rhat.name <- names(upper.Rhat[upper.Rhat == max.Rhat])
mv.Rhat <- g.diag$mpsrf

# plot growth functions
dataFile <- paste(D$Dtomo, "turtles/Australia/data/BPKAPT_Master5_growth_01-Apr-2014.txt", sep = "")
g_data <- read.table(dataFile, sep = '\t', header = FALSE)
tme <- as.matrix(g_data[, 6:dim(g_data)[2]])
sze <- as.matrix(g_data[, 2:5])

stats.A <- matrix(data = NA, nrow = dim(g_data)[1], ncol = 5)

for (k in 1:dim(g_data)[1]){
  stats.A[k,] <- statMat$quantiles[rownames(statMat$quantiles) == paste("A[", k, "]", sep = "")]

}
tme2 <- stats.A[,3] + cbind(matrix(data = 0, ncol=1,nrow=109), tme)

#stats.LinfMu <- statMat$quantiles[rownames(statMat$quantiles) == "Linf",]
stats.LinfMu <- statMat$quantiles[rownames(statMat$quantiles) == "LinfMu",]
stats.kAlpha <- statMat$quantiles[rownames(statMat$quantiles) == "kAlpha",]
stats.kBeta <- statMat$quantiles[rownames(statMat$quantiles) == "kBeta",]
median.k <- qbeta(0.5, stats.kAlpha[3], stats.kBeta[3])
low.k  <- qbeta(0.5, stats.kAlpha[1], stats.kBeta[1])
high.k <- qbeta(0.5, stats.kAlpha[5], stats.kBeta[5])

linTime <- seq(from = 0, to = 40, by = 0.1)
LtEstMin <- stats.LinfMu[1] * (1 - exp(-median.k * linTime))
LtEstMax <- stats.LinfMu[5] * (1 - exp(-median.k * linTime))

op <- par
par(bty = "l")
plot(linTime, LtEstMax, "n", 
     main = "Growths",
     xlab = "Time (yrs)",
     ylab = "CCL (cm)",
     ylim = c(0, 120))
lines(linTime, LtEstMin)
lines(linTime, LtEstMax)

for (k in 1:dim(g_data)[1]){  
  points(tme2[k,], sze[k, ], pch = 1, "b")

}

lines(linTime, rep(stats.LinfMu[1], times = length(linTime)))
lines(linTime, rep(stats.LinfMu[5], times = length(linTime)))
par(op)