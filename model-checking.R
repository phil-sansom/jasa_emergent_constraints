####################
## Model checking ##
####################

## Load libraries
library(mcmcse)

## Load data
load("samples.RData")

## Monte Carlo Standard Errors
mcses = lapply(samples, function(x) mcse.mat(x$samples))
mcses = simplify2array(mcses)
mcses = mcses[c(1:14,67:69),2,]

## Extract samples
samples = lapply(samples, function(x) x$samples[,c(1:14,67:69)])
samples = simplify2array(samples)
gc()
means = apply(samples, c(2,3), mean)
sds   = apply(samples, c(2,3), sd  )
quantile(as.numeric(mcses/abs(means)), 0.95)
quantile(as.numeric(mcses/sds       ), 0.95)

## Degrees-of-freedom
mode = function(x)
    which.max(hist(x, breaks = seq(-0.5,200.5,1), plot = FALSE)$counts) - 1
summary(apply(samples[,"nuh",], 2, mode))
summary(apply(samples[,"nuf",], 2, mode))

## Reanalysis standard error
summary(apply(1/sqrt(samples[,"tauw",]), 2, mean))

## Correlation matrices
correlations = lapply(1:dim(samples)[3], function(x) cor(samples[,,x]))
correlations = simplify2array(correlations)
corr.fit = apply(correlations, c(1,2), mean)
corr.lwr = apply(correlations, c(1,2), quantile, probs = 0.25)
corr.upr = apply(correlations, c(1,2), quantile, probs = 0.75)
