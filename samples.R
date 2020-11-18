###############################
## Gibbs-Metropolis sampling ##
###############################

## Load libraries
library(parallel)
library(coda)
options(mc.cores = 4)

## Load source
source("gibbs-sampler.R")

## Load data
load("nh_djf_rcp45_tas.RData")

## Sampling parameters
n.chains  = 4
n.samples = 10000  ## Total number of samples
thin      = 40     ## Thinning
burn.in   = 10000  ## Burn-in period

## Extract lon and lat
lon = sort(unique(cmip5$lon))
lat = sort(unique(cmip5$lat))

## Create grid for parallel execution
grid = expand.grid(lon = lon, lat = lat)

## Uncertainty multipliers
kappa = list()
kappa$kappad = 1.2
kappa$kappan = 1.2
kappa$kappaw = 1.2

## Parameters
parameters = c("muh","muf","muw","beta","tauh","tauf",
               "tauw","psisq","thetasq","nuh","nuf")

## Hyperparameters
hyper.parameters           = list()
hyper.parameters$amuh      = 0     ## Expectation of historical mean climate
hyper.parameters$bmuh      = 1e-6  ## Precision of historical mean climate
# hyper.parameters$amuf      = 0     ## Expectation of future mean climate
hyper.parameters$bmuf      = 1e-6  ## Precision of future mean climate
hyper.parameters$abeta     = 1     ## Expectation of emergent constraint
hyper.parameters$bbeta     = 1e-6  ## Precision of emergent constraint
hyper.parameters$atauh     = 1e-3  ## Shape of historical structural uncertainty
hyper.parameters$btauh     = 1e-3  ## Rate of historical structural uncertainty
hyper.parameters$atauf     = 1e-3  ## Shape of future structural uncertainty
hyper.parameters$btauf     = 1e-3  ## Rate of future structural uncertainty
hyper.parameters$apsisq    = 1e-3  ## Shape of historical internal variability
hyper.parameters$bpsisq    = 1e-3  ## Rate of historical internal variability
hyper.parameters$athetasq  = 1e-3  ## Shape of future internal variability
hyper.parameters$bthetasq  = 1e-3  ## Rate of future internal variability
hyper.parameters$anuh      = 1     ## Shape of historical degrees-of-freedom
hyper.parameters$bnuh      = 1/13  ## Rate of historical degrees-of-freedom
hyper.parameters$anuf      = 1     ## Shape of future degrees-of-freedom
hyper.parameters$bnuf      = 1/13  ## Rate of future degrees-of-freedom
hyper.parameters$atauw     = 1e-3  ## Shape of reanalysis structural uncertainty
hyper.parameters$btauw     = 1e-3  ## Rate of reanalysis structural uncertainty

## Proposal distributions
proposal.parameters         = list()
proposal.parameters$lambdah = 0.10  ## Historical proposal rate
proposal.parameters$lambdaf = 0.10  ## Future     proposal rate

## Function for parallel computation
fun = function(i, n.chains, n.samples, thin, burn.in, cmip5, reanalysis, grid,
               kappa, hyper.parameters, proposal.parameters) {

    ## Set seed
    set.seed(i)

    ## Extract data
    buffer  = cmip5[cmip5$lon == grid$lon[i] & cmip5$lat == grid$lat[i],]
    rbuffer = reanalysis[reanalysis$lon == grid$lon[i] &
                             reanalysis$lat == grid$lat[i],]

    ## Compute summary statistics
    m       = nlevels(buffer$model)
    nxsmr   = with(buffer, tapply(tas, list(model,experiment), length))
    mxsmr   = with(buffer, tapply(tas, list(model,experiment), mean))
    r       = nlevels(rbuffer$reanalysis)
    mwhr    = mean(rbuffer$tas)
    xsmr    = with(buffer, tapply(tas, list(model,experiment), function(x) x))

    ## Summary statistics
    summary.statistics         = list()
    summary.statistics$nxhm    = m
    summary.statistics$nxhmr   = nxsmr[,"historical"]
    summary.statistics$mxhmr   = mxsmr[,"historical"]
    summary.statistics$nxfmr   = nxsmr[,"rcp45"]
    summary.statistics$mxfmr   = mxsmr[,"rcp45"]
    summary.statistics$nwhr    = r
    summary.statistics$mwhr    = mean(rbuffer$tas)
    summary.statistics$whr     = rbuffer$tas
    summary.statistics$xhmr    = xsmr[,"historical"]
    summary.statistics$xfmr    = xsmr[,"rcp45"]

    ## Sampling loop
    samples = list()
    for (j in 1:n.chains) {

        statusj = 1
        while(statusj > 0) {

            ## Initial values
            initial.values = list()
            initial.values$muh     = rnorm(1, 0, 20)
            initial.values$muf     = rnorm(1, initial.values$muh, 20)
            initial.values$muw     = rnorm(1, 0, 20)
            initial.values$beta    = rnorm(1, 1, 2)
            initial.values$tauh    = rgamma(1, 0.457, 0.0239)
            initial.values$tauf    = rgamma(1, 0.457, 0.0239)
            initial.values$tauw    = rgamma(1, 0.457, 0.0239)
            initial.values$psisq   = rgamma(1, 0.457, 0.0239)
            initial.values$thetasq = rgamma(1, 0.457, 0.0239)
            initial.values$nuh     = rexp(1, 1/summary.statistics$nxhm)
            initial.values$nuf     = rexp(1, 1/summary.statistics$nxhm)

            ## Sampling
            samples[[j]] = gibbs.sampler(burn.in + n.samples,
                                         summary.statistics, kappa,
                                         initial.values, hyper.parameters,
                                         proposal.parameters)

            ## Update status
            statusj = samples[[j]]$status

        } ## statusj

    } ## j

    ## Extract samples
    mcmc.buffer = lapply(samples, function(x) mcmc(x$samples))
    mcmc.buffer = as.mcmc.list(mcmc.buffer)
    mcmc.buffer = window(mcmc.buffer, start = burn.in + 1)

    ## Extract acceptance probabilities
    acceptance = lapply(samples, function(x) mcmc(x$acceptance))
    acceptance = as.mcmc.list(acceptance)
    acceptance = window(acceptance, start = burn.in + 1)

    ## Gelman diagnostics
    gelman = gelman.diag(mcmc.buffer,
                         autoburnin = FALSE, multivariate = FALSE)
    if (any(is.na(gelman$psrf[,1]))) {
        not.converged = TRUE
    } else {
        not.converged = any(gelman$psrf[,1] > 1.10)
    }

    ## Initialise counter
    rho = burn.in + n.samples

    ## Check convergence and continue sampling if necessary
    while(not.converged) {

        ## Sampling loop
        samples = list()
        for (j in 1:n.chains) {

            ## Initial values
            initial.values = as.list(mcmc.buffer[[j]][n.samples, parameters])

            statusj = 1
            while(statusj > 0) {

                ## Sampling
                samples[[j]] = gibbs.sampler(n.samples, summary.statistics,
                                             kappa, initial.values,
                                             hyper.parameters,
                                             proposal.parameters)

                ## Update status
                statusj = samples[[j]]$status

             } ## statusj

        } ## j

        ## Extract samples
        mcmc.buffer = lapply(samples, function(x) mcmc(x$samples))
        mcmc.buffer = as.mcmc.list(mcmc.buffer)

        ## Extract acceptance probabilities
        acceptance = lapply(samples, function(x) mcmc(x$acceptance))
        acceptance = as.mcmc.list(acceptance)

        ## Gelman diagnostics
        gelman = gelman.diag(mcmc.buffer,
                             autoburnin = FALSE, multivariate = FALSE)
        if (any(is.na(gelman$psrf[,1]))) {
            not.converged = TRUE
        } else {
            not.converged = any(gelman$psrf[,1] > 1.10)
        }

        ## Increment counter
        rho = rho + n.samples

    } ## status

    ## Thin samples
    mcmc.buffer = window(mcmc.buffer,
                         start = mcpar(mcmc.buffer[[1]])[1] + thin - 1,
                         thin = thin)
    mcmc.buffer = as.matrix(mcmc.buffer)

    ## Thin acceptance probabilities
    acceptance = window(acceptance,
                        start = mcpar(acceptance[[1]])[1] + thin - 1,
                        thin = thin)
    acceptance = as.matrix(acceptance)

    ## Prediction
    yf   = rnorm(nrow(mcmc.buffer),
                 mcmc.buffer[,"muf"] +
                     mcmc.buffer[,"beta"]*(mcmc.buffer[,"yh"] - mcmc.buffer[,"muh"]),
                 kappa$kappad/sqrt(mcmc.buffer[,"tauf"]))
    phia = rgamma(nrow(mcmc.buffer),
                  0.5*mcmc.buffer[,"nuf"]/kappa$kappan^2,
                  0.5*mcmc.buffer[,"nuf"]/kappa$kappan^2*mcmc.buffer[,"thetasq"])
    yfa  = rnorm(nrow(mcmc.buffer), yf, 1/sqrt(mcmc.buffer[,"taua"]*phia))
    mcmc.buffer = cbind(mcmc.buffer,yf,phia,yfa)

    ## Return samples
    return(list(samples = mcmc.buffer, acceptance = acceptance,
                psrf = gelman$psrf, rho = rho))

}

## Error handling
wrapper = function(i, n.chains, n.samples, thin, burn.in, cmip5, reanalysis,
                   grid, kappa, hyper.parameters, proposal.parameters)
    try(fun(i, n.chains, n.samples, thin, burn.in, cmip5, reanalysis,
            grid, kappa, hyper.parameters, proposal.parameters), silent = TRUE)

## Parallel sampling
samples = mclapply(X = 1:nrow(grid), FUN = wrapper,
                   n.chains = n.chains, n.samples = n.samples, thin = thin,
                   burn.in = burn.in, cmip5 = cmip5, reanalysis = reanalysis,
                   grid = grid, kappa = kappa,
                   hyper.parameters = hyper.parameters,
                   proposal.parameters = proposal.parameters)

## Save samples
save(samples, n.chains, n.samples, thin, burn.in,
     cmip5, reanalysis, grid, kappa, hyper.parameters, proposal.parameters,
     file = "samples.RData")
