## Load libraries
library(compiler)

## Gibb's sampler
gibbs.sampler = function(n, summary.statistics, kappa, initial.values,
                         hyper.parameters, proposal.parameters) {

    ## Summary statistics
    nxhm     = summary.statistics$nxhm     ## Number of models
    nxhmr    = summary.statistics$nxhmr    ## Number of historical runs
    mxhmr    = summary.statistics$mxhmr    ## Mean of historical runs
    nxfmr    = summary.statistics$nxfmr    ## Number of future runs
    mxfmr    = summary.statistics$mxfmr    ## Mean of future runs
    nwhr     = summary.statistics$nwhr     ## Number of reanalyses
    mwhr     = summary.statistics$mwhr     ## Mean of reanalyses
    whr      = summary.statistics$whr      ## Reanalyses
    xhmr     = summary.statistics$xhmr     ## Historical runs
    xfmr     = summary.statistics$xfmr     ## Future runs

    ## Uncertainty multipliers
    kdsq     = kappa$kappad^2
    knsq     = kappa$kappan^2
    kwsq     = kappa$kappaw^2

    ## Hyper-parameters
    amuh     = hyper.parameters$amuh       ## Expectation of muh
    bmuh     = hyper.parameters$bmuh       ## Precision   of muh
    bmuf     = hyper.parameters$bmuf       ## Precision   of muf
    abeta    = hyper.parameters$abeta      ## Expectation of beta
    bbeta    = hyper.parameters$bbeta      ## Expectation of beta
    atauh    = hyper.parameters$atauh      ## Shape of tauh
    btauh    = hyper.parameters$btauh      ## Rate  of tauh
    atauf    = hyper.parameters$atauf      ## Shape of tauf
    btauf    = hyper.parameters$btauf      ## Rate  of tauf
    apsisq   = hyper.parameters$apsisq     ## Shape of psisq
    bpsisq   = hyper.parameters$bpsisq     ## Rate  of psisq
    athetasq = hyper.parameters$athetasq   ## Shape of thetasq
    bthetasq = hyper.parameters$bthetasq   ## Rate  of thetasq
    anuh     = hyper.parameters$anuh       ## Shape of nuh
    bnuh     = hyper.parameters$bnuh       ## Rate  of nuh
    anuf     = hyper.parameters$anuf       ## Shape of nuf
    bnuf     = hyper.parameters$bnuf       ## Rate  of nuf
    atauw    = hyper.parameters$atauw      ## Shape of tauw
    btauw    = hyper.parameters$btauw      ## Rate  of tauw

    ## Initial values
    muh      = initial.values$muh          ## Historical expected climate
    muf      = initial.values$muf          ## Future expected climate
    beta     = initial.values$beta         ## Emergent constraint
    tauh     = initial.values$tauh         ## Historical structural uncertainty
    tauf     = initial.values$tauf         ## Future structural uncertainty
    psisq    = initial.values$psisq        ## Historical internal variability
    thetasq  = initial.values$thetasq      ## Future internal variability
    nuh      = initial.values$nuh          ## Historical degrees-of-freedom
    nuf      = initial.values$nuf          ## Future degrees-of-freedom
    muw      = initial.values$muw          ## Reanalysis expectation
    tauw     = initial.values$tauw         ## Reanalysis structural uncertainty

    taudh    = tauh/kdsq
    tauz     = tauw/kwsq
    nuha     = nuh /knsq

    ## Initial values
    xhm  = rnorm(nxhm, muh                   , sqrt(1/tauh))
    xfm  = rnorm(nxhm, muf + beta*(xhm - muh), sqrt(1/tauf))
    taum = rgamma(nxhm, nuh /2, nuh *psisq  /2)
    phim = rgamma(nxhm, nuf /2, nuf *thetasq/2)
    taua = rgamma(1   , nuha/2, nuha*psisq  /2)

    yh   = rnorm(1, muh                  , sqrt(1/taudh))

    ## Proposal parameters
    lambdah  = proposal.parameters$lambdah ## Historical proposal rate
    lambdaf  = proposal.parameters$lambdaf ## Future     proposal rate

    ## Derived parameters
    m        = nxhm                        ## Number of models
    r        = nwhr                        ## Number of reanalyses

    ## Create objects to hold results
    samples = matrix(data = NA, nrow = n, ncol = 14 + 4*m,
                     dimnames = list(NULL, c("muh","muf","beta","tauh","tauf",
                                     "psisq","thetasq","nuh","nuf","muw","tauw",
                                     "yha","yh","taua",
                                     paste("xhm" , 1:m, sep = "."),
                                     paste("xfm" , 1:m, sep = "."),
                                     paste("taum", 1:m, sep = "."),
                                     paste("phim", 1:m, sep = "."))))
    ixhm  = 14 + 1:m
    ixfm  = 14 + m + 1:m
    itaum = 14 + 2*m + 1:m
    iphim = 14 + 3*m + 1:m
    acceptance = matrix(data = NA, nrow = n, ncol = 2,
                        dimnames = list(NULL, c("apnuh","apnuf")))

    ## Error status
    status = 0

    ## Sampling loop
    for (i in 1:n) {

        ## Update actual historical climate
        mean    = taua*yh + tauz*muw
        var     = 1/(taua + tauz)
        yha     = rnorm(1, mean*var, sqrt(var))

        ## Update historical structural uncertainty
        shape   = atauh + 0.5*(m + 1)
        rate    = btauh + 0.5*(sum((xhm - muh)^2) + (yh - muh)^2/kdsq)
        tauh    = rgamma(1, shape, rate)

        taudh   = tauh/kdsq

        ## Update historical model climates
        mean    = tauh*muh + tauf*beta*(xfm - muf + beta*muh) + taum*nxhmr*mxhmr
        var     = 1/(tauh + tauf*beta^2 + taum*nxhmr)
        xhm     = rnorm(m, mean*var, sqrt(var))

        ## Update future model climates
        mean    = tauf*(muf + beta*(xhm - muh)) + phim*taum*nxfmr*mxfmr
        var     = 1/(tauf + phim*taum*nxfmr)
        xfm     = rnorm(m, mean*var, sqrt(var))

        ## Update historical model internal variability
        shape   = nuh + nxhmr + nxfmr
        rate    = rep(nuh*psisq, m)
        for (j in 1:m)
            rate[j] = rate[j] + sum((xhmr[[j]] - xhm[j])^2, na.rm = TRUE) +
                phim[j]*sum((xfmr[[j]] - xfm[j])^2, na.rm = TRUE)
        taum    = rgamma(m, 0.5*shape, 0.5*rate)

        ## Update future model internal variability
        shape   = nuf + nxfmr
        rate    = rep(nuf*thetasq, m)
        for (j in 1:m)
            rate[j] = rate[j] + taum[j]*sum((xfmr[[j]] - xfm[j])^2, na.rm = TRUE)
        phim    = rgamma(m, 0.5*shape, 0.5*rate)

        ## Update expected historical climate
        mean    = bmuh*amuh + bmuf*muf + tauh*sum(xhm) -
            tauf*beta*sum(xfm - muf - beta*xhm) + taudh*yh
        var     = 1/(bmuh + bmuf + m*tauh + m*tauf*beta^2 + taudh)
        muh     = rnorm(1, mean*var, sqrt(var))

        ## Update expected future climate
        mean    = bmuf*muh + tauf*sum(xfm - beta*(xhm - muh))
        var     = 1/(bmuf + m*tauf)
        muf     = rnorm(1, mean*var, sqrt(var))

        ## Update emergent constraint
        mean    = bbeta*abeta + tauf*sum((xhm - muh)*(xfm - muf))
        var     = 1/(bbeta + tauf*sum((xhm - muh)^2))
        beta    = rnorm(1, mean*var, sqrt(var))

        ## Update future structural uncertainty
        shape   = atauf + 0.5*m
        rate    = btauf + 0.5*sum((xfm - muf - beta*(xhm - muh))^2)
        tauf    = rgamma(1, shape, rate)

        ## Update historical internal variability
        shape   = apsisq + 0.5*(m*nuh + nuha)
        rate    = bpsisq + 0.5*(nuh*sum(taum) + nuha*taua)
        psisq   = rgamma(1, shape, rate)

        ## Update future internal variability
        shape   = athetasq + 0.5*m*nuf
        rate    = bthetasq + 0.5*nuf*sum(phim)
        thetasq = rgamma(1, shape, rate)

        ## Update reanalysis mean ##
        mean    = tauw*r*mwhr + tauz*yha
        var     = 1/(tauw*r + tauz)
        muw     = rnorm(1, mean*var, sqrt(var))

        ## Update reanalysis uncertainty ##
        shape   = atauw + 0.5*(1 + r)
        rate    = btauw + 0.5*((yha - muw)^2/kwsq + sum((whr - muw)^2))
        tauw    = rgamma(1, shape, rate)

        tauz    = tauw/kwsq

        ## Update historical climate
        mean    = taudh*muh + taua*yha
        var     = 1/(taudh + taua)
        yh      = rnorm(1, mean*var, sqrt(var))

        ## Update natural variability
        shape   = 0.5*(nuha + 1)
        rate    = 0.5*(nuha*psisq + (yha - yh)^2)
        taua    = rgamma(1, shape, rate)

        ## Update historical degrees-of-freedom
        nuhs  = rgamma(1, nuh*lambdah, lambdah)  ## Proposal distribution
        p     = dgamma(nuh , anuh, bnuh, log = TRUE) +
            sum(dgamma(taum, 0.5*nuh , 0.5*nuh *psisq, log = TRUE)) +
            dgamma(taua, 0.5*nuh/knsq, 0.5*nuh *psisq/knsq, log = TRUE)
        ps    = dgamma(nuhs, anuh, bnuh, log = TRUE) +
            sum(dgamma(taum, 0.5*nuhs , 0.5*nuhs*psisq, log = TRUE)) +
            dgamma(taua, 0.5*nuhs/knsq, 0.5*nuhs*psisq/knsq, log = TRUE)
        q     = dgamma(nuh , nuhs*lambdah, lambdah, log = TRUE)
        qs    = dgamma(nuhs, nuh *lambdah, lambdah, log = TRUE)
        apnuh = min(exp(ps - p + q - qs),1)  ## Acceptance probability
        nuh   = ifelse(runif(1) < apnuh, nuhs, nuh)

        nuha  = nuh/knsq

        ## Update future degrees-of-freedom
        nufs  = rgamma(1, nuf*lambdaf, lambdaf)  ## Proposal distribution
        p     = dgamma(nuf , anuf, bnuf, log = TRUE) +
            sum(dgamma(phim, 0.5*nuf , 0.5*nuf *thetasq, log = TRUE))
        ps    = dgamma(nufs, anuf, bnuf, log = TRUE) +
            sum(dgamma(phim, 0.5*nufs, 0.5*nufs*thetasq, log = TRUE))
        q     = dgamma(nuf , nufs*lambdaf, lambdaf, log = TRUE)
        qs    = dgamma(nufs, nuf *lambdaf, lambdaf, log = TRUE)
        apnuf = min(exp(ps - p + q - qs),1)  ## Acceptance probability
        nuf   = ifelse(runif(1) < apnuf, nufs, nuf)

        ## Store samples
        samples[i,"muh"]     = muh      ## Historical expected climate
        samples[i,"muf"]     = muf      ## Future     expected climate
        samples[i,"beta"]    = beta     ## Emergent constraint
        samples[i,"tauh"]    = tauh     ## Historical structural uncertainty
        samples[i,"tauf"]    = tauf     ## Future     structural uncertainty
        samples[i,"psisq"]   = psisq    ## Historical internal variability
        samples[i,"thetasq"] = thetasq  ## Future     internal varibility
        samples[i,"nuh"]     = nuh      ## Historical degrees-of-freedom
        samples[i,"nuf"]     = nuf      ## Future     degrees-of-freedom
        samples[i,"yha"]     = yha      ## Historical actualised climate
        samples[i,"yh"]      = yh       ## Historical climate
        samples[i,"taua"]    = taua     ## Historical natural variability
        samples[i,ixhm]      = xhm      ## Historical model climates
        samples[i,ixfm]      = xfm      ## Future     model climates
        samples[i,itaum]     = taum     ## Historical model variability
        samples[i,iphim]     = phim     ## Future     model variability
        samples[i,"muw"]     = muw      ## Reanalysis expectation
        samples[i,"tauw"]    = tauw     ## Reanalysis variance

        ## Store acceptance probabilities
        acceptance[i,"apnuh"] = apnuh   ## Historical acceptance probability
        acceptance[i,"apnuf"] = apnuf   ## Future     acceptance probability

        ## Check for numerical errors
        if (any(is.na(samples[i,]))) {
            status = 1
            break
        }

    } ## Sampling loop

    ## Return samples
    return(list(samples = samples, acceptance = acceptance, status = status))

}
gibbs.sampler = cmpfun(gibbs.sampler)
