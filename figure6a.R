###############
## Figure 6a ##
###############

## Load MCMC output - no emergent constraint
load("samples-nobeta.RData")

## Extract useful quantities
lon      = sort(unique(cmip5$lon))
lat      = sort(unique(cmip5$lat))
n.lon    = length(lon)
n.lat    = length(lat)

## Create grid
grid = expand.grid(lon = lon, lat = lat)

## Select grid box - Svalbard
x = 38.75
y = 81.25
j = which(grid$lon == x & grid$lat == y)

## Extract data
samples0 = samples[[j]]$samples
rm(samples); gc()

## Load MCMC output - with emergent constraint
load("samples.RData")
obs     = reanalysis$tas[reanalysis$lon == x & reanalysis$lat == y]
cmip5   = cmip5[cmip5$lon == x & cmip5$lat == y, -c(1,2)]
samples = samples[[j]]$samples
gc()

## Quoted statistics
round(mean(samples[,"muf"] - samples[,"muh"]), 1)
round(quantile(samples[,"muf"] - samples[,"muh"], probs = c(0.05,0.95)), 1)
round(mean(samples[,"beta"]), 1)
round(quantile(samples[,"beta"], probs = c(0.05,0.95)), 1)
round(mean(samples[,"yh"] - samples[,"muh"]), 1)
round(quantile(samples[,"yh"] - samples[,"muh"], probs = c(0.05,0.95)), 1)
round(mean((samples[,"beta"] - 1)*(samples[,"yh"] - samples[,"muh"])), 1)
round(quantile((samples[,"beta"] - 1)*(samples[,"yh"] - samples[,"muh"]),
               probs = c(0.05,0.95)), 1)

## Extract data
models   = levels(cmip5$model)
n.models = nlevels(cmip5$model)
means = with(cmip5, tapply(tas, list(model,experiment), mean))
maxs  = with(cmip5, tapply(tas, list(model,experiment), max))
mins  = with(cmip5, tapply(tas, list(model,experiment), min))
means = as.data.frame(means)
maxs  = as.data.frame(maxs)
mins  = as.data.frame(mins)
means$response = means$rcp45 - means$historical
mins$response = mins$rcp45 - maxs$historical
maxs$response = maxs$rcp45 - mins$historical

## Compute plotting limits
xlim = c(floor(min(mins$historical)),ceiling(max(maxs$historical)))
#ylim = c(floor(min(mins$response)),ceiling(max(maxs$response)))
xx   = seq(xlim[1],xlim[2],0.5)
xat  = pretty(means$historical)
xlab = xat

## Sample projections along xx points
yr = array(NA, c(1e3,length(xx)))
for (i in 1:length(xx))
    yr[,i] = samples[,"muf"] - xx[i] +
      samples[,"beta"]*(xx[i] - samples[,"muh"]) +
      rnorm(1e3, 0, sqrt(1/samples[,"tauf"] +
                             1/(samples[,"phia"]*samples[,"taua"])))

#fit = apply(yr, 2, mean)
fit = mean(samples[,"muf"]) - xx +
    mean(samples[,"beta"])*(xx - mean(samples[,"muh"]))
lwr = apply(yr, 2, quantile, probs = 0.05)
upr = apply(yr, 2, quantile, probs = 0.95)

ylim = c(floor(min(lwr)),ceiling(max(upr)))

## Plotting options
width     = 6.5
pointsize = 10

## Open PDF device
pdf(file = "figure6a.pdf", width = 6, height = 3.7)

layout(matrix(c(1,2), 1, 2), c(0.8,0.2))

par(mar=c(3,3,1,1)+0.1, mgp = c(2,1,0), las = 1, cex = 1.0,
    xaxs = 'i', yaxs = 'i')

## Density multiplier
rho = 5

## Plot data
plot(means$historical, means$response,
     xlab = expression(paste("Historical climate (",phantom()^degree,"C)")),
     ylab = expression(paste("Climate response (",phantom()^degree,"C)")),
     main = "",
     xlim = xlim, ylim = ylim, pch = 3, xaxt = "n", yaxt = "n")

axis(1, at = xat, labels = xlab)
Axis(side = 2)

## Add ranges
for (i in 1:n.models) {

    if(maxs$historical[i] > mins$historical[i])
        arrows(x0 = mins$historical[i],
               x1 = maxs$historical[i],
               y0 = means$response[i],
               length = 0.025, angle = 90, code = 3)

    current = cmip5[cmip5$model == models[i],]
    responses = expand.grid(current$tas[current$experiment == 'historical'],
                            current$tas[current$experiment == 'rcp45'])
    responses = responses$Var2 - responses$Var1
    if(max(responses) > min(responses))
        arrows(x0 = means$historical[i],
               y0 = min(responses),
               y1 = max(responses),
               length = 0.025, angle = 90, code = 3)

} ## i

## Add regression
lines(xx, fit, col = "red", lwd = 2)
lines(xx, lwr, col = "red", lty = "dotted", lwd = 2)
lines(xx, upr, col = "red", lty = "dotted", lwd = 2)

# ## Add ensemble mean projection
# muh = mean(means$historical)
# muf = mean(means$rcp45)
# segments(x0 = muh, y0 = ylim[1], y1 = muf - muh      , col = "red", lty = "dashed")
# segments(x0 = xlim[1], x1 = muh, y0 = muf - muh, col = "red", lty = "dashed")

## Add observations
z = rnorm(1e6,samples[,"muw"],1/sqrt(samples[,"tauw"]))
zd = density(z)
lines(zd$x, ylim[1] + rho*zd$y, col = "blue", lwd = 2)
muz = mean(z)
beta = samples[,"muf"] - muz +
  samples[,"beta"]*(muz - samples[,"muh"]) +
    rnorm(1e3,0,sqrt(1/samples[,"tauf"] + 1/(samples[,"phia"]*samples[,"taua"])))
beta = mean(beta)
segments(x0 = muz, y0 = ylim[1], y1 = beta,
         col = "blue", lty = "dashed", lwd = 2)
segments(x0 = xlim[1], x1 = muz, y0 = beta,
         col = "blue", lty = "dashed", lwd = 2)

## Add densities
yh = samples[,"yha"]
yhd = density(yh)
# lines(yhd$x, ylim[1] + rho*yhd$y, col = "red")
yr  = samples[,"yfa"] - samples[,"yha"]
yrd = density(yr)
# lines(xlim[1] + rho*yrd$y, yrd$x, col = "red")
segments(x0 = mean(yh), y0 = ylim[1], y1 = mean(yr),
         col = "red", lty = "dashed", lwd = 2)
segments(x0 = xlim[1], x1 = mean(yh), y0 = mean(yr),
         col = "red", lty = "dashed", lwd = 2)

# ## Add densities
# yh0 = samples0$yha
# yhd0 = density(yh)
# # lines(yhd$x, ylim[1] + rho*yhd$y, col = "black")
yr0  = samples0[,"yfa"] - samples0[,"yha"]
yrd0 = density(yr0)
# # lines(xlim[1] + rho*yrd$y, yrd$x, col = "black")
# segments(x0 = mean(yh), y0 = ylim[1], y1 = mean(yr), col = "black", lty = "dashed")
# segments(x0 = xlim[1], x1 = mean(yh), y0 = mean(yr), col = "black", lty = "dashed")

## Add densities
muh  = samples[,"muh"]
muhd = density(muh)
# lines(yhd$x, ylim[1] + rho*yhd$y, col = "black")
mur  = samples[,"muf"] - samples[,"muh"]
murd = density(mur)
# lines(xlim[1] + rho*yrd$y, yrd$x, col = "black")
segments(x0 = mean(muh), y0 = ylim[1], y1 = mean(mur),
         col = "black", lty = "dashed", lwd = 2)
segments(x0 = xlim[1], x1 = mean(muh), y0 = mean(mur),
         col = "black", lty = "dashed", lwd = 2)

# abline(v = xlim[1],col=2)
# abline(v = xlim[2],col=2)

mtext("(a)", side = 3, line = 0, adj = -0.15)

par(mar=c(3.1,0,1.1,1.1), mgp = c(2,1,0), las = 1, cex = 1.0,
    xaxs = 'i', yaxs = 'i')

## Plot data
plot(yrd$y, yrd$x, type = "l", col = "red", xlim = c(0,max(yrd$y,yrd0$y)),
     ylim = ylim, xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "")
Axis(side = 2, labels = FALSE)
lines(yrd0$y, yrd0$x, col = "black")

## Close plotting device
dev.off()
