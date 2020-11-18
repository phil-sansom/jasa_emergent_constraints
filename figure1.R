#############################
## Figure 1 - Introduction ##
#############################

## Load data
load("nh_djf_rcp45_tas.RData")

## Convert to Celsius from Kelvin
cmip5$tas      = cmip5$tas
reanalysis$tas = reanalysis$tas

## Extract lon and lat
lon = sort(unique(cmip5$lon))
lat = sort(unique(cmip5$lat))

## Select point for plotting
x = lon[28]; y = lat[15]

## Extract data for plotting
z      = mean(reanalysis$tas[reanalysis$lon == x & reanalysis$lat == y])
buffer = cmip5[cmip5$lon == x & cmip5$lat == y,]
xs     = with(buffer, tapply(tas, list(model, experiment), mean))
vs     = with(buffer, tapply(tas, list(model, experiment), var ))
xs[,"rcp45"] = xs[,"rcp45"] - xs[,"historical"]
vs[,"rcp45"] = vs[,"rcp45"] + vs[,"historical"]

zz = z
xx = xs[,"historical"]
yy = xs[,"rcp45"]
sh = sqrt(vs[,"historical"])
sr = sqrt(vs[,"rcp45"])

## Approximate range of internal variability for each model
xmin = xx + qnorm(0.025)*sh
xmax = xx + qnorm(0.975)*sh
ymin = yy + qnorm(0.025)*sr
ymax = yy + qnorm(0.975)*sr

## Compute plotting parameters
xlim = c(floor(min(xmin,na.rm=TRUE)),ceiling(max(xmax,na.rm=TRUE)))
ylim = c(floor(min(ymin,na.rm=TRUE)),ceiling(max(ymax,na.rm=TRUE)))
xat = seq(xlim[1],xlim[2],2)
yat = seq(ylim[1],ylim[2],2)

width  = 6.5*2/3
height = width/1.30

## Open PDF device
pdf(file = "figure1.pdf", width = width, height = height)

par(ann = FALSE, las = 1, mar = c(2,2,0.25,0.5)+0.1,
    ps = 10, tcl = -1/3, xaxs = "i", yaxs = "i")

## Plot simulator outputs
plot(xx, yy, pch = 3, xlim = xlim, ylim = ylim, xaxt = "n", yaxt = "n", asp = 1)
title(xlab = expression(paste("Historical temperature (",degree,"C)")), line = 1)
title(ylab = expression(paste("Temperature response (",degree,"C)")), line = 1)
axis (side = 1, at = xat, mgp = c(2.5,1/4,0))
axis (side = 2, at = yat, mgp = c(2.5,1/2,0))

## Add ensemble mean
segments(x0 = mean(xx), y0 = ylim[1], y1 = mean(yy), lty = "dashed", lwd = 2)
segments(x0 = xlim[1], x1 = mean(xx), y0 = mean(yy), lty = "dashed", lwd = 2)

## Add regression line and projection
model = lm(yy ~ xx)
abline(model, lwd = 2)
zr = predict(model, data.frame(xx = zz))
segments(x0 = zz, y0 = ylim[1], y1 = zr, lty = "dotted", lwd = 2)
segments(x0 = xlim[1], x1 = zz, y0 = zr, lty = "dotted", lwd = 2)

## Add ranges
for (i in 1:length(xx)) {

    if(!is.na(sh[i]) & !is.na(sr[i])) {

        arrows(x0 = xmin[i],
               x1 = xmax[i],
               y0 = yy  [i],
               length = 0.025, angle = 90, code = 3)

        arrows(x0 = xx  [i],
               y0 = ymin[i],
               y1 = ymax[i],
               length = 0.025, angle = 90, code = 3)
    }

}

## Close plotting device
dev.off()
