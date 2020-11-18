##############
## Figure 5 ##
##############

## Load libraries
library(RColorBrewer)

## Load source
source("contourplot.R")

## Load MCMC output - no emergent constraint
load("samples-nobeta.RData")

## Extract lon and lat
lon = sort(unique(cmip5$lon))
lat = sort(unique(cmip5$lat))
n.lon = length(lon)
n.lat = length(lat)

## Extract projections
yrs        = sapply(samples, function(x) x$samples[,"yf"] - x$samples[,"yh"])
yr0        = apply(yrs, 2, mean)
dim(yr0)   = c(n.lon,n.lat)
yr0.v      = apply(yrs, 2, var)
dim(yr0.v) = c(n.lon,n.lat)
rm(samples); gc()

## Load MCMC output - with emergent constraint
load("samples.RData")

## Extract projections
yrs         = sapply(samples, function(x) x$samples[,"yf"] - x$samples[,"yh"])
yr          = apply(yrs, 2, mean)
dim(yr)     = c(n.lon,n.lat)
yr.v        = apply(yrs, 2, var)
dim(yr.v)   = c(n.lon,n.lat)
deltah      = sapply(samples, function(x) x$samples[,"yh"] - x$samples[,"muh"])
deltah      = apply(deltah, 2, mean)
dim(deltah) = c(n.lon,n.lat)
beta        = sapply(samples, function(x) x$samples[,"beta"])
beta        = apply(beta, 2, mean)
dim(beta)   = c(n.lon,n.lat)
rm(samples); gc()

## Plotting options
white = '#FFFFFF'
x = lon
y = lat
xat = seq(-180,+180,+30)
yat = seq(+ 45,+ 90,+15)
proj4string =
    CRS('+proj=laea +lat_0=+90 +lon_0=+0 +ellps=GRS80 +units=m +no_defs')
width  = 6.5
pointsize = 15

## Historical discrepancy
zat = seq(-3,+10,+1)
col = c(rev(brewer.pal(9,"Blues")[1:2]),"white","white",brewer.pal(9, "Reds"))

contour.plot(x, y, deltah, xat, yat, zat,
             col, proj4string, width, pointsize,
             draw.box = FALSE,
             ztitle = expression(paste("Historical discrepancy ",Y[H] - mu[H],
                                       " (",phantom()^degree,"C)")),
             key.position = "right",
             label = "(a)", label.pos = -0.05,
             dev = "pdf", file = "figure5a.pdf")

## Emergent constraint
zat = seq(0.3,1.5,0.1)
col = c(rev(brewer.pal(6,"Blues")),"white","white",brewer.pal(6, "Reds")[1:4])

contour.plot(x, y, beta, xat, yat, zat,
             col, proj4string, width, pointsize,
             draw.box = FALSE,
             ztitle = expression(paste("Emergent constraint ",beta)),
             key.position = "right",
             label = "(b)", label.pos = -0.05,
             dev = "pdf", file = "figure5b.pdf")

## Response mean difference
zat = seq(-3.0,+2.5,0.5)
col = c(rev(brewer.pal(5, "Blues")),"white","white",brewer.pal(5, "Reds")[1:4])

contour.plot(x, y, yr-yr0, xat, yat, zat,
             col, proj4string, width, pointsize,
             draw.box = FALSE,
             ztitle = expression(paste("Difference in mean of response ",Y[F]-Y[H],
                                       " (",phantom()^degree,"C)")),
             key.position = "right",
             label = "(c)", label.pos = -0.05,
             dev = "pdf", file = "figure5c.pdf")

## Response var ratio
zat = seq(0.4,1.3,0.1)
col = c(rev(brewer.pal(5, "Blues")),"white","white",brewer.pal(5, "Reds")[1:2])

contour.plot(x, y, sqrt(yr.v/yr0.v), xat, yat, zat,
             col, proj4string, width, pointsize,
             draw.box = FALSE,
             ztitle = expression(paste("Ratio of SD of response ",Y[F]-Y[H])),
             key.position = "right",
             label = "(d)", label.pos = -0.05,
             dev = "pdf", file = "figure5d.pdf")
