############################
## Supplementary Material ##
############################

## Load libraries
library(RColorBrewer)

## Load source
source("contourplot.R")

## Load MCMC output
load("samples.RData")

## Extract lon and lat
lon = sort(unique(cmip5$lon))
lat = sort(unique(cmip5$lat))

## Plotting options
white = '#FFFFFF'
x = lon
y = lat
xat = seq(-180,+180,+30)
yat = seq(+ 45,+ 90,+15)
proj4string =
    CRS('+proj=laea +lat_0=+90 +lon_0=+0 +ellps=GRS80 +units=m +no_defs')
width  = 137/25.4
pointsize = 10


## Historical consensus
muh      = sapply(samples, function(x) x$samples[,"muh"])
muh      = apply(muh, 2, mean)
dim(muh) = c(length(lon),length(lat))

zat = seq(-40,+15,+5)
col = c(rev(brewer.pal(8,"Blues")),brewer.pal(8, "Reds")[1:3])

contour.plot(x, y, muh, xat, yat, zat,
             col, proj4string, width, pointsize,
             draw.box = FALSE,
             ztitle = expression(paste("Representative historical climate ",mu[H],
                                       " (",phantom()^degree,"C)")),
             key.position = "bottom",
             dev = "pdf", file = "figures1.pdf")

## Future consensus
muf      = sapply(samples, function(x) x$samples[,"muf"])
muf      = apply(muf, 2, mean)
dim(muf) = c(length(lon),length(lat))

contour.plot(x, y, muf, xat, yat, zat,
             col, proj4string, width, pointsize,
             draw.box = FALSE,
             ztitle = expression(paste("Representative future climate ",mu[F],
                                       " (",phantom()^degree,"C)")),
             key.position = "bottom",
             dev = "pdf", file = "figures2.pdf")

## Reanalysis consensus
muw      = sapply(samples, function(x) x$samples[,"muw"])
muw      = apply(muw, 2, mean)
dim(muw) = c(length(lon),length(lat))

contour.plot(x, y, muw, xat, yat, zat,
             col, proj4string, width, pointsize,
             draw.box = FALSE,
             ztitle = expression(paste("Representative reanalysis ",mu[W],
                                       " (",phantom()^degree,"C)")),
             key.position = "bottom",
             dev = "pdf", file = "figures3.pdf")


## Historical model spread
sigmah      = sapply(samples, function(x) sqrt(1/x$samples[,"tauh"]))
sigmah      = apply(sigmah, 2, mean)
dim(sigmah) = c(length(lon),length(lat))

col = c("white",brewer.pal(9,"Reds"))
zat = seq(0,10,1)

contour.plot(x, y, sigmah, xat, yat, zat,
             col, proj4string, width, pointsize, draw.box = FALSE,
             ztitle = expression(paste("Historical model uncertainty ",sigma[H],
                                       " (",phantom()^degree,"C)")),
             key.position = "bottom",
             dev = "pdf", file = "figures4.pdf")

## Response model spread
sigmaf      = sapply(samples, function(x) sqrt(1/x$samples[,"tauf"]))
sigmaf      = apply(sigmaf, 2, mean)
dim(sigmaf) = c(length(lon),length(lat))

col = c("white",brewer.pal(6,"Reds"))
zat = seq(0,7,1)

contour.plot(x, y, sigmaf, xat, yat, zat,
             col, proj4string, width, pointsize, draw.box = FALSE,
             ztitle = expression(paste("Response model uncertainty ",sigma["F|H"],
                                       " (",phantom()^degree,"C)")),
             key.position = "bottom",
             dev = "pdf", file = "figures5.pdf")

## Reanalysis spread
sigmaw      = sapply(samples, function(x) 1/sqrt(x$samples[,"tauw"]))
sigmaw      = apply(sigmaw, 2, mean)
dim(sigmaw) = c(length(lon),length(lat))

col = c("white",brewer.pal(9,"Reds"))
sigmaw[sigmaw > 4.5] = 4.99
zat = seq(0.0,5.0,0.5)
zlabels = format(c(seq(0.0,4.5,0.5),7.2))

contour.plot(x, y, sigmaw, xat, yat, zat, col, proj4string, width, pointsize,
             draw.box = FALSE, zlabels = zlabels,
             ztitle = expression(paste("Reanalysis uncertainty ",sigma[W],
                                       " (",phantom()^degree,"C)")),
             key.position = "bottom",
             dev = "pdf", file = "figures6.pdf")

## Representative historical variability
psi      = sapply(samples, function(x) sqrt(x$samples[,"psisq"]))
psi      = apply(psi, 2, mean)
dim(psi) = c(length(lon),length(lat))

col = c("white",brewer.pal(9,"Reds"))
psi[psi > 0.9] = 0.999
zat = seq(0.0,1.0,0.1)
zlabels = c(seq(0.0,0.9,0.1),1.2)

contour.plot(x, y, psi, xat, yat, zat,
             col, proj4string, width, pointsize,
             draw.box = FALSE, zlabels = zlabels,
             ztitle = expression(paste("Representative internal variability ",psi,
                                       " (",phantom()^degree,"C)")),
             key.position = "bottom",
             dev = "pdf", file = "figures7.pdf")

## Representative future variability
theta      = sapply(samples, function(x) sqrt(x$samples[,"thetasq"]))
theta      = apply(theta, 2, mean)
dim(theta) = c(length(lon),length(lat))

col = c(rev(brewer.pal(6,"Blues")[1:4]),"white","white",brewer.pal(6,"Reds"))
zat = seq(0.0,2.4,0.2)

contour.plot(x, y, theta, xat, yat, zat,
             col, proj4string, width, pointsize, draw.box = FALSE,
             ztitle = expression(paste("Representative future variability ratio ",theta)),
             key.position = "bottom",
             dev = "pdf", file = "figures8.pdf")
