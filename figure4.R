############################
## Figure 4 - Projections ##
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

## Extract historical & future climate
yh = sapply(samples, function(x) x$samples[,"yh"  ])
yf = sapply(samples, function(x) x$samples[,"yf"  ])
yr = yf - yh
yh = apply(yh, 2, mean)
yf = apply(yf, 2, mean)
yr = apply(yr, 2, mean)
dim(yh) = c(length(lon),length(lat))
dim(yf) = c(length(lon),length(lat))
dim(yr) = c(length(lon),length(lat))

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

## Historical climate
zat = seq(-40,+15,+5)
col = c(rev(brewer.pal(8,"Blues")),brewer.pal(8, "Reds")[1:3])

contour.plot(x, y, yh, xat, yat, zat,
             col, proj4string, width, pointsize,
             draw.box = FALSE,
             ztitle = expression(paste("Historical climate ",Y[H],
                                       " (",phantom()^degree,"C)")),
             key.position = "bottom",
             label = "(a)", label.pos = -0.03,
             dev = "pdf", file = "figure4a.pdf")

## Future climate
zat = seq(-40,+15,+5)
col = c(rev(brewer.pal(8,"Blues")),brewer.pal(8, "Reds")[1:3])

contour.plot(x, y, yf, xat, yat, zat,
             col, proj4string, width, pointsize,
             draw.box = FALSE,
             ztitle = expression(paste("Future climate ",Y[F],
                                       " (",phantom()^degree,"C)")),
             key.position = "bottom",
             label = "(b)", label.pos = -0.03,
             dev = "pdf", file = "figure4b.pdf")

## Climate response
zat = c(0,2,4,6,8,10,12,14)
col = c("white",brewer.pal(6, "Reds"))
# zat = c(-4,-2,0,2,4,6,8,10,12,14)
# col = c(rev(brewer.pal(6,"Blues")[1]),"white","white",brewer.pal(6, "Reds"))

contour.plot(x, y, yr, xat, yat, zat,
             col, proj4string, width, pointsize,
             draw.box = FALSE,
             ztitle = expression(paste("Climate response ",Y[F] - Y[H],
                                       " (",phantom()^degree,"C)")),
             key.position = "bottom",
             label = "(c)", label.pos = -0.03,
             dev = "pdf", file = "figure4c.pdf")
