## Load libraries
library(sp)
library(maptools)
library(rgdal)
library(rgeos)
library(raster)

## Function to make contour plot, overlay map and write to file
contour.plot = function(x, y, z, xat, yat, zat, col, proj4string = NULL,
    width = 7, pointsize = 12, resolution = 2.5, key.position = 'right',
    draw.box = TRUE, ztitle = NULL, grid.labels = TRUE, title.col = 1,
    label = NULL, label.pos = 0,
    xlabels = format(xat, trim = TRUE), ylabels = format(yat), title = NULL,
    zlabels = format(zat), dev = c('x11','eps','pdf'), file = 'Rplots') {

    ############################
    ## Prepare background map ##
    ############################
    
    ## Load world map
    data(wrld_simpl)

    ## Trim world map

    ## This creates a curves on the parallels, if we use a simple 4 cornered
    ## box then the top and bottom become straight lines not curves when
    ## transformed
    parallel = seq(from = min(xat), to = max(xat), by = diff(range(xat))/360)
    xcoords = c(min(xat), parallel, rev(parallel))
    ycoords = c(min(yat), rep(max(yat), times = length(parallel)),
        rep(min(yat), times = length(parallel)))
    coords = cbind(xcoords,ycoords)
    colnames(coords) = c('x','y')
    extent = Polygon(coords = coords, hole = FALSE)
    extent = Polygons(srl = list(extent), ID = 'extent')
    extent = SpatialPolygons(Srl = list(extent), pO = 1:1,
        proj4string = attr(x = wrld_simpl, which = 'proj4string'))
    world.map = gIntersection(wrld_simpl, extent, byid = TRUE, 
                              checkValidity = 2L)


    #######################
    ## Prepare gridlines ##
    #######################

    ## Make meridians
    meridian.lines = list()
    for(i in 1:length(xat)) {
        xmer = rep(xat[i], times = 2)
        ymer = range(yat)
        meridian.lines[[i]] = Line(cbind(x = xmer, y = ymer))
    }
    meridian.lines = Lines(meridian.lines, ID = 'meridians')

    ## Make meridians
    parallel.lines = list()
    for(i in 1:length(yat)) {
        xpar = seq(from = min(xat), to = max(xat), length.out = 360)
        ypar = rep(yat[i], times = 360)
        parallel.lines[[i]] = Line(cbind(x = xpar, y = ypar))
    }
    parallel.lines = Lines(parallel.lines, ID = 'parallels')

    ## Make gridlines
    grid.lines = SpatialLines(list(meridian.lines, parallel.lines),
        proj4string = CRS("+proj=longlat +datum=WGS84"))
   

    ####################
    ## Prepare raster ##
    ####################
        
    ## Create raster
    xres = x[2] - x[1]
    yres = y[2] - y[1]
    z = raster(x = t(z)[length(y):1,],
        xmn = min(x) - xres/2, xmx = max(x) + xres/2,
        ymn = min(y) - yres/2, ymx = max(y) + yres/2,
        crs = CRS("+proj=longlat +datum=WGS84"))

    ## Resample to reduce holes when joining ends
    ## Enlarging sampling area fixes rectangular and segment plots
    ## but buggers up hemispheres.
    template = raster(nrows = diff(range(yat))/resolution,
        ncols = diff(range(xat))/resolution,
#        xmn = min(xat) - resolution, xmx = max(xat) + resolution,
        xmn = min(xat), xmx = max(xat),
#        ymn = min(yat) - resolution, ymx = min(90,max(yat) + resolution))
        ymn = min(yat), ymx = max(yat))
    z = resample(z,template)
    
#    ## Trim raster
#    extent = extent(min(xat), max(xat), min(yat), max(yat))
#    z      = crop(z, extent)

#    ## Experimental
#    extent = extent(min(xat) - xres, max(xat) + xres, min(yat) - yres,
#        max(yat) + yres)
#    z = crop(z, extent)

    
    #################################
    ## Transform plotting elements ##
    #################################
    
    if(!is.null(proj4string)) {
        world.map = spTransform(world.map, proj4string)
        grid.lines = spTransform(grid.lines, proj4string)
        z = projectRaster(z, crs = proj4string, over = FALSE)
        z = trim(z, padding = 1)
    }

    
    ## This section follows the code in raster::filledContour which is a
    ## wrapper for the standard plotting function filled.contour
    
    ## Resample raster for contouring
#    z = sampleRegular(z, size = 1e+05, asRaster = TRUE, useGDAL = TRUE)

    ## Extract extent
    x = xFromCol(z, 1:ncol(z))
    y = yFromRow(z, nrow(z):1)

    ## Transform raster to matrix for plotting
    z = t(matrix(data = getValues(z),
        ncol = ncol(z), byrow = TRUE)[nrow(z):1,])

    ## Experimental
    ## Need to extract 
    
    
    #################################
    ## Set up plotting environment ##
    #################################

    ## Set output device options
    ps.options(onefile = FALSE, horizontal = FALSE, paper = 'special')
    pdf.options(onefile = FALSE)
    
    ## Calculate line height and character width
    line.height = 1.2 * pointsize / 72
    char.width  = line.height * 3 / 4
    
    ## Calculate max label widths
    dev = match.arg(dev)
    switch(dev,
           x11 = x11(pointsize = pointsize),
           eps = postscript(width = 7, height = 7, pointsize = pointsize),
           pdf = pdf(pointsize = pointsize))
    xwidth = max(strwidth(s = xlabels, units = 'inches'))/line.height
    ywidth = max(strwidth(s = ylabels, units = 'inches'))/line.height
    zwidth = max(strwidth(s = zlabels, units = 'inches'))/line.height
    dev.off()

    if(key.position == 'right') {
        ## Calculate margins
        top.mar = ifelse(test = is.null(title), yes = 1, no = 2)
        map.mar = c(1,1,top.mar,1)
        ## If no projection then leave room for labels
        if(is.null(proj4string) & grid.labels) {
            map.mar = map.mar + c(1,ywidth,0,0)
        }
        key.mar = c(map.mar[1], 0, map.mar[3], 1+zwidth)
        if(!is.null(ztitle)) key.mar[4] = key.mar[4] + 1
        map.mar = map.mar + 0.1
        key.mar = key.mar + 0.1
        map.mai = map.mar * line.height
        key.mai = key.mar * line.height
        
        ## Width of map and key regions
        key.width = 1 ## Set key width in lines
        key.width = key.width * line.height + sum(key.mai[c(2,4)])
        map.width = width - key.width
        
        ## Height of page
        height = map.width - sum(map.mai[c(2,4)])
        ## If the map is not projected then use xat to fit map exactly
        if(is.null(proj4string)) {
            height = height * diff(range(yat)) / diff(range(xat))
        } else {
            height = height * diff(range(y)) / diff(range(x))
        }
        height = height + sum(map.mai[c(1,3)])
    }
    if(key.position == 'bottom') {
        ## Calculate margins
        top.mar = ifelse(test = is.null(title), yes = 1, no = 2)
        map.mar = c(1,1,top.mar,1)
        ## If no projection then leave room for labels
        if(is.null(proj4string) & grid.labels) {
            map.mar = map.mar + c(1,ywidth,0,0)
        }
        key.mar = c(2, map.mar[1], 0, 1)     #
        if(!is.null(ztitle)) key.mar[1] = key.mar[1] + 1
        map.mar = map.mar + 0.1
        key.mar = key.mar + 0.1
        map.mai = map.mar * line.height
        key.mai = key.mar * line.height
        
        ## Height of key region
        key.height = 1 ## Set key height in lines #
        key.height = key.height * line.height + sum(key.mai[c(1,3)]) #
##        map.height =  - key.width #
       
        ## Height of map section and page
        map.height = width - sum(map.mai[c(2,4)]) #
        ## If the map is not projected then use xat to fit map exactly
        if(is.null(proj4string)) {
            map.height = map.height * diff(range(yat)) / diff(range(xat)) #
        } else {
            map.height = map.height * diff(range(y)) / diff(range(x)) #
        }
        map.height = map.height + sum(map.mai[c(1,3)])# 
        height = map.height + key.height #
    }
    if(key.position == 'none') {
        ## Calculate margins
        top.mar = ifelse(test = is.null(title), yes = 1, no = 2)
        map.mar = c(1,1,top.mar,1)
        ## If no projection then leave room for labels
        if(is.null(proj4string) & grid.labels) {
            map.mar = map.mar + c(1,ywidth,0,0)
        }
        map.mar = map.mar + 0.1
        map.mai = map.mar * line.height
        
        ## Height of map section and page
        map.height = width - sum(map.mai[c(2,4)]) #
        ## If the map is not projected then use xat to fit map exactly
        if(is.null(proj4string)) {
            map.height = map.height * diff(range(yat)) / diff(range(xat)) #
        } else {
            map.height = map.height * diff(range(y)) / diff(range(x)) #
        }
        map.height = map.height + sum(map.mai[c(1,3)])# 
        height = map.height
    }


    ###############
    ## Plot data ##
    ###############
    
    ## Open plotting device
    switch(dev,
           x11 = x11(width = width, height = height, pointsize = pointsize),
           eps = postscript(file = file, width = width, height = height,
               pointsize = pointsize),#, title = title),
           pdf = pdf(file = file, width = width, height = height,
               pointsize = pointsize))#, title = title))
    if(key.position == 'right') {
        layout(matrix(1:2,1,2), widths = c(map.width,key.width))
    }
    if(key.position == 'bottom') {
        layout(matrix(1:2,2,1), heights = c(map.height,key.height))
    }
    
    ## Plot contour map  
    par(cex = 1, mar = map.mar, xaxs = 'i', yaxs = 'i', las = 1)
    plot.new()
    if(is.null(proj4string)) {
        plot.window(xlim = range(xat), ylim = range(yat), asp = 1)
    } else {
        plot.window(xlim = range(x), ylim = range(y), asp = 1)
    }
    .filled.contour(x = x, y = y, z = z, levels = zat, col = col)
    if(is.null(proj4string)) {
        if(grid.labels) {
            axis(side = 1, at = xat, labels = xlabels, hadj = 0.5)
            axis(side = 2, at = yat, labels = ylabels, hadj = 1)
        } else {
            axis(side = 1, at = xat, labels = FALSE)
            axis(side = 2, at = yat, labels = FALSE)
        }
    }
    plot(world.map, add = TRUE)
    plot(grid.lines, lty = 'dashed', add = TRUE)
    if(draw.box) box()
    if(!is.null(title)) title(main = title, line = 1, col.main = title.col)
    if(!is.null(label)) mtext(text = label, side = 3, line = 0, adj = label.pos)
    
    ## Add key
    if(key.position == 'right') {
        par(cex = 1, mar = key.mar, xaxs = 'i', yaxs = 'i', las = 1,
            mgp = c(3,1+zwidth,0))
        plot.new()
        plot.window(xlim = c(0,1), ylim = range(zat))
        rect(xleft = 0, ybottom = zat[-length(zat)] ,
             xright = 1 , ytop = zat[-1], col = col)
        axis(side = 4, at = zat, labels = zlabels, hadj = 1)
        if(!is.null(ztitle))
            mtext(text = ztitle, side = 4, line = key.mar[4] - 1.1, las = 3)
    }
    if(key.position == 'bottom') {
        par(cex = 1, mar = key.mar, xaxs = 'i', yaxs = 'i', las = 1,
            mgp = c(3,1,0))
        plot.new()
        plot.window(xlim = range(zat), ylim = c(0,1))
        rect(xleft = zat[-length(zat)], ybottom= 0,
             xright = zat[-1], ytop = 1, col = col)
        axis(side = 1, at = zat, labels = zlabels)#, hadj = 1)
        if(!is.null(ztitle))
            mtext(text = ztitle, side = 1, line = key.mar[1] - 1.1)
    }
    
    ## Close device if writing to file
    if(dev != 'x11') dev.off()
}
