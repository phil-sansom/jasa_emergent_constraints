###################
## Master script ##
###################

## Run sampler with emergent constraint and save output
## Package dependencies: coda, compiler, parallel
source("samples.R")

## Run sampler without emergent constraint and save output
## Package dependencies: coda, compiler, parallel
source("samples-nobeta.R")

## Perform model checking in Section 4.2
## Package dependencies: mcmcse
source("model-checking.R")

## Make figures 1,4,5 & 6
## Package dependencies: sp,maptools,rgdal,rgeos,raster
source("figure1.R")
source("figure4.R")
source("figure5.R")
source("figure6a.R")
source("figure6b.R")

## Make figures for supplementary material
## Package dependencies: sp,maptools,rgdal,rgeos,raster
source("supplementary-material.R")
