library(raster)
library(sp)
library(rgdal)
library(RANN)

setwd("/home/nronnei/gis/class/spatial_analysis/final_project/r/")
raw_dem <- readGDAL('../data/gdem_wgs84_raw_clipped.tif')
clipper <- readOGR(dsn="../data/clipper.shp", "clipper")

dem <- raster(raw_dem)
dem.coords <- coordinates(dem)
datum <- read.csv(file="../data/usa_111k_datum_correction.csv", header = T, sep = ",")

## Attempt KD-Tree lookup
results <- nn2(dem.coords, datum, 4)

## Note:
# This approach will not work if we attempt to implement the interpolation using a 
# SpatialPointDataFrame because it loads all ~360m points into memory. It crashes
# the computer every time. Instead, we should load matrices, use the KD tree to do 
# IDW just like we do in the Python code. It'll be a pain to implemement, but a good 
# learning experience.

write.csv(grid, file="../data/gdem_points.csv", row.names=FALSE, na="")
