library(rgdal)
library(raster)


## GLOBALS
PROJECT_PATH = "/home/nronnei/gis/class/spatial_analysis/final_project/"
RASTER_IN = "data/viewsheds/srtm/"
# This should be a number from 0 - 1 respresenting the minimum (decimal) precent of the 
# realizations in which a point must appear in to be counted as visible.
VISIBLE_THRESHOLD = 1
setwd(PROJECT_PATH)

## Load Data
poi <- readOGR(paste(PROJECT_PATH, "data/points/yellowstone_poi/yellowstone_poi.shp", sep = ""), "yellowstone_poi")

poi$VISIBLE <- 0


input_files <- list.files(paste(PROJECT_PATH, RASTER_IN, sep = ""), include.dirs = F)

for (i in 1:length(input_files)) {
  raster_path <- paste(RASTER_IN, input_files[i], sep = "")
  mc_raster <- raster(raster_path)
  
  poi$VISIBLE <- extract(mc_raster, poi)
  print(paste("POINT ", i, ": ", length(which(poi$VISIBLE >= VISIBLE_THRESHOLD)), sep = ""))
}