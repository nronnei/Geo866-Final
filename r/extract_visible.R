library(rgdal)
library(raster)


## GLOBALS
PROJECT_PATH = "/home/nronnei/gis/class/spatial_analysis/final_project/"
RASTER_IN = "data/viewsheds/monte-carlo/"
setwd(PROJECT_PATH)

## Load Data
poi <- readOGR(paste(PROJECT_PATH, "data/points/yellowstone_poi/yellowstone_poi.shp", sep = ""), "yellowstone_poi")

poi$VISIBLE <- 0


input_files <- list.files(paste(PROJECT_PATH, "data/viewsheds/ned/", sep = ""), include.dirs = F)

for (i in 1:length(input_files)) {
  raster_path <- paste(PROJECT_PATH, "data/viewsheds/ned/", input_files[i], sep = "")
  mc_raster <- raster(raster_path)
  
  poi$VISIBLE <- extract(mc_raster, poi)
  
  print(paste("POINT ", i, ": ", length(which(poi$VISIBLE==1)), sep = ""))
}