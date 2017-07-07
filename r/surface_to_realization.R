library(raster)

## GLOBALS
PROJECT_PATH = "/home/nronnei/gis/class/spatial_analysis/final_project/"
RASTER_IN = "data/srtm/error_surfaces/ked/"
RASTER_OUT = "data/srtm/error_realizations/ked/"
setwd(PROJECT_PATH)


## Load Data

srtm <- raster("./data/srtm/srtm_utm_clipped.tif")
error_surfaces <- list.files(RASTER_IN, full.names = T)
file_count <- length(error_surfaces)
i <- 1

for (surface in error_surfaces) {
  if (endsWith(basename(surface), ".aux.xml")) {
    
    print("Skipping .aux.xml metadata")
    file_count <- file_count - 1
    
  } else {
    
    print(paste("Processing Raster", i, "of", file_count, sep = " "))
    
    sim <- raster(surface)
    realization <- srtm - sim
    path <- paste(RASTER_OUT, "realization_", i, ".tif", sep = "")
    writeRaster(realization, path)
    
    ## Increment the file count
    i <- i + 1
  }
}

