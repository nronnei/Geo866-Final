library(raster)

## GLOBALS
PROJECT_PATH = "/home/nronnei/gis/class/spatial_analysis/final_project/"
VIEWSHED_DIR = "data/viewsheds/"
RASTER_OUT = "data/viewsheds/monte-carlo/"
prj.utm <- CRS("+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
setwd(PROJECT_PATH)

## Load Data
view.dirs <- list.dirs(VIEWSHED_DIR, recursive = F)

for (n in 1:length(view.dirs)) {
  if (substr(basename(view.dirs[n]), 1, 6) == "point_") {
    viewsheds <- list.files(view.dirs[n], include.dirs = F, full.names = T)
    print(paste("Processing", view.dirs[n], sep = " "))
    cv <- 0
    for (i in 1:length(viewsheds)) {
      
      ## Sum raster values
      if (i == 1) {
        cv <- raster(viewsheds[i], crs = prj.utm)
      } else {
        cv <- cv + raster(viewsheds[i], crs = prj.utm)
      }
    }
    # Get Average
    mean_viewshed <- cv / length(viewsheds)
    writeRaster(mean_viewshed, paste(RASTER_OUT, basename(view.dirs[n]), ".tif", sep = ""), format = "GTiff")
  }
}

