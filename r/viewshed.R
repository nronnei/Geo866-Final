library(rgdal)
library(maptools)
library(spatstat)
library(rgrass7)

#####################
##     GLOBALS     ##
#####################
PROJECT_PATH = "/home/nronnei/gis/class/spatial_analysis/final_project/"
REALIZATION_PATH = "data/srtm/error_realizations/ked/"
REFERENCE_PATH = "data/ned/dem_clipped.tif"
VIEWSHED_PATH = "data/viewsheds/test/"
GRASS_PATH = "/usr/lib/grass70"

loc <- initGRASS(gisBase = GRASS_PATH, home = PROJECT_PATH, gisDbase = tempdir(), override = T)
setwd(PROJECT_PATH)

#####################
##   READING DATA  ##
#####################

pts <- readOGR("./data/points/eer_sample_points/eer_sample_points.shp", "eer_sample_points")

################################
### GRASS REGION CONFIGURING ###
################################

# CONFIGURING
execGRASS("r.in.gdal", flags = "o", parameters = list(input = "./data/srtm/error_realizations/ked/realization_0.tif", output = "dem"))
execGRASS("g.region", flags = "p", parameters = list(raster = "dem"))

#######################################
### GRASS LINE-OF-SIGHT CALCULATION ###
#######################################

# LOOP
for (i in 1:length(pts@coords)) {
  # CONFIGURE PATH
  out_path <- paste(PROJECT_PATH, VIEWSHED_PATH, "viewshed_", i, ".tif", sep = "")
  print(paste("Calculating Viewshed at ", pts@coords[i,1], ", ", pts@coords[i,2], sep = ""))
  
  # GRASS LOS CALCULATION
  execGRASS("r.viewshed", parameters = list(input = "dem", output = "los", coordinates = pts@coords[i,], observer_elevation = 1.8, max_distance = 16093), flags = c("overwrite", "b"))   

  # WRITE VIEWSHED OUT
  execGRASS("r.out.gdal", 
            parameters = list(input = "los", output = out_path, format = "GTiff"))
}

