library(rgdal)
library(sp)
library(rgrass7)


#####################
##     GLOBALS     ##
#####################


PROJECT_PATH = "/home/nronnei/gis/class/spatial_analysis/final_project/"

REALIZATION_PATH = "data/srtm/error_realizations/ked/"

VIEWSHED_PATH = "data/viewsheds/"

# This is the path to your GRASS executable. You can determine your 
# GRASS_PATH by typing `type grass7x` on *nix systems or `where grass7x`
# (replacing x with your minor version number).
GRASS_PATH = "/usr/lib/grass70"

setwd(PROJECT_PATH)


###################
##   FUNCTIONS   ##
###################

process_realization <- function(viewshed.points, pid.field = "id", viewshed.dem, out.directory, overwrite = TRUE) {
  
  # Generate viewshed for several points and put them in folders for each point.
  # 
  # viewshed.points: object of class SpatialPointsDataFrame containing points to generate viewsheds for
  # pid: character vector that's the name of the point ID field, used for directory creation/selection
  # viewshed.dem: character vector that's the path to the DEM to use in r.viewshed
  # out.directory: character vectore that's the path of the output directory for all the point directories
  # overwrite: UNDER DEVELOPMENT - currently does nothing
  
  loc <- initGRASS(gisBase = GRASS_PATH, home = PROJECT_PATH, gisDbase = tempdir(), override = T)
  
  
  ## CONFIGURE GRASS
  execGRASS("g.proj",
            flags = "j",
            parameters = list(georef = viewshed.dem))
  execGRASS("r.in.gdal", 
            flags = "o", 
            parameters = list(input = viewshed.dem, output = "dem"))
  execGRASS("g.region", 
            flags = "p", 
            parameters = list(raster = "dem"))

  
  point.dirs <- list.dirs(out.directory, recursive = F)

  ## PROCESS POINTS
  for (i in 1:length(viewshed.points)) {
    
    ## CONFIGURE
    id <- viewshed.points[pid.field]@data[i,]
    view.pt <- viewshed.points@coords[i,]
    pt.dir <- paste(out.directory, "point_", id, "/", sep = "")
    out.path <- paste(pt.dir, basename(viewshed.dem), sep = "")
    
    
    ## MANAGE DIRECTORY
    if (dir.exists(pt.dir)) {
      real.dir <- list.files(pt.dir)
      if (file.exists(out.path)) {
        # IF EXISTS, overwrite file
        file.remove(out_path)
      }
    } else {
      # If they don't exist, create the directories
      dir.create(pt.dir)
    }
    
    
    # GRASS LOS CALCULATION
    cat("\n\n\n\n")
    print(paste("Calculating Viewshed of Point ", id, sep = ""))
    print(out.path)
    execGRASS("r.viewshed", 
             flags = c("overwrite", "b"),
             parameters = list(input = "dem", output = "viewshed", coordinates = view.pt, observer_elevation = 1.8, max_distance = 16093))   
    
    
    # WRITE VIEWSHED OUT
    execGRASS("r.out.gdal", 
             parameters = list(input = "viewshed", output = out.path, format = "GTiff", createopt="PROFILE=GDALGeoTIFF"))
    
  }
}


#####################
##   READING DATA  ##
#####################

## Read viewshed points
pts <- readOGR("./data/points/eer_sample_points/eer_sample_points.shp", "eer_sample_points")
realizations <- list.files(REALIZATION_PATH)

for (i in 1:length(realizations)) {
  realization <- paste(REALIZATION_PATH, realizations[i], sep = "")
  process_realization(pts, pid.field = "pid", realization, VIEWSHED_PATH)
}
