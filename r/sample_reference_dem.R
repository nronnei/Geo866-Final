library(rgdal)
library(raster)
library(sp)
library(gstat)
library(RANN)

PROJECT_PATH = "/home/nronnei/gis/class/spatial_analysis/final_project/"

great_circle_dist <- function(long1, lat1, long2, lat2) {
  # Source: https://www.r-bloggers.com/great-circle-distance-calculations-in-r/
  R <- 6371.01 # Earth mean radius [km]
  d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
  return(d) # Distance in km
}

idw.dc <- function(study_points, kd_results, dc_points, k = 2) {
  
  # Performs IDW interpolation to a coordinate (Lon, Lat dictionary) using
  # coordinates stored in a 2D point tree structure. Distances are
  # spherical, as this function assumes coordinates are lat-lon.

  idw.z <- rep(0, length(study_points[,1]))

  for (n in 1:length(study_points)) {
    
    # Get query point lon/lat
    q.lon <- slot(study_points, "coords")[n, 1]
    q.lat <- slot(study_points, "coords")[n, 2]
    
    # Get info for each neighbor to use in IDW
    neighbors <- kd_results$nn.idx[n,]
    n.lon <- slot(dc_points, "coords")[neighbors, 1]
    n.lat <- slot(dc_points, "coords")[neighbors, 2]
    n.vals <- dc_points$datumCorrection[neighbors]
    
    # Calculate weights
    weights <- 1 / (great_circle_dist(q.lon, q.lat, n.lon, n.lat) ^ k)
    
    # Calculate weighted values
    w.vals <- n.vals * weights

    # Do IDW
    if (sum(weights) > 0) {
      
      idw.z[n] <- sum(w.vals) / sum(weights)
      
    } else {
      
      idw.z[n] <- NA
      
    }
  }
  return(idw.z)
}


setwd(PROJECT_PATH)
src_prj <- CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
dst_prj <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## Loading Data
# Load in study area boundary
sa <- readOGR("./data/study_area/study_area_nad83.shp/study_area_nad83.shp", "study_area_nad83")
bb <- bbox(spTransform(readOGR("./data/osm/eer_48km_wgs84.shp", "eer_48km_wgs84"), src_prj)) 
# Load raster
dem <- raster("./data/ned/ned_nad83.tif")
# Load datum correction
dc <- read.csv("./data/points/usa_111k_datum_correction.csv", header = T, sep = ",")
coordinates(dc) <- ~Lon+Lat

## Sampling
# Get 300 random points from within our study area
sa.sam <- spsample(sa, 300, "random", bb = bb) # As far as I can tell, bbox doesn't really do much...
plot(sa, main = "Sample Points")
points(sa.sam, pch = 1, cex = 1)
summary(sa.sam)

elev <- extract(dem, sa.sam, method = 'bilinear')
# Convert to SpatialPointsDataFrame to keep everything in one place
sa.sam <- SpatialPointsDataFrame(sa.sam, data = data.frame(elev = elev))
# Add placeholder for datum correction values so it will play nice with Python tools.
sa.sam$dc <- 0

# Write out table for comparison to Python implementation of the datum correction
sa.sam.wgs <- spTransform(sa.sam, dst_prj)
write.csv(sa.sam.wgs, "./data/points/study_points_raw.csv")


kd_results <- nn2(slot(dc, "coords"), slot(sa.sam, "coords"), k = 6)
print(head(kd_results$nn.idx))

sa.sam$dc <- idw.dc(sa.sam, kd_results, dc)

sa.sam.wgs <- spTransform(sa.sam, dst_prj)
write.csv(sa.sam.wgs, "./data/points/study_points_R.csv")

