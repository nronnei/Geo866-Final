library(rgdal)
library(raster)
library(sp)
library(gstat)
library(RColorBrewer)
library(classInt)
library(RANN)

## GLOBALS
PROJECT_PATH <- "/home/nronnei/gis/class/spatial_analysis/final_project/"
POINTS_PATH_OUT <- "r/sample_data_utm.rds"
SAMPLE_SIZE <- 300

##
#####
## Custom Functions
#####
## 


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


polygonFromBbox <- function(bbox, s_srs, t_srs = NULL) {
  
  if (class(bbox) == "matrix") {
    coords <- matrix(c(
      bbox[1,1], bbox[2,1],
      bbox[1,1], bbox[2,2],
      bbox[1,2], bbox[2,2],
      bbox[1,2], bbox[2,1],
      bbox[1,1], bbox[2,1]
    ), ncol = 2, byrow = T)
    
    basePoly <- Polygon(coords)
    extent <- SpatialPolygons(list(Polygons(list(basePoly), id = "a")), proj4string = s_srs)
    
    if (is.null(t_srs)) {
      
      return(extent)
      
    } else {
      
      extent <- spTransform(extent, t_srs)
      return(extent)
      
    }
    
  } else if (class(bbox) == "Extent") {
    
    coords <- matrix(c(
      bbox[1,1], bbox[2,1],
      bbox[1,1], bbox[2,2],
      bbox[1,2], bbox[2,2],
      bbox[1,2], bbox[2,1],
      bbox[1,1], bbox[2,1]
    ), ncol = 2, byrow = T)
    
    basePoly <- Polygon(coords)
    extent <- SpatialPolygons(list(Polygons(list(basePoly), id = "a")), proj4string = s_srs)
    
    if (is.null(t_srs)) {
      
      return(extent)
      
    } else {
      
      extent <- spTransform(extent, t_srs)
      return(extent)
      
    }
  }
  
}


##
#####
## Data Prep
#####
##

setwd(PROJECT_PATH)

## Projections
prj.nad83 <- CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
prj.wgs84 <- CRS("+proj=longlat +datum=WGS84 +no_defs")
prj.utm <- CRS("+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")


## Prep study area boundary, sample area boundary
study_area.nad83 <- readOGR("./data/study_area/study_area_nad83.shp/study_area_nad83.shp", "study_area_nad83")
study_area.utm <- spTransform(study_area.nad83, prj.utm)
sam.45.bb <- bbox(spTransform(readOGR("./data/osm/eer_45km_wgs84.shp", "eer_45km_wgs84"), prj.nad83))
sam.48.bb <- bbox(spTransform(readOGR("./data/osm/eer_48km_wgs84.shp", "eer_48km_wgs84"), prj.nad83)) 


## Prep road data
eer.utm <- readOGR("./data/osm/eer_utm.shp", "eer_utm")


## Load DEMs
srtm <- raster("./data/srtm/srtm_utm_clipped.tif")
srtm.extent <- extent(srtm)
ned <- raster("./data/ned/ned_nad83.tif")



## Sample Reference DEM
# Load datum correction points
dc <- read.csv("./data/points/usa_111k_datum_correction.csv", header = T, sep = ",")
coordinates(dc) <- ~Lon+Lat

# Get 300 random points from within our study area
sa.sam <- spsample(study_area.nad83, SAMPLE_SIZE, "random", bb = sam.45.bb) # slightly smaller bb ensures all sample points fall w/in study area
plot(study_area.nad83, main = "Sample Points")
points(sa.sam, pch = 1, cex = 1)
summary(sa.sam)

# Sample NED Values
elev <- extract(ned, sa.sam, method = 'bilinear')
# Convert to SpatialPointsDataFrame to keep everything in one place
sa.sam <- SpatialPointsDataFrame(sa.sam, data = data.frame(n_dem = elev))

# Use KD-Tree / IDW to get Datum Correction Values
kd_results <- nn2(slot(dc, "coords"), slot(sa.sam, "coords"), k = 6)
print(head(kd_results$nn.idx))
sa.sam$dc <- idw.dc(sa.sam, kd_results, dc)
sa.sam$c_elev <- sa.sam$n_dem + sa.sam$dc

# Convert to UTM to make life easier and maps more rectangular
sam.utm <- spTransform(sa.sam, prj.utm)

# Compare srtm to NED
sam.utm$s_dem <- extract(srtm, sam.utm)
sam.utm$diff <- sam.utm$s_dem - sam.utm$c_elev



##
#####
## Background Investigation
#####
##


## Set up some color palletes
col.cont.greens <- rev(colorRampPalette(brewer.pal(5,"Greys"))(50))
col.cont.spect <-  rev(colorRampPalette(brewer.pal(7,"Spectral"))(50))
col.disc.spect <- rev(brewer.pal(7,"Spectral"))


## Plot DEM 
image(srtm, main = "SRTM (3as Resolution) Yellowstone National Park, WY",
      col = col.cont.greens, xlab = "Easting", ylab = "Northing")
mtext("East Entrance Road in Neon Blue")
lines(eer.utm, col = "#82fffc", lwd = 2)
dev.print(png, "./deliverables/img/srtm-with-roads.png", height=600, width=800)


## Plot the study points with UTM 12N Projection
# Set up breaks, colors
dc.brks.jenks <- classIntervals(sam.utm$dc, 7, "jenks")
diff.col <-  findColours(dc.brks.jenks, col.disc.spect)
# Plot it
plot(study_area.utm, main = "Datum Corrections by Z Score", axes = T)
mtext(paste("Values Range Between ", round(min(sam.utm$dc), digits = 3), "m and ", round(max(sam.utm$dc), digits = 3), "m", sep = "" ))
points(sam.utm["dc"], pch = 3, col = diff.col, cex = abs((sam.utm$dc - mean(sam.utm$dc)) / sd(sam.utm$dc)))
points(sam.utm["dc"], pch = 20, col = diff.col, cex = .75) 
dev.print(png, "./deliverables/img/dc-by-size.png", height=600, width=800) 


## Check out differences in DEM values
sam.hist <- hist(sam.utm$diff, main = "Difference: SRTM - NED", breaks = 30) 
dev.print(png, "./deliverables/img/diff-hist.png", height=600, width=800) 
sam.qqn <- qqnorm(sam.utm$diff, main = "QQ Norm of Difference: SRTM - NED", ylab = "Difference(m)")
dev.print(png, "./deliverables/img/diff-qqn.png", height=600, width=800)


## Plot DEM Difference
# Set up breaks, colors
diff.brks.sd <- classIntervals((sam.utm$diff - mean(sam.utm$diff)) / sd(sam.utm$diff), 7, "sd")
diff.col <-  findColours(diff.brks.sd, col.disc.spect)
# Plot it
plot(study_area.utm, main = "DEM Differences by Z-Score", axes = T)
mtext(paste("Values Range Between ", round(min(sam.utm$diff), digits = 3), "m and ", round(max(sam.utm$diff), digits = 3), "m", sep = "" ))
points(sam.utm["diff"], pch = 3, col = diff.col, cex = abs((sam.utm$diff - mean(sam.utm$diff)) / sd(sam.utm$diff)))
points(sam.utm["diff"], pch = 20, col = diff.col, cex = .75)
dev.print(png, "./deliverables/img/diff-by-zscore.png", height=600, width=800)


## Creating Slope Raster
srtm.slope <- terrain(srtm, opt = "slope", unit = "degrees")
sam.utm$s_slope <- extract(srtm.slope, sam.utm)

## Creating Aspect Raster
srtm.aspect <- terrain(srtm, opt = "aspect", unit = "degrees")
sam.utm$s_aspect <- extract(srtm.aspect, sam.utm)

## Store sample data so we don't mess up our sample in the TSA phase
saveRDS(sam.utm, POINTS_PATH_OUT)

##
#####
## Trend Surface Analysis
#####
##




## Exploratory Analysis



# spplot(sam, layout = list("sp.polygon", study_area))