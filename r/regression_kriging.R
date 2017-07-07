library(rgdal)
library(raster)
library(sp)
library(gstat)
PROJECT_PATH = "/home/nronnei/gis/class/spatial_analysis/final_project/"
setwd(PROJECT_PATH)


## Projections
prj.nad83 <- CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
prj.wgs84 <- CRS("+proj=longlat +datum=WGS84 +no_defs")
prj.utm <- CRS("+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")


##
#####
## Data Exploration - Regression Kriging (KED)
#####
##


## Begin exploratory analysis
# Read in data
yel <- readRDS("./r/sample_data_utm.rds")

# Compare data and basic transformations
# RAW
hist(yel$c_elev, breaks = 15, main = "Histogram of Elevation Values")
qqnorm(yel$c_elev, main = "QQ Normal - Elevation Values")
qqline(yel$c_elev)

# SQRT
hist(sqrt(yel$c_elev), breaks = 15, main = "Histogram of sqrt(Elevation Values)")
qqnorm(sqrt(yel$c_elev), main = "QQ Normal - sqrt(Elevation Values)")
qqline(sqrt(yel$c_elev))

# LOG
hist(log(yel$c_elev), breaks = 15, main = "Histogram of log(Elevation Values)")
qqnorm(log(yel$c_elev), main = "QQ Normal - log(Elevation Values)")
qqline(log(yel$c_elev))

# Judging by the histograms and QQ Plots, the log transformation
# appear to be the most normal. However, it adds a much higher 
# degree of complexity for limited gain. So, we'll stick with 
# the raw values.

# Local Regression
local.east <- loess(diff ~ x, data = yel)
plot(yel$x, yel$diff)
points(yel$x, local.east$fitted, col="orange", pch=4)
title("Elevation: Local Regression in Easting")

local.north <- loess(diff ~ y, data = yel)
plot(yel$y, yel$diff)
points(yel$y, local.north$fitted, col="orange", pch=4)
title("Elevation: Local Regression in Northing")

local.slope <- loess(diff ~ s_slope, data = yel)
plot(yel$s_slope, yel$diff)
points(yel$s_slope, local.slope$fitted, col="orange", pch=4)
title("Elevation: Local Regression in Slope")

local.aspect <- loess(diff ~ s_aspect, data = yel)
plot(yel$s_aspect, yel$diff)
points(yel$s_aspect, local.aspect$fitted, col="orange", pch=4)
title("Elevation: Local Regression in Aspect")

local.gdem <- loess(diff ~ s_dem, data = yel)
plot(yel$s_dem, yel$diff)
points(yel$s_dem, local.gdem$fitted, col = "orange", pch=4)


##
#####
## TSA - Regression Kriging (KED)
#####
##


# What we're doing here is trying to predict the reference DEM's values
# based on the available DEM's values, slope, aspect, and location. Some
# of these variables may not be significant, and they will be discarded.

## OLS
dem.ols <- lm(c_elev ~ s_slope + I(x^2) + x + y + x*y , yel)
summary(dem.ols)

yel$ols_res <- dem.ols$residuals
ols.resid.evar <- variogram(ols_res ~ 1, data = yel)
plot(ols.resid.evar)

ols.vgm <- vgm(psill = 58000, nugget = 5000, range = 7500, model = "Ste")
ols.vgm.fit <- fit.variogram(ols.resid.evar, ols.vgm)
plot(ols.resid.evar, model = ols.vgm.fit)

## GLS - ROUND 1
g.gls <- gstat(id = "dem", formula = c_elev ~ s_slope + I(x^2) + x + y + x*y, data = yel, model = ols.vgm.fit)
dem.gls0 <- predict(g.gls, newdata = yel, BLUE = T)
dem.gls0$gls.res <- yel$c_elev - dem.gls0$dem.pred
summary(dem.gls0$gls.res)

# Empirical variogram of GLS residuals
gls.resid.evar <- variogram(gls.res ~ 1, data = dem.gls0)
# plot(gls.resid.evar, main = "Variogram of GLS Residuals")
# plot(ols.resid.evar, pch = 1, main = "Variogram of OLS Residuals")

# Variogram model of GLS residuals
gls.vgm.fit <- fit.variogram(gls.resid.evar, model = ols.vgm.fit)
plot(gls.resid.evar, model = gls.vgm.fit, main = "Variogram Model of GLS Residuals - Round 1")


## GLS - ROUND 2
g.gls <- gstat(id = "dem", formula = c_elev ~ s_slope + I(x^2) + x + y + x*y, data = yel, model = gls.vgm.fit)
dem.gls1 <- predict(g.gls, newdata = yel, BLUE = T)
dem.gls1$gls.res <- yel$c_elev - dem.gls1$dem.pred
summary(dem.gls1$gls.res)

# Empirical variogram of GLS residuals
gls.resid.evar <- variogram(gls.res ~ 1, data = dem.gls1)
# plot(gls.resid.evar, main = "Variogram of GLS Residuals")
# plot(ols.resid.evar, pch = 1, main = "Variogram of OLS Residuals")

# Variogram model of GLS residuals
gls.vgm.fit <- fit.variogram(gls.resid.evar, model = gls.vgm.fit)
plot(gls.resid.evar, model = gls.vgm.fit, main = "Variogram Model of GLS Residuals - Round 2")


#
#####
## Kriging - Regression Kriging (KED)
#####
##


## Get the grid we'll interpolate to
srtm.sgdf <- readGDAL("./data/srtm/srtm_utm_clipped.tif")
pt.mat <- coordinates(srtm.sgdf)
dem.spdf <- as.data.frame(pt.mat)
coordinates(dem.spdf) <- ~x+y
sp::proj4string(dem.spdf) <- prj.utm

# This... is... REGRESSION KRIGING!!! Admittedly, a really bad example of it...

dem.rk <- krige(formula = gls.res ~ 1, dem.gls1, newdata = dem.spdf, 
                beta = mean(dem.gls1$gls.res), model = gls.vgm.fit, nmax = 10, nsim = 1)

# spplot(dem.rk['var1.pred'])

for (i in 1:length(slot(dem.rk, "data"))) {
  
  sim <- as.data.frame(slot(dem.rk, "coords"))
  sim$z <- slot(dem.rk, "data")[[i]]
  path <- paste("./data/srtm/error_surfaces/ked/sim_a_", i, ".tif", sep = "")
  
  sim.r <- rasterFromXYZ(sim, crs = prj.utm)
  writeRaster(sim.r, filename = path, format = "GTiff")
  
}

