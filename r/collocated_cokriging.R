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
## Data Exploration - Collocated Cokriging
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
local.srtm <- loess(c_elev ~ s_dem, data = yel)
plot(yel$s_dem, yel$c_elev)
points(yel$s_dem, local.srtm$fitted, col="orange", pch=4)
title("Elevation: Local Regression on SRTM Values")


##
#####
## TSA - Collocated Cokriging
#####
##

# What we're doing here is trying to predict the reference DEM's values
# based solely on the available srtm values. Collocated Cokriging relies
# on a single, highly correlated variable and long range spatial 
# autocorrelation to make predictions. In our case we will use it to 
# generate realizations of the SRTM error surface


dem.ols <- lm(c_elev ~ s_dem, yel)
summary(dem.ols)
hist(dem.ols$residuals, breaks = 15, main = "OLS Residuals - Collocated Cokriging")
yel$ols_res <- dem.ols$residuals

# Plot diagnostics
# par(mfrow=c(2,2))
# plot(dem.ols)   # Diagnostic Plots
# dev.print(png, "./deliverables/img/co_krige/ols-dem-diagnostics.png", height=600, width=800) 
# dev.off()


##
#####
## Variography - Regression Kriging (KED)
#####
##


## Empirical variogram
# Isotropic
dem.ols.evar <- variogram(ols_res ~ 1, data = yel)
plot(dem.ols.evar, main = "Variogram Model OLS Residuals")
# Anisotropic
a1 <- c(0, 30, 60, 90, 120, 150)
a2 <- c(0, 45, 90, 135)
a3 <- c(300, 150)
dem.ols.evar.anis <- variogram(ols_res ~ 1, data = yel, alpha = a1)
plot(dem.ols.evar.anis, main = "Anisotropic Variogram of OLS Residuals - Collocated Cokriging")

# Variogram model
# dem.ols.vgm <- fit.variogram(dem.ols.evar, vgm(model = "Wav"))
dem.ols.vgm <- vgm(model = "Wav", nugget = 10, psill = 50, range = 5000)
dem.ols.vgma <- vgm(model = "Gau", nugget = 30, psill = 80, range = 2000, add.to = dem.ols.vgm)

# Plot it
plot(dem.ols.evar, model = dem.ols.vgma, main = "Variogram Model of OLS Residuals")

#####################################################
# Use GLS to generate a trend prediction using the  #
# variogram model developed on the OLS residuals.   #
#####################################################

## ROUND 1 - GLS
g.gls <- gstat(id = "dem", formula = c_elev ~ s_dem, data = yel, model = dem.ols.vgma)
dem.gls0 <- predict(g.gls, newdata = yel, BLUE = T)
dem.gls0$gls.res <- yel$c_elev - dem.gls0$dem.pred
summary(dem.gls0$gls.res)

# Empirical variogram of GLS residuals
dem.gls.evar <- variogram(gls.res ~ 1, data = dem.gls0)
plot(dem.gls.evar, main = "Variogram of GLS Residuals")
plot(dem.ols.evar, pch = 1, main = "Variogram of OLS Residuals")

# Variogram model of GLS residuals
dem.gls.vgm0 <- vgm(model = "Wav", nugget = 20, psill = 40, range = 5500)
dem.gls.vgm0a <- vgm(model = "Gau", nugget = 52, psill = 60, range = 2000, add.to = dem.gls.vgm0)
plot(dem.gls.evar, model = dem.gls.vgm0a, main = "Variogram Model of GLS Residuals - Round 1")


## ROUND 2 - GLS
g.gls <- gstat(id = "dem", formula = c_elev ~ s_dem , data = yel, model = dem.gls.vgm0a)
dem.gls1 <- predict(g.gls, newdata = yel, BLUE = T)
dem.gls1$gls.res <- yel$c_elev - dem.gls1$dem.pred
summary(dem.gls1$gls.res)

# Empirical variogram of GLS residuals
dem.gls.evar1 <- variogram(gls.res ~ 1, data = dem.gls1)
plot(dem.gls.evar1, main = "Variogram of GLS Residuals")
# Variogram model of GLS residuals
dem.gls.vgm1 <- vgm(model = "Wav", nugget = 20, psill = 43, range = 5200)
dem.gls.vgm1a <- vgm(model = "Gau", nugget = 52, psill = 60, range = 1800, add.to = dem.gls.vgm1)
plot(dem.gls.evar1, model = dem.gls.vgm1a, main = "Variogram Model of GLS Residuals - Round 2")


## ROUND 3 - GLS
g.gls <- gstat(id = "dem", formula = c_elev ~ s_dem , data = yel, model = dem.gls.vgm1a)
dem.gls2 <- predict(g.gls, newdata = yel, BLUE = T)
dem.gls2$gls.res <- yel$c_elev - dem.gls2$dem.pred

# Values haven't changed enough to need refitting, so we're done
summary(dem.gls2$gls.res)
summary(dem.gls1$gls.res)

dem.gls.vgm <- dem.gls.vgm1a


##
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

# This... is... COLLOCATED COKRIGING!!! Admittedly, a really bad example of it...

dem.co <- krige(formula = gls.res ~ 1, dem.gls2, newdata = dem.spdf, 
                beta = mean(dem.gls2$gls.res), model = dem.gls.vgm, nmax = 8, nsim = 5)

for (i in 1:length(slot(dem.co, "data"))) {
  
  sim <- as.data.frame(slot(dem.co, "coords"))
  sim$z <- slot(dem.co, "data")[[i]]
  path <- paste("./data/srtm/error_surfaces/sim_0", i, ".tif", sep = "")
  
  sim.r <- rasterFromXYZ(sim, crs = prj.utm)
  writeRaster(sim.r, filename = path, format = "GTiff")
  
}
