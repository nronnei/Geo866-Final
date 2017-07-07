library(rgdal)
library(raster)
library(sp)
library(gstat)

PROJECT_PATH <- "/home/nronnei/gis/class/spatial_analysis/final_project/"
SIM_OUT_PATH <- "./data/srtm/error_surfaces/ked/sim_d_"
N_SIM <- 10

setwd(PROJECT_PATH)


## Projections
prj.nad83 <- CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
prj.wgs84 <- CRS("+proj=longlat +datum=WGS84 +no_defs")
prj.utm <- CRS("+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")



## Read in sample data
yel <- readRDS("./r/sample_data_utm.rds")


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

# This... is... REGRESSION KRIGING!!! Admittedly, a really bad example of it...

dem.rk <- krige(formula = gls.res ~ 1, dem.gls1, newdata = dem.spdf, 
                beta = mean(dem.gls1$gls.res), model = gls.vgm.fit, nmax = 10, nsim = N_SIM)


for (i in 1:length(slot(dem.rk, "data"))) {
  
  sim <- as.data.frame(slot(dem.rk, "coords"))
  sim$z <- slot(dem.rk, "data")[[i]]
  path <- paste(SIM_OUT_PATH, i, ".tif", sep = "")
  
  sim.r <- rasterFromXYZ(sim, crs = prj.utm)
  writeRaster(sim.r, filename = path, format = "GTiff")
  
}

