library(rgdal)
library(raster)
library(sp)
library(gstat)
library(RColorBrewer)
library(classInt)
library(RANN)
library(gridExtra)
PROJECT_PATH = "/home/nronnei/gis/class/spatial_analysis/final_project/"
setwd(PROJECT_PATH)


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
hist(yel$n_dem, breaks = 15, main = "Histogram of Elevation Values")
qqnorm(yel$n_dem, main = "QQ Normal - Elevation Values")
qqline(yel$n_dem)

# SQRT
hist(sqrt(yel$n_dem), breaks = 15, main = "Histogram of sqrt(Elevation Values)")
qqnorm(sqrt(yel$n_dem), main = "QQ Normal - sqrt(Elevation Values)")
qqline(sqrt(yel$n_dem))

# LOG
hist(log(yel$n_dem), breaks = 15, main = "Histogram of log(Elevation Values)")
qqnorm(log(yel$n_dem), main = "QQ Normal - log(Elevation Values)")
qqline(log(yel$n_dem))

# Judging by the histograms and QQ Plots, the log transformation
# appear to be the most normal. However, it adds a much higher 
# degree of complexity for limited gain. So, we'll stick with 
# the raw values.

# Local Regression
local.east <- loess(n_dem ~ x, data = yel)
plot(yel$x, yel$n_dem)
points(yel$x, local.east$fitted, col="orange", pch=4)
title("Elevation: Local Regression in Easting")

local.north <- loess(n_dem ~ y, data = yel)
plot(yel$y, yel$n_dem)
points(yel$y, local.north$fitted, col="orange", pch=4)
title("Elevation: Local Regression in Northing")

local.slope <- loess(n_dem ~ g_slope, data = yel)
plot(yel$g_slope, yel$n_dem)
points(yel$g_slope, local.slope$fitted, col="orange", pch=4)
title("Elevation: Local Regression in Slope")

local.aspect <- loess(n_dem ~ g_aspect, data = yel)
plot(yel$g_aspect, yel$n_dem)
points(yel$g_aspect, local.aspect$fitted, col="orange", pch=4)
title("Elevation: Local Regression in Aspect")

# This doesn't look promising, maybe its only useful in a particular direction.
# Lets look at some variograms
asp.e_var <- variogram(n_dem ~ g_aspect, yel)
plot(asp.e_var)
asp.de_var <- variogram(n_dem ~ g_aspect, yel, alpha = c(0, 45, 90, 135))
plot(asp.de_var)
# Nothin' much. Directional variograms look pretty much all the same, more or less.


##
#####
## TSA - Regression Kriging (KED)
#####
##

# What we're doing here is trying to predict the reference DEM's values
# based on the available DEM's values, slope, aspect, and location. Some
# of these variables may not be significant, and they will be discarded.

## Let's explore some models
# Start basic
dem.lm0 <- lm(n_dem ~ x + y + g_slope + g_aspect + g_dem, yel)
summary(dem.lm0)
# This definitely confirms that aspect doesn't much matter here.

# Get a little more complex
dem.lm1 <- lm(log(n_dem) ~ x + y + I(x^2) + I(y^2) + I(x^3) + I(y^3) + g_slope + g_aspect, yel)
summary(dem.lm1)

# Do some piecewise removal
dem.lm2 <- lm(n_dem ~ x + y + I(x^2) + I(y^2) + I(x^3) + g_slope + g_dem, yel)
summary(dem.lm2)
dem.lm3 <- lm(n_dem ~ x + y + I(x^2) + I(x^3) + g_slope + g_dem, yel)
summary(dem.lm3)
dem.lm4 <- lm(n_dem ~ y + I(x^2) + I(x^3) + g_slope + g_dem, yel)
summary(dem.lm4)
dem.lm5 <- lm(n_dem ~ I(y^3) + I(x^3) + g_slope + g_dem, yel)
summary(dem.lm5)
dem.lm6 <- lm(n_dem ~ y + x + g_dem + g_slope, yel)
summary(dem.lm6)
# None of these models are any better than the original model

# This is what we came up with:
dem.ols <- lm(n_dem ~ y + x + g_dem + g_slope, yel)
summary(dem.ols)
# Plot diagnostics
par(mfrow=c(2,2))
plot(dem.ols)   # Diagnostic Plots
dev.print(png, "./deliverables/img/ols-dem-diagnostics.png", height=600, width=800) 
dev.off()

yel$ols_res <- dem.ols$residuals

##
#####
## Variography - Regression Kriging (KED)
#####
##

##########################
## Investigate with OLS ##
##########################

# Empirical variogram
dem.ols.evar <- variogram(ols_res ~ 1, data = yel)
# Variogram model
dem.ols.vgm <- fit.variogram(dem.ols.evar, vgm(model = "Sph"))
# Plot it
plot(dem.ols.evar, model = dem.ols.vgm, main = "Variogram of OLS Residuals",
     sub = paste("psill = ", dem.ols.vgm$psill, ", nugget = 0, range = ", dem.ols.vgm$range, sep = "" ))

#####################################################
# Use GLS to generate a trend prediction using the  #
# variogram model developed on the OLS residuals.   #
#####################################################
g <- gstat("dem", n_dem ~ y + x + g_dem + g_slope, data = yel, model = dem.)

##
#####
## Kriging
#####
##


## Get the grid we'll interpolate to
# gdem.sgdf <- readGDAL("./data/dem/dem_utm_clipped.tif")
# pt.mat <- coordinates(gdem.sgdf)
# dem.spdf <- as.data.frame(pt.mat)
# coordinates(dem.spdf) <- ~x+y


# dem.uk <- krige(g_dem ~ y + x + yel$g_slope, yel, dem.sgdf, model = dem.res.vgm)






