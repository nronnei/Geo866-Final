# water_tsa.R
# Very skeletal R code to work with 
# A water quality dataset for Tucson (Mosbacher MA thesis data) from
# http://www.u.arizona.edu/ic/math574/574datasets.html (dead link)
# DO NOT USE THIS DATASET FOR ANYTHING SERIOUS! IT HAS BEEN 
# MODIFIED AND SOME VALUES ARE NOT THE MEASUREMENTS MADE BY MOSBACHER!

# Written/modified October, 2010, by A. Shortridge
# Last modified 10/2014


#par(ask=T)
## Importing Raster Data
# Raster data can be exported from (most) GIS as GeoTiffs. Here's how to read a GeoTiff
library(rgdal)

dem<-readGDAL("dem.tif")

#dem<-dem[1:150,1:150]  # Just grab the first 150 rows and 150 columns
names(dem) <- "z"
print(summary(dem))
spplot(dem['z'], col.regions=terrain.colors(20), scales=list(draw=TRUE))
# Check out the raster library on your own computer if interested - 
# there's tons of basic GIS raster tools in it.

## Now to the Trend Surface Analysis and Residual Variography!
# if the dataset is located in your workspace, 
# the following line should read it in.
wat <- read.table("tucson.txt", header=T, na.strings="NA")
names(wat)
wat[1:10,] # look at the  first few lines 

# Some visualizations - make more!
plot(wat$x, wat$y, asp=1)

# Create an sp object and a gridded set of points to interpolate to
library(sp)
coordinates(wat) <- ~x+y

# Make a grid over the wat region
wat.grid <- makegrid(wat, n = 2400)
names(wat.grid) <- c('x','y')
coordinates(wat.grid) <- ~x+y   # convert to SpatialPoints
wat.grid <- SpatialPointsDataFrame(wat.grid,as.data.frame(x=coordinates(wat.grid)[,1]))
gridded(wat.grid) <- TRUE

# Do an IDW interpolation using the gstat command
library(gstat)

idw.2 <- idw(Sodium~1, wat, newdata=wat.grid, idp=2, nmax=12)
wat.grid$idw <- idw.2$var1.pred

# Plot the Image
spplot(wat.grid, 'idw', col.regions = terrain.colors(20), 
       sp.layout = list('sp.points', wat, pch=3),
       main = 'Tucson Sodium Concentration (IDW Interpolation)')

# Histograms, etc
hist(wat$Sodium, nclass=12, main = "Histogram of Raw Sodium")
hist(sqrt(wat$Sodium), nclass=12, main = "Histogram of Sqrt(Sodium)")
hist(log(wat$Sodium), nclass=12, main = "Histogram of Log(Sodium)")

# Exploratory Trend Techniques
local.east <- loess(Sodium ~ x, data=wat)
local.north <- loess(Sodium ~ y, data=wat)

plot(wat$x, wat$Sodium)
points(wat$x, local.east$fitted, col="orange", pch=4)
title("Sodium: Local Regression in Easting")

plot(wat$y, wat$Sodium)
points(wat$y, local.north$fitted, col="orange", pch=4)
title("Sodium: Local Regression in Northing")

## Look at some summary stats
wat.tsm0 <- lm(Sodium ~ x + y, data=wat)
summary(wat.tsm0)
par(mfrow=c(2,2))
plot(wat.tsm0)
dev.off


## Try local regressions on sqrt and log transformations! What looks best? ##
## SQRT Transformation
local.s.east <- loess(sqrt(Sodium) ~ x, data=wat)
local.s.north <- loess(sqrt(Sodium) ~ y, data=wat)

plot(wat$x, sqrt(wat$Sodium))
points(wat$x, local.s.east$fitted, col="orange", pch=4)
title("Sqrt(Sodium): Local Regression in Easting")

plot(wat$y, sqrt(wat$Sodium))
points(wat$y, local.s.north$fitted, col="orange", pch=4)
title("Sqrt(Sodium): Local Regression in Northing")


## Log Transformation
local.l.east <- loess(log(Sodium) ~ x, data=wat)
local.l.north <- loess(log(Sodium) ~ y, data=wat)

plot(wat$x, log(wat$Sodium))
points(wat$x, local.l.east$fitted, col="orange", pch=4)
title("Log(Sodium): Local Regression in Easting")

plot(wat$y, log(wat$Sodium))
points(wat$y, local.l.north$fitted, col="orange", pch=4)
title("Log(Sodium): Local Regression in Northing")





## Are there any outliers?
## Can you determine what order of polynomial might be best for this data? ##






# Various Trend Surface Models (one is here as an example)
wat.tsm1 <- lm(log(Sodium) ~ x + y, data=wat)
summary(wat.tsm1)  # Linear: R2 = xx
par(mfrow=c(2,2))
plot(wat.tsm1)   # Diagnostic Plots
dev.off
#####
##
# My trend model
##
#####
wat.tsm2 <- lm(log(Sodium) ~ x + y + I(y^2) + I(y^3), data=wat)
summary(wat.tsm2)
par(mfrow=c(2,2))
plot(wat.tsm2)   # Diagnostic Plots
dev.off
###############################################################


# use step() to automatically do stepwise removal, or consider removing terms manually!

# Plot the model and residuals. Feel free to create and look at other plots!
wat.grid$tsm1 <- predict(wat.tsm1, wat.grid)
# Nick's model
wat.grid$tsm2 <- predict(wat.tsm2, wat.grid)

spplot(wat.grid, 'tsm1', col.regions = terrain.colors(20), 
       sp.layout = list('sp.points', wat, pch=3),
       main= 'Tucson Sodium Concentration (Trend Model Surface)')

# Plot Nick's model
spplot(wat.grid, 'tsm2', col.regions = terrain.colors(20), 
       sp.layout = list('sp.points', wat, pch=3),
       main="Nick's Tucson Sodium Concentration (Trend Model Surface)")

# extract residuals, then interpolate across the grid
wat$resids1 <- wat.tsm1$resid
idw.2 <- idw(resids1~1, wat, newdata=wat.grid, idp=2, nmax=12)
wat.grid$resids1 <- idw.2$var1.pred
# Plot it
bluered <- colorRampPalette(c('blue','gray','red'))
spplot(wat.grid, 'resids1', col.regions = bluered(20), 
       sp.layout = list('sp.points', wat, col='black', pch=3),
       main = 'Tucson Sodium TSM Residuals - Just for visualization!')

# Nick's model
# extract residuals, then interpolate across the grid
wat$resids2 <- wat.tsm2$resid
idw.2.2 <- idw(resids2~1, wat, newdata=wat.grid, idp=2, nmax=12)
wat.grid$resids2 <- idw.2.2$var1.pred
# Plot it 
bluered <- colorRampPalette(c('blue','gray','red'))
spplot(wat.grid, 'resids2', col.regions = bluered(20), 
       sp.layout = list('sp.points', wat, col='black', pch=3),
       main = "Nick's Tucson Sodium TSM Residuals")


#### Variograms for Model Residuals
library(gstat)  # other libraries also have variogram functions

# Generate an omnidirectional variogram
rvar <- variogram(resids1~1, wat)
plot(rvar, main="Omnidirectional Variogram for Sodium")

# Now generate a robust variogram
rvar2 <- variogram(resids1~1, wat, cressie=T)
plot(rvar2, main="Robust Estimator")

# Now construct and generate a plot of 4 directional variograms
rvar3 <- variogram(resids1~1, wat, alpha=c(0,45,90,135))
plot(rvar3, main="Directional Variogram for Sodium")


## For Nick's Model
# Generate an omnidirectional variogram
nvar <- variogram(resids2~1, wat)
plot(nvar, main="Omnidirectional Variogram for Sodium")

# Now generate a robust variogram
nvar2 <- variogram(resids2~1, wat, cressie=T)
plot(nvar2, main="Robust Estimator")

# Now construct and generate a plot of 4 directional variograms
nvar3 <- variogram(resids2~1, wat, alpha=c(0,45,90,135))
plot(nvar3, main="Directional Variogram for Sodium")


