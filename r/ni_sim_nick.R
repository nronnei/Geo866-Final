# ni_simulate.R
# This R script performs conditional and 
# unconditional simulations on an area equivalent 
# to the Jura dataset.
# written by A. Shortridge, 3/2004, 3/2008, 11/2009, 11/2013

# build.convex.grid
# A function to develop a gridded set of xyz locations within
# the convex hull of a dataset of points. npts is the 
# approximate number of points to generate
build.convex.grid <- function (x, y, npts) {
  library(splancs) # for gridding and inout functions
  # First make a convex hull border (splancs poly)
  ch <- chull(x, y) # index for pts on convex hull 
  ch <- c(ch, ch[1])
  border <- cbind(x[ch], y[ch])  # This works as a splancs poly
  
  # Now fill it with grid points
  xy.grid <- gridpts(border, npts)
  
  return(xy.grid)
}

# Get the library and data loaded
library(gstat)
library(sp)

data(jura)
jura <- jura.pred
names(jura)
jura[1:10,]   # list the first few rows
length(jura$Ni)  # 259 observations
coordinates(jura) <- c('Xloc','Yloc')  # Turn jura into a Spatial object

########################################
#### Create a surface for prediction ###
cj <- coordinates(jura)
jg <- data.frame(build.convex.grid(cj[,1], cj[,2], 10000))
names(jg) <- c('Xloc', 'Yloc')
jg$x <- jg$Xloc
jg$y <- jg$Yloc

coordinates(jg) <- c('Xloc','Yloc') # Make the thing a spatial object
jg <- as(jg, "SpatialPixelsDataFrame")


################################
## Calculate Variogram for Ni ##

ni.vg1 <- variogram(Ni~1, jura, width=0.08, cutoff=2.5)
plot(ni.vg1)
ni.vmod <- fit.variogram(ni.vg1, vgm("Sph", nugget=10, psill=50, range=1.2)) 
nug <- round(ni.vmod$psill[1], 2)
sill <-round(ni.vmod$psill[2], 2)
range <- round(ni.vmod$range[2], 2)
plot(ni.vg1, model=ni.vmod, 
     main=paste('Spherical, nug=', nug, '; sill=', sill, '; range=', range))
#dev.print(png, "ni_vario.png", height=400, width=500)

#####################
## Simple Kriging #
ni.sk <- krige(Ni~1, jura, jg, model=ni.vmod, nmax=16)
spplot(ni.sk, zcol='var1.pred', col.regions=heat.colors(20), main="Nickel: Simple Kriging Predictions")
#dev.print(png, "ni_predmap.png", height=500, width=500)

#############################
## Unconditional Simulation #

# First set up the 'dummy' gstat object with the desired model
dummy <- gstat(formula=Ni~1, locations=~Xloc+Yloc, dummy=TRUE, 
               model=ni.vmod, beta=mean(jura$Ni), nmax=16)

# Then use predict.gstat to simulate data from dummy to the newdata locations
sk.uncon <- predict(dummy, newdata=jg, nsim=4)

names(sk.uncon)

spplot(sk.uncon, col.regions=heat.colors(20), main='Jura Ni: Unconditional Simulations')
#dev.print(png, "ni_uncon.png", height=680, width=650)

###########################
## Conditional Simulation #

sk.con <- krige(Ni ~ 1, jura, jg, model = ni.vmod, beta=mean(jura$Ni), nmax=16, nsim=4)

spplot(sk.con, col.regions=heat.colors(20), main='Jura Ni: Conditional Simulations (SK)')
#dev.print(png, "ni_con.png", height=680, width=650)

summary(sk.con$sim1)
summary(sk.con$sim2)
var(sk.con$sim1)
summary(ni.sk$var1.pred)
var(ni.sk$var1.pred)

## End of script