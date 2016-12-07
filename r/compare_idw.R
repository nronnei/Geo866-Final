setwd("/home/nronnei/gis/class/spatial_analysis/final_project/")
dc.r <- read.csv("./data/points/study_points_R.csv", header = T, sep = ",")
dc.p <- read.csv("./data/points/study_points_Python.csv", header = T, sep = ",")


com.x <- dc.r$x[dc.r$x == dc.p$x]
com.y <- dc.r$y[dc.r$y == dc.p$y]
com.dc <- dc.r$dc[dc.r$dc == dc.p$dc]

hist(dc.p$dc)
hist(dc.r$dc)

# Doesn't look great. Values are close, but not identical. 
# Shape of the distributions is similar as well, but still not great.