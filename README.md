# Monte Carlo Mountain
*A little experiment with Geostats to demonstrate the propagation of errors in
DEMs and an example of how their effects can be accounted for.*

## Purpose
Yellowstone National Park is one of the most scenic attraction sites in America.
However, finding an optimum overlook site along the US-14 east highway entrance
of Yellowstone National Park is a challenge that needed to be addressed. There
are many variables that determine whether or not a certain point on the highway
can be characterized to be a scenic spot. This paper will try to add value by
addressing the issue through utilization of geospatial techniques.

## Data
We employ the following datasets:
  - ASTER GDEM v2 (1 arcsecond resolution)
  - USGS NED 2013 (1/3 arcsecond resolution)
  - OpenStreetMap Road and Feature Data (Higways, Mountain Peaks, etc.)
  - Yellowstone National Park Boundary

We will use NED as a "ground truth", allowing us to determine the errors at
about 150 "sample points"

## Members of the Group
  - Nick Ronnei
  - Yusri Jamaluddin
  - Ameen Khadim
  - Yachen Xie
