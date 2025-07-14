## trait space & island biogeography
## Corey Bradshaw & John Llewelyn
## Global Ecology, Flinders University
## globalecologyflinders.com
## July 2025

library(dismo)
library(dplyr)
library(galah)
library(gawdis)
library(gbm)
library(ggplot2)
library(gridExtra)
library(mFD)
library(readxl)
library(sf)
library(terra)

############################
## all Australian islands
############################

## important all island data
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/data")
ausisldat <- read.csv("ausIslDist2mainlAlbers.csv", header=TRUE)
head(ausisldat)
dim(ausisldat)

## histograms
hist(log10(ausisldat$area), xlab="log10(island area)", ylab="frequency", main="") # area in metres squared
hist(log10(ausisldat$distance), xlab="log10(distance to mainland)", ylab="frequency", main="") # distance in metres

## plot spatial data
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/gis")
# AUS islands
ausisles <- vect("ausislands.shp")
ausisles
plot(ausisles, col="red", border=NULL, lwd=0.5)

## mainlands
mainl <- vect("mainlandsLL.shp")
mainl
plot(mainl, col="tan", alpha=0.5, border=NULL, lwd=0.5, add=TRUE)

## import corrected island data (removing Tasmania, its islands, islands in Torres Strait closer to PNG,
## some islands off Timor, and a few south of Kangaroo Island) — remove 'stepping stone' effect where possible
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
ausisldatcorr <- read.csv("d2ausislandscorr.csv", header=TRUE)
head(ausisldatcorr)
dim(ausisldatcorr)
range(ausisldatcorr$area_2/1e6)
range(ausisldatcorr$distance/1000)
sort(ausisldatcorr$distance/1000) # area)


## histograms
hist(log10(ausisldatcorr$area_2/1e6), xlab="log10 island area (km2)", ylab="frequency", main="") # area in metres squared
range(ausisldatcorr$area_2/1e6)
sort(ausisldatcorr$area_2/1e6, decreasing=TRUE)[1:10] # area in km2
hist(log10(ausisldatcorr$distance/1000), xlab="log10 distance to mainland (km))", ylab="frequency", main="") # distance in metres
sort(ausisldatcorr$distance/1000, decreasing=TRUE)[1:10] # area in km2

## AUS islands corrected as above
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/gis")
ausislescorr <- vect("ausislandsAlbersGeom_corrected.shp")
ausislescorr
plot(ausislescorr, col="red", border=NULL, lwd=0.5)

## mainlands
mainl <- vect("mainlandsLL.shp")
mainl
crsAlbersGeom <- "epsg:9473"
mainAlbersGeom <- terra::project(mainl, crsAlbersGeom)
plot(mainAlbersGeom, col="tan", alpha=0.5, border=NULL, lwd=0.5, add=TRUE)

## set bounding box for Aus islands
ausislescorrLL <- terra::project(ausislescorr, "epsg:4326")
ausisles.ext <- ext(ausislescorrLL)
minlon <- ausisles.ext[1]
maxlon <- ausisles.ext[2]
minlat <- ausisles.ext[3]
maxlat <- ausisles.ext[4]

# split in half for smaller packet downloads
minlon1 <- minlon
maxlon1 <- mean(c(ausisles.ext[1],ausisles.ext[2]))
minlon2 <- maxlon1
maxlon2 <- mean(c(minlon2,maxlon))
minlon3 <- maxlon2
maxlon3 <- maxlon

minlat1 <- minlat
maxlat1 <- -30
minlat2 <- -30
maxlat2 <- maxlat

sf_ausbbox <- c(minlon,minlat,maxlon,maxlat)
crs(ausisles)

sf_ausbbox <- st_bbox(sf_ausbbox, crs = 4326)
print(sf_ausbbox)

sf_ausbbox_polygon <- st_as_sfc(sf_ausbbox)

sf_ausbbox_df <- st_sf(geometry = sf_ausbbox_polygon)
print(sf_ausbbox_df)

sf_ausbbox1 <- c(minlon1,minlat,maxlon1,maxlat)
attr(sf_ausbbox1, "names")[3] <- "xmax"
sf_ausbbox1 <- st_bbox(sf_ausbbox1, crs = 4326)
print(sf_ausbbox1)
sf_ausbbox_polygon1 <- st_as_sfc(sf_ausbbox1)
sf_ausbbox_df1 <- st_sf(geometry = sf_ausbbox_polygon1)
print(sf_ausbbox_df1)

sf_ausbbox2 <- c(minlon2,minlat,maxlon2,maxlat)
attr(sf_ausbbox2, "names")[1] <- "xmin"
attr(sf_ausbbox2, "names")[3] <- "xmax"
sf_ausbbox2 <- st_bbox(sf_ausbbox2, crs = 4326)
print(sf_ausbbox2)
sf_ausbbox_polygon2 <- st_as_sfc(sf_ausbbox2)
sf_ausbbox_df2 <- st_sf(geometry = sf_ausbbox_polygon2)
print(sf_ausbbox_df2)

sf_ausbbox3 <- c(minlon3,minlat1,maxlon3,maxlat1)
attr(sf_ausbbox3, "names")[1] <- "xmin"
attr(sf_ausbbox3, "names")[4] <- "ymax"
sf_ausbbox3 <- st_bbox(sf_ausbbox3, crs = 4326)
print(sf_ausbbox3)
sf_ausbbox_polygon3 <- st_as_sfc(sf_ausbbox3)
sf_ausbbox_df3 <- st_sf(geometry = sf_ausbbox_polygon3)
print(sf_ausbbox_df3)

sf_ausbbox4 <- c(minlon3,minlat2,maxlon3,maxlat2)
attr(sf_ausbbox4, "names")[1] <- "xmin"
attr(sf_ausbbox4, "names")[2] <- "ymin"
sf_ausbbox4 <- st_bbox(sf_ausbbox4, crs = 4326)
print(sf_ausbbox4)
sf_ausbbox_polygon4 <- st_as_sfc(sf_ausbbox4)
sf_ausbbox_df4 <- st_sf(geometry = sf_ausbbox_polygon4)
print(sf_ausbbox_df4)


# get ALA data
# set configuration
galah_config(email = "cjabradshaw@gmail.com")
galah_config(caching = TRUE)
show_all(reasons)
galah_config(download_reason_id = 4)
galah_config(verbose = TRUE)

# mammals
ausmamala.dat1 <- galah_call(method = "data",
                         type = "occurrences") |>
  identify("Mammalia") |>
  geolocate(sf_ausbbox_polygon, type = "bbox") |>
  atlas_occurrences()
head(ausmamala.dat)

# save to RDS
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/rds")
saveRDS(ausmamala.dat1, file = "ausmamala1.rds")

# **********************************************
# skip step above for data downloaded 09.07.2025
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/rds")
ausmamala.dat1 <- readRDS("ausmamala1.rds")
# **********************************************

# birds
ausavala1.dat <- galah_call(method = "data",
                        type = "occurrences") |>
  identify("Aves") |>
  geolocate(sf_ausbbox_polygon1, type = "bbox") |>
  atlas_occurrences()
head(ausavala1.dat)

setwd("~/Documents/Papers/Biogeography/Aus Isl traits/rds")
saveRDS(ausavala1.dat, file = "ausavala1.rds")

ausavala2.dat <- galah_call(method = "data",
                            type = "occurrences") |>
  identify("Aves") |>
  geolocate(sf_ausbbox_polygon2, type = "bbox") |>
  atlas_occurrences()
head(ausavala2.dat)

setwd("~/Documents/Papers/Biogeography/Aus Isl traits/rds")
saveRDS(ausavala2.dat, file = "ausavala2.rds")

ausavala3.dat <- galah_call(method = "data",
                            type = "occurrences") |>
  identify("Aves") |>
  geolocate(sf_ausbbox_polygon3, type = "bbox") |>
  atlas_occurrences()

setwd("~/Documents/Papers/Biogeography/Aus Isl traits/rds")
saveRDS(ausavala3.dat, file = "ausavala3.rds")

ausavala4.dat <- galah_call(method = "data",
                            type = "occurrences") |>
  identify("Aves") |>
  geolocate(sf_ausbbox_polygon4, type = "bbox") |>
  atlas_occurrences()

setwd("~/Documents/Papers/Biogeography/Aus Isl traits/rds")
saveRDS(ausavala4.dat, file = "ausavala4.rds")

# **********************************************
# skip step above for data downloaded 09.07.2025
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/rds")
ausavala1.dat <- readRDS("ausavala1.rds")
ausavala2.dat <- readRDS("ausavala2.rds")
ausavala3.dat <- readRDS("ausavala3.rds")
ausavala4.dat <- readRDS("ausavala4.rds")
# **********************************************

# combine Aves datasets
ausavala.dat <- rbind(ausavala1.dat, ausavala2.dat, ausavala3.dat, ausavala4.dat)
dim(ausavala.dat)
rm(ausavala1.dat, ausavala2.dat, ausavala3.dat, ausavala4.dat)

## reptiles
ausreptala.dat <- galah_call(method = "data",
                          type = "occurrences") |>
  identify("Reptilia") |>
  geolocate(sf_ausbbox_polygon, type = "bbox") |>
  atlas_occurrences()
head(ausreptala.dat)

setwd("~/Documents/Papers/Biogeography/Aus Isl traits/rds")
saveRDS(ausreptala.dat, file = "ausreptala.rds")

# **********************************************
# skip step above for data downloaded 09.07.2025
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/rds")
ausreptala.dat <- readRDS("ausreptala.rds")
# **********************************************


## amphibians
ausamphala.dat <- galah_call(method = "data",
                             type = "occurrences") |>
  identify("Amphibia") |>
  geolocate(sf_ausbbox_polygon, type = "bbox") |>
  atlas_occurrences()
head(ausamphala.dat)

setwd("~/Documents/Papers/Biogeography/Aus Isl traits/rds")
saveRDS(ausamphala.dat, file = "ausamphala.rds")

# **********************************************
# skip step above for data downloaded 09.07.2025
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/rds")
ausamphala.dat <- readRDS("ausamphala.rds")
# **********************************************


## remove entries for most of mainland Australia (5 overlapping polygons)
rem.poly1.minlat <- -31.2
rem.poly1.maxlat <- -22.2
rem.poly1.minlon <- 116.7
rem.poly1.maxlon <- 148.6

rem.poly2.minlat <- -35.2
rem.poly2.maxlat <- -23.1
rem.poly2.minlon <- 138.8
rem.poly2.maxlon <- 150.3

rem.poly3.minlat <- -31.3
rem.poly3.maxlat <- -16.0
rem.poly3.minlon <- 125.2
rem.poly3.maxlon <- 136.7

rem.poly4.minlat <- -33.4
rem.poly4.maxlat <- -21.6
rem.poly4.minlon <- 116.0
rem.poly4.maxlon <- 123.7

rem.poly5.minlat <- -21.8
rem.poly5.maxlat <- -13.0
rem.poly5.minlon <- 130.6
rem.poly5.maxlon <- 135.2

rem.poly6.minlat <- -23.3
rem.poly6.maxlat <- -14.8
rem.poly6.minlon <- 141.8
rem.poly6.maxlon <- 144.8

rem.poly7.minlat <- -25.8
rem.poly7.maxlat <- -17.9
rem.poly7.minlon <- 122.5
rem.poly7.maxlon <- 145.9

rem.poly8.minlat <- -65.0
rem.poly8.maxlat <- -50.0
rem.poly8.minlon <- 90.0
rem.poly8.maxlon <- 150.0

rem.poly9.minlat <- -37.5
rem.poly9.maxlat <- -33.6
rem.poly9.minlon <- 140.2
rem.poly9.maxlon <- 149.8

# Mammalia
ausmamala.dat2 <- ausmamala.dat1[!(ausmamala.dat1$decimalLatitude >= rem.poly1.minlat & 
                                    ausmamala.dat1$decimalLatitude <= rem.poly1.maxlat &
                                    ausmamala.dat1$decimalLongitude >= rem.poly1.minlon &
                                    ausmamala.dat1$decimalLongitude <= rem.poly1.maxlon),]
ausmamala.dat3 <- ausmamala.dat2[!(ausmamala.dat2$decimalLatitude >= rem.poly2.minlat &
                                     ausmamala.dat2$decimalLatitude <= rem.poly2.maxlat &
                                    ausmamala.dat2$decimalLongitude >= rem.poly2.minlon &
                                    ausmamala.dat2$decimalLongitude <= rem.poly2.maxlon),]
ausmamala.dat4 <- ausmamala.dat3[!(ausmamala.dat3$decimalLatitude >= rem.poly3.minlat &
                                     ausmamala.dat3$decimalLatitude <= rem.poly3.maxlat &
                                     ausmamala.dat3$decimalLongitude >= rem.poly3.minlon &
                                     ausmamala.dat3$decimalLongitude <= rem.poly3.maxlon),]
ausmamala.dat5 <- ausmamala.dat4[!(ausmamala.dat4$decimalLatitude >= rem.poly4.minlat &
                                     ausmamala.dat4$decimalLatitude <= rem.poly4.maxlat &
                                     ausmamala.dat4$decimalLongitude >= rem.poly4.minlon &
                                     ausmamala.dat4$decimalLongitude <= rem.poly4.maxlon),]
ausmamala.dat6 <- ausmamala.dat5[!(ausmamala.dat5$decimalLatitude >= rem.poly5.minlat &
                                    ausmamala.dat5$decimalLatitude <= rem.poly5.maxlat &
                                    ausmamala.dat5$decimalLongitude >= rem.poly5.minlon &
                                    ausmamala.dat5$decimalLongitude <= rem.poly5.maxlon),]
ausmamala.dat7 <- ausmamala.dat6[!(ausmamala.dat6$decimalLatitude >= rem.poly6.minlat &
                                     ausmamala.dat6$decimalLatitude <= rem.poly6.maxlat &
                                     ausmamala.dat6$decimalLongitude >= rem.poly6.minlon &
                                     ausmamala.dat6$decimalLongitude <= rem.poly6.maxlon),]
ausmamala.dat8 <- ausmamala.dat7[!(ausmamala.dat7$decimalLatitude >= rem.poly7.minlat &
                                     ausmamala.dat7$decimalLatitude <= rem.poly7.maxlat &
                                     ausmamala.dat7$decimalLongitude >= rem.poly7.minlon &
                                     ausmamala.dat7$decimalLongitude <= rem.poly7.maxlon),]
ausmamala.dat9 <- ausmamala.dat8[!(ausmamala.dat8$decimalLatitude >= rem.poly8.minlat & 
                                    ausmamala.dat8$decimalLatitude <= rem.poly8.maxlat &
                                    ausmamala.dat8$decimalLongitude >= rem.poly8.minlon &
                                    ausmamala.dat8$decimalLongitude <= rem.poly8.maxlon),]
ausmamala.dat10 <- as.data.frame(ausmamala.dat9[!(ausmamala.dat9$decimalLatitude >= rem.poly9.minlat &
                                      ausmamala.dat9$decimalLatitude <= rem.poly9.maxlat &
                                     ausmamala.dat9$decimalLongitude >= rem.poly9.minlon &
                                     ausmamala.dat9$decimalLongitude <= rem.poly9.maxlon),])
dim(ausmamala.dat1)
dim(ausmamala.dat10)
paste("-",round((dim(ausmamala.dat1)[1] - dim(ausmamala.dat10)[1]) / dim(ausmamala.dat1)[1] * 100, 1), "%", sep="")
head(ausmamala.dat10)

# remove temporary datasets
rm(ausmamala.dat2, ausmamala.dat3, ausmamala.dat4, ausmamala.dat5, 
   ausmamala.dat6, ausmamala.dat7, ausmamala.dat8, ausmamala.dat9)

# Aves
ausavala.dat2 <- ausavala.dat[!(ausavala.dat$decimalLatitude >= rem.poly1.minlat & 
                                    ausavala.dat$decimalLatitude <= rem.poly1.maxlat &
                                    ausavala.dat$decimalLongitude >= rem.poly1.minlon &
                                    ausavala.dat$decimalLongitude <= rem.poly1.maxlon),]
ausavala.dat3 <- ausavala.dat2[!(ausavala.dat2$decimalLatitude >= rem.poly2.minlat &
                                     ausavala.dat2$decimalLatitude <= rem.poly2.maxlat &
                                     ausavala.dat2$decimalLongitude >= rem.poly2.minlon &
                                     ausavala.dat2$decimalLongitude <= rem.poly2.maxlon),]
ausavala.dat4 <- ausavala.dat3[!(ausavala.dat3$decimalLatitude >= rem.poly3.minlat &
                                     ausavala.dat3$decimalLatitude <= rem.poly3.maxlat &
                                     ausavala.dat3$decimalLongitude >= rem.poly3.minlon &
                                     ausavala.dat3$decimalLongitude <= rem.poly3.maxlon),]
ausavala.dat5 <- ausavala.dat4[!(ausavala.dat4$decimalLatitude >= rem.poly4.minlat &
                                    ausavala.dat4$decimalLatitude <= rem.poly4.maxlat &
                                    ausavala.dat4$decimalLongitude >= rem.poly4.minlon &
                                    ausavala.dat4$decimalLongitude <= rem.poly4.maxlon),]
ausavala.dat6 <- ausavala.dat5[!(ausavala.dat5$decimalLatitude >= rem.poly5.minlat &
                                    ausavala.dat5$decimalLatitude <= rem.poly5.maxlat &
                                    ausavala.dat5$decimalLongitude >= rem.poly5.minlon &
                                    ausavala.dat5$decimalLongitude <= rem.poly5.maxlon),]
ausavala.dat7 <- ausavala.dat6[!(ausavala.dat6$decimalLatitude >= rem.poly6.minlat &
                                     ausavala.dat6$decimalLatitude <= rem.poly6.maxlat &
                                     ausavala.dat6$decimalLongitude >= rem.poly6.minlon &
                                     ausavala.dat6$decimalLongitude <= rem.poly6.maxlon),]
ausavala.dat8 <- ausavala.dat7[!(ausavala.dat7$decimalLatitude >= rem.poly7.minlat &
                                     ausavala.dat7$decimalLatitude <= rem.poly7.maxlat &
                                     ausavala.dat7$decimalLongitude >= rem.poly7.minlon &
                                     ausavala.dat7$decimalLongitude <= rem.poly7.maxlon),]
ausavala.dat9 <- ausavala.dat8[!(ausavala.dat8$decimalLatitude >= rem.poly8.minlat &
                                    ausavala.dat8$decimalLatitude <= rem.poly8.maxlat &
                                    ausavala.dat8$decimalLongitude >= rem.poly8.minlon &
                                    ausavala.dat8$decimalLongitude <= rem.poly8.maxlon),]
ausavala.dat10 <- as.data.frame(ausavala.dat9[!(ausavala.dat9$decimalLatitude >= rem.poly9.minlat &
                                      ausavala.dat9$decimalLatitude <= rem.poly9.maxlat &
                                     ausavala.dat9$decimalLongitude >= rem.poly9.minlon &
                                     ausavala.dat9$decimalLongitude <= rem.poly9.maxlon),])
dim(ausavala.dat)
dim(ausavala.dat10)
paste("-",round((dim(ausavala.dat)[1] - dim(ausavala.dat10)[1]) / dim(ausavala.dat)[1] * 100, 1), "%", sep="")
head(ausavala.dat10)

# remove temporary datasets
rm(ausavala.dat2, ausavala.dat3, ausavala.dat4, ausavala.dat5,
   ausavala.dat6, ausavala.dat7, ausavala.dat8, ausavala.dat9)


# reptiles
ausreptala.dat2 <- ausreptala.dat[!(ausreptala.dat$decimalLatitude >= rem.poly1.minlat & 
                                  ausreptala.dat$decimalLatitude <= rem.poly1.maxlat &
                                  ausreptala.dat$decimalLongitude >= rem.poly1.minlon &
                                  ausreptala.dat$decimalLongitude <= rem.poly1.maxlon),]
ausreptala.dat3 <- ausreptala.dat2[!(ausreptala.dat2$decimalLatitude >= rem.poly2.minlat &
                                   ausreptala.dat2$decimalLatitude <= rem.poly2.maxlat &
                                   ausreptala.dat2$decimalLongitude >= rem.poly2.minlon &
                                   ausreptala.dat2$decimalLongitude <= rem.poly2.maxlon),]
ausreptala.dat4 <- ausreptala.dat3[!(ausreptala.dat3$decimalLatitude >= rem.poly3.minlat &
                                   ausreptala.dat3$decimalLatitude <= rem.poly3.maxlat &
                                   ausreptala.dat3$decimalLongitude >= rem.poly3.minlon &
                                   ausreptala.dat3$decimalLongitude <= rem.poly3.maxlon),]
ausreptala.dat5 <- ausreptala.dat4[!(ausreptala.dat4$decimalLatitude >= rem.poly4.minlat &
                                   ausreptala.dat4$decimalLatitude <= rem.poly4.maxlat &
                                   ausreptala.dat4$decimalLongitude >= rem.poly4.minlon &
                                   ausreptala.dat4$decimalLongitude <= rem.poly4.maxlon),]
ausreptala.dat6 <- ausreptala.dat5[!(ausreptala.dat5$decimalLatitude >= rem.poly5.minlat &
                                   ausreptala.dat5$decimalLatitude <= rem.poly5.maxlat &
                                   ausreptala.dat5$decimalLongitude >= rem.poly5.minlon &
                                   ausreptala.dat5$decimalLongitude <= rem.poly5.maxlon),]
ausreptala.dat7 <- ausreptala.dat6[!(ausreptala.dat6$decimalLatitude >= rem.poly6.minlat &
                                   ausreptala.dat6$decimalLatitude <= rem.poly6.maxlat &
                                   ausreptala.dat6$decimalLongitude >= rem.poly6.minlon &
                                   ausreptala.dat6$decimalLongitude <= rem.poly6.maxlon),]
ausreptala.dat8 <- ausreptala.dat7[!(ausreptala.dat7$decimalLatitude >= rem.poly7.minlat &
                                   ausreptala.dat7$decimalLatitude <= rem.poly7.maxlat &
                                   ausreptala.dat7$decimalLongitude >= rem.poly7.minlon &
                                   ausreptala.dat7$decimalLongitude <= rem.poly7.maxlon),]
ausreptala.dat9 <- ausreptala.dat8[!(ausreptala.dat8$decimalLatitude >= rem.poly8.minlat &
                                   ausreptala.dat8$decimalLatitude <= rem.poly8.maxlat &
                                   ausreptala.dat8$decimalLongitude >= rem.poly8.minlon &
                                   ausreptala.dat8$decimalLongitude <= rem.poly8.maxlon),]
ausreptala.dat10 <- as.data.frame(ausreptala.dat9[!(ausreptala.dat9$decimalLatitude >= rem.poly9.minlat &
                                                  ausreptala.dat9$decimalLatitude <= rem.poly9.maxlat &
                                                  ausreptala.dat9$decimalLongitude >= rem.poly9.minlon &
                                                  ausreptala.dat9$decimalLongitude <= rem.poly9.maxlon),])
dim(ausreptala.dat)
dim(ausreptala.dat10)
paste("-",round((dim(ausreptala.dat)[1] - dim(ausreptala.dat10)[1]) / dim(ausreptala.dat)[1] * 100, 1), "%", sep="")
head(ausreptala.dat10)

# remove temporary datasets
rm(ausreptala.dat2, ausreptala.dat3, ausreptala.dat4, ausreptala.dat5,
   ausreptala.dat6, ausreptala.dat7, ausreptala.dat8, ausreptala.dat9)



# amphibians
ausamphala.dat2 <- ausamphala.dat[!(ausamphala.dat$decimalLatitude >= rem.poly1.minlat & 
                                      ausamphala.dat$decimalLatitude <= rem.poly1.maxlat &
                                      ausamphala.dat$decimalLongitude >= rem.poly1.minlon &
                                      ausamphala.dat$decimalLongitude <= rem.poly1.maxlon),]
ausamphala.dat3 <- ausamphala.dat2[!(ausamphala.dat2$decimalLatitude >= rem.poly2.minlat &
                                       ausamphala.dat2$decimalLatitude <= rem.poly2.maxlat &
                                       ausamphala.dat2$decimalLongitude >= rem.poly2.minlon &
                                       ausamphala.dat2$decimalLongitude <= rem.poly2.maxlon),]
ausamphala.dat4 <- ausamphala.dat3[!(ausamphala.dat3$decimalLatitude >= rem.poly3.minlat &
                                       ausamphala.dat3$decimalLatitude <= rem.poly3.maxlat &
                                       ausamphala.dat3$decimalLongitude >= rem.poly3.minlon &
                                       ausamphala.dat3$decimalLongitude <= rem.poly3.maxlon),]
ausamphala.dat5 <- ausamphala.dat4[!(ausamphala.dat4$decimalLatitude >= rem.poly4.minlat &
                                       ausamphala.dat4$decimalLatitude <= rem.poly4.maxlat &
                                       ausamphala.dat4$decimalLongitude >= rem.poly4.minlon &
                                       ausamphala.dat4$decimalLongitude <= rem.poly4.maxlon),]
ausamphala.dat6 <- ausamphala.dat5[!(ausamphala.dat5$decimalLatitude >= rem.poly5.minlat &
                                       ausamphala.dat5$decimalLatitude <= rem.poly5.maxlat &
                                       ausamphala.dat5$decimalLongitude >= rem.poly5.minlon &
                                       ausamphala.dat5$decimalLongitude <= rem.poly5.maxlon),]
ausamphala.dat7 <- ausamphala.dat6[!(ausamphala.dat6$decimalLatitude >= rem.poly6.minlat &
                                       ausamphala.dat6$decimalLatitude <= rem.poly6.maxlat &
                                       ausamphala.dat6$decimalLongitude >= rem.poly6.minlon &
                                       ausamphala.dat6$decimalLongitude <= rem.poly6.maxlon),]
ausamphala.dat8 <- ausamphala.dat7[!(ausamphala.dat7$decimalLatitude >= rem.poly7.minlat &
                                       ausamphala.dat7$decimalLatitude <= rem.poly7.maxlat &
                                       ausamphala.dat7$decimalLongitude >= rem.poly7.minlon &
                                       ausamphala.dat7$decimalLongitude <= rem.poly7.maxlon),]
ausamphala.dat9 <- ausamphala.dat8[!(ausamphala.dat8$decimalLatitude >= rem.poly8.minlat &
                                       ausamphala.dat8$decimalLatitude <= rem.poly8.maxlat &
                                       ausamphala.dat8$decimalLongitude >= rem.poly8.minlon &
                                       ausamphala.dat8$decimalLongitude <= rem.poly8.maxlon),]
ausamphala.dat10 <- as.data.frame(ausamphala.dat9[!(ausamphala.dat9$decimalLatitude >= rem.poly9.minlat &
                                                      ausamphala.dat9$decimalLatitude <= rem.poly9.maxlat &
                                                      ausamphala.dat9$decimalLongitude >= rem.poly9.minlon &
                                                      ausamphala.dat9$decimalLongitude <= rem.poly9.maxlon),])
dim(ausamphala.dat)
dim(ausamphala.dat10)
paste("-",round((dim(ausamphala.dat)[1] - dim(ausamphala.dat10)[1]) / dim(ausamphala.dat)[1] * 100, 1), "%", sep="")
head(ausamphala.dat10)

# remove temporary datasets
rm(ausamphala.dat2, ausamphala.dat3, ausamphala.dat4, ausamphala.dat5,
   ausamphala.dat6, ausamphala.dat7, ausamphala.dat8, ausamphala.dat9)


# convert to terra vector point object
ausmamala.vect <- terra::vect(ausmamala.dat10, geom=c("decimalLongitude", "decimalLatitude"), crs="epsg:4326")
# keep only 'recordID', 'scientificName', 'eventDate' columns
ausmamala.vect2 <- ausmamala.vect[, -c(3,5,6)]
ausmamala.vect2
plot(ausmamala.vect2, col="blue", border=NULL, cex=0.3)
plot(mainl, col="tan", alpha=0.5, border=NULL, lwd=0.5, add=TRUE)

ausavala.vect <- vect(ausavala.dat10, geom=c("decimalLongitude", "decimalLatitude"), crs="epsg:4326")
# keep only 'recordID', 'scientificName', 'eventDate' columns
ausavala.vect2 <- ausavala.vect[, -c(3,5,6)]
ausavala.vect2
plot(ausavala.vect2, col="blue", border=NULL, cex=0.3)
plot(mainl, col="tan", alpha=0.5, border=NULL, lwd=0.5, add=TRUE)

ausreptala.vect <- vect(ausreptala.dat10, geom=c("decimalLongitude", "decimalLatitude"), crs="epsg:4326")
# keep only 'recordID', 'scientificName', 'eventDate' columns
ausreptala.vect2 <- ausreptala.vect[, -c(3,5,6)]
ausreptala.vect2
plot(ausreptala.vect2, col="blue", border=NULL, cex=0.3)
plot(mainl, col="tan", alpha=0.5, border=NULL, lwd=0.5, add=TRUE)

ausamphala.vect <- vect(ausamphala.dat10, geom=c("decimalLongitude", "decimalLatitude"), crs="epsg:4326")
# keep only 'recordID', 'scientificName', 'eventDate' columns
ausamphala.vect2 <- ausamphala.vect[, -c(3,5,6)]
ausamphala.vect2
plot(ausamphala.vect2, col="blue", border=NULL, cex=0.3)
plot(mainl, col="tan", alpha=0.5, border=NULL, lwd=0.5, add=TRUE)

# save vectors to .rds files
# save rasters to RDS
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/rds")
saveRDS(ausmamala.vect2, file = "ausmamalavect2.rds")
saveRDS(ausavala.vect2, file = "ausavalavect2.rds")
saveRDS(ausreptala.vect2, file = "ausreptalavect2.rds")
saveRDS(ausamphala.vect2, file = "ausreptalavect2.rds")
saveRDS(ausislescorr, file = "ausislescorr.rds")

# write to shapefiles for faster processing in QGIS
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/gis")
writeVector(ausmamala.vect2, filename = "ausmamala.shp", overwrite = TRUE)# export to shapefile
writeVector(ausavala.vect2, filename = "ausavala.shp", overwrite = TRUE)# export to shapefile
writeVector(ausreptala.vect2, filename = "ausreptala.shp", overwrite = TRUE)# export to shapefile
writeVector(ausamphala.vect2, filename = "ausamphala.shp", overwrite = TRUE)# export to shapefile

#***************************************
## re-read shapefiles for new session
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/gis")
ausmamala.vect2 <- vect("ausmamala.shp")
ausavala.vect2 <- vect("ausavala.shp")
ausreptala.vect2 <- vect("ausreptala.shp")
ausamphala.vect2 <- vect("ausamphala.shp")
#***************************************


# once processed in QGIS, import intersected shapefiles
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/gis")
ausisles.mamala <- vect("ausislmamala.shp")
ausisles.avala <- vect("ausislavala.shp")
ausisles.reptala <- vect("ausislreptala.shp")
ausisles.amphala <- vect("ausislamphala.shp")


## summarise species richness by island
# Mammalia
ausisles.mamala.df <- as.data.frame(ausisles.mamala)
head(ausisles.mamala.df)
table(ausisles.mamala.df$scientific)

# clean up
ausisles.mamala.df1 <- ausisles.mamala.df # if including non-native mammals
## remove invasives
ausisles.mamala.df1 <- ausisles.mamala.df[!(ausisles.mamala.df$scientific %in% c("Axis axis",
                                                             "ARTIODACTYLA",
                                                             "Axis porcinus",
                                                             "Bos",
                                                             "Bos (Bos) taurus",
                                                             "BOVIDAE",
                                                             "Bubalus bubalis",
                                                             "CAMELIDAE",
                                                             "Camelus",
                                                             "Camelus dromedarius",
                                                             "CANIDAE",
                                                             "Canis",
                                                             "Canis familiaris",
                                                             "Capra",
                                                             "Capra hircus",
                                                             "Cavia porcellus",
                                                             "CERVIDAE",
                                                             "Cervus",
                                                             "Cervus elaphus",
                                                             "Cervus timorensis",
                                                             "Cervus unicolor",
                                                             "Dama",
                                                             "Dama dama",
                                                             "Dama dama dama",
                                                             "Equus",
                                                             "Equus (Asinus) asinus",
                                                             "Equus (Equus) caballus",
                                                             "FELIDAE",
                                                             "Felis",
                                                             "Felis catus",
                                                             "LAGOMORPHA",
                                                             "Lepus",
                                                             "Lepus capensis",
                                                             "Lepus capensis occidentalis",
                                                             "Mus",
                                                             "Mus musculus",
                                                             "Mus musculus domesticus",
                                                             "Mustela",
                                                             "Mustela putorius",
                                                             "Oryctolagus",
                                                             "Oryctolagus cuniculus",
                                                             "Oryctolagus cuniculus cuniculus",
                                                             "Rattus norvegicus",
                                                             "Rattus rattus",
                                                             "Sus",
                                                             "Sus scrofa",
                                                             "Vulpes",
                                                             "Vulpes vulpes")),]

## remove marine mammals
ausisles.mamala.df2 <- ausisles.mamala.df1[!(ausisles.mamala.df1$scientific %in% c("Arctocephalus",
                                                               "Arctocephalus forsteri",
                                                               "Arctocephalus gazella",
                                                               "Arctocephalus pusillus",
                                                               "Arctocephalus pusillus doriferus",
                                                               "Arctocephalus tropicalis",
                                                               "Arctophoca",
                                                               "Arctophoca forsteri",
                                                               "BALAENIDAE",
                                                               "Balaenoptera",
                                                               "Balaenoptera acutorostrata",
                                                               "Balaenoptera bonaerensis",
                                                               "Balaenoptera edeni",
                                                               "Balaenoptera musculus",
                                                               "Balaenoptera musculus brevicauda",
                                                               "Balaenoptera musculus intermedia",
                                                               "Balaenoptera omurai",
                                                               "Balaenoptera physalus",
                                                               "BALAENOPTERIDAE",
                                                               "CETACEA",
                                                               "DELPHINIDAE",
                                                               "Delphinus",
                                                               "Delphinus delphis",
                                                               "Delphinus delphis delphis",
                                                               "Dugong",
                                                               "Dugong dugon",
                                                               "DUGONGIDAE",
                                                               "Eubalaena",
                                                               "Eubalaena australis",
                                                               "Globicephala",
                                                               "Globicephala macrorhynchus",
                                                               "Globicephala melas",
                                                               "Globicephala melas edwardii",
                                                               "Grampus griseus",
                                                               "Hydrurga",
                                                               "Hydrurga leptonyx",
                                                               "Hyperoodon planifrons",
                                                               "Kogia",
                                                               "Kogia breviceps",
                                                               "Kogia sima",
                                                               "Lagenodelphis hosei",
                                                               "Lagenorhynchus",
                                                               "Lobodon carcinophaga",
                                                               "Lagenorhynchus obscurus",
                                                               "Leptonychotes weddelli",
                                                               "Lissodelphis peronii",
                                                               "Megaptera",
                                                               "Megaptera novaeangliae",
                                                               "Mesoplodon bowdoini",
                                                               "Mesoplodon grayi",
                                                               "Mesoplodon hectori",
                                                               "Mesoplodon layardii",
                                                               "Mesoplodon",
                                                               "Mirounga",
                                                               "Mirounga leonina",
                                                               "Mesoplodon densirostris",
                                                               "MYSTICETI",
                                                               "Neophoca",
                                                               "Neophoca cinerea",
                                                               "ODONTOCETI",
                                                               "Orcinus orca",
                                                               "Ommatophoca rossi",
                                                               "Orcinus orca",
                                                               "OTARIIDAE",
                                                               "Phocoena phocoena",
                                                               "Phocoena dioptrica",
                                                               "Phocaena dioptrica",
                                                               "Physeter macrocephalus",
                                                               "Pseudorca crassidens",
                                                               "Physeter",
                                                               "Sousa",
                                                               "Sousa sahulensis",
                                                               "Stenella attenuata",
                                                               "Stenella coeruleoalba",
                                                               "Stenella longirostris",
                                                               "Stenella longirostris roseiventris",
                                                               "Stenella",
                                                               "Steno",
                                                               "Steno bredanensis",
                                                               "Tasmacetus shepherdi",
                                                               "Tursiops",
                                                               "Tursiops aduncus",
                                                               "Tursiops truncatus",
                                                               "Tursiops australis",
                                                               "ZIPHIIDAE",
                                                               "Ziphius cavirostris")),]

# remove higher-taxonomic level redundancies
ausisles.mamala.df3 <- ausisles.mamala.df2[!(ausisles.mamala.df2$scientific %in% c("BURRAMYIDAE",
                                                               "CARNIVORA",
                                                               "CHIROPTERA",
                                                               "DASYURIDAE",
                                                               "DIPROTODONTIA",
                                                               "LEPORIDAE",
                                                               "MACROPODIDAE",
                                                               "MACROPODINAE",
                                                               "MAMMALIA",
                                                               "MARSUPIALIA",
                                                               "MICROCHIROPTERA",
                                                               "MOLOSSIDAE",
                                                               "MONOTREMATA",
                                                               "MURIDAE",
                                                               "PERAMELEMORPHIA",
                                                               "PERAMELINAE",
                                                               "PERAMELIDAE",
                                                               "PETAURIDAE",
                                                               "PETAUROIDEA",
                                                               "PHALANGERIDAE",
                                                               "POTOROIDAE",
                                                               "PROTOTHERIA",
                                                               "PSEUDOCHEIRIDAE",
                                                               "PTEROPODIDAE",
                                                               "PTEROPODINAE",
                                                               "RATTUS GROUP",
                                                               "RHINOLOPHOIDEA",
                                                               "RODENTIA",
                                                               "VESPERTILIONIDAE",
                                                               "VESPERTILIONINAE",
                                                               "YANGOCHIROPTERA")),]

# remove single genera that cannot be definitively attributed to a single species
ausisles.mamala.df4 <- ausisles.mamala.df3[!(ausisles.mamala.df3$scientific %in% c("Antechinus",
                                                               "Bettongia",
                                                               "Cercartetus",
                                                               "Chalinolobus",
                                                               "Dasyurus",
                                                               "Hipposideros",
                                                               "Isoodon",
                                                               "Melomys",
                                                               "Macropus",
                                                               "Miniopterus",
                                                               "Myotis",
                                                               "Notomys",
                                                               "Notamacropus",
                                                               "Nyctophilus",
                                                               "Osphranter",
                                                               "Ozimops",
                                                               "Perameles",
                                                               "Petaurus",
                                                               "Petrogale",
                                                               "Phascogale",
                                                               "Pipistrellus",
                                                               "Planigale",
                                                               "Pseudomys",
                                                               "Potorous",
                                                               "Pteropus",
                                                               "Rattus",
                                                               "Scotorepens",
                                                               "Sminthopsis",
                                                               "Taphozous",
                                                               "Thylogale",
                                                               "Vespadelus",
                                                               "Zyzomys")),]

# remove extinct species
ausisles.mamala.df5 <- ausisles.mamala.df4[!(ausisles.mamala.df4$scientific %in% c("Thylacinus cynocephalus")),]

# combine taxonomic names (including subspecies) and rename into single species
ausisles.mamala.df6 <- ausisles.mamala.df5 %>%
  mutate(scientific = recode(scientific,
                                 "Antechinus flavipes flavipes" = "Antechinus flavipes",
                                 "Antechinus flavipes leucogaster" = "Antechinus flavipes",
                                 "Antechinus minimus maritimus" = "Antechinus minimus",
                                 "Antechinus mimetes mimetes" = "Antechinus mimetes",
                                 "Antechinus minimus minimus" = "Antechinus minimus",
                                 "Bettongia lesueur graii" = "Bettongia lesueur",
                                 "Bettongia lesueur lesueur" = "Bettongia lesueur",
                                 "Bettongia penicillata penicillata" = "Bettongia penicillata",
                                 "Bettongia gaimardi cuniculus" = "Bettongia gaimardi",
                                 "Cercartetus nanus nanus" = "Cercartetus nanus",
                                 "Chalinolobus nigrogriseus rogersi" = "Chalinolobus nigrogriseus",
                                 "Dasyurus maculatus maculatus" = "Dasyurus maculatus",
                                 "Hipposideros ater aruensis" = "Hipposideros ater",
                                 "Hipposideros ater gilberti" = "Hipposideros ater",
                                 "Hipposideros cervinus cervinus" = "Hipposideros cervinus",
                                 "Hipposideros diadema reginae" = "Hipposideros diadema",
                                 "Hydromys" = "Hydromys chrysogaster",
                                 "Isoodon auratus auratus" = "Isoodon auratus",
                                 "Isoodon auratus barrowensis" = "Isoodon auratus", 
                                 "Isoodon obesulus nauticus" = "Isoodon obesulus",
                                 "Isoodon obesulus obesulus" = "Isoodon obesulus",
                                 "Isoodon obesulus affinis" = "Isoodon obesulus",
                                 "Isoodon macrourus torosus" = "Isoodon macrourus",
                                 "Isoodon macrourus macrourus" = "Isoodon macrourus",
                                 "Lagorchestes" = "Lagorchestes hirsutus",
                                 "Lagorchestes hirsutus bernieri" = "Lagorchestes hirsutus",
                                 "Lagostrophus fasciatus fasciatus" = "Lagostrophus fasciatus",
                                 "Lagorchestes conspicillatus conspicillatus" = "Lagorchestes conspicillatus",
                                 "Lasiorhinus" = "Lasiorhinus latifrons",
                                 "Leporillus" = "Leporillus conditor",
                                 "Macropus fuliginosus fuliginosus" = "Macropus fuliginosus",
                                 "Macropus fuliginosus melanops" = "Macropus fuliginosus",
                                 "Macrotis" = "Macrotis lagotis",
                                 "Mastacomys" = "Mastacomys fuscus",
                                 "Mastacomys fuscus fuscus" = "Mastacomys fuscus",
                                 "Mastacomys fuscus mordicus" = "Mastacomys fuscus",
                                 "Mesembriomys gouldii melvillensis" = "Mesembriomys gouldii",
                                 "Miniopterus orianae orianae" = "Miniopterus orianae",
                                 "Miniopterus orianae oceanensis" = "Miniopterus orianae",
                                 "Myrmecobius fasciatus rufus" = "Myrmecobius fasciatus",
                                 "Ningaui" = "Ningaui yvonneae",
                                 "Notamacropus agilis agilis" = "Notamacropus agilis",
                                 "Notamacropus agilis jardinii" = "Notamacropus agilis",
                                 "Notamacropus agilis nigrescens" = "Notamacropus agilis",
                                 "Notamacropus eugenii eugenii" = "Notamacropus eugenii",
                                 "Notamacropus eugenii derbianus" = "Notamacropus eugenii",
                                 "Notamacropus rufogriseus banksianus" = "Notamacropus rufogriseus",
                                 "Notamacropus rufogriseus rufogriseus" = "Notamacropus rufogriseus",
                                 "Nyctimene" = "Nyctimene robinsoni",
                                 "Nyctophilus major tor" = "Nyctophilus major",
                                 "Nyctophilus geoffroyi pacificus" = "Nyctophilus geoffroyi",
                                 "Osphranter robustus erubescens" = "Osphranter robustus",
                                 "Osphranter robustus isabellinus" = "Osphranter robustus",
                                 "Osphranter robustus robustus" = "Osphranter robustus",
                                 "Perameles gunnii gunnii" = "Perameles gunnii", 
                                 "Petaurus breviceps breviceps" = "Petaurus breviceps",
                                 "Petrogale lateralis pearsoni" = "Petrogale lateralis",
                                 "Petrogale lateralis lateralis" = "Petrogale lateralis",
                                 "Petrogale concinna monastria" = "Petrogale concinna",
                                 "Petrogale xanthopus xanthopus" = "Petrogale xanthopus",
                                 "Phascogale tapoatafa tapoatafa" = "Phascogale tapoatafa",
                                 "PHASCOLARCTIDAE" = "Phascolarctos cinereus",
                                 "Phascolarctos" = "Phascolarctos cinereus",
                                 "Planigale maculata maculata" = "Planigale maculata",
                                 "Potorous tridactylus apicalis" = "Potorous tridactylus",
                                 "Potorous tridactylus tridactylus" = "Potorous tridactylus",
                                 "Potorous tridactylus trisulcatus" = "Potorous tridactylus",
                                 "Pseudocheirus" = "Pseudocheirus peregrinus",
                                 "Pseudocheirus peregrinus peregrinus" = "Pseudocheirus peregrinus",
                                 "Pseudocheirus peregrinus convolutor" = "Pseudocheirus peregrinus",
                                 "Pseudocheirus peregrinus cooki" = "Pseudocheirus peregrinus",
                                 "Pseudochirops" = "Pseudochirops archeri",
                                 "Pseudomys albocinereus squalorum" = "Pseudomys albocinereus",
                                 "Pteropus alecto gouldii" = "Pteropus alecto",
                                 "Rattus fuscipes assimilis" = "Rattus fuscipes",
                                 "Rattus fuscipes greyi" = "Rattus fuscipes",
                                 "Rattus fuscipes fuscipes" = "Rattus fuscipes",
                                 "Rattus lutreolus lutreolus" = "Rattus lutreolus",
                                 "Rattus lutreolus velutinus" = "Rattus lutreolus",
                                 "Rattus tunneyi culmorum" = "Rattus tunneyi",
                                 "Rattus tunneyi tunneyi" = "Rattus tunneyi",
                                 "Saccolaimus saccolaimus nudicluniatus" = "Saccolaimus saccolaimus",
                                 "Sarcophilus" = "Sarcophilus harrisii",
                                 "Sminthopsis crassicaudata crassicaudata" = "Sminthopsis crassicaudata",
                                 "Sminthopsis fuliginosa aitkeni" = "Sminthopsis fuliginosa",
                                 "Sminthopsis fuliginosa fuliginosa" = "Sminthopsis fuliginosa",
                                 "Sminthopsis fuliginosus aitkeni" = "Sminthopsis fuliginosa",
                                 "Sminthopsis leucopus leucopus" = "Sminthopsis leucopus",
                                 "Sminthopsis murina murina" = "Sminthopsis murina",
                                 "Syconycteris australis australis" = "Syconycteris australis",
                                 "TACHYGLOSSIDAE" = "Tachyglossus aculeatus",
                                 "Tachyglossus" = "Tachyglossus aculeatus",
                                 "Tachyglossus aculeatus acanthion" = "Tachyglossus aculeatus",
                                 "Tachyglossus aculeatus setosus" = "Tachyglossus aculeatus",
                                 "Tachyglossus aculeatus acanthion" = "Tachyglossus aculeatus",
                                 "Tachyglossus aculeatus aculeatus" = "Tachyglossus aculeatus",
                                 "Tachyglossus aculeatus multiaculeatus" = "Tachyglossus aculeatus",
                                 "Thylogale stigmatica wilcoxi" = "Thylogale stigmatica",
                                 "Trichosurus" = "Trichosurus vulpecula",
                                 "Trichosurus vulpecula vulpecula" = "Trichosurus vulpecula",
                                 "Trichosurus vulpecula arnhemensis" = "Trichosurus vulpecula",
                                 "Trichosurus vulpecula fuliginosus" = "Trichosurus vulpecula",
                                 "Trichosurus vulpecula hypoleucus" = "Trichosurus vulpecula",
                                 "Trichosurus vulpecula johnstoni" = "Trichosurus vulpecula",
                                 "Vombatus ursinus hirsutus" = "Vombatus ursinus",
                                 "Vombatus" = "Vombatus ursinus",
                                 "Vombatus ursinus ursinus" = "Vombatus ursinus",
                                 "Vombatus ursinus tasmaniensis" = "Vombatus ursinus",
                                 "VOMBATIDAE" = "Vombatus ursinus",
                                 "Wallabia" = "Wallabia bicolor",
                                 "Wallabia bicolor bicolor" = "Wallabia bicolor",
                                 "Wyulda" = "Wyulda squamicaudata"))

# check species
table(ausisles.mamala.df6$scientific)

head(ausisles.mamala.df6)
ausisles.mamala.spp <- ausisles.mamala.df6[,c("scientific", "FID")]
head(ausisles.mamala.spp)
ausisles.mamala.spp.unique <- unique(ausisles.mamala.spp)
head(ausisles.mamala.spp.unique)
ausisles.mamala.spp.unique.tab <- as.data.frame(table(ausisles.mamala.spp.unique$FID)) # number of unique species per island)
ausisles.mamala.spp.unique.tab
head(ausisles.mamala.spp.unique.tab)
colnames(ausisles.mamala.spp.unique.tab) <- c("FID", "mamSR")
head(ausisles.mamala.spp.unique.tab)

# Aves
ausisles.avala.df <- as.data.frame(ausisles.avala)
head(ausisles.avala.df)
table(ausisles.avala.df$scientific)

## species list
avspp <- table(ausisles.avala.df$scientific)
avspp

avala.dat1 <- ausisles.avala.df# if keeping invasives
## remove invasives
avala.dat1 <- ausisles.avala.df[!(ausisles.avala.df$scientific %in% c("Agapornis",
                                                          "Agapornis roseicollis",
                                                          "Alauda",
                                                          "Alauda arvensis",
                                                          "Alauda arvensis arvensis",
                                                          "Acridotheres tristis",
                                                          "Acridotheres tristis tristis",
                                                          "Anas (Anas) platyrhynchos",
                                                          "Anas superciliosa X platyrhynchos",
                                                          "Branta canadensis",
                                                          "Branta canadensis maxima",
                                                          "Phasianus colchicus",
                                                          "Carduelis",
                                                          "Carduelis carduelis",
                                                          "Carduelis carduelis britannica",
                                                          "Columba (Columba) livia",
                                                          "Columba (Columba) livia livia",
                                                          "Gallus gallus",
                                                          "Passer (Passer) domesticus",
                                                          "Passer (Passer) domesticus domesticus",
                                                          "Passer (Passer) montanus",
                                                          "Streptopelia risoria",
                                                          "Streptopelia",
                                                          "Sturnus",
                                                          "Sturnus (Sturnus) vulgaris",
                                                          "Sturnus (Sturnus) vulgaris vulgaris")),]

## remove vagrants
avala.dat2 <- avala.dat1[!(avala.dat1$scientific %in% c("Eudyptes",
                                                            "Eudyptes chrysocome",
                                                            "Eudyptes chrysolophus",
                                                            "Eudyptes chrysolophus schlegeli",
                                                            "Eudyptes pachyrhynchus")),]

## remove higher-taxonomic level ambiguiies
avala.dat3 <- avala.dat2[!(avala.dat2$scientific %in% c("ACCIPITRIDAE",
                                                            "ACCIPITRIFORMES",
                                                            "Accipitrinae",
                                                            "ACANTHIZIDAE",
                                                            "ACROCEPHALIDAE",
                                                            "ALCEDINIDAE",
                                                            "ANATIDAE",
                                                            "ANSERIFORMES",
                                                            "APODIDAE",
                                                            "ARDEIDAE",
                                                            "Ardeinae",
                                                            "ARTAMIDAE",
                                                            "AVES",
                                                            "CACATUIDAE",
                                                            "CAMPEPHAGIDAE",
                                                            "CAPRIMULGIDAE",
                                                            "CHARADRIIDAE",
                                                            "CHARADRIIFORMES",
                                                            "CLIMACTERIDAE",
                                                            "COLUMBIDAE",
                                                            "COLUMBIFORMES",
                                                            "CORACIIFORMES",
                                                            "COLUMBIFORMES",
                                                            "COLUMBIDAE",
                                                            "CORCORACIDAE",
                                                            "CORVIDAE",
                                                            "CUCULIDAE",
                                                            "Cuculinae",
                                                            "CUCULIFORMES",
                                                            "DIOMEDEIDAE",
                                                            "ESTRILDIDAE",
                                                            "FALCONIDAE",
                                                            "FALCONIFORMES",
                                                            "FREGATIDAE",
                                                            "FRINGILLIDAE",
                                                            "GALLIFORMES",
                                                            "GRUIFORMES",
                                                            "HAEMATOPODIDAE",
                                                            "HIRUNDINIDAE",
                                                            "HYDROBATIDAE",
                                                            "LARIDAE",
                                                            "Larinae",
                                                            "LOCUSTELLIDAE",
                                                            "MALURIDAE",
                                                            "MELIPHAGIDAE",
                                                            "MONARCHIDAE",
                                                            "MOTACILLIDAE",
                                                            "NECTARINIIDAE",
                                                            "OCEANITIDAE",
                                                            "ORIOLIDAE",
                                                            "PARADISAEIDAE",
                                                            "PARDALOTIDAE",
                                                            "PASSERIFORMES",
                                                            "PELECANIDAE",
                                                            "PELECANIFORMES",
                                                            "PETROICIDAE",
                                                            "PHALACROCORACIDAE",
                                                            "PHASIANIDAE",
                                                            "PHYLLOSCOPIDAE",
                                                            "Platycercinae",
                                                            "Ploceidae",
                                                            "PODICIPEDIDAE",
                                                            "PODICIPEDIFORMES",
                                                            "POMATOSTOMIDAE",
                                                            "PROCELLARIIDAE",
                                                            "PROCELLARIIFORMES",
                                                            "PSITTACIDAE",
                                                            "PSITTACIFORMES",
                                                            "RALLIDAE",
                                                            "Rallinae",
                                                            "RHIPIDURIDAE",
                                                            "SCOLOPACIDAE",
                                                            "SPHENISCIDAE",
                                                            "STERCORARIIDAE",
                                                            "Sterninae",
                                                            "STRIGIFORMES",
                                                            "STURNIDAE",
                                                            "SULIDAE",
                                                            "THRESKIORNITHIDAE",
                                                            "Threskiornithinae",
                                                            "TURNICIDAE",
                                                            "TYTONIDAE",
                                                            "ZOSTEROPIDAE")),]

## simplify species names (& treat subspecies as species)
avala.dat4 <- avala.dat3 %>%
  mutate(scientific = recode(scientific,
                                 "Acanthiza (Acanthiza) apicalis" = "Acanthiza apicalis",
                                 "Acanthiza (Acanthiza) apicalis albiventris" = "Acanthiza apicalis",
                                 "Acanthiza (Acanthiza) apicalis apicalis" = "Acanthiza apicalis",
                                 "Acanthiza (Acanthiza) apicalis whitlocki" = "Acanthiza apicalis",
                                 "Acanthiza (Acanthiza) katherina" = "Acanthiza katherina",
                                 "Acanthiza (Acanthiza) pusilla" = "Acanthiza pusilla",
                                 "Acanthiza (Acanthiza) pusilla pusilla" = "Acanthiza pusilla",
                                 "Acanthiza (Acanthiza) pusilla samueli" = "Acanthiza pusilla",
                                 "Acanthiza (Acanthiza) pusilla zietzi" = "Acanthiza pusilla",
                                 "Acanthiza (Geobasileus) chrysorrhoa" = "Acanthiza chrysorrhoa",
                                 "Acanthiza (Geobasileus) chrysorrhoa leighi" = "Acanthiza chrysorrhoa",
                                 "Acanthiza (Geobasileus) inornata" = "Acanthiza inornata",
                                 "Acanthiza (Geobasileus) iredalei" = "Acanthiza iredalei",
                                 "Acanthiza (Geobasileus) iredalei hedleyi" = "Acanthiza iredalei",
                                 "Acanthiza (Geobasileus) iredalei iredalei" = "Acanthiza iredalei",
                                 "Acanthiza (Geobasileus) iredalei rosinae" = "Acanthiza iredalei",
                                 "Acanthiza (Geobasileus) reguloides" = "Acanthiza reguloides",
                                 "Acanthiza (Geobasileus) reguloides australis" = "Acanthiza reguloides",
                                 "Acanthiza (Geobasileus) reguloides reguloides" = "Acanthiza reguloides",
                                 "Acanthiza (Geobasileus) uropygialis" = "Acanthiza uropygialis",
                                 "Acanthiza (Subacanthiza) lineata" = "Acanthiza lineata",
                                 "Acanthiza (Subacanthiza) lineata alberti" = "Acanthiza lineata",
                                 "Acanthiza (Subacanthiza) lineata clelandi" = "Acanthiza lineata",
                                 "Acanthiza (Subacanthiza) lineata lineata" = "Acanthiza lineata",
                                 "Acanthiza (Subacanthiza) nana" = "Acanthiza nana",
                                 "Acanthiza (Subacanthiza) nana modesta" = "Acanthiza nana",
                                 "Acanthiza (Subacanthiza) nana nana" = "Acanthiza nana",
                                 'Acanthiza (Acanthiza) pusilla dawsonensis' = 'Acanthiza pusilla',
                                 "Acanthorhynchus tenuirostris halmaturinus" = "Acanthorhynchus tenuirostris",
                                 "Acanthorhynchus tenuirostris tenuirostris" = "Acanthorhynchus tenuirostris",
                                 "Accipiter (Leucospiza) fasciatus" = "Accipiter fasciatus",
                                 "Accipiter (Leucospiza) fasciatus didimus" = "Accipiter fasciatus",
                                 "Accipiter (Leucospiza) fasciatus fasciatus" = "Accipiter fasciatus",
                                 "Accipiter (Leucospiza) novaehollandiae" = "Accipiter novaehollandiae",
                                 "Accipiter (Paraspizias) cirrocephalus" = "Accipiter cirrocephalus",
                                 "Accipiter (Paraspizias) cirrocephalus cirrocephalus" = "Accipiter cirrocephalus",
                                 "Acrocephalus (Acrocephalus) australis" = "Acrocephalus australis",
                                 "Acrocephalus (Acrocephalus) australis australis" = "Acrocephalus australis",
                                 'Acrocephalus (Acrocephalus) orientalis' = 'Acrocephalus orientalis',
                                 "Aegotheles" = "Aegotheles cristatus",
                                 "Aegotheles (Aegotheles) cristatus" = "Aegotheles cristatus",
                                 "Aegotheles (Aegotheles) cristatus cristatus" = "Aegotheles cristatus",
                                 'Aerodramus' = 'Aerodramus terraereginae',
                                 'Aerodramus terraereginae terraereginae' = 'Aerodramus terraereginae',
                                 "Alectura lathami lathami" = "Alectura lathami",
                                 'Alectura lathami purpureicollis' = 'Alectura lathami',
                                 'Alisterus' = 'Alisterus scapularis',
                                 'Amaurornis moluccana ruficrissa' = 'Amaurornis ruficrissa',
                                 'Amytornis (Amytornis) housei' = 'Amytornis housei',
                                 "Amytornis (Amytornis) modestus" = "Amytornis modestus",
                                 "Amytornis (Amytornis) modestus curnamona" = "Amytornis modestus",
                                 "Amytornis (Amytornis) purnelli" = "Amytornis purnelli",
                                 "Amytornis (Amytornis) textilis" = "Amytornis textilis",
                                 "Amytornis (Amytornis) textilis myall" = "Amytornis textilis",
                                 "Amytornis (Amytornis) textilis textilis" = "Amytornis textilis",
                                 "Amytornis (Cryptamytis) merrotsyi" = "Amytornis merrotsyi",
                                 "Amytornis (Cryptamytis) merrotsyi merrotsyi" = "Amytornis merrotsyi",
                                 "Amytornis (Cryptamytis) merrotsyi pedleri" = "Amytornis merrotsyi",
                                 "Amytornis (Magnamytis) striatus" = "Amytornis striatus",
                                 "Amytornis (Magnamytis) striatus howei" = "Amytornis striatus",
                                 "Amytornis (Magnamytis) striatus striatus" = "Amytornis striatus",
                                 "Amytornis (Magnamytis) whitei aenigma" = "Amytornis whitei",
                                 'Anas (Dafila) eatoni' = 'Anas eatoni',
                                 "Anas (Anas) superciliosa" = "Anas superciliosa",
                                 "Anas (Anas) superciliosa superciliosa" = "Anas superciliosa",
                                 'Anas superciliosa X Anas platyrhynchos' = 'Anas superciliosa',
                                 "Anas (Nettion) castanea" = "Anas castanea",
                                 "Anas (Nettion) gracilis" = "Anas gracilis",
                                 "Anas (Nettion) gracilis gracilis" = "Anas gracilis",
                                 "Anhinga" = "Anhinga novaehollandiae",
                                 "Anhinga novaehollandiae novaehollandiae" = "Anhinga novaehollandiae",
                                 'Anous minutus minutus' = 'Anous minutus',
                                 'Anous stolidus pileatus' = 'Anous stolidus',
                                 "Anseranas" = "Anseranas semipalmata",
                                 "Anthochaera (Anellobia) chrysoptera" = "Anthochaera chrysoptera",
                                 "Anthochaera (Anellobia) chrysoptera chrysoptera" = "Anthochaera chrysoptera",
                                 "Anthochaera (Anellobia) chrysoptera halmaturina" = "Anthochaera chrysoptera",
                                 "Anthochaera (Anellobia) chrysoptera tasmanica" = "Anthochaera chrysoptera",
                                 "Anthochaera (Anellobia) lunulata" = "Anthochaera lunulata",
                                 "Anthochaera (Anthochaera) carunculata" = "Anthochaera carunculata",
                                 "Anthochaera (Anthochaera) carunculata carunculata" = "Anthochaera carunculata",
                                 "Anthochaera (Anthochaera) carunculata clelandi" = "Anthochaera carunculata",
                                 "Anthochaera (Anthochaera) carunculata woodwardi" = "Anthochaera carunculata",
                                 "Anthochaera (Anthochaera) paradoxa" = "Anthochaera paradoxa",
                                 "Anthochaera (Xanthomyza) phrygia" = "Anthochaera phrygia",
                                 "Anthus" = "Anthus novaeseelandiae",
                                 "Anthus (Anthus) novaeseelandiae" = "Anthus novaeseelandiae",
                                 "Anthus (Anthus) novaeseelandiae bilbali" = "Anthus novaeseelandiae",
                                 "Anthus (Anthus) novaeseelandiae novaeseelandiae" = "Anthus novaeseelandiae",
                                 "Aphelocephala leucopsis leucopsis" = "Aphelocephala leucopsis",
                                 'Aplonis (Aplonis) fusca hulliana' = 'Aplonis fusca',
                                 'Aplonis (Lamprocorax) metallica metallica' = 'Aplonis metallica',
                                 'Aplonis (Lamprocorax) metallica' = 'Aplonis metallica',
                                 'Aprosmictus erythropterus coccineopterus' = 'Aprosmictus erythropterus',
                                 "Apus" = "Apus pacificus",
                                 "Apus (Apus) pacificus" = "Apus pacificus",
                                 "Apus (Apus) pacificus pacificus" = "Apus pacificus",
                                 "Aquila" = "Aquila audax",
                                 "Aquila (Uroaetus) audax" = "Aquila audax",
                                 "Aquila (Uroaetus) audax audax" = "Aquila audax",
                                 "Ardea alba modesta" = "Ardea alba",
                                 "Ardea intermedia plumifera" = "Ardea intermedia",
                                 "Arenaria" = "Arenaria interpres",
                                 "Arenaria interpres interpres" = "Arenaria interpres",
                                 "Artamus (Angroyan) cinereus" = "Artamus cinereus",
                                 "Artamus (Angroyan) cinereus cinereus" = "Artamus cinereus",
                                 "Artamus (Angroyan) cinereus melanops" = "Artamus cinereus",
                                 "Artamus (Angroyan) cyanopterus" = "Artamus cyanopterus",
                                 "Artamus (Angroyan) cyanopterus cyanopterus" = "Artamus cyanopterus",
                                 "Artamus (Angroyan) cyanopterus perthi" = "Artamus cyanopterus",
                                 "Artamus (Angroyan) minor" = "Artamus minor",
                                 'Artamus (Angroyan) minor derbyi' = 'Artamus minor',
                                 'Artamus (Angroyan) minor minor' = 'Artamus minor',
                                 "Artamus (Artamus) leucorynchus" = "Artamus leucorynchus",
                                 'Artamus (Artamus) leucorynchus leucopygialis' = 'Artamus leucorynchus',
                                 "Artamus (Campbellornis) personatus" = "Artamus personatus",
                                 "Artamus (Campbellornis) superciliosus" = "Artamus superciliosus",
                                 'Aviceda (Aviceda) subcristata' = 'Aviceda subcristata',
                                 "Aythya (Nyroca) australis" = "Aythya australis",
                                 "Barnardius" = "Barnardius zonarius",
                                 "Barnardius zonarius barnardi" = "Barnardius zonarius",
                                 "Barnardius zonarius semitorquatus" = "Barnardius zonarius",
                                 "Barnardius zonarius zonarius" = "Barnardius zonarius",
                                 "Biziura lobata menziesi" = "Biziura lobata",
                                 "Botaurus" = "Botaurus poiciloptilus",
                                 "Bubulcus" = "Bubulcus ibis",
                                 "Bubulcus ibis coromandus" = "Bubulcus ibis",
                                 "Burhinus" = "Burhinus grallarius",
                                 "Burhinus (Burhinus)" = "Burhinus grallarius",
                                 "Burhinus (Burhinus) grallarius" = "Burhinus grallarius",
                                 'Butorides striata macrorhyncha' = 'Butorides striata',
                                 'Butorides striata rogersi' = 'Butorides striata',
                                 'Butorides striata stagnatilis' = 'Butorides striata',
                                 "Cacatua (Cacatua) galerita" = "Cacatua galerita",
                                 'Cacatua (Cacatua) galerita fitzroyi' = 'Cacatua galerita',
                                 "Cacatua (Cacatua) galerita galerita" = "Cacatua galerita",
                                 "Cacatua (Licmetis) sanguinea" = "Cacatua sanguinea",
                                 "Cacatua (Licmetis) sanguinea gymnopis" = "Cacatua sanguinea",
                                 "Cacatua (Licmetis) sanguinea sanguinea" = "Cacatua sanguinea",
                                 "Cacatua (Licmetis) sanguinea westralensis" = "Cacatua sanguinea",
                                 "Cacatua (Licmetis) tenuirostris" = "Cacatua tenuirostris",
                                 "Cacomantis (Cacomantis) variolosus" = "Cacomantis variolosus",
                                 "Cacomantis (Cacomantis) variolosus variolosus" = "Cacomantis variolosus",
                                 "Cacomantis (Vidgenia) castaneiventris" = "Cacomantis castaneiventris",
                                 "Cacomantis (Vidgenia) flabelliformis" = "Cacomantis flabelliformis",
                                 "Cacomantis (Vidgenia) flabelliformis flabelliformis" = "Cacomantis flabelliformis",
                                 'Cacomantis (Cacomantis) variolosus dumetorum' = 'Cacomantis variolosus',
                                 "Calamanthus campestris campestris" = "Calamanthus campestris",
                                 'Calamanthus campestris dorrie' = 'Calamanthus campestris',
                                 "Calamanthus campestris isabellinus" = "Calamanthus campestris",
                                 'Calamanthus campestris hartogi' = 'Calamanthus campestris',
                                 "Calamanthus campestris winiam" = "Calamanthus campestris",
                                 "Calamanthus fuliginosus bourneorum" = "Calamanthus fuliginosus",
                                 "Calidris (Calidris) canutus" = "Calidris canutus",
                                 "Calidris (Calidris) canutus rogersi" = "Calidris canutus",
                                 "Calidris (Calidris) falcinellus" = "Calidris falcinellus",
                                 "Calidris (Calidris) falcinellus sibirica" = "Calidris falcinellus",
                                 "Calidris (Calidris) tenuirostris" = "Calidris tenuirostris",
                                 "Calidris (Crocethia) alba" = "Calidris alba",
                                 "Calidris (Crocethia) alba alba" = "Calidris alba",
                                 "Calidris (Ereunetes) minuta" = "Calidris minuta",
                                 "Calidris (Ereunetes) subminuta" = "Calidris subminuta",
                                 "Calidris (Ereunetes) ruficollis" = "Calidris ruficollis",
                                 "Calidris (Erolia) acuminata" = "Calidris acuminata",
                                 "Calidris (Erolia) ferruginea" = "Calidris ferruginea",
                                 "Calidris (Erolia) melanotos" = "Calidris melanotos",
                                 "Caligavis chrysops chrysops" = "Caligavis chrysops",
                                 "Caligavis chrysops samueli" = "Caligavis chrysops",
                                 "Calyptorhynchus (Calyptorhynchus) banksii" = "Calyptorhynchus banksii",
                                 'Calyptorhynchus (Calyptorhynchus) banksii banksii' = 'Calyptorhynchus banksii',
                                 "Calyptorhynchus (Calyptorhynchus) banksii graptogyne" = "Calyptorhynchus banksii",
                                 "Calyptorhynchus (Calyptorhynchus) banksii naso" = "Calyptorhynchus banksii",
                                 "Calyptorhynchus (Calyptorhynchus) lathami" = "Calyptorhynchus lathami",
                                 'Calyptorhynchus (Calyptorhynchus) lathami erebus' = 'Calyptorhynchus lathami',
                                 "Calyptorhynchus (Calyptorhynchus) lathami halmaturinus" = "Calyptorhynchus lathami",
                                 'Calyptorhynchus (Calyptorhynchus) lathami lathami' = 'Calyptorhynchus lathami',
                                 "Casuarius" = "Casuarius casuarius",
                                 'Casuarius casuarius johnsonii (southern population)' = 'Casuarius casuarius',
                                 'Centropus' = 'Centropus phasianinus',
                                 'Centropus phasianinus melanurus' = 'Centropus phasianinus',
                                 'Centropus phasianinus phasianinus' = 'Centropus phasianinus',
                                 "Cereopsis" = "Cereopsis novaehollandiae",
                                 "Cereopsis novaehollandiae novaehollandiae" = "Cereopsis novaehollandiae",
                                 "Certhionyx" = "Certhionyx variegatus",
                                 "Certhionyx (Certhionyx) variegatus" = "Certhionyx variegatus",
                                 'Ceyx azureus azureus' = 'Ceyx azureus',
                                 'Ceyx azureus ruficollaris' = 'Ceyx azureus',
                                 "Chalcites lucidus plagosus" = "Chalcites lucidus",
                                 'Chalcites lucidus lucidus' = 'Chalcites lucidus',
                                 'Chalcites minutillus barnardi' = 'Chalcites minutillus',
                                 'Chalcites minutillus minutillus' = 'Chalcites minutillus',
                                 'Chalcites minutillus russatus' = 'Chalcites minutillus',
                                 'Chalcophaps longirostris longirostris' = 'Chalcophaps longirostris',
                                 'Chalcophaps longirostris rogersi' = 'Chalcophaps longirostris',
                                 "Charadrius (Charadrius) bicinctus" = "Charadrius bicinctus",
                                 "Charadrius (Charadrius) bicinctus bicinctus" = "Charadrius bicinctus",
                                 "Charadrius (Charadrius) dubius" = "Charadrius dubius",
                                 "Charadrius (Charadrius) dubius dubius" = "Charadrius dubius",
                                 "Charadrius (Charadrius) leschenaultii" = "Charadrius leschenaultii",
                                 "Charadrius (Charadrius) leschenaultii leschenaultii" = "Charadrius leschenaultii",
                                 "Charadrius (Charadrius) mongolus" = "Charadrius mongolus",
                                 "Charadrius (Charadrius) mongolus mongolus" = "Charadrius mongolus",
                                 "Charadrius (Charadrius) ruficapillus" = "Charadrius ruficapillus",
                                 "Charadrius (Eupoda) veredus" = "Charadrius veredus",
                                 "Chenonetta" = "Chenonetta jubata",
                                 "Chlamydera guttata guttata" = "Chlamydera guttata",
                                 'Chlamydera nuchalis nuchalis' = 'Chlamydera nuchalis',
                                 "Chlidonias (Chlidonias) leucopterus" = "Chlidonias leucopterus",
                                 "Chlidonias (Pelodes) hybrida" = "Chlidonias hybrida",
                                 "Chlidonias (Pelodes) hybrida javanicus" = "Chlidonias hybrida",
                                 "Chroicocephalus" = "Chroicocephalus novaehollandiae",
                                 "Chroicocephalus novaehollandiae novaehollandiae" = "Chroicocephalus novaehollandiae",
                                 "Cincloramphus (Cincloramphus) cruralis" = "Cincloramphus cruralis",
                                 "Cincloramphus (Maclennania) mathewsi" = "Cincloramphus mathewsi",
                                 "Cinclosoma (Cinclosoma) punctatum" = "Cinclosoma punctatum",
                                 "Cinclosoma (Cinclosoma) punctatum anachoreta" = "Cinclosoma punctatum",
                                 "Cinclosoma (Cinclosoma) punctatum punctatum" = "Cinclosoma punctatum",
                                 "Cinclosoma (Malleeavis) castanotum" = "Cinclosoma castanotum",
                                 "Cinclosoma (Malleeavis) castanotum castanotum" = "Cinclosoma castanotum",
                                 "Cinclosoma (Malleeavis) clarum" = "Cinclosoma clarum",
                                 "Cinclosoma (Samuela) cinnamomeum" = "Cinclosoma cinnamomeum",
                                 "Cinclosoma (Samuela) cinnamomeum cinnamomeum" = "Cinclosoma cinnamomeum",
                                 "Cisticola (Cisticola) exilis" = "Cisticola exilis",
                                 "Cisticola (Cisticola) exilis exilis" = "Cisticola exilis",
                                 'Cisticola (Cisticola) exilis diminuta' = 'Cisticola exilis',
                                 'Cisticola (Cisticola) exilis lineocapilla' = 'Cisticola exilis',
                                 "Cisticola (Cisticola) juncidis" = "Cisticola juncidis",
                                 'Cisticola (Cisticola) juncidis laveryi' = 'Cisticola juncidis',
                                 "Cladorhynchus" = "Cladorhynchus leucocephalus",
                                 'Climacteris (Climacterobates) erythrops' = 'Climacteris erythrops',
                                 'Climacteris (Climacteris) melanurus' = 'Climacteris melanurus',
                                 "Climacteris (Climacteris) melanurus melanurus" = "Climacteris melanurus",
                                 "Climacteris (Climacteris) picumnus" = "Climacteris picumnus",
                                 "Climacteris (Climacteris) picumnus picumnus" = "Climacteris picumnus",
                                 "Climacteris (Climacteris) rufus" = "Climacteris rufus",
                                 "Climacteris (Climacterobates) affinis" = "Climacteris affinis",
                                 "Climacteris (Climacterobates) affinis affinis" = "Climacteris affinis",
                                 "Climacteris (Climacterobates) affinis superciliosus" = "Climacteris affinis",
                                 'Collocalia' = 'Collocalia esculenta',
                                 'Colluricincla (Myiolestes) boweri' = 'Colluricincla boweri',
                                 "Colluricincla" = "Colluricincla harmonica",
                                 "Colluricincla (Colluricincla) harmonica" = "Colluricincla harmonica",
                                 'Colluricincla (Colluricincla) harmonica brunnea' = 'Colluricincla harmonica',
                                 "Colluricincla (Colluricincla) harmonica harmonica" = "Colluricincla harmonica",
                                 "Colluricincla (Colluricincla) harmonica rufiventris" = "Colluricincla harmonica",
                                 'Colluricincla (Colluricincla) harmonica superciliosa' = 'Colluricincla harmonica',
                                 'Colluricincla (Myiolestes) megarhyncha' = 'Colluricincla megarhyncha',
                                 'Colluricincla (Myiolestes) megarhyncha parvula' = 'Colluricincla megarhyncha',
                                 'Colluricincla (Myiolestes) rufogaster' = 'Colluricincla rufogaster',
                                 'Colluricincla (Myiolestes) rufogaster rufogaster' = 'Colluricincla rufogaster',
                                 'Colluricincla (Colluricincla) woodwardi' = 'Colluricincla woodwardi',
                                 "Columba (Janthoenas) leucomela" = "Columba leucomela",
                                 'Columba (Janthoenas) vitiensis godmanae' = 'Columba vitiensis',
                                 'Conopophila (Conopophila) albogularis' = 'Conopophila albogularis',
                                 'Conopophila (Conopophila) rufogularis' = 'Conopophila rufogularis',
                                 'Coracina (Paragraucalus) lineata' = 'Coracina lineata',
                                 "Coracina (Coracina) novaehollandiae" = "Coracina novaehollandiae",
                                 "Coracina (Coracina) novaehollandiae melanops" = "Coracina novaehollandiae",
                                 "Coracina (Coracina) novaehollandiae subpallida" = "Coracina novaehollandiae",
                                 'Coracina (Coracina) novaehollandiae novaehollandiae' = 'Coracina novaehollandiae',
                                 "Coracina (Coracina) papuensis" = "Coracina papuensis",
                                 'Coracina (Coracina) papuensis apsleyi' = 'Coracina papuensis',
                                 'Coracina (Coracina) papuensis hypoleuca' = 'Coracina papuensis',
                                 'Coracina (Coracina) papuensis oriomo' = 'Coracina papuensis',
                                 "Coracina (Coracina) papuensis robusta" = "Coracina papuensis",
                                 "Coracina (Pteropodocys) maxima" = "Coracina maxima",
                                 "Corcorax melanorhamphos melanorhamphos" = "Corcorax melanorhamphos",
                                 "Corcorax melanorhamphos whitei" = "Corcorax melanorhamphos",
                                 "Cormobates" = "Cormobates leucophaea",
                                 "Cormobates leucophaea grisescens" = "Cormobates leucophaea",
                                 "Cormobates leucophaea leucophaea" = "Cormobates leucophaea",
                                 'Cormobates leucophaea metastasis' = 'Cormobates leucophaea',
                                 "Corvus coronoides coronoides" = "Corvus coronoides",
                                 'Corvus coronoides perplexus' = 'Corvus coronoides',
                                 'Corvus orru cecilae' = 'Corvus orru',
                                 'Corvus orru orru' = 'Corvus orru',
                                 "Corvus tasmanicus tasmanicus" = "Corvus tasmanicus",
                                 "Coturnix (Coturnix) pectoralis" = "Coturnix pectoralis",
                                 "Coturnix ypsilophora australis" = "Coturnix ypsilophora", 
                                 "Cracticus nigrogularis nigrogularis" = "Cracticus nigrogularis",
                                 'Cracticus nigrogularis picatus' = 'Cracticus nigrogularis',
                                 "Cracticus torquatus leucopterus" = "Cracticus torquatus",
                                 "Cracticus torquatus torquatus" = "Cracticus torquatus",
                                 'Cyanoramphus novaezelandiae subflavescens' = 'Cyanoramphus novaezelandiae',
                                 'Cyclopsitta diophthalma macleayana' = 'Cyclopsitta diophthalma',
                                 "Cygnus" = "Cygnus atratus",
                                 "Cygnus (Chenopis) atratus" = "Cygnus atratus",
                                 "Dacelo (Dacelo) leachii" = "Dacelo leachii",
                                 "Dacelo (Dacelo) leachii leachii" = "Dacelo leachii",
                                 'Dacelo (Dacelo) leachii cervina' = 'Dacelo leachii',
                                 "Dacelo (Dacelo) novaeguineae" = "Dacelo novaeguineae",
                                 "Dacelo (Dacelo) novaeguineae novaeguineae" = "Dacelo novaeguineae",
                                 "Daphoenositta" = "Daphoenositta chrysoptera",
                                 "Daphoenositta (Neositta)" = "Daphoenositta chrysoptera",
                                 "Daphoenositta (Neositta) chrysoptera" = "Daphoenositta chrysoptera",
                                 "Daphoenositta (Neositta) chrysoptera pileata" = "Daphoenositta chrysoptera",
                                 'Daphoenositta (Neositta) chrysoptera leucoptera' = 'Daphoenositta chrysoptera',
                                 "Daption capense australe" = "Daption capense",
                                 "Daption capense capense" = "Daption capense",
                                 "Dasyornis (Maccoyornis) broadbenti" = "Dasyornis broadbenti",
                                 "Dasyornis (Maccoyornis) broadbenti broadbenti" = "Dasyornis broadbenti",
                                 'Dendrocygna arcuata australis' = 'Dendrocygna arcuata',
                                 'Dicaeum (Dicaeum) geelvinkianum' = 'Dicaeum geelvinkianum',
                                 'Dicaeum (Dicaeum) hirundinaceum' = 'Dicaeum hirundinaceum',
                                 'Dicaeum (Dicaeum) hirundinaceum hirundinaceum' = 'Dicaeum hirundinaceum',
                                 'Dicrurus bracteatus atrabectus' = 'Dicrurus bracteatus',
                                 'Dicrurus bracteatus baileyi' = 'Dicrurus bracteatus',
                                 'Dicrurus bracteatus bracteatus' = 'Dicrurus bracteatus',
                                 "Diomedea antipodensis antipodensis" = "Diomedea antipodensis",
                                 "Diomedea antipodensis gibsoni" = "Diomedea antipodensis",
                                 "Dromaius" = "Dromaius novaehollandiae",
                                 "Dromaius novaehollandiae baudinianus" = "Dromaius novaehollandiae",
                                 "Dromaius novaehollandiae diemenianus" = "Dromaius novaehollandiae",
                                 "Dromaius novaehollandiae novaehollandiae" = "Dromaius novaehollandiae",
                                 'Ducula (Myristicivora) bicolor' = 'Ducula bicolor',
                                 'Ducula (Myristicivora) spilorrhoa' = 'Ducula spilorrhoa',
                                 'Eclectus' = 'Eclectus polychloros',
                                 'Edolisoma' = 'Edolisoma tenuirostre',
                                 "Edolisoma tenuirostre tenuirostre" = "Edolisoma tenuirostre",
                                 'Edolisoma tenuirostre melvillensis' = 'Edolisoma tenuirostre',
                                 "Egretta garzetta nigripes" = "Egretta garzetta",
                                 "Egretta sacra sacra" = "Egretta sacra",
                                 "Elseyornis" = "Elseyornis melanops",
                                 "Emblema" = "Emblema pictum",
                                 "Entomyzon cyanotis cyanotis" = "Entomyzon cyanotis",
                                 "Eolophus roseicapilla albiceps" = "Eolophus roseicapilla",
                                 "Eolophus roseicapilla roseicapilla" = "Eolophus roseicapilla",
                                 'Eolophus roseicapilla kuhli' = 'Eolophus roseicapilla',
                                 "Eopsaltria (Eopsaltria) australis" = "Eopsaltria australis",
                                 'Eopsaltria (Eopsaltria) australis chrysorrhoa' = 'Eopsaltria australis',
                                 "Eopsaltria (Eopsaltria) australis australis" = "Eopsaltria australis",
                                 "Eopsaltria (Eopsaltria) griseogularis" = "Eopsaltria griseogularis",
                                 "Eopsaltria (Eopsaltria) griseogularis griseogularis" = "Eopsaltria griseogularis",
                                 "Eopsaltria (Eopsaltria) griseogularis rosinae" = "Eopsaltria griseogularis",
                                 "Ephippiorhynchus (Ephippiorhynchus) asiaticus" = "Ephippiorhynchus asiaticus",
                                 "Epthianura (Aurepthianura) aurifrons" = "Epthianura aurifrons",
                                 'Epthianura (Aurepthianura) crocea' = 'Epthianura crocea',
                                 "Epthianura (Aurepthianura) crocea crocea" = "Epthianura crocea",
                                 "Epthianura (Aurepthianura) crocea tunneyi" = "Epthianura crocea",
                                 "Epthianura (Epthianura) albifrons" = "Epthianura albifrons",
                                 "Epthianura (Parepthianura) tricolor" = "Epthianura tricolor",
                                 'Esacus' = 'Esacus magnirostris',
                                 "Eudynamys orientalis cyanocephalus" = "Eudynamys orientalis",
                                 'Eudynamys' = 'Eudynamys orientalis',
                                 "Eudyptula" = "Eudyptula minor",
                                 "Eudyptula minor novaehollandiae" = "Eudyptula minor",
                                 'Eulabeornis castaneoventris castaneoventris' = 'Eulabeornis castaneoventris',
                                 "Eurostopodus (Eurostopodus) argus" = "Eurostopodus argus",
                                 "Eurostopodus (Eurostopodus) mystacalis" = "Eurostopodus mystacalis",
                                 'Eurostopodus (Eurostopodus) mystacalis mystacalis' = 'Eurostopodus mystacalis',
                                 'Eurystomus orientalis pacificus' = 'Eurystomus orientalis',
                                 "Excalfactoria chinensis australis" = "Excalfactoria chinensis",
                                 'Falco (Ieracidea) berigora novaeguineae' = 'Falco berigora',
                                 "Falco (Falco) longipennis" = "Falco longipennis",
                                 "Falco (Falco) longipennis longipennis" = "Falco longipennis",
                                 "Falco (Falco) longipennis murchisonianus" = "Falco longipennis",
                                 "Falco (Hierofalco) hypoleucos" = "Falco hypoleucos",
                                 "Falco (Hierofalco) peregrinus" = "Falco peregrinus",
                                 "Falco (Hierofalco) peregrinus macropus" = "Falco peregrinus",
                                 "Falco (Hierofalco) subniger" = "Falco subniger",
                                 "Falco (Ieracidea) berigora" = "Falco berigora",
                                 "Falco (Ieracidea) berigora berigora" = "Falco berigora",
                                 "Falco (Ieracidea) berigora occidentalis" = "Falco berigora",
                                 "Falco (Tinnunculus) cenchroides" = "Falco cenchroides",
                                 "Falco (Tinnunculus) cenchroides cenchroides" = "Falco cenchroides",
                                 "Falcunculus frontatus frontatus" = "Falcunculus frontatus",
                                 'Fregata ariel ariel' = 'Fregata ariel',
                                 'Fregata minor palmerstoni' = 'Fregata minor',
                                 'Fregetta grallaria grallaria' = 'Fregetta grallaria', 
                                 "Fregetta tropica tropica" = "Fregetta tropica",
                                 "Fulica" = "Fulica atra",
                                 'Fulica atra australis' = 'Fulica atra',
                                 "Fulica (Fulica) australias" = "Fulica australis",
                                 "Gallinago (Gallinago) hardwickii" = "Gallinago hardwickii",
                                 'Gallinago (Gallinago) megala' = 'Gallinago megala',
                                 'Gallinago (Gallinago) stenura' = 'Gallinago stenura',
                                 "Gallinula" = "Gallinula tenebrosa",
                                 "Gallinula (Gallinula) tenebrosa tenebrosa" = "Gallinula tenebricosa",
                                 "Gallinula (Gallinula) tenebrosa" = "Gallinula tenebricosa",
                                 'Gavicalis versicolor versicolor' = 'Gavicalis versicolor',
                                 "Gavicalis virescens forresti" = "Gavicalis virescens",
                                 "Gavicalis virescens sonorus" = "Gavicalis virescens",
                                 "Gavicalis virescens virescens" = "Gavicalis virescens",
                                 "Gelochelidon nilotica affinis" = "Gelochelidon nilotica",
                                 "Gelochelidon nilotica macrotarsa" = "Gelochelidon nilotica",
                                 "Geopelia placida placida" = "Geopelia placida",
                                 'Geopelia humeralis headlandi' = 'Geopelia humeralis',
                                 'Geopelia humeralis humeralis' = 'Geopelia humeralis',
                                 'Geopelia humeralis inexpectata' = 'Geopelia humeralis',
                                 "Geophaps (Geophaps) scripta" = "Geophaps scripta",
                                 "Geophaps (Lophophaps) plumifera" = "Geophaps plumifera",
                                 "Geophaps (Lophophaps) plumifera ferruginea" = "Geophaps plumifera",
                                 "Geophaps (Lophophaps) plumifera leucogaster" = "Geophaps plumifera",
                                 'Geophaps (Geophaps) scripta peninsulae' = 'Geophaps scripta',
                                 'Geophaps (Geophaps) smithii' = 'Geophaps smithii',
                                 'Geophaps (Geophaps) smithii blaauwi' = 'Geophaps smithii',
                                 'Geophaps (Geophaps) smithii smithii' = 'Geophaps smithii',
                                 'Gerygone chloronota chloronota' = 'Gerygone chloronota',
                                 "Gerygone fusca fusca" = "Gerygone fusca",
                                 "Gerygone fusca mungi" = "Gerygone fusca",
                                 'Gerygone magnirostris brunneipectus' = 'Gerygone magnirostris',
                                 'Gerygone magnirostris cairnsensis' = 'Gerygone magnirostris',
                                 'Gerygone magnirostris magnirostris' = 'Gerygone magnirostris',
                                 "Gerygone olivacea olivacea" = "Gerygone olivacea",
                                 'Gerygone olivacea rogersi' = 'Gerygone olivacea',
                                 'Gerygone palpebrosa flavida' = 'Gerygone palpebrosa',
                                 "Glareola (Glareola) maldivarum" = "Glareola maldivarum",
                                 "Gliciphila melanops melanops" = "Gliciphila melanops",
                                 "Glossopsitta" = "Glossopsitta concinna",
                                 "Glossopsitta concinna concinna" = "Glossopsitta concinna",
                                 "Grallina cyanoleuca cyanoleuca" = "Grallina cyanoleuca",
                                 'Gygis alba candida' = 'Gygis alba',
                                 "Gymnorhina" = "Gymnorhina tibicen",
                                 "Gymnorhina tibicen dorsalis" = "Gymnorhina tibicen",
                                 "Gymnorhina tibicen hypoleuca" = "Gymnorhina tibicen",
                                 "Gymnorhina tibicen telonocua" = "Gymnorhina tibicen",
                                 "Gymnorhina tibicen tibicen" = "Gymnorhina tibicen",
                                 "Gymnorhina tibicen tyrannica" = "Gymnorhina tibicen",
                                 'Gymnorhina tibicen eylandtensis' = 'Gymnorhina tibicen',
                                 "Haematopus fuliginosus fuliginosus" = "Haematopus fuliginosus",
                                 'Haematopus fuliginosus opthalmicus' = 'Haematopus fuliginosus',
                                 "Haliaeetus (Pontoaetus) leucogaster" = "Haliaeetus leucoryphus",
                                 'Haliastur indus girrenera' = 'Haliastur indus',
                                 "Hieraaetus (Hieraaetus) morphnoides morphnoides" = "Hieraaetus morphnoides",
                                 "Hieraaetus (Hieraaetus) morphnoides" = "Hieraaetus morphnoides",
                                 "Himantopus himantopus leucocephalus" = "Himantopus himantopus",
                                 "Himantopus himantopus leucocephalus" = "Himantopus himantopus",
                                 "Himantopus" = "Himantopus himantopus",
                                 "Hirundapus caudacutus caudacutus" = "Hirundapus caudacutus",
                                 "Hirundo (Hirundo) neoxena" = "Hirundo neoxena",
                                 "Hirundo (Hirundo) neoxena neoxena" = "Hirundo neoxena",
                                 'Hirundo (Hirundo) neoxena carteri' = 'Hirundo neoxena',
                                 "Hirundo (Hirundo) rustica" = "Hirundo rustica",
                                 "Hirundo (Hirundo) tahitica" = "Hirundo tahitica",
                                 "Hylacola cauta cauta" = "Hylacola cauta",
                                 "Hylacola cauta halmaturina" = "Hylacola cauta",
                                 "Hylacola pyrrhopygia parkeri" = "Hylacola pyrrhopygia",
                                 "Hylacola pyrrhopygia pedleri" = "Hylacola pyrrhopygia",
                                 "Hypotaenidia philippensis mellori" = "Hypotaenidia philippensis",
                                 'Hypotaenidia philippensis tounelieri' = 'Hypotaenidia philippensis',
                                 'Irediparra gallinacea novaehollandiae' = 'Irediparra gallinacea',
                                 'Ixobrychus flavicollis australis' = 'Ixobrychus flavicollis',
                                 'Lalage (Karua) leucomela' = 'Lalage leucomela',
                                 'Lalage (Karua) leucomela yorki' = 'Lalage leucomela',
                                 'Lalage (Karua) leucomela leucomela' = 'Lalage leucomela',
                                 'Lalage (Karua) leucomela rufiventris' = 'Lalage leucomela',
                                 "Lalage (Lalage) tricolor" = "Lalage tricolor",
                                 "Larus dominicanus dominicanus" = "Larus dominicanus",
                                 "Larus pacificus georgii" = "Larus pacificus",
                                 "Larus pacificus pacificus" = "Larus pacificus",
                                 "Lewinia pectoralis pectoralis" = "Lewinia pectoralis",
                                 "Lichenostomus cratitius cratitius" = "Lichenostomus cratitius",
                                 "Lichenostomus cratitius occidentalis" = "Lichenostomus cratitius",
                                 'Lichmera' = 'Lichmera indistincta',
                                 "Lichmera (Lichmera) indistincta" = "Lichmera indistincta",
                                 "Lichmera (Lichmera) indistincta indistincta" = "Lichmera indistincta",
                                 'Lichmera (Lichmera) indistincta ocularis' = 'Lichmera indistincta',
                                 "Limosa lapponica baueri" = "Limosa lapponica",
                                 'Limosa lapponica menzbieri' = 'Limosa lapponica',
                                 "Limosa limosa melanuroides" = "Limosa limosa",
                                 'Lonchura (Munia) castaneothorax castaneothorax' = 'Lonchura castaneothorax',
                                 "Lonchura (Lonchura) punctulata" = "Lonchura punctulata",
                                 "Lonchura (Munia) castaneothorax" = "Lonchura castaneothorax",
                                 "Lophochroa leadbeateri leadbeateri" = "Lophochroa leadbeateri",
                                 "Lophochroa leadbeateri mollis" = "Lophochroa leadbeateri",
                                 "Lugensa" = "Lugensa brevirostris",
                                 'Macropygia' = 'Macropygia phasianella',
                                 'Macropygia (Macropygia) phasianella' = 'Macropygia phasianella',
                                 'Macropygia (Macropygia) phasianella phasianella' = 'Macropygia phasianella',
                                 "Malurus (Leggeornis) amabilis" = "Malurus amabilis",
                                 "Malurus (Leggeornis) assimilis" = "Malurus assimilis",
                                 'Malurus (Leggeornis) assimilis rogersi' = 'Malurus assimilis',
                                 "Malurus (Leggeornis) assimilis assimilis" = "Malurus assimilis",
                                 "Malurus (Leggeornis) elegans" = "Malurus elegans",
                                 "Malurus (Leggeornis) lamberti" = "Malurus lamberti",
                                 "Malurus (Leggeornis) pulcherrimus" = "Malurus pulcherrimus",
                                 "Malurus (Malurus) coronatus" = "Malurus coronatus",
                                 'Malurus (Malurus) coronatus coronatus' = 'Malurus coronatus',
                                 'Malurus (Malurus) coronatus macgillivrayi' = 'Malurus coronatus',
                                 "Malurus (Malurus) cyaneus" = "Malurus cyaneus",
                                 "Malurus (Malurus) cyaneus ashbyi" = "Malurus cyaneus",
                                 "Malurus (Malurus) cyaneus cyaneus" = "Malurus cyaneus",
                                 "Malurus (Malurus) cyaneus cyanochlamys" = "Malurus cyaneus",
                                 "Malurus (Malurus) cyaneus leggei" = "Malurus cyaneus",
                                 "Malurus (Malurus) splendens" = "Malurus splendens",
                                 "Malurus (Malurus) splendens melanotus" = "Malurus splendens",
                                 "Malurus (Malurus) splendens musgravi" = "Malurus splendens",
                                 "Malurus (Malurus) splendens splendens" = "Malurus splendens",
                                 "Malurus (Musciparus) leucopterus" = "Malurus leucopterus",
                                 "Malurus (Musciparus) leucopterus leuconotus" = "Malurus leucopterus",
                                 'Malurus (Musciparus) leucopterus edouardi' = 'Malurus leucopterus',
                                 "Malurus (Musciparus) leucopterus leucopterus" = "Malurus leucopterus",
                                 'Malurus (Musciparus) melanocephalus' = 'Malurus melanocephalus',
                                 'Malurus (Musciparus) melanocephalus cruentatus' = 'Malurus melanocephalus',
                                 'Malurus (Musciparus) melanocephalus melanocephalus' = 'Malurus melanocephalus',
                                 "Manorina (Manorina) melanophrys" = "Manorina melanophrys",
                                 "Manorina (Myzantha) flavigula" = "Manorina flavigula",
                                 "Manorina (Myzantha) flavigula flavigula" = "Manorina flavigula",
                                 'Manorina (Myzantha) flavigula melvillensis' = 'Manorina flavigula',
                                 "Manorina (Myzantha) flavigula wayensis" = "Manorina flavigula",
                                 "Manorina (Myzantha) melanocephala" = "Manorina melanocephala",
                                 "Manorina (Myzantha) melanocephala melanocephala" = "Manorina melanocephala",
                                 "Manorina (Myzantha) melanotis" = "Manorina melanotis",
                                 "Menura (Menura) novaehollandiae" = "Menura novaehollandiae",
                                 'Megapodius (Megapodius) reinwardt' = 'Megapodius reinwardt',
                                 'Megapodius (Megapodius) reinwardt castanonotus' = 'Megapodius reinwardt',
                                 'Megapodius (Megapodius) reinwardt tumulus' = 'Megapodius reinwardt',
                                 'Megapodius (Megapodius) reinwardt reinwardt' = 'Megapodius reinwardt',
                                 'Megapodius (Megapodius) reinwardt yorki' = 'Megapodius reinwardt',
                                 "Melanodryas (Amaurodryas) vittata" = "Melanodryas vittata",
                                 "Melanodryas (Melanodryas) cucullata" = "Melanodryas cucullata",
                                 'Melanodryas (Melanodryas) cucullata melvillensis' = 'Melanodryas cucullata',
                                 "Melanodryas (Melanodryas) cucullata cucullata" = "Melanodryas cucullata",
                                 "Melanodryas (Melanodryas) cucullata westralensis" = "Melanodryas cucullata",
                                 'Meliphaga (Microptilotis) albilineata' = 'Meliphaga albilineata',
                                 'Meliphaga (Microptilotis) fordiana' = 'Meliphaga fordiana',
                                 'Meliphaga (Microptilotis) gracilis' = 'Meliphaga gracilis',
                                 'Meliphaga (Microptilotis) gracilis gracilis' = 'Meliphaga gracilis',
                                 'Meliphaga (Meliphaga) lewinii' = 'Meliphaga lewinii',
                                 'Meliphaga (Meliphaga) lewinii lewinii' = 'Meliphaga lewinii',
                                 'Meliphaga (Meliphaga) notata' = 'Meliphaga notata',
                                 'Meliphaga (Meliphaga) notata notata' = 'Meliphaga notata',
                                 'Melithreptus (Melithreptus) affinis' = 'Melithreptus affinis',
                                 'Melithreptus (Melithreptus) albogularis' = 'Melithreptus albogularis',
                                 'Melithreptus (Melithreptus) albogularis albogularis' = 'Melithreptus albogularis',
                                 'Melithreptus' = 'Melithreptus brevirostris',
                                 'Melithreptus (Eidopsarus) brevirostris' = 'Melithreptus brevirostris',
                                 'Melithreptus (Eidopsarus) brevirostris brevirostris' = 'Melithreptus brevirostris',
                                 'Melithreptus (Melithreptus) chloropsis' = 'Melithreptus chloropsis',
                                 'Melithreptus (Eidopsarus) brevirostris leucogenys' = 'Melithreptus brevirostris',
                                 'Melithreptus (Eidopsarus) brevirostris magnirostris' = 'Melithreptus brevirostris',
                                 'Melithreptus (Eidopsarus) brevirostris pallidiceps' = 'Melithreptus brevirostris',
                                 'Melithreptus (Eidopsarus) gularis' = 'Melithreptus gularis',
                                 'Melithreptus (Eidopsarus) gularis gularis' = 'Melithreptus gularis',
                                 'Melithreptus (Melithreptus) lunatus' = 'Melithreptus lunatus',
                                 'Merops (Merops) ornatus' = 'Merops ornatus',
                                 'Microcarbo' = 'Microcarbo melanoleucos',
                                 'Microcarbo melanoleucos melanoleucos' = 'Microcarbo melanoleucos',
                                 'Microeca (Microeca) fascinans' = 'Microeca fascinans',
                                 'Microeca (Microeca) fascinans assimilis' = 'Microeca fascinans',
                                 'Microeca (Microeca) fascinans fascinans' = 'Microeca fascinans',
                                 'Microeca (Microeca) flavigaster' = 'Microeca flavigaster',
                                 'Microeca (Microeca) flavigaster flavigaster' = 'Microeca flavigaster',
                                 'Microeca (Microeca) flavigaster tormenti' = 'Microeca flavigaster',
                                 'Microeca (Kempiella) griseoceps' = 'Microeca griseoceps',
                                 'Milvus migrans affinis' = 'Milvus migrans',
                                 'Mirafra (Mirafra) javanica' = 'Mirafra javanica',
                                 'Mirafra (Mirafra) javanica melvillensis' = 'Mirafra javanica',
                                 'Mirafra (Mirafra) javanica horsfieldii' = 'Mirafra javanica',
                                 'Mirafra (Mirafra) javanica secunda' = 'Mirafra javanica',
                                 'Mirafra (Mirafra) javanica woodwardi' = 'Mirafra javanica',
                                 'Monarcha (Monarcha) frater' = 'Monarcha frater',
                                 'Monarcha (Monarcha) melanopsis' = 'Monarcha melanopsis',
                                 'Motacilla (Motacilla) alba' = 'Motacilla alba',
                                 'Motacilla (Calobates) cinerea' = 'Motacilla cinerea',
                                 'Motacilla (Budytes) tschutschensis' = 'Motacilla tschutschensis',
                                 'Motacilla (Budytes) tschutschensis tschutschensis' = 'Motacilla tschutschensis',
                                 'Myiagra (Piezorhynchus) alecto' = 'Myiagra alecto',
                                 'Myiagra (Piezorhynchus) alecto melvillensis' = 'Myiagra alecto',
                                 'Myiagra (Piezorhynchus) alecto wardelli' = 'Myiagra alecto',
                                 'Myiagra (Myiagra) cyanoleuca' = 'Myiagra cyanoleuca',
                                 'Myiagra (Seisura) inquieta' = 'Myiagra inquieta',
                                 'Myiagra (Seisura) nana' = 'Myiagra nana',
                                 'Myiagra (Myiagra) rubecula' = 'Myiagra rubecula',
                                 'Myiagra (Myiagra) rubecula concinna' = 'Myiagra rubecula',
                                 'Myiagra (Myiagra) rubecula okyri' = 'Myiagra rubecula',
                                 'Myiagra (Myiagra) rubecula rubecula' = 'Myiagra rubecula',
                                 'Myiagra (Myiagra) rubecula papuana' = 'Myiagra rubecula',
                                 'Myiagra (Myiagra) rubecula yorki' = 'Myiagra rubecula',
                                 'Myiagra (Myiagra) ruficollis' = 'Myiagra ruficollis',
                                 'Myiagra (Myiagra) ruficollis mimikae' = 'Myiagra ruficollis', 
                                 'Myzomela (Myzomela) erythrocephala' = 'Myzomela erythrocephala',
                                 'Myzomela (Myzomela) erythrocephala erythrocephala' = 'Myzomela erythrocephala',
                                 'Myzomela (Myzomela) erythrocephala infuscata' = 'Myzomela erythrocephala',
                                 'Myzomela (Cosmeteira) obscura' = 'Myzomela obscura',
                                 'Myzomela (Cosmeteira) obscura harterti' = 'Myzomela obscura',
                                 'Myzomela (Myzomela) sanguinolenta' = 'Myzomela sanguinolenta',
                                 'Myzomela (Myzomela) sanguinolenta sanguinolenta' = 'Myzomela sanguinolenta',
                                 'Nectarinia (Cyrtostomus) jugularis' = 'Nectarinia jugularis',
                                 'Nectarinia (Cyrtostomus) jugularis frenata' = 'Nectarinia jugularis',
                                 'Neochmia (Neochmia) phaeton' = 'Neochmia phaeton',
                                 'Neochmia (Neochmia) phaeton phaeton' = 'Neochmia phaeton',
                                 'Neochmia (Neochmia) ruficauda clarescens' = 'Neochmia ruficauda',
                                 'Neochmia (Aegintha) temporalis' = 'Neochmia temporalis',
                                 'Neochmia (Aegintha) temporalis temporalis' = 'Neochmia temporalis',
                                 'Neophema (Neonanodes) chrysogaster' = 'Neophema chrysogaster',
                                 'Neophema (Neonanodes) chrysostoma' = 'Neophema chrysostoma',
                                 'Neophema (Neonanodes) elegans' = 'Neophema elegans',
                                 'Neophema (Neonanodes) elegans elegans' = 'Neophema elegans',
                                 'Neophema (Neonanodes) petrophila' = 'Neophema petrophila',
                                 'Neophema (Neonanodes) petrophila zietzi' = 'Neophema petrophila',
                                 'Neophema (Neonanodes) petrophila petrophila' = 'Neophema petrophila',
                                 'Neophema (Neophema) pulchella' = 'Neophema pulchella',
                                 'Neophema (Neophema) splendida' = 'Neophema splendida',
                                 'Neosericornis citreogularis citreogularis' = 'Neosericornis citreogularis',
                                 'Nesoptilotis leucotis depauperata' = 'Nesoptilotis leucotis',
                                 'Nesoptilotis leucotis leucotis' = 'Nesoptilotis leucotis',
                                 'Nesoptilotis leucotis novaenorciae' = 'Nesoptilotis leucotis',
                                 'Nesoptilotis leucotis schoddei' = 'Nesoptilotis leucotis',
                                 'Nesoptilotis leucotis thomasi' = 'Nesoptilotis leucotis',
                                 'Nettapus (Cheniscus) coromandelianus' = 'Nettapus coromandelianus',
                                 'Nettapus (Cheniscus) pulchellus' = 'Nettapus pulchellus',
                                 'Ninox (Ninox) boobook' = 'Ninox boobook',
                                 'Ninox (Ninox) boobook boobook' = 'Ninox boobook',
                                 'Ninox (Ninox) boobook halmaturina' = 'Ninox boobook',
                                 'Ninox (Ninox) boobook ocellata' = 'Ninox boobook',
                                 'Ninox (Hieracoglaux) connivens' = 'Ninox connivens',
                                 'Ninox (Hieracoglaux) connivens connivens' = 'Ninox connivens',
                                 'Ninox (Hieracoglaux) connivens peninsularis' = 'Ninox connivens',
                                 'Ninox (Ninox) novaeseelandiae' = 'Ninox novaeseelandiae',
                                 'Ninox (Ninox) novaeseelandiae albaria' = 'Ninox novaeseelandiae',
                                 'Ninox (Rhabdoglaux) rufa' = 'Ninox rufa',
                                 'Ninox (Rhabdoglaux) rufa queenslandica' = 'Ninox rufa',
                                 'Ninox (Rhabdoglaux) strenua' = 'Ninox strenua',
                                 'Numenius (Numenius) madagascariensis' = 'Numenius madagascariensis',
                                 'Numenius (Mesoscolopax) minutus' = 'Numenius minutus',
                                 'Numenius (Phaeopus) phaeopus' = 'Numenius phaeopus',
                                 'Numenius (Phaeopus) phaeopus variegatus' = 'Numenius phaeopus',
                                 'Nycticorax caledonicus australasiae' = 'Nycticorax caledonicus',
                                 'Oceanites' = 'Oceanites oceanicus',
                                 'Oceanites oceanicus oceanicus' = 'Oceanites oceanicus',
                                 'Ocyphaps' = 'Ocyphaps lophotes',
                                 'Ocyphaps lophotes lophotes' = 'Ocyphaps lophotes',
                                 'Onychoprion anaethetus anaethetus' = 'Onychoprion anaethetus',
                                 'Onychoprion fuscatus serrata' = 'Onychoprion fuscatus',
                                 'Oreoica' = 'Oreoica gutturalis',
                                 'Oreoica gutturalis pallescens' = 'Oreoica gutturalis',
                                 'Oreoica gutturalis gutturalis' = 'Oreoica gutturalis',
                                 'Oriolus (Mimeta) flavocinctus' = 'Oriolus flavocinctus',
                                 'Oriolus (Mimeta) flavocinctus flavocinctus' = 'Oriolus flavocinctus',
                                 'Oriolus (Mimeta) flavocinctus tiwi' = 'Oriolus flavocinctus',
                                 'Oriolus (Mimeta) sagittatus' = 'Oriolus sagittatus',
                                 'Oriolus (Mimeta) sagittatus affinis' = 'Oriolus sagittatus',
                                 'Oriolus (Mimeta) sagittatus sagittatus' = 'Oriolus sagittatus',
                                 'Oriolus (Mimeta) sagittatus grisescens' = 'Oriolus sagittatus',
                                 'Pachycephala (Timixos) inornata' = 'Pachycephala inornata',
                                 'Pachycephala (Alisterornis) lanioides' = 'Pachycephala lanioides',
                                 'Pachycephala (Pachycephala) melanura' = 'Pachycephala melanura',
                                 'Pachycephala (Pachycephala) melanura melanura' = 'Pachycephala melanura',
                                 'Pachycephala (Pachycephala) melanura robusta' = 'Pachycephala melanura',
                                 'Pachycephala (Timixos) olivacea' = 'Pachycephala olivacea',
                                 'Pachycephala (Timixos) olivacea hesperus' = 'Pachycephala olivacea',
                                 'Pachycephala (Timixos) olivacea olivacea' = 'Pachycephala olivacea',
                                 'Pachycephala (Pachycephala) pectoralis' = 'Pachycephala pectoralis',
                                 'Pachycephala (Pachycephala) pectoralis fuliginosa' = 'Pachycephala pectoralis',
                                 'Pachycephala (Pachycephala) pectoralis contempta' = 'Pachycephala pectoralis',
                                 'Pachycephala (Pachycephala) pectoralis pectoralis' = 'Pachycephala pectoralis',
                                 'Pachycephala (Pachycephala) pectoralis youngi' = 'Pachycephala pectoralis',
                                 'Pachycephala (Alisterornis) rufiventris' = 'Pachycephala rufiventris',
                                 'Pachycephala (Alisterornis) rufiventris falcata' = 'Pachycephala rufiventris',
                                 'Pachycephala (Alisterornis) rufiventris rufiventris' = 'Pachycephala rufiventris',
                                 'Pachycephala (Alisterornis) rufiventris minor' = 'Pachycephala rufiventris',
                                 'Pachycephala (Alisterornis) rufiventris pallida' = 'Pachycephala rufiventris',
                                 'Pachycephala (Timixos) rufogularis' = 'Pachycephala rufogularis',
                                 'Pachycephala (Mattingleya) simplex' = 'Pachycephala simplex',
                                 'Pachycephala (Mattingleya) simplex peninsulae' = 'Pachycephala simplex',
                                 'Pachycephala (Mattingleya) simplex simplex' = 'Pachycephala simplex',
                                 'Pachyptila turtur turtur' = 'Pachyptila turtur',
                                 'Pandion haliaetus cristatus' = 'Pandion haliaetus',
                                 'Pardalotus (Pardalotus) punctatus' = 'Pardalotus punctatus',
                                 'Pardalotus (Pardalotus) punctatus punctatus' = 'Pardalotus punctatus',
                                 'Pardalotus (Pardalotus) punctatus xanthopyge' = 'Pardalotus punctatus',
                                 'Pardalotus (Pardalotinus) rubricatus' = 'Pardalotus rubricatus',
                                 'Pardalotus (Pardalotinus) striatus' = 'Pardalotus striatus',
                                 'Pardalotus (Pardalotinus) striatus melvillensis' = 'Pardalotus striatus',
                                 'Pardalotus (Pardalotinus) striatus ornatus' = 'Pardalotus striatus',
                                 'Pardalotus (Pardalotinus) striatus melanocephalus' = 'Pardalotus striatus',
                                 'Pardalotus (Pardalotinus) striatus uropygialis' = 'Pardalotus striatus',
                                 'Pardalotus (Pardalotinus) striatus striatus' = 'Pardalotus striatus',
                                 'Pardalotus (Pardalotinus) striatus substriatus' = 'Pardalotus striatus',
                                 'Peneothello (Peneoenanthe)' = 'Peneothello pulverulenta',
                                 'Peneothello (Peneoenanthe) pulverulenta' = 'Peneothello pulverulenta',
                                 'Pelagodroma marina dulciae' = 'Pelagodroma marina',
                                 'Pelecanoides urinatrix urinatrix' = 'Pelecanoides urinatrix',
                                 'Petrochelidon (Petrochelidon) ariel' = 'Petrochelidon ariel',
                                 'Petrochelidon (Hylochelidon) nigricans' = 'Petrochelidon nigricans',
                                 'Petrochelidon (Hylochelidon) nigricans neglecta' = 'Petrochelidon nigricans',
                                 'Petrochelidon (Hylochelidon) nigricans nigricans' = 'Petrochelidon nigricans',
                                 'Petroica (Petroica) boodang' = 'Petroica boodang',
                                 'Petroica (Petroica) boodang boodang' = 'Petroica boodang',
                                 'Petroica (Petroica) boodang campbelli' = 'Petroica boodang',
                                 'Petroica (Petroica) boodang leggii' = 'Petroica boodang',
                                 'Petroica (Petroica) goodenovii' = 'Petroica goodenovii',
                                 'Petroica (Petroica) multicolor' = 'Petroica multicolor',
                                 'Petroica (Littlera) phoenicea' = 'Petroica phoenicea',
                                 'Petroica (Erythrodryas) rodinogaster' = 'Petroica rodinogaster',
                                 'Petroica (Erythrodryas) rodinogaster inexpectata' = 'Petroica rodinogaster',
                                 'Petroica (Erythrodryas) rosea' = 'Petroica rosea',
                                 'Petrophassa albipennis albipennis' = 'Petrophassa albipennis',
                                 'Pezoporus wallicus wallicus' = 'Pezoporus wallicus',
                                 'Phaethon lepturus dorotheae' = 'Phaethon lepturus',
                                 'Phaethon rubricauda roseotinctus' = 'Phaethon rubricauda',
                                 'Phaethon rubricauda westralis' = 'Phaethon rubricauda',
                                 'Phalacrocorax (Phalacrocorax) carbo' = 'Phalacrocorax carbo',
                                 'Phalacrocorax (Phalacrocorax) carbo novaehollandiae' = 'Phalacrocorax carbo',
                                 'Phalacrocorax (Anacarbo) fuscescens' = 'Phalacrocorax fuscescens',
                                 'Phalacrocorax (Phalacrocorax) sulcirostris' = 'Phalacrocorax sulcirostris',
                                 'Phalacrocorax (Phalacrocorax) varius' = 'Phalacrocorax varius',
                                 'Phalacrocorax (Phalacrocorax) varius hypoleucos' = 'Phalacrocorax varius',
                                 'Phaps (Phaps) chalcoptera' = 'Phaps chalcoptera',
                                 'Phaps (Phaps) elegans' = 'Phaps elegans',
                                 'Phaps (Phaps) elegans occidentalis' = 'Phaps elegans',
                                 'Phaps (Phaps) elegans elegans' = 'Phaps elegans',
                                 'Phaps (Histriophaps) histrionica' = 'Phaps histrionica',
                                 'Philemon (Philemon) argenticeps' = 'Philemon argenticeps',
                                 'Philemon (Philemon) argenticeps argenticeps' = 'Philemon argenticeps',
                                 'Philemon (Philemon) buceroides gordoni' = 'Philemon buceroides',
                                 'Philemon (Philemon) buceroides' = 'Philemon buceroides',
                                 'Philemon (Philemon) buceroides yorki' = 'Philemon buceroides',
                                 'Philemon (Microphilemon) citreogularis sordidus' = 'Philemon citreogularis',
                                 'Philemon (Microphilemon) citreogularis' = 'Philemon citreogularis',
                                 'Philemon (Microphilemon) citreogularis citreogularis' = 'Philemon citreogularis',
                                 'Philemon (Tropidorhynchus) corniculatus' = 'Philemon corniculatus',
                                 'Philemon (Tropidorhynchus) corniculatus monachus' = 'Philemon corniculatus',
                                 'Philemon (Tropidorhynchus) corniculatus corniculatus' = 'Philemon corniculatus',
                                 'Phylidonyris (Meliornis) niger' = 'Phylidonyris niger',
                                 'Phylidonyris (Meliornis) novaehollandiae' = 'Phylidonyris novaehollandiae',
                                 'Phylidonyris (Meliornis) novaehollandiae campbelli' = 'Phylidonyris novaehollandiae',
                                 'Phylidonyris (Meliornis) niger niger' = 'Phylidonyris niger',
                                 'Phylidonyris (Meliornis) novaehollandiae novaehollandiae' = 'Phylidonyris novaehollandiae',
                                 'Phylidonyris (Phylidonyris)' = 'Phylidonyris pyrrhopterus',
                                 'Phylidonyris (Phylidonyris) pyrrhopterus' = 'Phylidonyris pyrrhopterus',
                                 'Phylidonyris (Phylidonyris) pyrrhopterus halmaturinus' = 'Phylidonyris pyrrhopterus',
                                 'Phylidonyris (Phylidonyris) pyrrhopterus pyrrhopterus' = 'Phylidonyris pyrrhopterus',
                                 'Phylloscopus' = 'Phylloscopus borealis',
                                 'Phylloscopus (Acanthopneuste) borealis' = 'Phylloscopus borealis',
                                 'Pitta (Erythropitta) erythrogaster' = 'Pitta erythrogaster',
                                 'Pitta (Pitta) iris' = 'Pitta iris',
                                 'Pitta (Pitta) iris johnstoneiana' = 'Pitta iris',
                                 'Pitta (Pitta) iris iris' = 'Pitta iris',
                                 'Pitta (Pitta) versicolor' = 'Pitta versicolor',
                                 'Pitta (Pitta) versicolor simillima' = 'Pitta versicolor',
                                 'Pitta (Pitta) versicolor versicolor' = 'Pitta versicolor',
                                 'Platalea (Platibis) flavipes' = 'Platalea flavipes',
                                 'Platalea (Platalea) regia' = 'Platalea regia',
                                 'Platycercus (Violania) adscitus' = 'Platycercus adscitus',
                                 'Platycercus (Violania) adscitus palliceps' = 'Platycercus adscitus',
                                 'Platycercus (Platycercus) caledonicus' = 'Platycercus caledonicus',
                                 'Platycercus (Platycercus) elegans' = 'Platycercus elegans',
                                 'Platycercus (Platycercus) elegans elegans' = 'Platycercus elegans',
                                 'Platycercus (Platycercus) elegans flaveolus' = 'Platycercus elegans',
                                 'Platycercus (Platycercus) elegans fleurieuensis' = 'Platycercus elegans',
                                 'Platycercus (Platycercus) elegans melanopterus' = 'Platycercus elegans',
                                 'Platycercus (Platycercus) elegans subadelaidae' = 'Platycercus elegans',
                                 'Platycercus (Violania) eximius' = 'Platycercus eximius',
                                 'Platycercus (Violania) eximius eximius' = 'Platycercus eximius',
                                 'Platycercus (Violania) icterotis' = 'Platycercus icterotis',
                                 'Platycercus (Violania) venustus' = 'Platycercus venustus',
                                 'Platycercus (Violania) venustus venustus' = 'Platycercus venustus',
                                 'Pluvialis squatarola squatarola' = 'Pluvialis squatarola',
                                 'Podargus ocellatus plumiferus' = 'Podargus ocellatus',
                                 'Podargus strigoides brachypterus' = 'Podargus strigoides',
                                 'Podargus strigoides strigoides' = 'Podargus strigoides',
                                 'Podargus strigoides phalaenoides' = 'Podargus strigoides',
                                 'Podiceps' = 'Podiceps cristatus',
                                 'Podiceps cristatus australis' = 'Podiceps cristatus',
                                 'Poecilodryas (Poecilodryas) cerviniventris' = 'Poecilodryas cerviniventris',
                                 'Poecilodryas (Poecilodryas) superciliosa' = 'Poecilodryas superciliosa',
                                 'Poephila (Poephila) acuticauda' = 'Poephila acuticauda',
                                 'Poephila (Neopoephila) personata' = 'Poephila personata',
                                 'Poephila (Neopoephila) personata personata' = 'Poephila personata',
                                 'Polytelis anthopeplus monarchoides' = 'Polytelis anthopeplus',
                                 'Pomatostomus (Morganornis) ruficeps' = 'Pomatostomus ruficeps',
                                 'Pomatostomus (Morganornis) superciliosus' = 'Pomatostomus superciliosus',
                                 'Pomatostomus (Morganornis) superciliosus ashbyi' = 'Pomatostomus superciliosus',
                                 'Pomatostomus (Morganornis) superciliosus centralis' = 'Pomatostomus superciliosus',
                                 'Pomatostomus (Morganornis) superciliosus gilgandra' = 'Pomatostomus superciliosus',
                                 'Pomatostomus (Morganornis) superciliosus superciliosus' = 'Pomatostomus superciliosus',
                                 'Pomatostomus (Pomatostomus) temporalis' = 'Pomatostomus temporalis',
                                 'Pomatostomus (Pomatostomus) temporalis temporalis' = 'Pomatostomus temporalis',
                                 'Pomatostomus (Pomatostomus) temporalis rubeculus' = 'Pomatostomus temporalis',
                                 'Poodytes gramineus goulburni' = 'Poodytes gramineus',
                                 'Porphyrio (Notornis) albus' = 'Porphyrio albus',
                                 'Porphyrio (Porphyrio) porphyrio' = 'Porphyrio porphyrio',
                                 'Porphyrio (Porphyrio) porphyrio melanotus' = 'Porphyrio porphyrio',
                                 'Porzana' = 'Porzana fluminea',
                                 'Porzana (Porzana)' = 'Porzana fluminea',
                                 'Porzana (Porzana) fluminea' = 'Porzana fluminea',
                                 'Probosciger' = 'Probosciger aterrimus',
                                 'Probosciger aterrimus macgillivrayi' = 'Probosciger aterrimus',
                                 'Procellaria (Procellaria) aequinoctialis' = 'Procellaria aequinoctialis',
                                 'Procellaria (Procellaria) aequinoctialis aequinoctialis' = 'Procellaria aequinoctialis',
                                 'Procellaria (Adamastor) cinerea' = 'Procellaria cinerea',
                                 'Procellaria (Procellaria) parkinsoni' = 'Procellaria parkinsoni',
                                 'Procellaria (Procellaria) westlandica' = 'Procellaria westlandica',
                                 'Psephotus' = 'Psephotus haematonotus',
                                 'Pseudobulweria rostrata rostrata' = 'Pseudobulweria rostrata',
                                 'Psephotus haematonotus haematonotus' = 'Psephotus haematonotus',
                                 'Psitteuteles' = 'Psitteuteles versicolor',
                                 'Psophodes (Sphenostoma) cristatus' = 'Psophodes cristatus',
                                 'Psophodes (Phodopses) leucogaster' = 'Psophodes leucogaster',
                                 'Psophodes (Phodopses) leucogaster lashmari' = 'Psophodes leucogaster',
                                 'Psophodes (Phodopses) leucogaster leucogaster' = 'Psophodes leucogaster',
                                 'Psophodes (Phodopses) nigrogularis' = 'Psophodes nigrogularis',
                                 'Psophodes (Phodopses) nigrogularis nigrogularis' = 'Psophodes nigrogularis',
                                 'Psophodes (Sphenostoma) occidentalis' = 'Psophodes occidentalis',
                                 'Psophodes (Psophodes) olivaceus' = 'Psophodes olivaceus',
                                 'Psophodes (Psophodes) olivaceus lateralis' = 'Psophodes olivaceus',
                                 'Pterodroma (Aestrelata) cervicalis' = 'Pterodroma cervicalis',
                                 'Pterodroma (Pterodroma) gouldi' = 'Pterodroma gouldi',
                                 'Pterodroma (Hallstroma) heraldica' = 'Pterodroma heraldica',
                                 'Pterodroma (Pterodroma) lessonii' = 'Pterodroma lessonii',
                                 'Pterodroma (Cookilaria) leucoptera' = 'Pterodroma leucoptera',
                                 'Pterodroma (Pterodroma) macroptera' = 'Pterodroma macroptera',
                                 'Pterodroma (Hallstroma) neglecta' = 'Pterodroma neglecta',
                                 'Pterodroma (Hallstroma) neglecta neglecta' = 'Pterodroma neglecta',
                                 'Pterodroma (Cookilaria) nigripennis' = 'Pterodroma nigripennis',
                                 'Pterodroma (Pterodroma) solandri' = 'Pterodroma solandri',
                                 'Ptilinopus magnificus assimilis' = 'Ptilinopus magnificus',
                                 'Ptilinopus regina ewingii' = 'Ptilinopus regina',
                                 'Ptilinopus regina regina' = 'Ptilinopus regina',
                                 'Ptilinopus superbus superbus' = 'Ptilinopus superbus',
                                 'Ptiloris (Craspedophora) magnificus' = 'Ptiloris magnificus',
                                 'Ptiloris (Ptiloris) paradiseus' = 'Ptiloris paradiseus',
                                 'Ptiloris (Ptiloris) victoriae' = 'Ptiloris victoriae',
                                 'Ptilotula penicillata leilavalensis' = 'Ptilotula penicillata',
                                 'Ptilotula penicillata penicillata' = 'Ptilotula penicillata',
                                 'Ptilotula plumula graingeri' = 'Ptilotula plumula',
                                 'Ptilotula plumula plumula' = 'Ptilotula plumula',
                                 'Puffinus (Puffinus) assimilis' = 'Puffinus assimilis',
                                 'Puffinus (Puffinus) assimilis assimilis' = 'Puffinus assimilis',
                                 'Puffinus (Puffinus) assimilis tunneyi' = 'Puffinus assimilis',
                                 'Puffinus (Puffinus) gavia' = 'Puffinus gavia',
                                 'Puffinus (Puffinus) huttoni' = 'Puffinus huttoni',
                                 'Pycnonotus (Pycnonotus) jocosus' = 'Pycnonotus jocosus',
                                 'Radjah radjah rufitergum' = 'Radjah radjah',
                                 'Rallina (Rallina) tricolor' = 'Rallina tricolor',
                                 'Recurvirostra' = 'Recurvirostra novaehollandiae',
                                 'Rhipidura (Rhipidura) albiscapa' = 'Rhipidura albiscapa',
                                 'Rhipidura (Rhipidura) albiscapa albiscapa' = 'Rhipidura albiscapa',
                                 'Rhipidura (Rhipidura) albiscapa alisteri' = 'Rhipidura albiscapa',
                                 'Rhipidura (Rhipidura) albiscapa preissi' = 'Rhipidura albiscapa',
                                 'Rhipidura (Howeavis) dryas' = 'Rhipidura dryas',
                                 'Rhipidura (Howeavis) dryas dryas' = 'Rhipidura dryas',
                                 'Rhipidura (Rhipidura) fuliginosa' = 'Rhipidura fuliginosa',
                                 'Rhipidura (Rhipidura) fuliginosa cervina' = 'Rhipidura fuliginosa',
                                 'Rhipidura (Sauloprocta) leucophrys' = 'Rhipidura leucophrys',
                                 'Rhipidura (Sauloprocta) leucophrys leucophrys' = 'Rhipidura leucophrys',
                                 'Rhipidura (Rhipidura) phasiana' = 'Rhipidura phasiana',
                                 'Rhipidura (Howeavis) rufifrons' = 'Rhipidura rufifrons',
                                 'Rhipidura (Howeavis) rufifrons intermedia' = 'Rhipidura rufifrons',
                                 'Rhipidura (Howeavis) rufifrons rufifrons' = 'Rhipidura rufifrons',
                                 'Rhipidura (Setosura) rufiventris' = 'Rhipidura rufiventris',
                                 'Rhipidura (Setosura) rufiventris isura' = 'Rhipidura rufiventris',
                                 'Scythrops novaehollandiae novaehollandiae' = 'Scythrops novaehollandiae',
                                 'Sericornis (Arfakornis) beccarii' = 'Sericornis beccarii',
                                 'Sericornis (Sericornis) frontalis' = 'Sericornis frontalis',
                                 'Sericornis (Sericornis) frontalis ashbyi' = 'Sericornis frontalis',
                                 'Sericornis (Sericornis) frontalis balstoni' = 'Sericornis frontalis',
                                 'Sericornis (Sericornis) frontalis frontalis' = 'Sericornis frontalis',
                                 'Sericornis (Sericornis) frontalis maculatus' = 'Sericornis frontalis',
                                 'Sericornis (Sericornis) frontalis mellori' = 'Sericornis frontalis',
                                 'Sericornis (Sericornis) frontalis rosinae' = 'Sericornis frontalis',
                                 'Sericornis (Arfakornis) magnirostra' = 'Sericornis magnirostra',
                                 'Sericornis (Arfakornis) magnirostra magnirostra' = 'Sericornis magnirostra',
                                 'Sericulus (Sericulus) chrysocephalus' = 'Sericulus chrysocephalus',
                                 'Smicrornis brevirostris brevirostris' = 'Smicrornis brevirostris',
                                 'Smicrornis brevirostris occidentalis' = 'Smicrornis brevirostris',
                                 'Smicrornis' = 'Smicrornis brevirostris',
                                 'Sphecotheres' = 'Sphecotheres vieilloti',
                                 'Sphecotheres vieilloti ashbyi' = 'Sphecotheres vieilloti',
                                 'Sphecotheres vieilloti flaviventris' = 'Sphecotheres vieilloti',
                                 'Sphecotheres vieilloti vieilloti' = 'Sphecotheres vieilloti',
                                 'Spilopelia chinensis tigrina' = 'Spilopelia chinensis',
                                 'Stagonopleura (Zonaeginthus) bella' = 'Stagonopleura bella',
                                 'Stagonopleura (Zonaeginthus) bella interposita' = 'Stagonopleura bella',
                                 'Stagonopleura (Zonaeginthus) bella samueli' = 'Stagonopleura bella',
                                 'Stagonopleura (Stagonopleura) guttata' = 'Stagonopleura guttata',
                                 'Stagonopleura (Zonaeginthus) oculata' = 'Stagonopleura oculata',
                                 'Stercorarius antarcticus lonnbergi' = 'Stercorarius antarcticus',
                                 'Sterna (Sterna) dougallii' = 'Sterna dougallii',
                                 'Sterna (Sterna) dougallii bangsi' = 'Sterna dougallii',
                                 'Sterna (Sterna) dougallii gracilis' = 'Sterna dougallii',
                                 'Sterna (Sterna) hirundo' = 'Sterna hirundo',
                                 'Sterna (Sterna) hirundo longipennis' = 'Sterna hirundo',
                                 'Sterna (Sterna) paradisaea' = 'Sterna paradisaea',
                                 'Sterna (Sterna) striata' = 'Sterna striata',
                                 'Sterna (Gygisterna) sumatrana' = 'Sterna sumatrana',
                                 'Sterna (Gygisterna) sumatrana sumatrana' = 'Sterna sumatrana',
                                 'Sterna (Sterna) vittata' = 'Sterna vittata',
                                 'Sternula' = 'Sternula albifrons',
                                 'Sternula albifrons sinensis' = 'Sternula albifrons',
                                 'Sternula nereis nereis' = 'Sternula nereis',
                                 'Sternula nereis exsul' = 'Sternula nereis',
                                 'Stipiturus malachurus halmaturinus' = 'Stipiturus malachurus',
                                 'Stipiturus malachurus intermedius' = 'Stipiturus malachurus',
                                 'Stipiturus malachurus parimeda' = 'Stipiturus malachurus',
                                 'Stipiturus malachurus polionotum' = 'Stipiturus malachurus',
                                 'Stipiturus malachurus hartogi' = 'Stipiturus malachurus',
                                 'Stizoptera bichenovii annulosa' = 'Stizoptera bichenovii',
                                 'Strepera (Strepera) fuliginosa' = 'Strepera fuliginosa',
                                 'Strepera (Strepera) graculina crissalis' = 'Strepera graculina',
                                 'Strepera (Strepera) graculina graculina' = 'Strepera graculina',
                                 'Strepera (Strepera) graculina' = 'Strepera graculina',
                                 'Strepera (Strepera) graculina robinsoni' = 'Strepera graculina',
                                 'Strepera (Neostrepera) versicolor' = 'Strepera versicolor',
                                 'Strepera (Neostrepera) versicolor halmaturina' = 'Strepera versicolor',
                                 'Strepera (Neostrepera) versicolor intermedia' = 'Strepera versicolor',
                                 'Strepera (Neostrepera) versicolor melanoptera' = 'Strepera versicolor',
                                 'Strepera (Neostrepera) versicolor versicolor' = 'Strepera versicolor',
                                 'Struthidea cinerea cinerea' = 'Struthidea cinerea',
                                 'Sula dactylatra personata' = 'Sula dactylatra',
                                 'Sula dactylatra tasmani' = 'Sula dactylatra',
                                 'Sula leucogaster plotus' = 'Sula leucogaster',
                                 'Sula sula rubripes' = 'Sula sula',
                                 'Symposiachrus' = 'Symposiachrus trivirgatus',
                                 'Symposiachrus trivirgatus albiventris' = 'Symposiachrus trivirgatus',
                                 'Synoicus chinensis victoriae' = 'Synoicus chinensis',
                                 'Synoicus ypsilophora australis' = 'Synoicus ypsilophora',
                                 'Synoicus ypsilophora ypsilophora' = 'Synoicus ypsilophora',
                                 'Tachybaptus' = 'Tachybaptus novaehollandiae',
                                 'Tachybaptus novaehollandiae novaehollandiae' = 'Tachybaptus novaehollandiae',
                                 'Tadorna (Casarca) tadornoides' = 'Tadorna tadornoides',
                                 'Taeniopygia' = 'Taeniopygia guttata',
                                 'Taeniopygia guttata castanotis' = 'Taeniopygia guttata',
                                 'Tanysiptera (Uralcyon) sylvia' = 'Tanysiptera sylvia',
                                 'Tanysiptera (Uralcyon) sylvia sylvia' = 'Tanysiptera sylvia',
                                 'Thalassarche bulleri bulleri' = 'Thalassarche bulleri',
                                 'Thalassarche cauta cauta' = 'Thalassarche cauta',
                                 'Thalasseus bengalensis torresii' = 'Thalasseus bengalensis',
                                 'Thalasseus bergii cristatus' = 'Thalasseus bergii',
                                 'Thinornis' = 'Thinornis cucullatus',
                                 'Thinornis cucullatus cucullatus' = 'Thinornis cucullatus',
                                 'Threskiornis moluccus moluccus' = 'Threskiornis moluccus',
                                 'Todiramphus (Todiramphus) chloris sordidus' = 'Todiramphus chloris',
                                 'Todiramphus (Todiramphus) chloris' = 'Todiramphus chloris',
                                 'Todiramphus (Lazulena) macleayii' = 'Todiramphus macleayii',
                                 'Todiramphus (Lazulena) macleayii incinctus' = 'Todiramphus macleayii',
                                 'Todiramphus (Cyanalcyon) pyrrhopygius' = 'Todiramphus pyrrhopygius',
                                 'Todiramphus (Todiramphus) sanctus' = 'Todiramphus sanctus',
                                 'Todiramphus (Todiramphus) sanctus sanctus' = 'Todiramphus sanctus',
                                 'Todiramphus (Todiramphus) sanctus vagans' = 'Todiramphus sanctus',
                                 'Tregellasia capito capito' = 'Tregellasia capito',
                                 "Trichoglossus" = 'Trichoglossus haematodus',
                                 'Trichoglossus haematodus moluccanus' = 'Trichoglossus haematodus',
                                 'Trichoglossus haematodus rubritorquis' = 'Trichoglossus haematodus',
                                 'Tringa (Heteroscelus) brevipes' = 'Tringa brevipes',
                                 'Tringa (Rhyacophilus) glareola' = 'Tringa glareola',
                                 'Tringa (Heteroscelus) incana' = 'Tringa incana',
                                 'Tringa (Glottis) nebularia' = 'Tringa nebularia',
                                 'Tringa (Rhyacophilus) stagnatilis' = 'Tringa stagnatilis',
                                 'Tringa (Totanus) totanus' = 'Tringa totanus',
                                 'Turdus merula merula' = 'Turdus merula',
                                 'Turdus poliocephalus poliocephalus' = 'Turdus poliocephalus',
                                 'Turdus poliocephalus vinitinctus' = 'Turdus poliocephalus',
                                 'Turnix (Austroturnix) castanotus' = 'Turnix castanotus',
                                 'Turnix (Ortygodes) maculosus' = 'Turnix maculosus',
                                 'Turnix (Austroturnix) melanogaster' = 'Turnix melanogaster',
                                 'Turnix (Ortygodes) maculosus melanotus' = 'Turnix maculosus',
                                 'Turnix (Austroturnix) varius' = 'Turnix varius',
                                 'Turnix (Austroturnix) varius varius' = 'Turnix varius',
                                 'Turnix (Alphaturnia) velox' = 'Turnix velox',
                                 'Tyto javanica delicatula' = 'Tyto javanica',
                                 'Tyto novaehollandiae castanops' = 'Tyto novaehollandiae',
                                 'Tyto novaehollandiae kimberli' = 'Tyto novaehollandiae',
                                 'Tyto novaehollandiae melvillensis' = 'Tyto novaehollandiae',
                                 'Tyto novaehollandiae novaehollandiae' = 'Tyto novaehollandiae',
                                 'Tyto tenebricosa tenebricosa' = 'Tyto tenebricosa',
                                 'Vanellus (Lobipluvia) miles' = 'Vanellus miles',
                                 'Vanellus (Lobipluvia) miles miles' = 'Vanellus miles',
                                 'Vanellus (Lobipluvia) miles novaehollandiae' = 'Vanellus miles',
                                 'Vanellus (Lobivanellus) tricolor' = 'Vanellus tricolor',
                                 'Xenus' = 'Xenus cinereus',
                                 'Zanda funerea whiteae' = 'Zanda funerea',
                                 'Zanda funerea xanthanota' = 'Zanda funerea',
                                 'Zapornia pusilla palustris' = 'Zapornia pusilla',
                                 'Zoothera (Zoothera) heinei' = 'Zoothera heinei',
                                 'Zoothera (Zoothera) lunulata' = 'Zoothera lunulata',
                                 'Zoothera (Zoothera) lunulata halmaturina' = 'Zoothera lunulata',
                                 'Zoothera (Zoothera) lunulata lunulata' = 'Zoothera lunulata',
                                 'Zosterops lateralis chloronotus' = 'Zosterops lateralis',
                                 'Zosterops lateralis chlorocephalus' = 'Zosterops lateralis',
                                 'Zosterops lateralis cornwalli' = 'Zosterops lateralis',
                                 'Zosterops lateralis tephropleurus' = 'Zosterops lateralis',
                                 'Zosterops lateralis lateralis' = 'Zosterops lateralis',
                                 'Zosterops lateralis vegetus' = 'Zosterops lateralis',
                                 'Zosterops lateralis westernensis' = 'Zosterops lateralis',
                                 'Zosterops luteus balstoni' = 'Zosterops luteus',
                                 'Zosterops luteus luteus' = 'Zosterops luteus',
                                 'Zosterops lateralis pinarochrous' = 'Zosterops lateralis'))

# remove genera unidentifiable to species
avala.dat5 <- avala.dat4[!(avala.dat4$scientific %in% c("Acanthiza",
                                                            "Accipiter",
                                                            "Acrocephalus",
                                                            "Amaurornis",
                                                            "Anas",
                                                            "Anas (Anas)",
                                                            "Anous",
                                                            "Anthochaera",
                                                            "Aphelocephala",
                                                            "Aplonis (Aplonis)",
                                                            "Ardea",
                                                            "Ardenna",
                                                            "Artamus",
                                                            "Cacatua",
                                                            'Cacatua (Cacatua)',
                                                            "Cacomantis",
                                                            "Calamanthus",
                                                            "Calidris",
                                                            "Calyptorhynchus",
                                                            "Chalcites",
                                                            "Charadrius",
                                                            "Chlidonias",
                                                            "Cincloramphus",
                                                            "Cinclosoma",
                                                            "Circus",
                                                            "Cisticola",
                                                            "Climacteris",
                                                            "Conopophila",
                                                            "Coracina",
                                                            "Corvus",
                                                            "Coturnix",
                                                            "Coturnix (Coturnix)",
                                                            "Cracticus",
                                                            "Dacelo",
                                                            "Diomedea",
                                                            "Ducula",
                                                            "Egretta",
                                                            "Elanus",
                                                            "Eopsaltria",
                                                            "Epthianura",
                                                            "Eurostopodus",
                                                            "Falco",
                                                            "Falco (Falco)",
                                                            "Fregata",
                                                            "Gallinago",
                                                            "Geopelia",
                                                            "Geophaps (Lophophaps)",
                                                            "Gerygone",
                                                            "Haematopus",
                                                            "Haliastur",
                                                            "Hirundo",
                                                            "Hirundo (Hirundo)",
                                                            "Hydrobates",
                                                            "Hydroprogne",
                                                            "Hylacola",
                                                            "Lalage",
                                                            "Larus",
                                                            "Lichenostomus",
                                                            "Limosa",
                                                            "Lonchura (Padda)",
                                                            "Macronectes",
                                                            "Malurus",
                                                            "Manorina",
                                                            "Meliphaga",
                                                            "Meliphaga (Microptilotis)",
                                                            "Microeca",
                                                            "Monarcha",
                                                            "Myiagra",
                                                            "Neochmia",
                                                            "Neochmia (Neochmia)",
                                                            "Neophema",
                                                            "Ninox",
                                                            "Numenius",
                                                            "Onychoprion",
                                                            "Pachycephala",
                                                            "Pachyptila",
                                                            "Pardalotus",
                                                            "Petrochelidon",
                                                            "Petroica",
                                                            "Phalacrocorax",
                                                            "Phaps",
                                                            "Phaps (Phaps)",
                                                            "Philemon",
                                                            "Phylidonyris",
                                                            "Pitta (Erythropitta)",
                                                            "Ptilinopus",
                                                            "Ptilotula",
                                                            "Platalea",
                                                            "Platycercus",
                                                            "Pluvialis",
                                                            "Poecilodryas",
                                                            "Poephila",
                                                            "Polytelis",
                                                            "Pomatostomus",
                                                            "Pterodroma",
                                                            "Puffinus",
                                                            "Rhipidura",
                                                            "Sericornis",
                                                            "Spilopelia",
                                                            "Stagonopleura",
                                                            "Stercorarius",
                                                            "Sterna",
                                                            "Sterna (Sterna)",
                                                            "Stipiturus",
                                                            "Strepera",
                                                            "Synoicus",
                                                            "Thalassarche",
                                                            "Thalasseus",
                                                            "Threskiornis",
                                                            "Todiramphus",
                                                            "Tringa",
                                                            "Tringa (Heteroscelus)",
                                                            "Turdus",
                                                            "Turnix",
                                                            "Tyto",
                                                            "Zanda",
                                                            "Zoothera",
                                                            "Zosterops")),]


avspp <- table(avala.dat5$scientific)
avspp
avspp.names <- attr(avspp, "names")
avspp.names

ausisles.avala.df <- as.data.frame(avala.dat5)
head(ausisles.avala.df)
ausisles.avala.spp <- ausisles.avala.df[,c("scientific", "FID")]
head(ausisles.avala.spp)
ausisles.avala.spp.unique <- unique(ausisles.avala.spp)
head(ausisles.avala.spp.unique)
ausisles.avala.spp.unique.tab <- as.data.frame(table(ausisles.avala.spp.unique$FID)) # number of unique species per island)
ausisles.avala.spp.unique.tab
head(ausisles.avala.spp.unique.tab)
colnames(ausisles.avala.spp.unique.tab) <- c("FID", "avSR")
head(ausisles.avala.spp.unique.tab)



# reptiles
ausisles.reptala.df <- as.data.frame(ausisles.reptala)
head(ausisles.reptala.df)
table(ausisles.reptala.df$scientific)

## species list
repspp <- table(ausisles.reptala.df$scientific)
repspp

## remove invasive/escapee species
reptala.dat1 <- ausisles.reptala.df[!(ausisles.reptala.df$scientific %in% c("Boa",
                                                                "Boa constrictor",
                                                                "Calotes",
                                                                "Chamaeleo",
                                                                "Gekko",
                                                                "Gekko gecko",
                                                                "Iguana iguana",
                                                                "Iguanidae",
                                                                "Testudo")),]

## remove higher-taxonomic level ambiguities
reptala.dat2 <- reptala.dat1[!(reptala.dat1$scientific %in% c("ARCHOSAURIA",
                                                                "AGAMIDAE",
                                                                "Boidae",
                                                                "CARETTOCHELYIDAE",
                                                                "CARPHODACTYLIDAE",
                                                                "CHELIDAE",
                                                                "CHELONIIDAE",
                                                                "COLUBRIDAE",
                                                                "CROCODYLIA",
                                                                "CROCODYLIDAE",
                                                                "DERMOCHELYIDAE",
                                                                "DIPLODACTYLIDAE",
                                                                "ELAPIDAE",
                                                                "EMYDIDAE",
                                                                "GEKKONIDAE",
                                                                "GEKKOTA",
                                                                "Geoemydidae",
                                                                "HOMALOPSIDAE",
                                                                "Hydrophiinae",
                                                                "LACERTILIA",
                                                                "LEPIDOSAURIA",
                                                                "PYGOPODIDAE",
                                                                "PYTHONIDAE",
                                                                "REPTILIA",
                                                                "SCINCIDAE",
                                                                "SERPENTES",
                                                                "SQUAMATA",
                                                                "TESTUDINES",
                                                                "Testudinidae",
                                                                "TYPHLOPIDAE",
                                                                "VARANIDAE",
                                                                "Viperidae")),]

reptala.dat3 <- reptala.dat2 %>%
  mutate(scientific = recode(scientific,
                                 "Acanthophis" = "Acanthophis antarcticus",
                                 "Anepischetosia" = "Anepischetosia maccoyi",
                                 "Antaresia maculosa maculosa" = "Antaresia maculosa",
                                 "Antaresia maculosa pensinsularis" = "Antaresia maculosa",
                                 "Boiga" = "Boiga irregularis",
                                 "Brachyurophis fasciolatus fasciatus" = "Brachyurophis fasciolatus",
                                 "Brachyurophis fasciolatus fasciolatus" = "Brachyurophis fasciolatus",
                                 "Candoia" = "Candoia carinata",
                                 "Caretta" = "Caretta caretta",
                                 "Caretta caretta gigas" = "Caretta caretta",
                                 "Cerberus" = "Cerberus australis",
                                 "Chelodina (Chelydera) burrungandjii" = "Chelodina burrungandjii",
                                 "Chelodina (Chelodina) canni" = "Chelodina canni",
                                 "Chelodina (Chelydera) kuchlingi" = "Chelodina kuchlingi",
                                 "Chelodina (Chelodina) longicollis" = "Chelodina longicollis",
                                 "Chelodina (Chelydera) expansa" = "Chelodina expansa",
                                 "Chelodina (Macrochelodina) oblonga" = "Chelodina oblonga",
                                 "Chelodina (Chelydera) rugosa" = "Chelodina rugosa",
                                 "Chelodina (Chelodina) steindachneri" = "Chelodina steindachneri",
                                 "Chelonia" = "Chelonia mydas",
                                 "Chelonia mydas japonica" = "Chelonia mydas",
                                 "Chlamydosaurus" = "Chlamydosaurus kingii",
                                 "Christinus" = "Christinus marmoratus",
                                 "Cryptoblepharus litoralis horneri" = "Cryptoblepharus litoralis",
                                 "Cryptoblepharus litoralis litoralis" = "Cryptoblepharus litoralis",
                                 "Cryptoblepharus pulcher clarus" = "Cryptoblepharus pulcher",
                                 "Cryptoblepharus pulcher pulcher" = "Cryptoblepharus pulcher",
                                 "Ctenophorus isolepis gularis" = "Ctenophorus isolepis",
                                 "Ctenophorus isolepis isolepis" = "Ctenophorus isolepis",
                                 "Ctenophorus isolepis citrinus" = "Ctenophorus isolepis",
                                 "Ctenophorus maculatus griseus" = "Ctenophorus maculatus",
                                 "Ctenophorus maculatus dualis" = "Ctenophorus maculatus",
                                 "Ctenophorus maculatus maculatus" = "Ctenophorus maculatus",
                                 "Ctenotus decaneurus yampiensis" = "Ctenotus decaneurus",
                                 "Ctenotus decaneurus decaneurus" = "Ctenotus decaneurus",
                                 "Ctenotus gemmula (Swan Coastal Plain subpopulation)" = "Ctenotus gemmula",
                                 "Ctenotus grandis grandis" = "Ctenotus grandis",
                                 "Ctenotus grandis titan" = "Ctenotus grandis",
                                 "Ctenotus hebetior hebetior" = "Ctenotus hebetior",
                                 "Ctenotus pantherinus ocellifer" = "Ctenotus pantherinus",
                                 "Ctenotus pantherinus acripes" = "Ctenotus pantherinus",
                                 "Ctenotus pantherinus calx" = "Ctenotus pantherinus",
                                 "Ctenotus pantherinus pantherinus" = "Ctenotus pantherinus",
                                 "Ctenotus rimacola camptris" = "Ctenotus rimacola",
                                 "Ctenotus rimacola rimacola" = "Ctenotus rimacola",
                                 "Ctenotus strauchii strauchii" = "Ctenotus strauchii",
                                 "Ctenotus uber uber" = "Ctenotus uber",
                                 "Ctenotus uber johnstonei" = "Ctenotus uber",
                                 "Cyclodomorphus melanops elongatus" = "Cyclodomorphus melanops",
                                 "Cyclodomorphus melanops melanops" = "Cyclodomorphus melanops",
                                 "Cyclodomorphus melanops siticulosus" = "Cyclodomorphus melanops",
                                 "Delma concinna concinna" = "Delma concinna",
                                 "Demansia reticulata cupreiceps" = "Demansia reticulata",
                                 "Dermochelys" = "Dermochelys coriacea",
                                 "Diplodactylus granariensis granariensis" = "Diplodactylus granariensis",
                                 "Diplodactylus granariensis rex" = "Diplodactylus granariensis",
                                 "Egernia saxatilis saxatilis" = "Egernia saxatilis",
                                 "Egernia saxatilis intermedia" = "Egernia saxatilis",
                                 "Egernia stokesii badia" = "Egernia stokesii",
                                 "Egernia stokesii zellingi" = "Egernia stokesii",
                                 "Egernia stokesii stokesii" = "Egernia stokesii",
                                 "Elseya (Pelocomastes) albagula" = "Elseya albagula",
                                 "Elseya (Elseya) dentata" = "Elseya dentata",
                                 "Elseya (Elseya) flaviventralis" = "Elseya flaviventralis",
                                 "Elseya (Pelocomastes) irwini" = "Elseya irwini",
                                 "Elseya (Pelocomastes) lavarackorum" = "Elseya lavarackorum",
                                 "Emoia atrocostata atrocostata" = "Emoia atrocostata",
                                 "Emoia atrocostata australis" = "Emoia atrocostata",
                                 "Emydura subglobosa subglobosa" = "Emydura subglobosa",
                                 "Emydura subglobosa angkibaanya" = "Emydura subglobosa",
                                 "Emydura subglobosa worrelli" = "Emydura subglobosa",
                                 "Eremiascincus pardalis pardalis" = "Eremiascincus pardalis",
                                 "Eretmochelys" = "Eretmochelys imbricata",
                                 "Eretmochelys imbricata bissa" = "Eretmochelys imbricata", 
                                 "Eretmochelys imbricata squamata" = "Eretmochelys imbricata",
                                 "Eulamprus tympanum marnieae" = "Eulamprus tympanum",
                                 "Eulamprus tympanum tympanum" = "Eulamprus tympanum",
                                 "Gowidon" = "Gowidon longirostris",
                                 "Harrisoniascincus" = "Harrisoniascincus zia",
                                 "Hemidactylus" = "Hemidactylus frenatus",
                                 "Hemiergis decresiensis continentis" = "Hemiergis decresiensis",
                                 "Hemiergis decresiensis decresiensis" = "Hemiergis decresiensis",
                                 "Hemiergis initialis brookeri" = "Hemiergis initialis",
                                 "Hemiergis initialis initialis" = "Hemiergis initialis",
                                 "Hemiergis talbingoensis talbingoensis" = "Hemiergis talbingoensis",
                                 "Hemiergis talbingoensis davisi" = "Hemiergis talbingoensis",
                                 "Hydrophis platurus platurus" = "Hydrophis platurus",
                                 "Intellagama" = "Intellagama lesueurii",
                                 "Intellagama lesueurii lesueurii" = "Intellagama lesueurii",
                                 "Intellagama lesueurii howitti" = "Intellagama lesueurii",
                                 "Lerista microtis schwaneri" = "Lerista microtis",
                                 "Lerista macropisthopus fusciceps" = "Lerista macropisthopus",
                                 "Lerista macropisthopus galea" = "Lerista macropisthopus",
                                 "Lerista macropisthopus remota" = "Lerista macropisthopus",
                                 "Lerista macropisthopus macropisthopus" = "Lerista macropisthopus",
                                 "Lerista microtis intermedia" = "Lerista microtis",
                                 "Lerista microtis microtis" = "Lerista microtis",
                                 "Lerista planiventralis decora" = "Lerista planiventralis",
                                 "Lerista planiventralis maryani" = "Lerista planiventralis",
                                 "Lerista planiventralis planiventralis" = "Lerista planiventralis",
                                 "Liasis olivaceus olivaceus" = "Liasis olivaceus",
                                 "Liasis olivaceus barroni" = "Liasis olivaceus",
                                 "Liopholis margaretae margaretae" = "Liopholis margaretae",
                                 "Liopholis pulchra longicauda" = "Liopholis pulchra",
                                 "Liopholis pulchra pulchra" = "Liopholis pulchra",
                                 "Liopholis slateri virgata" = "Liopholis slateri",
                                 "Liopholis slateri slateri" = "Liopholis slateri",
                                 "Menetia" = "Menetia greyii",
                                 "Menetia surda cresswelli" = "Menetia surda",
                                 "Menetia surda surda" = "Menetia surda",
                                 "Morelia spilota variegata" = "Morelia spilota",
                                 "Morelia spilota spilota" = "Morelia spilota",
                                 "Morethia ruficauda exquisita" = "Morethia ruficauda",
                                 "Morethia ruficauda ruficauda" = "Morethia ruficauda",
                                 "Nebulifera" = "Nebulifera robusta",
                                 "Nephrurus levis levis" = "Nephrurus levis",
                                 "Nephrurus levis occidentalis" = "Nephrurus levis",
                                 "Nephrurus levis pilbarensis" = "Nephrurus levis",
                                 "Nephrurus wheeleri wheeleri" = "Nephrurus wheeleri",
                                 "Notechis" = "Notechis scutatus",
                                 "Notoscincus ornatus ornatus" = "Notoscincus ornatus",
                                 "Notoscincus ornatus wotjulum"  = "Notoscincus ornatus",
                                 "Oligosoma" = "Oligosoma lichenigerum",
                                 "Orraya" = "Orraya occultus",
                                 "Oxyuranus scutellatus scutellatus" = "Oxyuranus scutellatus",
                                 "Phyllurus ossa hobsoni" = "Phyllurus ossa",
                                 "Phyllurus ossa tamoya" = "Phyllurus ossa",
                                 "Pletholax gracilis gracilis" = "Pletholax gracilis",
                                 "Pogona minor minor" = "Pogona minor",
                                 "Pogona minor minima" = "Pogona minor",
                                 "Pogona minor mitchelli" = "Pogona minor",
                                 "Pseudemydura" = "Pseudemydura umbrina",
                                 "Rankinia" = "Rankinia diemensis",
                                 "Rankinia diemensis (Grampians)" = "Rankinia diemensis",
                                 "Rheodytes" = "Rheodytes leukops",
                                 "Rhinoplocephalus" = "Rhinoplocephalus bicolor",
                                 "Saiphos" = "Saiphos equalis",
                                 "Simoselaps" = "Simoselaps bertholdi",
                                 "Strophurus ciliaris aberrans" = "Strophurus ciliaris",
                                 "Strophurus ciliaris ciliaris" = "Strophurus ciliaris",
                                 "Strophurus taenicauda taenicauda" = "Strophurus taenicauda",
                                 "Strophurus taenicauda albiocularis" = "Strophurus taenicauda",
                                 "Strophurus taenicauda triaureus" = "Strophurus taenicauda",
                                 "Strophurus spinigerus inornatus" = "Strophurus spinigerus",
                                 "Strophurus spinigerus spinigerus" = "Strophurus spinigerus",
                                 "Tiliqua rugosa aspera" = "Tiliqua rugosa",
                                 "Tiliqua rugosa konowi" = "Tiliqua rugosa",
                                 "Tiliqua rugosa rugosa" = "Tiliqua rugosa",
                                 "Tiliqua rugosa palarra" = "Tiliqua rugosa",
                                 "Tiliqua scincoides intermedia" = "Tiliqua scincoides",
                                 "Tiliqua scincoides scincoides" = "Tiliqua scincoides",
                                 "Trachemys scripta elegans" = "Trachemys scripta",
                                 "Tropicagama" = "Tropicagama temporalis",
                                 "Tropidechis" = "Tropidechis carinatus",
                                 "Tropidonophis mairii mairii" = "Tropidonophis mairii",
                                 "Varanus panoptes rubidus" = "Varanus panoptes",
                                 "Varanus panoptes panoptes" = "Varanus panoptes"))

reptala.dat4 <- reptala.dat3[!(reptala.dat3$scientific %in% c("Acritoscincus",
                                                                  "Acrochordus",
                                                                  "Aipysurus",
                                                                  "Amalosia",
                                                                  "Amphibolurus",
                                                                  "Anomalopus",
                                                                  "Anilios",
                                                                  "Antaresia",
                                                                  "Aprasia",
                                                                  "Aspidites",
                                                                  "Austrelaps",
                                                                  "Bellatorias",
                                                                  "Brachyurophis",
                                                                  "Cacophis",
                                                                  "Calyptotis",
                                                                  "Carinascincus",
                                                                  "Carlia",
                                                                  "Chelodina",
                                                                  "Chelodina (Chelodina)",
                                                                  "Chelodina (Macrochelodina)",
                                                                  "Concinnia",
                                                                  "Crocodylus",
                                                                  "Cryptoblepharus",
                                                                  "Cryptophis",
                                                                  "Ctenophorus",
                                                                  "Ctenotus",
                                                                  "Cyclodomorphus",
                                                                  "Cyrtodactylus",
                                                                  "Delma",
                                                                  "Demansia",
                                                                  "Denisonia",
                                                                  "Dendrelaphis",
                                                                  "Diporiphora",
                                                                  "Drysdalia",
                                                                  "Diplodactylus",
                                                                  "Egernia",
                                                                  "Elseya (Elseya)",
                                                                  "Elseya",
                                                                  "Emoia",
                                                                  "Emydocephalus",
                                                                  "Eremiascincus",
                                                                  "Eugongylus",
                                                                  "Eulamprus",
                                                                  "Furina",
                                                                  "Gehyra",
                                                                  "Glaphyromorphus",
                                                                  "Hemiaspis",
                                                                  "Hemiergis",
                                                                  "Heteronotia",
                                                                  "Hoplocephalus",
                                                                  "Hydrophis",
                                                                  "Indotyphlops",
                                                                  "Lampropeltis",
                                                                  "Lampropholis",
                                                                  "Lepidodactylus",
                                                                  "Lerista",
                                                                  "Lialis",
                                                                  "Liasis",
                                                                  "Liburnascincus",
                                                                  "Liopholis",
                                                                  "Lophognathus",
                                                                  "Lophosaurus",
                                                                  "Lucasium",
                                                                  "Lycodon",
                                                                  "Lygisaurus",
                                                                  "Morelia",
                                                                  "Morethia",
                                                                  "Myron",
                                                                  "Myuchelys",
                                                                  "Nactus",
                                                                  "Nephrurus",
                                                                  "Notoscincus",
                                                                  "Oedura",
                                                                  "Oxyuranus",
                                                                  "Phyllurus",
                                                                  "Pogona",
                                                                  "Proablepharus",
                                                                  "Pseudechis",
                                                                  "Pseudemoia",
                                                                  "Pseudonaja",
                                                                  "Pseudothecadactylus",
                                                                  "Pygopus",
                                                                  "Python",
                                                                  "Ramphotyphlops",
                                                                  "Rhynchoedura",
                                                                  "Saltuarius",
                                                                  "Saproscincus",
                                                                  "Silvascincus",
                                                                  "Simalia",
                                                                  "Species Inquirenda",
                                                                  "Stegonotus",
                                                                  "Strophurus",
                                                                  "Suta",
                                                                  "Tiliqua",
                                                                  "Tropidonophis",
                                                                  "Tympanocryptis",
                                                                  "Underwoodisaurus",
                                                                  "Varanus",
                                                                  "Vermicella")),]

repspp <- table(reptala.dat4$scientific)
repspp
repspp.names <- attr(repspp, "names")
repspp.names
length(repspp.names)
repspp.names

ausisles.reptala.df <- as.data.frame(reptala.dat4)
head(ausisles.reptala.df)
ausisles.reptala.spp <- ausisles.reptala.df[,c("scientific", "FID")]
head(ausisles.reptala.spp)
ausisles.reptala.spp.unique <- unique(ausisles.reptala.spp)
head(ausisles.reptala.spp.unique)
ausisles.reptala.spp.unique.tab <- as.data.frame(table(ausisles.reptala.spp.unique$FID)) # number of unique species per island)
ausisles.reptala.spp.unique.tab
head(ausisles.reptala.spp.unique.tab)
colnames(ausisles.reptala.spp.unique.tab) <- c("FID", "reptSR")
head(ausisles.reptala.spp.unique.tab)


## amphibians
ausisles.amphala.df <- as.data.frame(ausisles.amphala)
head(ausisles.amphala.df)
table(ausisles.amphala.df$scientific)

## species list
amphspp <- table(ausisles.amphala.df$scientific)
amphspp

## remove higher-taxonomic level ambiguities & invasives
amphala.dat1 <- ausisles.amphala.df[!(ausisles.amphala.df$scientific %in% c("Rhinella marina",
                                                                "AMPHIBIA",
                                                                "ANURA",
                                                                "LIMNODYNASTIDAE",
                                                                "MICROHYLIDAE",
                                                                "MYOBATRACHIDAE",
                                                                "Rhinella")),]

amphala.dat2 <- amphala.dat1 %>%
  mutate(scientific = recode(scientific,
                                 "Limnodynastes dumerilii dumerilii" = "Limnodynastes dumerilii",
                                 "Limnodynastes dumerilii insularis" = "Limnodynastes dumerilii",
                                 "Limnodynastes dumerilii variegatus" = "Limnodynastes dumerilii",
                                 "Litoria verreauxii verreauxii" = "Litoria verreauxii",
                                 "Papurana" = "Papurana daemeli"))

amphala.dat3 <- amphala.dat2[!(amphala.dat2$scientific %in% c("Crinia",
                                                                  "Cophixalus",
                                                                  "Cyclorana",
                                                                  "Limnodynastes",
                                                                  "Litoria",
                                                                  "Litoria sp. cf. cooloolensis (North Stradbroke Is population)",
                                                                  "Neobatrachus",
                                                                  "Notaden",
                                                                  "Platyplectrum",
                                                                  "Pseudophryne",
                                                                  "Uperoleia")),]

amphspp <- table(amphala.dat3$scientific)
amphspp
amphspp.names <- attr(amphspp, "names")
amphspp.names
length(amphspp.names)
amphspp.names

ausisles.amphala.df <- as.data.frame(amphala.dat3)
head(ausisles.amphala.df)
ausisles.amphala.spp <- ausisles.amphala.df[,c("scientific", "FID")]
head(ausisles.amphala.spp)
ausisles.amphala.spp.unique <- unique(ausisles.amphala.spp)
head(ausisles.amphala.spp.unique)
ausisles.amphala.spp.unique.tab <- as.data.frame(table(ausisles.amphala.spp.unique$FID)) # number of unique species per island)
ausisles.amphala.spp.unique.tab
head(ausisles.amphala.spp.unique.tab)
colnames(ausisles.amphala.spp.unique.tab) <- c("FID", "amphSR")
head(ausisles.amphala.spp.unique.tab)




## merge with island data
mamausisl <- merge(ausisldatcorr, ausisles.mamala.spp.unique.tab, by="FID", all.x=TRUE)
head(mamausisl)
mamausisl$mamSR <- ifelse(is.na(mamausisl$mamSR), 0, mamausisl$mamSR) # replace NAs with 0
mamausislnoZero <- mamausisl[mamausisl$mamSR > 0,] # remove islands with no species
dim(mamausislnoZero)
mamausisl3spp <- mamausisl[mamausisl$mamSR > 2,] # remove islands with < 3 species
dim(mamausisl3spp)

median(mamausislnoZero$area_2/1e6, na.rm=T)
mad(mamausislnoZero$area_2/1e6, na.rm=T, constant=0.5)
sd(mamausislnoZero$area_2/1e6, na.rm=T)/sqrt(dim(mamausislnoZero)[1])

avausisl <- merge(ausisldatcorr, ausisles.avala.spp.unique.tab, by="FID", all.x=TRUE)
head(avausisl)
avausisl$avSR <- ifelse(is.na(avausisl$avSR), 0, avausisl$avSR) # replace NAs with 0
avausislnoZero <- avausisl[avausisl$avSR > 0,] # remove islands with no species
dim(avausislnoZero)
avausisl3spp <- avausisl[avausisl$avSR > 2,] # remove islands with < 3 species
dim(avausisl3spp)

median(avausislnoZero$area_2/1e6, na.rm=T)
mad(avausislnoZero$area_2/1e6, na.rm=T, constant=0.5)
sd(avausislnoZero$area_2/1e6, na.rm=T)/sqrt(dim(avausislnoZero)[1])

reptausisl <- merge(ausisldatcorr, ausisles.reptala.spp.unique.tab, by="FID", all.x=TRUE)
head(reptausisl)
reptausisl$reptSR <- ifelse(is.na(reptausisl$reptSR), 0, reptausisl$reptSR) # replace NAs with 0
reptausislnoZero <- reptausisl[reptausisl$reptSR > 0,] # remove islands with no species
dim(reptausislnoZero)
reptausisl3spp <- reptausisl[reptausisl$reptSR > 2,] # remove islands with < 3 species
dim(reptausisl3spp)

median(reptausislnoZero$area_2/1e6, na.rm=T)
mad(reptausislnoZero$area_2/1e6, na.rm=T, constant=0.5)
sd(reptausislnoZero$area_2/1e6, na.rm=T)/sqrt(dim(reptausislnoZero)[1])

amphausisl <- merge(ausisldatcorr, ausisles.amphala.spp.unique.tab, by="FID", all.x=TRUE)
head(amphausisl)
amphausisl$amphSR <- ifelse(is.na(amphausisl$amphSR), 0, amphausisl$amphSR) # replace NAs with 0
amphausislnoZero <- amphausisl[amphausisl$amphSR > 0,] # remove islands with no species
dim(amphausislnoZero)
amphausisl3spp <- amphausisl[amphausisl$amphSR > 2,] # remove islands with < 3 species
dim(amphausisl3spp)

median(amphausislnoZero$area_2/1e6, na.rm=T)
mad(amphausislnoZero$area_2/1e6, na.rm=T, constant=0.5)
sd(amphausislnoZero$area_2/1e6, na.rm=T)/sqrt(dim(amphausislnoZero)[1])


## plot species richness by area and distance with ggplot2
# Mammalia
mamausSARplot <- ggplot(mamausislnoZero, aes(x=area, y=mamSR)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("island area (m^2)") +
  ylab("species richness") +
  ggtitle("species richness vs island area") +
  theme_minimal()

mamausSDplot <- ggplot(mamausislnoZero, aes(x=distance, y=mamSR)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("distance to mainland (m)") +
  ylab("species richness") +
  ggtitle("species richness vs distance to mainland") +
  theme_minimal()

grid.arrange(mamausSARplot, mamausSDplot, ncol=2)

# z coefficients
mamausSARlm <- lm(log10(mamSR) ~ log10(area_2/1e6), data=mamausislnoZero)
summary(mamausSARlm)

mamausislnoZero$lmamSR <- log10(mamausislnoZero$mamSR)
mamausislnoZero$larea <- log10(mamausislnoZero$area_2/1e6)
mamausislSARlinfit <- lm(lmamSR ~ larea, data=mamausislnoZero)
mamausislSARlinfit.slope <- mamausislSARlinfit$coefficients[2]
mamausislSARlinfit.intercept <- mamausislSARlinfit$coefficients[1]
mamausislSARlinfit.slope # most frequent z = 0.2–0.4 https://www.pnas.org/doi/10.1073/pnas.0510605103
mamausislSARlinfit.slope.se <- summary(mamausislSARlinfit)$coefficients[4]
mamausislSARlinfit.slope.se
mamausislSARlinfit.intercept
10^mamausislSARlinfit.intercept
mamausislSARlinfit.intercept.se <- summary(mamausislSARlinfit)$coefficients[3]


# Aves
avausSARplot <- ggplot(avausislnoZero, aes(x=area_2, y=avSR)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("island area (m^2)") +
  ylab("species richness") +
  ggtitle("species richness vs island area") +
  theme_minimal()

avausSDplot <- ggplot(avausislnoZero, aes(x=distance, y=avSR)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("distance to mainland (m)") +
  ylab("species richness") +
  ggtitle("species richness vs distance to mainland") +
  theme_minimal()

grid.arrange(avausSARplot, avausSDplot, ncol=2)

# z coefficients
avausSARlm <- lm(log10(avSR) ~ log10(area_2/1e6), data=avausislnoZero)
summary(avausSARlm)

head(avausislnoZero)
avausislSARlinfit <- lm(log10(avSR) ~ log10(area_2/1e6), data=avausislnoZero)
summary(avausislSARlinfit)
avausislSARlinfit.slope <- avausislSARlinfit$coefficients[2]
avausislSARlinfit.intercept <- avausislSARlinfit$coefficients[1]
avausislSARlinfit.slope # most frequent z = 0.2–0.4 https://www.pnas.org/doi/10.1073/pnas.0510605103
avausislSARlinfit.slope.se <- summary(avausislSARlinfit)$coefficients[4]
avausislSARlinfit.slope.se
avausislSARlinfit.intercept
10^avausislSARlinfit.intercept
avausislSARlinfit.intercept.se <- summary(avausislSARlinfit)$coefficients[3]


# reptiles
reptausSARplot <- ggplot(reptausislnoZero, aes(x=area_2, y=reptSR)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("island area (m^2)") +
  ylab("species richness") +
  ggtitle("species richness vs island area") +
  theme_minimal()

reptausSDplot <- ggplot(reptausislnoZero, aes(x=distance, y=reptSR)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("distance to mainland (m)") +
  ylab("species richness") +
  ggtitle("species richness vs distance to mainland") +
  theme_minimal()

grid.arrange(reptausSARplot, reptausSDplot, ncol=2)

# z coefficients
reptausSARlm <- lm(log10(reptSR) ~ log10(area_2/1e6), data=reptausislnoZero)
summary(reptausSARlm)

head(reptausislnoZero)
# remove islands with < 3 spp
reptausislnoZero$larea <- log10(reptausislnoZero$area_2/1e6)
reptausislnoZero$ldistance <- log10(reptausislnoZero$distance/1000)
reptausislnoZero$lreptSR <- log10(reptausislnoZero$reptSR)
head(reptausislnoZero)

reptausislSARlinfit <- lm(lreptSR ~ larea, data=reptausislnoZero)
summary(reptausislSARlinfit)
reptausislSARlinfit.slope <- reptausislSARlinfit$coefficients[2]
reptausislSARlinfit.intercept <- reptausislSARlinfit$coefficients[1]
reptausislSARlinfit.slope # most frequent z = 0.2–0.4 https://www.pnas.org/doi/10.1073/pnas.0510605103
reptausislSARlinfit.slope.se <- summary(reptausislSARlinfit)$coefficients[4]
reptausislSARlinfit.slope.se
reptausislSARlinfit.intercept
10^reptausislSARlinfit.intercept
reptausislSARlinfit.intercept.se <- summary(reptausislSARlinfit)$coefficients[3]



# amphibians
amphausSARplot <- ggplot(amphausislnoZero, aes(x=area_2, y=amphSR)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("island area (m^2)") +
  ylab("species richness") +
  ggtitle("species richness vs island area") +
  theme_minimal()

amphausSDplot <- ggplot(amphausislnoZero, aes(x=distance, y=amphSR)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("distance to mainland (m)") +
  ylab("species richness") +
  ggtitle("species richness vs distance to mainland") +
  theme_minimal()

grid.arrange(amphausSARplot, amphausSDplot, ncol=2)

# z coefficients
amphausSARlm <- lm(log10(amphSR) ~ log10(area_2/1e6), data=amphausislnoZero)
summary(amphausSARlm)

head(amphausislnoZero)
amphausislnoZero$larea <- log10(amphausislnoZero$area_2/1e6)
amphausislnoZero$ldistance <- log10(amphausislnoZero$distance/1000)
amphausislnoZero$lamphSR <- log10(amphausislnoZero$amphSR)
head(amphausislnoZero)

amphausislSARlinfit <- lm(lamphSR ~ larea, data=amphausislnoZero)
summary(amphausislSARlinfit)
amphausislSARlinfit.slope <- amphausislSARlinfit$coefficients[2]
amphausislSARlinfit.intercept <- amphausislSARlinfit$coefficients[1]
amphausislSARlinfit.slope # most frequent z = 0.2–0.4 https://www.pnas.org/doi/10.1073/pnas.0510605103
amphausislSARlinfit.slope.se <- summary(amphausislSARlinfit)$coefficients[4]
amphausislSARlinfit.slope.se
amphausislSARlinfit.intercept
10^amphausislSARlinfit.intercept
amphausislSARlinfit.intercept.se <- summary(amphausislSARlinfit)$coefficients[3]


## boosted regression tree to test relative influence of area and distance to mainland
# Mammalia
head(mamausislnoZero)
mamausislnoZero$ldistance <- log10(mamausislnoZero$distance/1000)
mamausislnoZero$larea <- log10(mamausislnoZero$area_2/1e6)
head(mamausislnoZero)

which(colnames(mamausislnoZero)=="larea")
which(colnames(mamausislnoZero)=="ldistance")
which(colnames(mamausislnoZero)=="lmamSR")

ausmamIB.brt <- gbm.step(mamausislnoZero, gbm.x = attr(mamausislnoZero, "names")[c(12:13)],
                         gbm.y = attr(mamausislnoZero, "names")[11], family="gaussian", max.trees=100000,
                         tolerance = 0.0002, learning.rate = 0.0005, bag.fraction=0.75,
                         tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(ausmamIB.brt)
barplot(summary(ausmamIB.brt)$rel.inf, names.arg = summary(ausmamIB.brt)$var, xlab="relative influence", ylab="", col="blue")
ausmamIB.brt.summ <- summary(ausmamIB.brt)

ausmamIB.brt.plot <- ggplot(ausmamIB.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
ausmamIB.brt.plot.flip <- ausmamIB.brt.plot + coord_flip()
ausmamIB.brt.plot.flip

ausmamIB.brt.CV.cor <- 100 * ausmamIB.brt$cv.statistics$correlation.mean
ausmamIB.brt.CV.cor.se <- 100 * ausmamIB.brt$cv.statistics$correlation.se
print(c(ausmamIB.brt.CV.cor, ausmamIB.brt.CV.cor.se))

## plot partial dependence plots
gbm.plot(ausmamIB.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 mammal species richness", x.label="log10 area (km2)", plot.layout=c(1,1))
gbm.plot(ausmamIB.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T,
         y.label="log10 mammal species richness", x.label="log10 distance to mainland (km)", plot.layout=c(1,1))

## export partial dependency values
mam.brt.partial.deps <- list()
mam.pred.names <- ausmamIB.brt$var.names
eq.sp.pts <- 100

setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
for(i in seq_along(mam.pred.names)) {
  # Use gbm.plot to get partial dependency values
  # Set plot.layout = c(1,1) to prevent automatic plotting
  pd_data <- plot.gbm(ausmamIB.brt, i.var=i, continuous.resolution = eq.sp.pts, return.grid=T)
  mam.brt.partial.deps[[mam.pred.names[i]]] <- pd_data
  write.csv(mam.brt.partial.deps[[i]], file = paste0("mam_", mam.pred.names[i], "_partial_deps.csv"))
}

# Aves
head(avausislnoZero)
avausislnoZero$larea <- log10(avausislnoZero$area_2/1e6)
avausislnoZero$ldistance <- log10(avausislnoZero$distance/1000)
avausislnoZero$lavSR <- log10(avausislnoZero$avSR)
head(avausislnoZero)

which(colnames(avausislnoZero)=="larea")
which(colnames(avausislnoZero)=="ldistance")
which(colnames(avausislnoZero)=="lavSR")

ausavIB.brt <- gbm.step(avausislnoZero, gbm.x = attr(avausislnoZero, "names")[c(11:12)],
                        gbm.y = attr(avausislnoZero, "names")[13], family="gaussian", max.trees=100000,
                        tolerance = 0.0002, learning.rate = 0.0005, bag.fraction=0.75,
                        tree.complexity = 2, silent=F, tolerance.method = "auto")

summary(ausavIB.brt)
barplot(summary(ausavIB.brt)$rel.inf, names.arg = summary(ausavIB.brt)$var, xlab="relative influence", ylab="", col="blue")
ausavIB.brt.summ <- summary(ausavIB.brt)

ausavIB.brt.plot <- ggplot(ausavIB.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
ausavIB.brt.plot.flip <- ausavIB.brt.plot + coord_flip()
ausavIB.brt.plot.flip

ausavIB.brt.CV.cor <- 100 * ausavIB.brt$cv.statistics$correlation.mean
ausavIB.brt.CV.cor.se <- 100 * ausavIB.brt$cv.statistics$correlation.se
print(c(ausavIB.brt.CV.cor, ausavIB.brt.CV.cor.se))

## plot partial dependence plots
gbm.plot(ausavIB.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 bird species richness", x.label="log10 area (km2)", plot.layout=c(1,1))
gbm.plot(ausavIB.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T,
         y.label="log10 bird species richness", x.label="log10 distance to mainland (km)", plot.layout=c(1,1))

## export partial dependency values
av.brt.partial.deps <- list()
av.pred.names <- ausavIB.brt$var.names
eq.sp.pts <- 100

setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
for(i in seq_along(av.pred.names)) {
  # Use gbm.plot to get partial dependency values
  # Set plot.layout = c(1,1) to prevent automatic plotting
  pd_data <- plot.gbm(ausavIB.brt, i.var=i, continuous.resolution = eq.sp.pts, return.grid=T)
  av.brt.partial.deps[[av.pred.names[i]]] <- pd_data
  write.csv(av.brt.partial.deps[[i]], file = paste0("av_", av.pred.names[i], "_partial_deps.csv"))
}


# reptiles
head(reptausislnoZero)
reptausislnoZero$larea <- log10(reptausislnoZero$area_2/1e6)
reptausislnoZero$ldistance <- log10(reptausislnoZero$distance/1000)
reptausislnoZero$lreptSR <- log10(reptausislnoZero$reptSR)
head(reptausislnoZero)

which(colnames(reptausislnoZero)=="larea")
which(colnames(reptausislnoZero)=="ldistance")
which(colnames(reptausislnoZero)=="lreptSR")

ausreptIB.brt <- gbm.step(reptausislnoZero, gbm.x = attr(reptausislnoZero, "names")[c(11:12)],
                        gbm.y = attr(reptausislnoZero, "names")[13], family="gaussian", max.trees=100000,
                        tolerance = 0.0002, learning.rate = 0.0005, bag.fraction=0.75,
                        tree.complexity = 2, silent=F, tolerance.method = "auto")

summary(ausreptIB.brt)
barplot(summary(ausreptIB.brt)$rel.inf, names.arg = summary(ausreptIB.brt)$var, xlab="relative influence", ylab="", col="blue")
ausreptIB.brt.summ <- summary(ausreptIB.brt)

ausreptIB.brt.plot <- ggplot(ausreptIB.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
ausreptIB.brt.plot.flip <- ausreptIB.brt.plot + coord_flip()
ausreptIB.brt.plot.flip

ausreptIB.brt.CV.cor <- 100 * ausreptIB.brt$cv.statistics$correlation.mean
ausreptIB.brt.CV.cor.se <- 100 * ausreptIB.brt$cv.statistics$correlation.se
print(c(ausreptIB.brt.CV.cor, ausreptIB.brt.CV.cor.se))

## plot partial dependence plots
gbm.plot(ausreptIB.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 reptile species richness", x.label="log10 area (km2)", plot.layout=c(1,1))
gbm.plot(ausreptIB.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T,
         y.label="log10 reptile species richness", x.label="log10 distance to mainland (km)", plot.layout=c(1,1))

## export partial dependency values
rept.brt.partial.deps <- list()
rept.pred.names <- ausreptIB.brt$var.names
eq.sp.pts <- 100

setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
for(i in seq_along(rept.pred.names)) {
  # Use gbm.plot to get partial dependency values
  # Set plot.layout = c(1,1) to prevent automatic plotting
  pd_data <- plot.gbm(ausreptIB.brt, i.var=i, continuous.resolution = eq.sp.pts, return.grid=T)
  rept.brt.partial.deps[[rept.pred.names[i]]] <- pd_data
  write.csv(rept.brt.partial.deps[[i]], file = paste0("rept_", rept.pred.names[i], "_partial_deps.csv"))
}


# amphibians
head(amphausislnoZero)
amphausislnoZero$larea <- log10(amphausislnoZero$area_2/1e6)
amphausislnoZero$ldistance <- log10(amphausislnoZero$distance/1000)
amphausislnoZero$lamphSR <- log10(amphausislnoZero$amphSR)
head(amphausislnoZero)

which(colnames(amphausislnoZero)=="larea")
which(colnames(amphausislnoZero)=="ldistance")
which(colnames(amphausislnoZero)=="lamphSR")

ausamphIB.brt <- gbm.step(amphausislnoZero, gbm.x = attr(amphausislnoZero, "names")[c(11:12)],
                          gbm.y = attr(amphausislnoZero, "names")[13], family="gaussian", max.trees=100000,
                          tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.75,
                          tree.complexity = 2, silent=F, tolerance.method = "auto")

summary(ausamphIB.brt)
barplot(summary(ausamphIB.brt)$rel.inf, names.arg = summary(ausamphIB.brt)$var, xlab="relative influence", ylab="", col="blue")
ausamphIB.brt.summ <- summary(ausamphIB.brt)

ausamphIB.brt.plot <- ggplot(ausamphIB.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
ausamphIB.brt.plot.flip <- ausamphIB.brt.plot + coord_flip()
ausamphIB.brt.plot.flip

ausamphIB.brt.CV.cor <- 100 * ausamphIB.brt$cv.statistics$correlation.mean
ausamphIB.brt.CV.cor.se <- 100 * ausamphIB.brt$cv.statistics$correlation.se
print(c(ausamphIB.brt.CV.cor, ausamphIB.brt.CV.cor.se))

## plot partial dependence plots
gbm.plot(ausamphIB.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 amphibian species richness", x.label="log10 area (km2)", plot.layout=c(1,1))
gbm.plot(ausamphIB.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T,
         y.label="log10 amphibian species richness", x.label="log10 distance to mainland (km)", plot.layout=c(1,1))

## export partial dependency values
amph.brt.partial.deps <- list()
amph.pred.names <- ausamphIB.brt$var.names
eq.sp.pts <- 100

setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
for(i in seq_along(amph.pred.names)) {
  # Use gbm.plot to get partial dependency values
  # Set plot.layout = c(1,1) to prevent automatic plotting
  pd_data <- plot.gbm(ausamphIB.brt, i.var=i, continuous.resolution = eq.sp.pts, return.grid=T)
  amph.brt.partial.deps[[amph.pred.names[i]]] <- pd_data
  write.csv(amph.brt.partial.deps[[i]], file = paste0("amph_", amph.pred.names[i], "_partial_deps.csv"))
}



## combine taxonomic datasets
head(mamausisl)
ausislala.mrg1 <- merge(mamausisl, avausisl, by="FID", all.x=TRUE, all.y=TRUE)
head(ausislala.mrg1)
ausislala.mrg2 <- merge(ausislala.mrg1, reptausisl, by="FID", all.x=TRUE, all.y=TRUE)
head(ausislala.mrg2)
ausislala.mrg3 <- merge(ausislala.mrg2, amphausisl, by="FID", all.x=TRUE, all.y=TRUE)
head(ausislala.mrg3)

# remove duplicate columns
ausislala.mrg3 <- ausislala.mrg3[, !duplicated(as.list(ausislala.mrg3))]
head(ausislala.mrg3)



# combine species richness columns
ausislala.mrg3$totSR <- rowSums(ausislala.mrg3[,c("mamSR", "avSR", "reptSR", "amphSR")], na.rm=TRUE)
ausislala.mrg3$totSR <- ifelse(ausislala.mrg3$totSR > 0, ausislala.mrg3$totSR, NA) # replace 0 with NA
ausislalaNoZero <- ausislala.mrg3[ausislala.mrg3$totSR > 0,] # remove islands with no species
head(ausislalaNoZero)
median(ausislalaNoZero$area_2.x/1e6, na.rm=T)
range(ausislalaNoZero$area_2.x/1e6, na.rm=T)
quantile(ausislalaNoZero$area_2.x/1e6, probs=0.875, na.rm=T)
quantile(ausislalaNoZero$area_2.x/1e6, probs=0.125, na.rm=T)

median(ausislalaNoZero$distance/1000, na.rm=T)
range(ausislalaNoZero$distance/1000, na.rm=T)
quantile(ausislalaNoZero$distance/1000, probs=0.875, na.rm=T)
quantile(ausislalaNoZero$distance/1000, probs=0.125, na.rm=T)


## remove duplicate columns
ausislalaNoZero.NoDups <- ausislalaNoZero[, !duplicated(as.list(ausislalaNoZero))]
head(ausislalaNoZero.NoDups)

# export .csv
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(ausislalaNoZero.NoDups, "AusIslalaNoZero.NoDups.csv", row.names=FALSE)

# plot with ggplot2
ausSARplot <- ggplot(ausislalaNoZero.NoDups, aes(x=area_2.x/1e6, y=totSR)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("island area (km^2)") +
  ylab("total species richness") +
  ggtitle("species richness vs island area") +
  theme_minimal()

ausSDplot <- ggplot(ausislalaNoZero.NoDups, aes(x=distance.x/1000, y=totSR)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("distance to mainland (km)") +
  ylab("total species richness") +
  ggtitle("species richness vs distance to mainland") +
  theme_minimal()

grid.arrange(ausSARplot, ausSDplot, ncol=2)

# boosted regression tree
head(ausislalaNoZero.NoDups)
ausislalaNoZero.NoDups$larea <- log10(ausislalaNoZero.NoDups$area_2.x/1e6)
ausislalaNoZero.NoDups$ldistance <- log10(ausislalaNoZero.NoDups$distance.x/1000)
ausislalaNoZero.NoDups$ltotSR <- log10(ausislalaNoZero.NoDups$totSR)
head(ausislalaNoZero.NoDups)
colnames(ausislalaNoZero.NoDups)
ausislalaNoNA <- ausislalaNoZero.NoDups[!is.na(ausislalaNoZero.NoDups$totSR),] # remove NAs
ausisl3spp <- ausislalaNoZero.NoDups[ausislalaNoZero.NoDups$totSR > 2,] # remove islands with < 3 species
dim(ausisl3spp)
head(ausisl3spp)

# z coefficients
vertausSARlm <- lm(log10(totSR) ~ log10(area_2.x/1e6), data=ausislalaNoNA)
summary(vertausSARlm)

head(ausislalaNoNA)
vertausislSARlinfit <- lm(log10(totSR) ~ log10(area_2.x/1e6), data=ausislalaNoNA)
summary(vertausislSARlinfit)
vertausislSARlinfit.slope <- vertausislSARlinfit$coefficients[2]
vertausislSARlinfit.intercept <- vertausislSARlinfit$coefficients[1]
vertausislSARlinfit.slope # most frequent z = 0.2–0.4 https://www.pnas.org/doi/10.1073/pnas.0510605103
vertausislSARlinfit.slope.se <- summary(vertausislSARlinfit)$coefficients[4]
vertausislSARlinfit.slope.se
vertausislSARlinfit.intercept
10^vertausislSARlinfit.intercept
vertausislSARlinfit.intercept.se <- summary(vertausislSARlinfit)$coefficients[3]


which(colnames(ausislalaNoNA)=="larea")
which(colnames(ausislalaNoNA)=="ldistance")
which(colnames(ausislalaNoNA)=="ltotSR")

ausIBtot.brt <- gbm.step(ausislalaNoNA, gbm.x = attr(ausislalaNoNA, "names")[c(15:16)],
                         gbm.y = attr(ausislalaNoNA, "names")[17], family="gaussian", max.trees=100000,
                         tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.75,
                         tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(ausIBtot.brt)
barplot(summary(ausIBtot.brt)$rel.inf, names.arg = summary(ausIBtot.brt)$var, xlab="relative influence", ylab="", col="blue")
ausIBtot.brt.summ <- summary(ausIBtot.brt)
ausIBtot.brt.plot <- ggplot(ausIBtot.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
ausIBtot.brt.plot.flip <- ausIBtot.brt.plot + coord_flip()
ausIBtot.brt.plot.flip

ausIBtot.brt.CV.cor <- 100 * ausIBtot.brt$cv.statistics$correlation.mean
ausIBtot.brt.CV.cor.se <- 100 * ausIBtot.brt$cv.statistics$correlation.se
print(c(ausIBtot.brt.CV.cor, ausIBtot.brt.CV.cor.se))

## plot partial dependence plots
gbm.plot(ausIBtot.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 total species richness", x.label="log10 area (m2)", plot.layout=c(1,1))
gbm.plot(ausIBtot.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T,
         y.label="log10 total species richness", x.label="log10 distance to mainland (m)", plot.layout=c(1,1))



## generate fits using slope and intercept values + SEs
# mammals
lareakm2.vec <- as.numeric(na.omit(log10(ausislalaNoZero.NoDups$area_2.x/1e6)))

iter <- 10000
predMamlSR.mat <- matrix(data=NA, nrow=iter, ncol=length(lareakm2.vec)) # storage matrix
for (i in 1:iter) {
  # resample slope
  slope.rsmp <- rnorm(n=length(lareakm2.vec), mean=as.numeric(mamausislSARlinfit.slope), 
                      sd=as.numeric(mamausislSARlinfit.slope.se))
  # resample intercept
  int.rsmp <- rnorm(n=length(lareakm2.vec), mean=as.numeric(mamausislSARlinfit.intercept), 
                    sd=as.numeric(mamausislSARlinfit.intercept.se)) 
  # predict log SR
  predMamlSR.mat[i,] <- int.rsmp + lareakm2.vec*slope.rsmp
}
predMamlSR.mn <- apply(predMamlSR.mat, MARGIN=2, median, na.rm=T)
predMamlSR.up <- apply(predMamlSR.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
predMamlSR.lo <- apply(predMamlSR.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)

predMamSR.dat <- data.frame(areakm2=10^lareakm2.vec,
                              predMamSR.mn=10^predMamlSR.mn,
                              predMamSR.up=10^predMamlSR.up,
                              predMamSR.lo=10^predMamlSR.lo)
predMamSR.sort <- predMamSR.dat[order(predMamSR.dat$areakm2, decreasing=T),]
head(predMamSR.sort)

# export
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(predMamSR.sort, "predMamSRsort.csv", row.names=FALSE)


# birds
predAvlSR.mat <- matrix(data=NA, nrow=iter, ncol=length(lareakm2.vec)) # storage matrix
for (i in 1:iter) {
  # resample slope
  slope.rsmp <- rnorm(n=length(lareakm2.vec), mean=as.numeric(avausislSARlinfit.slope), 
                      sd=as.numeric(avausislSARlinfit.slope.se))
  # resample intercept
  int.rsmp <- rnorm(n=length(lareakm2.vec), mean=as.numeric(avausislSARlinfit.intercept), 
                    sd=as.numeric(avausislSARlinfit.intercept.se)) 
  # predict log SR
  predAvlSR.mat[i,] <- int.rsmp + lareakm2.vec*slope.rsmp
}
predAvlSR.mn <- apply(predAvlSR.mat, MARGIN=2, median, na.rm=T)
predAvlSR.up <- apply(predAvlSR.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
predAvlSR.lo <- apply(predAvlSR.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)

predAvSR.dat <- data.frame(areakm2=10^lareakm2.vec,
                            predAvSR.mn=10^predAvlSR.mn,
                            predAvSR.up=10^predAvlSR.up,
                            predAvSR.lo=10^predAvlSR.lo)
predAvSR.sort <- predAvSR.dat[order(predAvSR.dat$areakm2, decreasing=T),]
head(predAvSR.sort)

# export
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(predAvSR.sort, "predAvSRsort.csv", row.names=FALSE)

# reptiles
predReptlSR.mat <- matrix(data=NA, nrow=iter, ncol=length(lareakm2.vec)) # storage matrix
for (i in 1:iter) {
  # resample slope
  slope.rsmp <- rnorm(n=length(lareakm2.vec), mean=as.numeric(reptausislSARlinfit.slope), 
                      sd=as.numeric(reptausislSARlinfit.slope.se))
  # resample intercept
  int.rsmp <- rnorm(n=length(lareakm2.vec), mean=as.numeric(reptausislSARlinfit.intercept), 
                    sd=as.numeric(reptausislSARlinfit.intercept.se)) 
  # predict log SR
  predReptlSR.mat[i,] <- int.rsmp + lareakm2.vec*slope.rsmp
}
predReptlSR.mn <- apply(predReptlSR.mat, MARGIN=2, median, na.rm=T)
predReptlSR.up <- apply(predReptlSR.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
predReptlSR.lo <- apply(predReptlSR.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)

predReptSR.dat <- data.frame(areakm2=10^lareakm2.vec,
                           predReptSR.mn=10^predReptlSR.mn,
                           predReptSR.up=10^predReptlSR.up,
                           predReptSR.lo=10^predReptlSR.lo)
predReptSR.sort <- predReptSR.dat[order(predReptSR.dat$areakm2, decreasing=T),]
head(predReptSR.sort)

# export
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(predReptSR.sort, "predReptSRsort.csv", row.names=FALSE)

# amphibians
predAmphlSR.mat <- matrix(data=NA, nrow=iter, ncol=length(lareakm2.vec)) # storage matrix
for (i in 1:iter) {
  # resample slope
  slope.rsmp <- rnorm(n=length(lareakm2.vec), mean=as.numeric(amphausislSARlinfit.slope), 
                      sd=as.numeric(amphausislSARlinfit.slope.se))
  # resample intercept
  int.rsmp <- rnorm(n=length(lareakm2.vec), mean=as.numeric(amphausislSARlinfit.intercept), 
                    sd=as.numeric(amphausislSARlinfit.intercept.se)) 
  # predict log SR
  predAmphlSR.mat[i,] <- int.rsmp + lareakm2.vec*slope.rsmp
}
predAmphlSR.mn <- apply(predAmphlSR.mat, MARGIN=2, median, na.rm=T)
predAmphlSR.up <- apply(predAmphlSR.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
predAmphlSR.lo <- apply(predAmphlSR.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)

predAmphSR.dat <- data.frame(areakm2=10^lareakm2.vec,
                             predAmphSR.mn=10^predAmphlSR.mn,
                             predAmphSR.up=10^predAmphlSR.up,
                             predAmphSR.lo=10^predAmphlSR.lo)
predAmphSR.sort <- predAmphSR.dat[order(predAmphSR.dat$areakm2, decreasing=T),]
head(predAmphSR.sort)

# export
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(predAmphSR.sort, "predAmphSRsort.csv", row.names=FALSE)


# reptiles
predReptlSR.mat <- matrix(data=NA, nrow=iter, ncol=length(lareakm2.vec)) # storage matrix
for (i in 1:iter) {
  # resample slope
  slope.rsmp <- rnorm(n=length(lareakm2.vec), mean=as.numeric(reptausislSARlinfit.slope), 
                      sd=as.numeric(reptausislSARlinfit.slope.se))
  # resample intercept
  int.rsmp <- rnorm(n=length(lareakm2.vec), mean=as.numeric(reptausislSARlinfit.intercept), 
                    sd=as.numeric(reptausislSARlinfit.intercept.se)) 
  # predict log SR
  predReptlSR.mat[i,] <- int.rsmp + lareakm2.vec*slope.rsmp
}
predReptlSR.mn <- apply(predReptlSR.mat, MARGIN=2, median, na.rm=T)
predReptlSR.up <- apply(predReptlSR.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
predReptlSR.lo <- apply(predReptlSR.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)

predReptSR.dat <- data.frame(areakm2=10^lareakm2.vec,
                             predReptSR.mn=10^predReptlSR.mn,
                             predReptSR.up=10^predReptlSR.up,
                             predReptSR.lo=10^predReptlSR.lo)
predReptSR.sort <- predReptSR.dat[order(predReptSR.dat$areakm2, decreasing=T),]
head(predReptSR.sort)

# export
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(predReptSR.sort, "predReptSRsort.csv", row.names=FALSE)

# vertebrates (mammals + birds + amphibians + reptiles)
predVertlSR.mat <- matrix(data=NA, nrow=iter, ncol=length(lareakm2.vec)) # storage matrix
for (i in 1:iter) {
  # resample slope
  slope.rsmp <- rnorm(n=length(lareakm2.vec), mean=as.numeric(vertausislSARlinfit.slope), 
                      sd=as.numeric(vertausislSARlinfit.slope.se))
  # resample intercept
  int.rsmp <- rnorm(n=length(lareakm2.vec), mean=as.numeric(vertausislSARlinfit.intercept), 
                    sd=as.numeric(vertausislSARlinfit.intercept.se)) 
  # predict log SR
  predVertlSR.mat[i,] <- int.rsmp + lareakm2.vec*slope.rsmp
}
predVertlSR.mn <- apply(predVertlSR.mat, MARGIN=2, median, na.rm=T)
predVertlSR.up <- apply(predVertlSR.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
predVertlSR.lo <- apply(predVertlSR.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)

predVertSR.dat <- data.frame(areakm2=10^lareakm2.vec,
                             predVertSR.mn=10^predVertlSR.mn,
                             predVertSR.up=10^predVertlSR.up,
                             predVertSR.lo=10^predVertlSR.lo)
predVertSR.sort <- predVertSR.dat[order(predVertSR.dat$areakm2, decreasing=T),]
head(predVertSR.sort)

# export
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(predVertSR.sort, "predVertSRsort.csv", row.names=FALSE)


## determine variation in z with incrementally removing islands with the fewest species
head(ausislalaNoZero.NoDups)

## mammals
mamdat <- ausislalaNoZero.NoDups[,c("FID", "mamSR", "larea")]
mamdat <- subset(mamdat, mamSR > 0)
head(mamdat)

# max number of total species per island to remove
rem.max.spp <- 10
z.val <- z.se <- rep(NA, rem.max.spp)
dat.sim <- mamdat
colnames(dat.sim)[2] <- "SR"
for (r in 1:rem.max.spp) {
 dat.incr <- subset(dat.sim, SR >= r)
 fit.larea <- lm(log10(SR) ~ larea, data=dat.incr)
 z.val[r] <- fit.larea$coefficients[2]
 z.se[r] <- summary(fit.larea)$coefficients[4]
}
rem.spp.vec <- 0:(rem.max.spp-1)

## plot z and its se in ggplot2
mam.z.dat <- data.frame(rem.spp=rem.spp.vec, z=z.val, z.se=z.se)
mam.z.plot <- ggplot(mam.z.dat, aes(x=rem.spp, y=z)) +
  geom_point(col="blue") +
  geom_line(col="blue") +
  geom_errorbar(aes(ymin=z-z.se, ymax=z+z.se), width=0.2, col="blue") +
  xlab("number of minimum species removed") +
  ylab("z") +
  theme_minimal()
mam.z.plot


## birds
avdat <- ausislalaNoZero.NoDups[,c("FID", "avSR", "larea")]
avdat <- subset(avdat, avSR > 0)
head(avdat)

# max number of total species per island to remove
rem.max.spp <- 10
z.val <- z.se <- rep(NA, rem.max.spp)
dat.sim <- avdat
colnames(dat.sim)[2] <- "SR"
for (r in 1:rem.max.spp) {
  dat.incr <- subset(dat.sim, SR >= r)
  fit.larea <- lm(log10(SR) ~ larea, data=dat.incr)
  z.val[r] <- fit.larea$coefficients[2]
  z.se[r] <- summary(fit.larea)$coefficients[4]
}
rem.spp.vec <- 0:(rem.max.spp-1)

## plot z and its se in ggplot2
av.z.dat <- data.frame(rem.spp=rem.spp.vec, z=z.val, z.se=z.se)
av.z.plot <- ggplot(av.z.dat, aes(x=rem.spp, y=z)) +
  geom_point(col="red") +
  geom_line(col="red") +
  geom_errorbar(aes(ymin=z-z.se, ymax=z+z.se), width=0.2, col="red") +
  xlab("number of minimum species removed") +
  ylab("z") +
  theme_minimal()
av.z.plot


## reptiles
reptdat <- ausislalaNoZero.NoDups[,c("FID", "reptSR", "larea")]
reptdat <- subset(reptdat, reptSR > 0)
head(reptdat)

# max number of total species per island to remove
rem.max.spp <- 10
z.val <- z.se <- rep(NA, rem.max.spp)
dat.sim <- reptdat
colnames(dat.sim)[2] <- "SR"
for (r in 1:rem.max.spp) {
  dat.incr <- subset(dat.sim, SR >= r)
  fit.larea <- lm(log10(SR) ~ larea, data=dat.incr)
  z.val[r] <- fit.larea$coefficients[2]
  z.se[r] <- summary(fit.larea)$coefficients[4]
}
rem.spp.vec <- 0:(rem.max.spp-1)

## plot z and its se in ggplot2
rept.z.dat <- data.frame(rem.spp=rem.spp.vec, z=z.val, z.se=z.se)
rept.z.plot <- ggplot(rept.z.dat, aes(x=rem.spp, y=z)) +
  geom_point(col="green") +
  geom_line(col="green") +
  geom_errorbar(aes(ymin=z-z.se, ymax=z+z.se), width=0.2, col="green") +
  xlab("number of minimum species removed") +
  ylab("z") +
  theme_minimal()
rept.z.plot


## amphibians
amphdat <- ausislalaNoZero.NoDups[,c("FID", "amphSR", "larea")]
amphdat <- subset(amphdat, amphSR > 0)
head(amphdat)

# max number of total species per island to remove
rem.max.spp <- 10
z.val <- z.se <- rep(NA, rem.max.spp)
dat.sim <- amphdat
colnames(dat.sim)[2] <- "SR"
for (r in 1:rem.max.spp) {
  dat.incr <- subset(dat.sim, SR >= r)
  fit.larea <- lm(log10(SR) ~ larea, data=dat.incr)
  z.val[r] <- fit.larea$coefficients[2]
  z.se[r] <- summary(fit.larea)$coefficients[4]
}
rem.spp.vec <- 0:(rem.max.spp-1)

## plot z and its se in ggplot2
amph.z.dat <- data.frame(rem.spp=rem.spp.vec, z=z.val, z.se=z.se)
amph.z.plot <- ggplot(amph.z.dat, aes(x=rem.spp, y=z)) +
  geom_point(col="brown") +
  geom_line(col="brown") +
  geom_errorbar(aes(ymin=z-z.se, ymax=z+z.se), width=0.2, col="brown") +
  xlab("number of minimum species removed") +
  ylab("z") +
  theme_minimal()
amph.z.plot

# combine all z plots
z.plots <- grid.arrange(mam.z.plot, av.z.plot, rept.z.plot, amph.z.plot, ncol=2)

vert.z.dat <- data.frame(rem.spp=rem.spp.vec, 
                      mam.z=mam.z.dat$z, mam.z.se=mam.z.dat$z.se,
                      av.z=av.z.dat$z, av.z.se=av.z.dat$z.se,
                      rept.z=rept.z.dat$z, rept.z.se=rept.z.dat$z.se,
                      amph.z=amph.z.dat$z, amph.z.se=amph.z.dat$z.se)
# export
head(vert.z.dat)
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(vert.z.dat, "vertZdat.csv", row.names=FALSE)



###########################
###########################
## trait analyses
###########################
###########################

## SahulTraits
## mammals

# import trait data
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/data/SahulTraits")
mamtraits <- read.csv("SahulTraitsMam.csv", header=TRUE)
head(mamtraits)

# create concatenated 'scientificName' column to match ALA data
# there is whitespace after some names. Trim it first
mamtraits$Genus <- trimws(mamtraits$Genus)
mamtraits$Species <- trimws(mamtraits$Species)
mamtraits$scientific <- paste(mamtraits$Genus, mamtraits$Species, sep=" ")
head(mamtraits)

# fix one name to match island data
mamtraits$scientific[mamtraits$scientific == "Rattus lutreola"] <- "Rattus lutreolus"
mamtraits$scientific

# merge with island ALA dataset
head(ausisles.mamala.df6)

## unique species only
ausisles.mamala.df.unique <- unique(ausisles.mamala.df6[,c("FID", "scientific")])
dim(ausisles.mamala.df.unique)
ausisles.mamala.traits <- merge(ausisles.mamala.df.unique, mamtraits, by="scientific", all.x=F)
head(ausisles.mamala.traits)
dim(ausisles.mamala.traits)
table(ausisles.mamala.traits$FID)

# which don't have any matches?
unique(ausisles.mamala.traits$scientific[is.na(ausisles.mamala.traits$Body.mass_Mean..g.)]) # all good now

# make a list of unique species and traits
sp_traits <- ausisles.mamala.traits
sp_traits <- sp_traits[,!names(sp_traits) %in% c("FID","X","Family","Genus","Species","Introduced..y.n.")]
sp_traits <- unique(sp_traits)

# check distributions of the traits being used - transform if skewed
hist(sp_traits$Body.mass_Mean..g.) # needs to be log transformed (skewed) and renamed logBodyMass (or something like that)
sp_traits$Log_Body_mass <- log(sp_traits$Body.mass_Mean..g.)
hist(sp_traits$Log_Body_mass)
sp_traits$Body.mass_Mean..g. <- NULL

hist(sp_traits$Total.body.Length_Mean..mm.) # needs to be log transformed (skewed) and renamed log...
sp_traits$Log_TL <- log(sp_traits$Total.body.Length_Mean..mm.)
hist(sp_traits$Log_TL)
sp_traits$Total.body.Length_Mean..mm. <- NULL

hist(sp_traits$SVL_Mean..mm.)
sp_traits$Log_SVL <- log(sp_traits$SVL_Mean..mm.)
hist(sp_traits$Log_SVL)
sp_traits$SVL_Mean..mm. <- NULL

hist(sp_traits$Longevity..days.)
sp_traits$Log_longevity <- log((sp_traits$Longevity..days.))
hist(sp_traits$Log_longevity)
sp_traits$Longevity..days. <- NULL

hist(sp_traits$Gestation.length..days.)
sp_traits$Log_Gestation_length <- log(sp_traits$Gestation.length..days.)
hist(sp_traits$Log_Gestation_length)
sp_traits$Gestation.length..days. <- NULL

hist(sp_traits$Litter.clutch.size)
hist(log(sp_traits$Litter.clutch.size)) #still skewed but a slight improvement
sp_traits$Log_clutch_size <- log(sp_traits$Litter.clutch.size)
sp_traits$Litter.clutch.size <- NULL

hist(sp_traits$Litters.clutches.per.year)
hist(log(sp_traits$Litters.clutches.per.year))
sp_traits$Log_litters_yr <- log(sp_traits$Litters.clutches.per.year)
sp_traits$Litters.clutches.per.year <- NULL

hist(sp_traits$Weaning.age..days..)
hist(log(sp_traits$Weaning.age..days..))
sp_traits$Log_weaning_age <- log(sp_traits$Weaning.age..days..)
sp_traits$Weaning.age..days.. <- NULL

hist(sp_traits$trophic_level)

hist(sp_traits$Sexual.maturity..days..)

hist(sp_traits$brain_mass_g)
hist(log(sp_traits$brain_mass_g))
sp_traits$Log_brain_mass <- log(sp_traits$brain_mass_g)
sp_traits$brain_mass_g <- NULL

dim(sp_traits)

sp_traits$scientific

# check for NAs
na_counts <- sapply(sp_traits, function(x) sum(is.na(x)))
na_counts

## make list of mammal communities by island (FID)
ausisles.mamala.traits.list <- split(ausisles.mamala.traits, ausisles.mamala.traits$FID)
View(ausisles.mamala.traits.list[[which.max(sapply(ausisles.mamala.traits.list, nrow))]]) # view the island with the most records

# make a presence absence matrix with row = FID and colname = species
combined_unique <- unique(ausisles.mamala.traits[, c("FID", "scientific")])
island_pres <- table(combined_unique$FID, combined_unique$scientific)

# convert table to a data frame with species as columns
island_pres_df <- as.data.frame.matrix(island_pres)

# add FID as its own column (from row names)
island_pres_df$FID <- rownames(island_pres_df)

# move FID to the first column
island_pres_df <- island_pres_df[, c("FID", setdiff(names(island_pres_df), "FID"))]

# reset row names
rownames(island_pres_df) <- NULL

# trait type data frame
tt <- data.frame(trait_name = names(sp_traits), trait_type=NA)

# fix the trait_type column
tt$trait_type[tt$trait_name %in% c("Sociality")] <- "N"
tt$trait_type[tt$trait_name %in% c("Plant..0.1.","Terrestrial.vertebrate..0.1.","Fish..0.1.",
                                   "Invertebrate..0.1.","Carrion..0.1.",
                                 "DigitsWebbed","DigitsClaws","DigitsFlippers","DigitsHooves",
                                 "LocoAquatic","LocoArboreal","LocoFossorial","LocoGliding",
                                 "LocoCursorial","LocoSaltation","LocoTerrestrial",
                                 "day_active","night_active","crepuscular_active",
                                 "niche_aerial","niche_arboreal","niche_ground","niche_marine",
                                 "niche_scansorial")] <- "F"
tt$trait_type[tt$trait_name %in% c("Total.body.Length_Mean..mm.","Sexual.maturity..days..","trophic_level",
                                   "Log_Body_mass","Log_TL","Log_SVL","Log_longevity","Log_Gestation_length",
                                 "Log_clutch_size","Log_litters_yr","Log_weaning_age","Log_brain_mass")] <- "Q"
# add the fuzzy column
tt$fuzzy_name[tt$trait_name %in% c("Plant..0.1.","Terrestrial.vertebrate..0.1.","Fish..0.1.",
                                   "Invertebrate..0.1.","Carrion..0.1.")] <- "diet"
tt$fuzzy_name[tt$trait_name %in% c("DigitsWebbed","DigitsClaws","DigitsFlippers","DigitsHooves")] <- "digits"
tt$fuzzy_name[tt$trait_name %in% c("LocoAquatic","LocoArboreal","LocoFossorial","LocoGliding",
                                   "LocoCursorial","LocoSaltation","LocoTerrestrial")] <- "locomotion"
tt$fuzzy_name[tt$trait_name %in% c("day_active","night_active","crepuscular_active")] <- "activity_time"
tt$fuzzy_name[tt$trait_name %in% c("niche_aerial","niche_arboreal","niche_ground","niche_marine",
                                   "niche_scansorial")] <- "niche"
#tt$fuzzy_name <- ifelse(is.na(tt$fuzzy_name), tt$trait_name, tt$fuzzy_name)

# remove non-trait rows  
tt <- tt[!tt$trait_name %in% c("scientific", "Limbs"),]

# species as row name in sp_traits
row.names(sp_traits) <- sp_traits$scientific
sp_traits$scientific <- NULL

# remove limbs column - they all have four
sp_traits$Limbs <- NULL

# format data in sp_traits. mFD required particular formatting
# loop over each trait defined in tt
for (i in seq_len(nrow(tt))) {
  trait <- tt$trait_name[i]
  type <- tt$trait_type[i]
  
  if (trait %in% colnames(sp_traits)) {
    if (type == "N") {
      sp_traits[[trait]] <- as.factor(sp_traits[[trait]])
    } else if (type == "Q" || type == "F") {
      sp_traits[[trait]] <- as.numeric(sp_traits[[trait]])
    } else {
      warning(paste("unknown trait type for", trait))
    }
  } else {
    warning(paste("Trait", trait, "not found in sp_traits"))
  }
}

# species traits summary:
traits_summ <- mFD::sp.tr.summary(
  tr_cat     = tt,   
  sp_tr      = sp_traits, 
  stop_if_NA = TRUE)

#traits_summ$"tr_types" 
#traits_summ$"mod_list"

# Summary of the assemblages * species dataframe:
island_pres_mat <- island_pres_df
row.names(island_pres_mat) <- island_pres_mat$FID
island_pres_mat$FID <- NULL
island_pres_mat[] <- lapply(island_pres_mat, as.numeric)
island_pres_mat <- as.matrix(island_pres_mat)

isl_summ <- mFD::asb.sp.summary(asb_sp_w = island_pres_mat)
isl_summ$"asb_sp_richn" # sp richness by island (non-zero only)

# remove variables with no variation
sp_traits$DigitsFlippers <- NULL
sp_traits$DigitsHooves <- NULL  # can change this if invasives are added
sp_traits$LocoGliding <- NULL
sp_traits$LocoSaltation <- NULL
sp_traits$niche_marine <- NULL
sp_traits$niche_scansorial <- NULL

# no cut trait type df back to traits in sp_traits
tt <- tt[tt$trait_name%in%names(sp_traits),]

sp_dist <- mFD::funct.dist(  # might want to change w.type
  sp_tr         = sp_traits,
  tr_cat        = tt,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# above appears to work but check with gawdis
# calculate functional distances
gaw.groups <- gawdis(sp_traits, w.type = "equal", 
                   groups = c(1,2,3,3,3,3,3,4,5,5,6,6,6,6,6,7,7,7,8,8,8,9,10,11,12,13,14,15,16,17),
                   fuzzy=c(3,5,6,7,8))   # warning about unbalanced distribution of some traits 
                                         # (traits where most species belong to one category)   

# remove the unbalanced traits
traits_to_remove <- c(
  "Terrestrial.vertebrate..0.1.", "Fish..0.1.", "Carrion..0.1.",
  "DigitsClaws", "LocoAquatic", "LocoFossorial", "LocoCursorial",
  "day_active", "niche_aerial", "niche_arboreal"
)
sp_traits_clean <- sp_traits[, !(names(sp_traits) %in% traits_to_remove)]
tt_clean <- tt[!(tt$trait_name %in% traits_to_remove), ]

# rerun to calculate distances
# calculate functional distances
gdist <- gawdis(sp_traits_clean,w.type = "equal", 
              groups = c(1,2,3,3,4,5,6,6,7,7,8,9,10,11,12,13,14,15,16,17),
              fuzzy=c(3,6,7))      

attr(gdist,"correls") # some traits contribute more to variation in trait space
attr(gdist,"weights") # all have positive weights

# compute multidimensional functional spaces (PCoA) and assess their quality
qual <- quality.fspaces(sp_dist = gdist, fdendro = "average",maxdim_pcoa = 44,
                        deviation_weighting = c("absolute", "squared"),fdist_scaling = c(TRUE, FALSE)) 

# how much variation explained by each PCoA?
round(qual$"quality_fspaces", 3)

# position of species on PCoA axes
sp_coords <- qual$details_fspaces$sp_pc_coord

# calculate functional diversity indices
# need matrix of 1s and 0s where row = community and column = species -> the presence_matrix
# cut matrix to islands with 3 or more species
rs <- which(rowSums(island_pres_mat) > 2)
enough <- island_pres_mat[rs,]

dim(sp_coords)
dim(enough)

FD <- alpha.fd.multidim(sp_faxes_coord = sp_coords[,paste("PC",1:2, sep="")], asb_sp_w = as.matrix(enough)) # cut to islands with 3 or more species, can only use 2 PCoAs
inds <- FD$functional_diversity_indices  
inds$FID <- as.numeric(rownames(inds)) # add FID to indices

# merge with island data
head(inds)
head(ausisldatcorr)
ausisles.mamala.traitIndices <- merge(inds, ausisldatcorr, by="FID", all.x=TRUE)
head(ausisles.mamala.traitIndices)

# trait diversity by island area
plot(log10(ausisles.mamala.traitIndices$area), log10(ausisles.mamala.traitIndices$fric), pch=19)

head(ausisles.mamala.traitIndices)
ausisles.mamala.traitIndices$larea <- log10(ausisles.mamala.traitIndices$area_2/1e6)
ausisles.mamala.traitIndices$lfric <- log10(ausisles.mamala.traitIndices$fric) 
ausisles.mamala.traitIndices$lSR <- log10(ausisles.mamala.traitIndices$sp_richn) 

# SAR
## power-law values
SARlinfit2 <- lm(lSR ~ larea, data=ausisles.mamala.traitIndices) # with fewer islands in this dataset
summary(SARlinfit2)

SARlinfit.slope <- SARlinfit2$coefficients[2]
SARlinfit.intercept <- SARlinfit2$coefficients[1]
SARlinfit.slope # most frequent z = 0.2–0.4 https://www.pnas.org/doi/10.1073/pnas.0510605103
SARlinfit.slope.se <- summary(SARlinfit2)$coefficients[4]
SARlinfit.slope.se
SARlinfit.intercept
10^SARlinfit.intercept

SARlinpred <- predict(SARlinfit2, newdata=data.frame(larea=ausisles.mamala.traitIndices$larea), interval = "confidence", level = 0.95)
SARlinpred.out <- data.frame(larea = ausisles.mamala.traitIndices$larea, SARlinpred)
head(SARlinpred.out)
SARlinpred.out$sarpred <- 10^SARlinpred.out$fit
SARlinpred.out$sarpredlo <- 10^SARlinpred.out$lwr
SARlinpred.out$sarpredup <- 10^SARlinpred.out$upr
SARlinpred.out$areakm2 <- 10^(SARlinpred.out$larea)
head(SARlinpred.out)

SARlinpred.sort <- SARlinpred.out[order(SARlinpred.out$areakm2),]
head(SARlinpred.sort)
dim(SARlinpred.sort)
head(ausisles.mamala.traitIndices)

# output data for outside graphing
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(SARlinpred.sort, file="AUS_SARlinpredmam.csv", row.names=FALSE)

## functional richness
FRIClinfit <- lm(lfric ~ larea, data=ausisles.mamala.traitIndices)
summary(FRIClinfit)
FRIClinfit.slope <- FRIClinfit$coefficients[2]
FRIClinfit.intercept <- FRIClinfit$coefficients[1]
FRIClinfit.slope # most frequent z = 0.2–0.4 https://www.pnas.org/doi/10.1073/pnas.0510605103
FRIClinfit.slope.se <- summary(FRIClinfit)$coefficients[4]
FRIClinfit.slope.se
FRIClinfit.intercept
10^FRIClinfit.intercept

FRIClinpred <- predict(FRIClinfit, newdata=data.frame(larea=ausisles.mamala.traitIndices$larea), interval = "confidence", level = 0.95)
FRIClinpred.out <- data.frame(larea = ausisles.mamala.traitIndices$larea, FRIClinpred)
head(FRIClinpred.out)
FRIClinpred.out$fricpred <- 10^FRIClinpred.out$fit
FRIClinpred.out$fricpredlo <- 10^FRIClinpred.out$lwr
FRIClinpred.out$fricpredup <- 10^FRIClinpred.out$upr
FRIClinpred.out$areakm2 <- 10^(FRIClinpred.out$larea)
head(FRIClinpred.out)

FRIClinpred.sort <- FRIClinpred.out[order(FRIClinpred.out$areakm2),]
head(FRIClinpred.sort)
dim(FRIClinpred.sort)
head(ausisles.mamala.traitIndices)

# output data for outside graphing
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(FRIClinpred.sort, file="AUS_FRIClinpredmam.csv", row.names=FALSE)

## functional nearest-neighbour distance
head(ausisles.mamala.traitIndices)
ausisles.mamala.traitIndices$lfnnd <- log10(ausisles.mamala.traitIndices$fnnd) # log transform
head(ausisles.mamala.traitIndices)

FNNDlinfit <- lm(lfnnd ~ larea, data=ausisles.mamala.traitIndices)
summary(FNNDlinfit)
FNNDlinfit.slope <- FNNDlinfit$coefficients[2]
FNNDlinfit.intercept <- FNNDlinfit$coefficients[1]
FNNDlinfit.slope # most frequent z = 0.2–0.4 https://www.pnas.org/doi/10.1073/pnas.0510605103
FNNDlinfit.slope.se <- summary(FNNDlinfit)$coefficients[4]
FNNDlinfit.slope.se
FNNDlinfit.intercept
10^FNNDlinfit.intercept

FNNDlinpred <- predict(FNNDlinfit, newdata=data.frame(larea=ausisles.mamala.traitIndices$larea), interval = "confidence", level = 0.95)
FNNDlinpred.out <- data.frame(larea = ausisles.mamala.traitIndices$larea, FNNDlinpred)
head(FNNDlinpred.out)
FNNDlinpred.out$fricpred <- 10^FNNDlinpred.out$fit
FNNDlinpred.out$fricpredlo <- 10^FNNDlinpred.out$lwr
FNNDlinpred.out$fricpredup <- 10^FNNDlinpred.out$upr
FNNDlinpred.out$areakm2 <- 10^(FNNDlinpred.out$larea)
head(FNNDlinpred.out)

FNNDlinpred.sort <- FNNDlinpred.out[order(FNNDlinpred.out$areakm2),]
head(FNNDlinpred.sort)
dim(FNNDlinpred.sort)
head(ausisles.mamala.traitIndices)

# output data for outside graphing
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(FNNDlinpred.sort, file="AUS_FNNDlinpredmam.csv", row.names=FALSE)



# island area versus species richness
FD.SR.area.plot <- ggplot(ausisles.mamala.traitIndices, aes(x=area, y=sp_richn)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("island area (m^2)") +
  ylab("mammal species richness") +
  ggtitle("species richness vs. island area") +
  theme_minimal()
FD.SR.area.plot

# plot functional trait index relationships to island area
FD.fric.plot <- ggplot(ausisles.mamala.traitIndices, aes(x=area, y=fric)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("island area (m^2)") +
  ylab("mammal species functional richness") +
  ggtitle("functional richness vs. island area") +
  theme_minimal()
FD.fric.plot

FD.fnnd.plot <- ggplot(ausisles.mamala.traitIndices, aes(x=area, y=fnnd)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("island area (m^2)") +
  ylab("mammal functional nearest neighbour distance") +
  ggtitle("functional nearest neighbour distance vs. island area") +
  theme_minimal()
FD.fnnd.plot

grid.arrange(FD.SR.area.plot, FD.fric.plot, FD.fnnd.plot, ncol=3, nrow=1)




# species richness vs. functional diversity
FD.SR.fdis.plot <- ggplot(ausisles.mamala.traitIndices, aes(x=sp_richn, y=fdis)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("mammal species richness") +
  ylab("functional dispersion") +
  ggtitle("functional dispersion vs. species richness") +
  theme_minimal()
FD.SR.fdis.plot


FD.SR.fnnd.plot <- ggplot(ausisles.mamala.traitIndices, aes(x=sp_richn, y=fnnd)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("mammal species richness") +
  ylab("functional nearest neighbour distance") +
  ggtitle("functional nearest neighbour distance vs. species richness") +
  theme_minimal()
FD.SR.fnnd.plot


## merge islSRnoZero & for export
mamTraitSRout <- merge(mamausislnoZero, ausisles.mamala.traitIndices, by="FID", all.x=F)
head(mamTraitSRout)
dim(mamTraitSRout)
setwd("/Users/brad0317/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(mamTraitSRout, file="AUSislmamTraitSR.csv", row.names=FALSE)


# scale SR & FR data to examine z on same scale
head(mamTraitSRout)
mamTraitSRout$lSR.sc <- scale(mamTraitSRout$lSR, center = TRUE, scale = TRUE)
mamTraitSRout$lfric.sc <- scale(logit(mamTraitSRout$fric), center = TRUE, scale = TRUE)

SARlinfit.sc <- lm(lSR.sc ~ log10(area_2.x/1e6), data=mamTraitSRout)
summary(SARlinfit.sc)

FRIClinfit.sc <- lm(lfric.sc ~ log10(area_2.x/1e6), data=mamTraitSRout)
summary(FRIClinfit.sc)
plot(log10(mamTraitSRout$area_2.x/1e6), mamTraitSRout$lfric.sc, xlab="island area (km2)",
     ylab="scaled logit functional richness", cex=0.8, col="lightblue")
abline(FRIClinfit.sc, col="blue", lwd=2, lty=2)



#############################
#############################
## bird traits (AVONET)
#############################
#############################

setwd("~/Documents/Papers/Biogeography/Aus Isl traits/data/bird traits/")

# read in 3 versions (different taxonomy used in each)
bl <- read_excel ("TraitData/AVONET1_Birdlife.xlsx", sheet = "AVONET1_BirdLife")
eb <- read_excel ("TraitData/AVONET2_eBird.xlsx", sheet = "AVONET2_eBird")
bt <- read_excel ("TraitData/AVONET3_BirdTree.xlsx", sheet = "AVONET3_BirdTree")

# missing data appears to have been imputed
dim(bl)
colnames(bl)
table(bl$Beak.Length_Culmen > 0)
table(bl$Mass > 0)

# colnames differ among datasets - seems to just be taxonomic and distribution information
#View(cbind(names(bl),(names(bl) %in% names(eb))))
#View(cbind(names(bl),(names(bl) %in% names(bt))))
#View(cbind(names(bt),(names(bt) %in% names(bl))))

names(bl)
keep <- c("Species1","Beak.Length_Culmen", "Beak.Length_Nares", "Beak.Width", "Beak.Depth", "Tarsus.Length","Wing.Length",
          "Kipps.Distance","Secondary1", "Hand-Wing.Index", "Tail.Length", "Mass","Mass.Source", "Mass.Refs.Other",
          "Inference", "Traits.inferred", "Reference.species",  "Habitat", "Habitat.Density", "Migration", "Trophic.Level",
          "Trophic.Niche", "Primary.Lifestyle", "Min.Latitude", "Max.Latitude", "Centroid.Latitude", "Centroid.Longitude",
          "Range.Size")

rm_keep <- c("Beak.Length_Culmen", "Secondary1","Kipps.Distance", "Mass.Source", "Mass.Refs.Other","Inference", 
            "Traits.inferred", "Reference.species", "Habitat", "Trophic.Niche", "Min.Latitude", "Max.Latitude",
            "Centroid.Latitude",  "Centroid.Longitude", "Range.Size") 
# removed "Habitat" and"Trophic.Niche" because there are similar variables with broader categories
# (less problematic when calculating Gower distances)
# removed distribution data to keep the island analyses independent of species distribution

# make the selection adjustment
keep <- keep[!keep%in%rm_keep]

bl_traits <- bl[,names(bl) %in% keep]
colnames(bl_traits)

# check for correlations with mass
# select only numeric columns
numeric_data <- bl_traits[sapply(bl_traits, is.numeric)]

# remove Mass column from set of predictors
other_vars <- setdiff(names(numeric_data), "Mass")

# remove rows with NA or non-positive values (log requires positive numbers)
filtered_data <- numeric_data[complete.cases(numeric_data[, c("Mass", other_vars)]), ]
filtered_data <- filtered_data[apply(filtered_data[, c("Mass", other_vars)], 1, function(x) all(x > 0)), ]

# log-transform data
log_data <- log10(filtered_data)

# compute correlations between log(Mass) and log(other variables)
sapply(other_vars, function(col) cor(log_data[[col]], log_data$Mass))

# Some are correlated with mass. Consider using residuals against mass for beak length, beak width, 
# beak depth, tarsus length, wing length, and tail length

# examine data
colnames(bl_traits)
bl_traits$Species1

# any NAs?
na_counts <- sapply(bl_traits, function(x) sum(is.na(x)))
na_counts

# merge with island ALA dataset
head(ausisles.avala.df)

# unique species only
ausisles.avala.df.unique <- unique(ausisles.avala.df[,c("FID", "scientific")])
head(ausisles.avala.df.unique)
colnames(ausisles.avala.df.unique)[2] <- 'Species1' # rename to match AVONET

# any NAs?
na_counts2 <- sapply(ausisles.avala.df.unique, function(x) sum(is.na(x)))
na_counts2

ausisles.avala.traits <- merge(ausisles.avala.df.unique, bl_traits, by="Species1", all.x=F)
head(ausisles.avala.traits)
dim(ausisles.avala.traits)
table(ausisles.avala.traits$FID)

# any NAs?
na_counts3 <- sapply(ausisles.avala.traits, function(x) sum(is.na(x)))
na_counts3

# which don't have any matches?
colnames(bl_traits)
unique(ausisles.avala.traits$Species1[is.na(ausisles.avala.traits$Mass)]) # all good now

# make a list of unique species and traits
avsp_traits <- ausisles.avala.traits
head(avsp_traits)
avsp_traits <- avsp_traits[,!names(avsp_traits) %in% c("FID")]
avsp_traits <- unique(avsp_traits)

# any NAs?
na_counts4 <- sapply(avsp_traits, function(x) sum(is.na(x)))
na_counts4

# check distributions of the traits being used - transform if skewed
hist(avsp_traits$Mass) # needs to be log transformed (skewed) and renamed logBodyMass (or something like that)
avsp_traits$Log_Body_mass <- log(avsp_traits$Mass)
hist(avsp_traits$Log_Body_mass)
avsp_traits$Mass <- NULL

hist(avsp_traits$Wing.Length) # needs to be log transformed (skewed) and renamed log...
avsp_traits$Log_WL <- log(avsp_traits$Wing.Length)
hist(avsp_traits$Log_WL)
avsp_traits$Wing.Length <- NULL

hist(avsp_traits$Tail.Length)
avsp_traits$Log_TailL <- log(avsp_traits$Tail.Length)
hist(avsp_traits$Log_TailL)
avsp_traits$Tail.Length <- NULL

hist(avsp_traits$Tarsus.Length)
avsp_traits$Log_TarsL <- log(avsp_traits$Tarsus.Length)
hist(avsp_traits$Log_TarsL)
avsp_traits$Tarsus.Length <- NULL

hist(avsp_traits$Beak.Length_Nares)
avsp_traits$Log_BeakLNares <- log(avsp_traits$Beak.Length_Nares)
hist(avsp_traits$Log_BeakLNares)
avsp_traits$Beak.Length_Nares <- NULL

hist(avsp_traits$Beak.Width)
avsp_traits$Log_BeakW <- log(avsp_traits$Beak.Width)
hist(avsp_traits$Log_BeakW)
avsp_traits$Beak.Width <- NULL

hist(avsp_traits$'Hand-Wing.Index')
avsp_traits$Log_HWI <- log(avsp_traits$'Hand-Wing.Index')
hist(avsp_traits$Log_HWI)
avsp_traits$'Hand-Wing.Index' <- NULL

table(avsp_traits$Migration)
# str(avsp_traits$Migration)
# avsp_traits$Migration <- as.integer(avsp_traits$Migration) # change to integer
# str(avsp_traits$Migration)
# hist(avsp_traits$Migration)

str(avsp_traits$Trophic.Level)
table(avsp_traits$Trophic.Level)

# check trophic level from mammals for recoding birds
mamtraits$scientific
mamtraits[which(mamtraits$scientific == "Sarcophilus harrisii"),"trophic_level"] # carnivore trophic level = 3
mamtraits[which(mamtraits$scientific == "Wallabia bicolor"),"trophic_level"] # herbivore trophic level = 1

avsp_traits$trophic_level <- as.integer(ifelse(avsp_traits$Trophic.Level == "Carnivore", 3,
                                      ifelse(avsp_traits$Trophic.Level == "Herbivore", 1,
                                        ifelse(avsp_traits$Trophic.Level == "Omnivore", 2, NA))))
table(avsp_traits$trophic_level)
avsp_traits$Trophic.Level <- NULL

table(avsp_traits$Habitat.Density)
 
table(avsp_traits$Primary.Lifestyle)

# any NAs?
na_counts4 <- sapply(avsp_traits, function(x) sum(is.na(x)))
na_counts4

dim(avsp_traits)
avsp_traits$Species1


# check for NAs
na_counts <- sapply(avsp_traits, function(x) sum(is.na(x)))
na_counts

## make list of mammal communities by island (FID)
ausisles.avala.traits.list <- split(ausisles.avala.traits, ausisles.avala.traits$FID)
#View(ausisles.avala.traits.list[[which.max(sapply(ausisles.avala.traits.list, nrow))]]) # view the island with the most records

# make a presence absence matrix with row = FID and colname = species
combined_unique <- unique(ausisles.avala.traits[, c("FID", "Species1")])
island_pres <- table(combined_unique$FID, combined_unique$Species1)

# convert table to a data frame with species as columns
island_pres_df <- as.data.frame.matrix(island_pres)

# any NAs?
na_counts5 <- sapply(island_pres_df, function(x) sum(is.na(x)))
na_counts5

# add FID as its own column (from row names)
island_pres_df$FID <- rownames(island_pres_df)

# move FID to the first column
island_pres_df <- island_pres_df[, c("FID", setdiff(names(island_pres_df), "FID"))]

# reset row names
rownames(island_pres_df) <- NULL
head(island_pres_df)

# trait type data frame
avtt <- data.frame(trait_name = names(avsp_traits), trait_type=NA)

# fix the trait_type column
avtt$trait_name
avtt$trait_type[avtt$trait_name %in% c("Primary.Lifestyle")] <- "N"
avtt$trait_type[avtt$trait_name %in% c("Habitat.Density","Migration")] <- "O"
avtt$trait_type[avtt$trait_name %in% c("Beak.Depth","Log_HWI","trophic_level",
                                   "Log_Body_mass","Log_WL","Log_TailL","Log_BeakLNares","Log_TarsL",
                                   "Log_BeakW","Log_litters_yr",
                                   "Log_weaning_age","Log_brain_mass")] <- "Q"
avtt$trait_type

# add the fuzzy column
# tt$fuzzy_name[tt$trait_name %in% c("Plant..0.1.","Terrestrial.vertebrate..0.1.","Fish..0.1.",
#                                    "Invertebrate..0.1.","Carrion..0.1.")] <- "diet"
# tt$fuzzy_name[tt$trait_name %in% c("DigitsWebbed","DigitsClaws","DigitsFlippers","DigitsHooves")] <- "digits"
# tt$fuzzy_name[tt$trait_name %in% c("LocoAquatic","LocoArboreal","LocoFossorial","LocoGliding",
#                                    "LocoCursorial","LocoSaltation","LocoTerrestrial")] <- "locomotion"
# tt$fuzzy_name[tt$trait_name %in% c("day_active","night_active","crepuscular_active")] <- "activity_time"
# tt$fuzzy_name[tt$trait_name%in%c("niche_aerial","niche_arboreal","niche_ground","niche_marine",
#                                  "niche_scansorial")] <- "niche"
#tt$fuzzy_name <- ifelse(is.na(tt$fuzzy_name), tt$trait_name, tt$fuzzy_name)

# remove non-trait rows  
avtt <- avtt[!avtt$trait_name %in% c("Species1","FID"),]

# species as row name in avsp_traits
row.names(avsp_traits) <- avsp_traits$Species1
avsp_traits$Species1 <- NULL
head(avsp_traits)

# format data in avsp_traits. mFD required particular formatting
# loop over each trait defined in tt
for (i in seq_len(nrow(avtt))) {
  trait <- avtt$trait_name[i]
  type <- avtt$trait_type[i]
  
  if (trait %in% colnames(avsp_traits)) {
    if (type == "N") {
      avsp_traits[[trait]] <- as.factor(avsp_traits[[trait]])
    } else if (type == "Q" || type == "F") {
      avsp_traits[[trait]] <- as.numeric(avsp_traits[[trait]])
    } else if (type == "O") {
        avsp_traits[[trait]] <- as.ordered(avsp_traits[[trait]])
    } else {
      warning(paste("unknown trait type for", trait))
    }
  } else {
    warning(paste("trait", trait, "not found in avsp_traits"))
  }
}

# any NAs?
na_counts6 <- sapply(avsp_traits, function(x) sum(is.na(x)))
na_counts6

# species traits summary:
traits_summ <- mFD::sp.tr.summary(
  tr_cat     = avtt,   
  sp_tr      = avsp_traits, 
  stop_if_NA = TRUE)

traits_summ$"tr_types" 
traits_summ$"mod_list"

# Summary of the assemblages * species dataframe:
island_pres_mat <- island_pres_df
row.names(island_pres_mat) <- island_pres_mat$FID
island_pres_mat$FID <- NULL
island_pres_mat[] <- lapply(island_pres_mat, as.numeric)
island_pres_mat <- as.matrix(island_pres_mat)

isl_summ <- mFD::asb.sp.summary(asb_sp_w = island_pres_mat)
isl_summ$"asb_sp_richn" # sp richness by island (non-zero only)

# remove variables with no variation
#avsp_traits$DigitsFlippers <- NULL
#avsp_traits$DigitsHooves <- NULL  # can change this if invasives are added
#avsp_traits$LocoGliding <- NULL
#avsp_traits$LocoSaltation <- NULL
#avsp_traits$niche_marine <- NULL
#avsp_traits$niche_scansorial <- NULL

# no cut trait type df back to traits in sp_traits
avtt <- avtt[avtt$trait_name %in% names(avsp_traits),]

sp_dist <- mFD::funct.dist(  # might want to change w.type
  sp_tr         = avsp_traits,
  tr_cat        = avtt,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# above appears to work but check with gawdis
# calculate functional distances
colnames(avsp_traits)
gdist <- gawdis(avsp_traits, w.type = "equal", 
                     groups = c(1,2,3,4,1,1,1,1,1,1,1,5)) 

# # remove the unbalanced traits
# traits_to_remove <- c(
#   "Terrestrial.vertebrate..0.1.", "Fish..0.1.", "Carrion..0.1.",
#   "DigitsClaws", "LocoAquatic", "LocoFossorial", "LocoCursorial",
#   "day_active", "niche_aerial", "niche_arboreal"
# )
# sp_traits_clean <- sp_traits[, !(names(sp_traits) %in% traits_to_remove)]
# tt_clean <- tt[!(tt$trait_name %in% traits_to_remove), ]

# rerun to calculate distances
# calculate functional distances
# gdist <- gawdis(sp_traits_clean,w.type = "equal", 
#                 groups = c(1,2,3,3,4,5,6,6,7,7,8,9,10,11,12,13,14,15,16,17),
#                 fuzzy=c(3,6,7))      

attr(gdist,"correls") # some traits contribute more to variation in trait space
attr(gdist,"weights") # all have positive weights

# compute multidimensional functional spaces (PCoA) and assess their quality
qual <- quality.fspaces(sp_dist = gdist, fdendro = "average", maxdim_pcoa = 44,
                        deviation_weighting = c("absolute", "squared"), fdist_scaling = c(TRUE, FALSE)) 

# how much variation explained by each PCoA?
round(qual$"quality_fspaces", 3)

# position of species on PCoA axes
sp_coords <- qual$details_fspaces$sp_pc_coord

# calculate functional diversity indices
# need matrix of 1s and 0s where row = community and column = species -> the presence_matrix
# cut matrix to islands with 3 or more species
rs <- which(rowSums(island_pres_mat) > 2)
enough <- island_pres_mat[rs,]

dim(sp_coords)
dim(enough)

FD <- alpha.fd.multidim(sp_faxes_coord = sp_coords[,paste("PC",1:2, sep="")],
                        asb_sp_w = as.matrix(enough)) # cut to islands with 3 or more species,
                                                      # can only use 2 PCoAs
avinds <- FD$functional_diversity_indices  
avinds$FID <- as.numeric(rownames(avinds)) # add FID to indices

# merge with island data
head(avinds)
head(ausisldatcorr)
ausisles.avala.traitIndices <- merge(avinds, ausisldatcorr, by="FID", all.x=TRUE)
head(ausisles.avala.traitIndices)

# trait diversity by island area
plot(log10(ausisles.avala.traitIndices$area), log10(ausisles.avala.traitIndices$fric), pch=19)

head(ausisles.avala.traitIndices)
ausisles.avala.traitIndices$larea <- log10(ausisles.avala.traitIndices$area_2/1e6)
ausisles.avala.traitIndices$lfric <- log10(ausisles.avala.traitIndices$fric) 
ausisles.avala.traitIndices$lSR <- log10(ausisles.avala.traitIndices$sp_richn) 

# SAR
## power-law values
SARlinfit2 <- lm(lSR ~ larea, data=ausisles.avala.traitIndices) # with fewer islands in this dataset
summary(SARlinfit2)

SARlinfit.slope <- SARlinfit2$coefficients[2]
SARlinfit.intercept <- SARlinfit2$coefficients[1]
SARlinfit.slope # most frequent z = 0.2–0.4 https://www.pnas.org/doi/10.1073/pnas.0510605103
SARlinfit.slope.se <- summary(SARlinfit2)$coefficients[4]
SARlinfit.slope.se
SARlinfit.intercept
10^SARlinfit.intercept

SARlinpred <- predict(SARlinfit2, newdata=data.frame(larea=ausisles.avala.traitIndices$larea), 
                      interval = "confidence", level = 0.95)
SARlinpred.out <- data.frame(larea = ausisles.avala.traitIndices$larea, SARlinpred)
head(SARlinpred.out)
SARlinpred.out$sarpred <- 10^SARlinpred.out$fit
SARlinpred.out$sarpredlo <- 10^SARlinpred.out$lwr
SARlinpred.out$sarpredup <- 10^SARlinpred.out$upr
SARlinpred.out$areakm2 <- 10^(SARlinpred.out$larea)
head(SARlinpred.out)

SARlinpred.sort <- SARlinpred.out[order(SARlinpred.out$areakm2),]
head(SARlinpred.sort)
dim(SARlinpred.sort)
head(ausisles.avala.traitIndices)

# output data for outside graphing
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(SARlinpred.sort, file="AUS_SARlinpredav.csv", row.names=FALSE)

## functional richness
FRIClinfit <- lm(lfric ~ larea, data=ausisles.avala.traitIndices)
summary(FRIClinfit)
FRIClinfit.slope <- FRIClinfit$coefficients[2]
FRIClinfit.intercept <- FRIClinfit$coefficients[1]
FRIClinfit.slope # most frequent z = 0.2–0.4 https://www.pnas.org/doi/10.1073/pnas.0510605103
FRIClinfit.slope.se <- summary(FRIClinfit)$coefficients[4]
FRIClinfit.slope.se
FRIClinfit.intercept
10^FRIClinfit.intercept

FRIClinpred <- predict(FRIClinfit, newdata=data.frame(larea=ausisles.avala.traitIndices$larea), interval = "confidence", level = 0.95)
FRIClinpred.out <- data.frame(larea = ausisles.avala.traitIndices$larea, FRIClinpred)
head(FRIClinpred.out)
FRIClinpred.out$fricpred <- 10^FRIClinpred.out$fit
FRIClinpred.out$fricpredlo <- 10^FRIClinpred.out$lwr
FRIClinpred.out$fricpredup <- 10^FRIClinpred.out$upr
FRIClinpred.out$areakm2 <- 10^(FRIClinpred.out$larea)
head(FRIClinpred.out)

FRIClinpred.sort <- FRIClinpred.out[order(FRIClinpred.out$areakm2),]
head(FRIClinpred.sort)
dim(FRIClinpred.sort)
head(ausisles.mamala.traitIndices)

# output data for outside graphing
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(FRIClinpred.sort, file="AUS_FRIClinpredav.csv", row.names=FALSE)

## functional nearest-neighbour distance
head(ausisles.avala.traitIndices)
ausisles.avala.traitIndices$lfnnd <- log10(ausisles.avala.traitIndices$fnnd) # log transform
head(ausisles.avala.traitIndices)

FNNDlinfit <- lm(lfnnd ~ larea, data=ausisles.avala.traitIndices)
summary(FNNDlinfit)
FNNDlinfit.slope <- FNNDlinfit$coefficients[2]
FNNDlinfit.intercept <- FNNDlinfit$coefficients[1]
FNNDlinfit.slope # most frequent z = 0.2–0.4 https://www.pnas.org/doi/10.1073/pnas.0510605103
FNNDlinfit.slope.se <- summary(FNNDlinfit)$coefficients[4]
FNNDlinfit.slope.se
FNNDlinfit.intercept
10^FNNDlinfit.intercept

FNNDlinpred <- predict(FNNDlinfit, newdata=data.frame(larea=ausisles.avala.traitIndices$larea), interval = "confidence", level = 0.95)
FNNDlinpred.out <- data.frame(larea = ausisles.avala.traitIndices$larea, FNNDlinpred)
head(FNNDlinpred.out)
FNNDlinpred.out$fricpred <- 10^FNNDlinpred.out$fit
FNNDlinpred.out$fricpredlo <- 10^FNNDlinpred.out$lwr
FNNDlinpred.out$fricpredup <- 10^FNNDlinpred.out$upr
FNNDlinpred.out$areakm2 <- 10^(FNNDlinpred.out$larea)
head(FNNDlinpred.out)

FNNDlinpred.sort <- FNNDlinpred.out[order(FNNDlinpred.out$areakm2),]
head(FNNDlinpred.sort)
dim(FNNDlinpred.sort)
head(ausisles.mamala.traitIndices)

# output data for outside graphing
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(FNNDlinpred.sort, file="AUS_FNNDlinpredav.csv", row.names=FALSE)



# island area versus species richness
FD.SR.area.plot <- ggplot(ausisles.avala.traitIndices, aes(x=area, y=sp_richn)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("island area (m^2)") +
  ylab("avian species richness") +
  ggtitle("species richness vs. island area") +
  theme_minimal()
FD.SR.area.plot

# plot functional trait index relationships to island area
FD.fric.plot <- ggplot(ausisles.avala.traitIndices, aes(x=area, y=fric)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("island area (m^2)") +
  ylab("avian species functional richness") +
  ggtitle("functional richness vs. island area") +
  theme_minimal()
FD.fric.plot

FD.fric.logit.plot <- ggplot(ausisles.avala.traitIndices, aes(x=area, y=logit(fric))) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  #scale_y_log10() +
  xlab("island area (m^2)") +
  ylab("logit avian species functional richness") +
  ggtitle("logit functional richness vs. island area") +
  theme_minimal()
FD.fric.logit.plot
FRIClogitfit <- lm(logit(fric) ~ larea, data=ausisles.avala.traitIndices)
summary(FRIClogitfit)


FD.fnnd.plot <- ggplot(ausisles.avala.traitIndices, aes(x=area, y=fnnd)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("island area (m^2)") +
  ylab("avian functional nearest neighbour distance") +
  ggtitle("functional nearest neighbour distance vs. island area") +
  theme_minimal()
FD.fnnd.plot

grid.arrange(FD.SR.area.plot, FD.fric.plot, FD.fnnd.plot, ncol=3, nrow=1)




# species richness vs. functional diversity
FD.SR.fdis.plot <- ggplot(ausisles.avala.traitIndices, aes(x=sp_richn, y=fdis)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("avian species richness") +
  ylab("functional dispersion") +
  ggtitle("functional dispersion vs. species richness") +
  theme_minimal()
FD.SR.fdis.plot


FD.SR.fnnd.plot <- ggplot(ausisles.avala.traitIndices, aes(x=sp_richn, y=fnnd)) +
  geom_point() +
  geom_smooth(method="lm", se=T, color="blue") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("avian species richness") +
  ylab("functional nearest neighbour distance") +
  ggtitle("functional nearest neighbour distance vs. species richness") +
  theme_minimal()
FD.SR.fnnd.plot


## merge islSRnoZero & for export
avTraitSRout <- merge(avausislnoZero, ausisles.avala.traitIndices, by="FID", all.x=F)
head(avTraitSRout)
dim(avTraitSRout)
setwd("/Users/brad0317/Documents/Papers/Biogeography/Aus Isl traits/out")
write.csv(avTraitSRout, file="AUSislavTraitSR.csv", row.names=FALSE)

hist(logit(avTraitSRout$fric))


# scale SR & FR data to examine z on same scale
head(avTraitSRout)
avTraitSRout$lSR.sc <- scale(avTraitSRout$lSR, center = TRUE, scale = TRUE)
avTraitSRout$lfric.sc <- scale(logit(avTraitSRout$fric), center = TRUE, scale = TRUE)

FRIClinfit.sc <- lm(lfric.sc ~ log10(area_2.x/1e6), data=avTraitSRout)
summary(FRIClinfit.sc)
plot(log10(avTraitSRout$area_2.x/1e6), avTraitSRout$lfric.sc, xlab="island area (km2)",
     ylab="scaled logit functional richness", pch=19, cex=0.8, col="pink")
abline(FRIClinfit.sc, col="red", lwd=2, lty=2)


## save workspace
setwd("~/Documents/Papers/Biogeography/Aus Isl traits/R images")
save.image(file="AUS_Isl_traits.RData")
