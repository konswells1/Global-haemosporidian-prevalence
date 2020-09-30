

########################
#
# Avian prevalence and lineage distribution
#
################



##########
# Settings Kons-device

Drive <- "D"
dir.library <- paste(Drive, ":/Kons_backup/R Library", sep = "")
.libPaths(dir.library)

dir.stem <- paste(Drive,  ":/Kons_backup/Projekte_kons/2020.05_AvMalaria.Glob_prevalence/", sep = "")
dir.data <- paste(dir.stem,  "AvMal.prev_Data", sep = "")
dir.analysis <- paste(dir.stem,  "avmprev_analysis/", sep = "")
dir.env <- paste(dir.stem,  "AvMal.prev_EnvData", sep = "")

dir.geodata <- paste(Drive,  ":/GeoData/", sep = "")

setwd(dir.stem)


##########
# Load R packages

library(Hmisc)
library(raster)
library(ape)
library(dplyr)
library(ape)
library(phytools)
library(maptools)
library(rgdal)
library(rgeos)
library(sp)
library(sf)
library(ggplot2)
library(ggthemes)
library(ggtree)



WGS84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"


##############
# Geodata



library(rnaturalearth)
worldmap <- ne_countries(scale = "medium", type = "countries", returnclass = "sf")
ggplot(worldmap) + geom_sf()


# WorldClim Data
##Worldclim <- getData("worldclim",var="bio",res=5)

# Zoogeographical regions
setwd(dir.data)
poly_Realms <- st_read("./zoo.regions_Holt2013/newRealms.shp")


realm_name <- sort(unique(poly_Realms$Realm))
nrealm <- length(realm_name)


##########
# Load Bird tree
setwd(dir.data)
BirdTree100 <- read.tree("BirdTree100.tre")
##BirdTreeC <- averageTree(BirdTree100[1:3])
BirdTree1 <- BirdTree100[[1]]
BirdTree1$tip.label <- sub("_"," ",BirdTree1$tip.label)
BirdTree1$tip.label <- tolower(BirdTree1$tip.label)
rm(BirdTree100)

#############
# Avian Malaria Data

setwd(dir.data)
# Data set global from Alan Fecchio 04/05/20
###Dat.avm.Global <- read.csv("GlobalDataSet_prepared_04052020.csv", header=TRUE, stringsAsFactors = FALSE)  
Dat.avm.Global <- read.csv("GlobalDataSet_23062020.csv", header=TRUE, stringsAsFactors = FALSE)  

## Omit record without coordinates
#delete <- which(is.na(Dat.avm.Global$Latitude) | is.na(Dat.avm.Global$Longitude))
#Dat.avm.Global <- Dat.avm.Global[-delete,]

# Data set Ecol Lett
Dat.avm.AusAsia <- read.csv("Australasia_200623.csv", header=TRUE, stringsAsFactors = FALSE) 

## Prepare infection data from global data set
d_scren.result <- Dat.avm.Global$Screen.Result
d_lineage.1.genus <- Dat.avm.Global$Lineage.1.Genus
d_lineage.1.name <- Dat.avm.Global$Lineage.1.Name
d_lineage.1.access.no <- Dat.avm.Global$Lineage.1.Accession..
d_lineage.1.cytB <- Dat.avm.Global$Lineage.1.Cyt.B
d_lineage.2.genus <- Dat.avm.Global$Lineage.2.Genus
d_lineage.2.name <- Dat.avm.Global$Lineage.2.Name
d_lineage.2.access.no <- Dat.avm.Global$Lineage.2.Accession..
d_lineage.2.cytB <- Dat.avm.Global$Lineage.2.Cyt.B
d_lineage.3.genus <- Dat.avm.Global$Lineage.3.Genus
d_lineage.3.name <- Dat.avm.Global$Lineage.3.Name
d_lineage.3.access.no <- Dat.avm.Global$Lineage.3.Accession..
d_lineage.3.cytB <- Dat.avm.Global$Lineage.3.Cyt.B
d_lineage.4.genus <- Dat.avm.Global$Lineage.4.Genus
d_lineage.4.name <- Dat.avm.Global$Lineage.4.Name
d_lineage.4.access.no <- Dat.avm.Global$Lineage.4.Accession..
d_lineage.4.cytB <- Dat.avm.Global$Lineage.4.Cyt.B
d_lineage.5.genus <- Dat.avm.Global$Lineage.5.Genus
d_lineage.5.name <- Dat.avm.Global$Lineage.5.Name
d_lineage.5.access.no <- Dat.avm.Global$Lineage.5.Accession..
d_lineage.5.cytB <- Dat.avm.Global$Lineage.5.Cyt.B
d_Screen.Leuco <- Dat.avm.Global$Screen.Leuco
d_Screen.Leuco <- ifelse((d_Screen.Leuco=='Yes' | d_Screen.Leuco=='YES'),1,0) 

d_infect_PL.1 <- ifelse(d_lineage.1.genus=="PL", 1,0)
d_infect_PL.2 <- ifelse(d_lineage.2.genus=="PL", 1,0)
d_infect_PL.3 <- ifelse(d_lineage.3.genus=="PL", 1,0)
d_infect_PL.4 <- ifelse(d_lineage.4.genus=="PL", 1,0)
d_infect_PL.5 <- ifelse(d_lineage.5.genus=="PL", 1,0)
d_infect_PL.N <- d_infect_PL.1 + d_infect_PL.2 + d_infect_PL.3 + d_infect_PL.4 + d_infect_PL.5
d_infect_PL.PA <- ifelse(d_infect_PL.N>0,1,0)

d_infect_PA.1 <- ifelse(d_lineage.1.genus=="PA", 1,0)
d_infect_PA.2 <- ifelse(d_lineage.2.genus=="PA", 1,0)
d_infect_PA.3 <- ifelse(d_lineage.3.genus=="PA", 1,0)
d_infect_PA.4 <- ifelse(d_lineage.4.genus=="PA", 1,0)
d_infect_PA.5 <- ifelse(d_lineage.5.genus=="PA", 1,0)
d_infect_PA.N <- d_infect_PA.1 + d_infect_PA.2 + d_infect_PA.3 + d_infect_PA.4 + d_infect_PA.5
d_infect_PA.PA <- ifelse(d_infect_PA.N>0,1,0)

d_infect_HA.1 <- ifelse(d_lineage.1.genus=="HA", 1,0)
d_infect_HA.2 <- ifelse(d_lineage.2.genus=="HA", 1,0)
d_infect_HA.3 <- ifelse(d_lineage.3.genus=="HA", 1,0)
d_infect_HA.4 <- ifelse(d_lineage.4.genus=="HA", 1,0)
d_infect_HA.5 <- ifelse(d_lineage.5.genus=="HA", 1,0)
d_infect_HA.N <- d_infect_HA.1 + d_infect_HA.2 + d_infect_HA.3 + d_infect_HA.4 + d_infect_HA.5
d_infect_HA.PA <- ifelse(d_infect_HA.N>0,1,0)

d_infect_LE.1 <- ifelse(d_lineage.1.genus=="LE", 1,0)
d_infect_LE.1[which(d_Screen.Leuco==0)] <- NA
d_infect_LE.2 <- ifelse(d_lineage.2.genus=="LE", 1,0)
d_infect_LE.2[which(d_Screen.Leuco==0)] <- NA
d_infect_LE.3 <- ifelse(d_lineage.3.genus=="LE", 1,0)
d_infect_LE.3[which(d_Screen.Leuco==0)] <- NA
d_infect_LE.4 <- ifelse(d_lineage.4.genus=="LE", 1,0)
d_infect_LE.4[which(d_Screen.Leuco==0)] <- NA
d_infect_LE.5 <- ifelse(d_lineage.5.genus=="LE", 1,0)
d_infect_LE.5[which(d_Screen.Leuco==0)] <- NA
d_infect_LE.N <- d_infect_LE.1 + d_infect_LE.2 + d_infect_LE.3 + d_infect_LE.4 + d_infect_LE.5
d_infect_LE.N[which(d_Screen.Leuco==0)] <- NA
d_infect_LE.PA <- ifelse(d_infect_LE.N>0,1,0)

Dat.AVM.All <- data.frame(Host = c(Dat.avm.Global$Host.Latin.Name, Dat.avm.AusAsia$Host))
Dat.AVM.All$Latitude <- c(Dat.avm.Global$Latitude, Dat.avm.AusAsia$Latitude)
Dat.AVM.All$Longitude <- c(Dat.avm.Global$Longitude , Dat.avm.AusAsia$Longitude)

Dat.AVM.All$Plas_pa <- c(d_infect_PL.PA, ifelse(Dat.avm.AusAsia$Plas==1,1,0))
Dat.AVM.All$ParaHaem_pa <- c(d_infect_PA.PA, ifelse(Dat.avm.AusAsia$ParaHaem==1,1,0))
Dat.AVM.All$Haem_pa <- c(d_infect_HA.PA, ifelse(Dat.avm.AusAsia$ParaHaem==1,1,0))
Dat.AVM.All$Leuc_pa <- c(d_infect_LE.PA, rep(0, dim(Dat.avm.AusAsia)[1]))

Dat.AVM.All$Plas_nsp <- c(d_infect_PL.N , rep(0, dim(Dat.avm.AusAsia)[1]))
Dat.AVM.All$ParaHaem_nsp <- c(d_infect_PA.N , rep(0, dim(Dat.avm.AusAsia)[1]))
Dat.AVM.All$Haem_nsp <- c(d_infect_HA.N , rep(0, dim(Dat.avm.AusAsia)[1]))
Dat.AVM.All$Leuc_nsp <- c(d_infect_LE.N, rep(0, dim(Dat.avm.AusAsia)[1]))  

Dat.AVM.All$Location <- c(Dat.avm.Global$Latitude, Dat.avm.AusAsia$Location)
Dat.AVM.All$Year <- c(Dat.avm.Global$Year, rep(NA, dim(Dat.avm.AusAsia)[1]))

setwd(dir.data)
write.csv(Dat.AVM.All, file= "Data.AVM.Global_Host.Infection_200717.csv")




###############
# Collate environmental data for location data

# Number of unique coordinates and vector for spatial points
NPoints <- dim(unique(Dat.AVM.All[, c("Longitude", "Latitude")]))[1]
pts_long <- unique(Dat.AVM.All[, c("Longitude", "Latitude")])$Longitude
pts_lat <- unique(Dat.AVM.All[, c("Longitude", "Latitude")])$Latitude
pts_pointID <- 1:NPoints 

# Spatial points
sp.pts <- sf::st_multipoint(matrix(c(pts_long,pts_lat), length(pts_long), 2))
st_crs(sp.pts) <- WGS84 

# Buffer
buff10.pts <- st_buffer(sp.pts, d=0.1)
##buff10.pts <- circles(sp.pts, d=10000, lonlat=TRUE)

PTS <- data.frame(ID.pts = 1:NPoints, long = pts_long, lat = pts_lat)


######
# Worldclim data 

rast_worldclim <- raster::getData("worldclim",var="bio",res=10)
rast_worldclim <- rast_worldclim[[c("bio1", "bio12", "bio14", "bio15")]] 

sp.pts <- SpatialPoints(data.frame(pts_long, pts_lat), proj4string = rast_worldclim@crs)
worldclim.values <- extract(rast_worldclim,sp.pts)   #st_extract(rast_worldclim,sp.pts)    #

PTS$bio1 <- worldclim.values[, "bio1"]
PTS$bio12 <- worldclim.values[, "bio12"]
PTS$bio14 <- worldclim.values[, "bio14"] 
PTS$bio15 <- worldclim.values[, "bio15"] 

######
# Altitude data


PTS$elevation <- rep(NA, NPoints)

elevation_world <- getData('worldclim', var='alt', res=2.5)

PTS$elevation <- extract(elevation_world,sp.pts)
PTS$elevation <- ifelse(PTS$elevation<1,1,PTS$elevation)



######
# Copernicus landcover data

##LCCS codes found here:  https://maps.elie.ucl.ac.be/CCI/viewer/download/CCI-LC_Maps_Legend.pdf

setwd(paste(Drive,  ":/COPERNICUS_Landcover/", sep=""))
Copernicus <- raster("ESACCI-LC-L4-LCCS-Map-300m-P1Y-2010-v2.0.7cds.nc")

coord <- PTS
coordinates(coord) <- ~long + lat
proj4string(coord) <- crs(Copernicus)

coord_vals <- raster::extract(Copernicus, coord, buffer = 10000)

PTS$prop_cropland <- rep(NA, NPoints)
PTS$prop_water <- rep(NA, NPoints)
PTS$prop_urban <- rep(NA, NPoints)
PTS$prop_grassland <- rep(NA, NPoints)
PTS$prop_snow_ice <- rep(NA, NPoints)
PTS$prop_treecover <- rep(NA, NPoints)

for(p in 1:NPoints){
	coord_vals <- raster::extract(Copernicus, coord[p,], buffer = 10000)
	coord_vals <- unlist(coord_vals) 
	PTS$prop_cropland[p] <- length(which(coord_vals %in% c(10, 20, 30))) / length(coord_vals)
	PTS$prop_water[p] <- length(which(coord_vals == 210)) / length(coord_vals)
	PTS$prop_urban[p] <- length(which(coord_vals == 190)) / length(coord_vals)
	PTS$prop_grassland[p] <- length(which(coord_vals == 130)) / length(coord_vals)
	PTS$prop_snow_ice[p] <- length(which(coord_vals == 220)) / length(coord_vals)
	PTS$prop_shrubland[p] <- length(which(coord_vals %in% c(120, 121, 122))) / length(coord_vals)
	PTS$prop_treecover[p] <- length(which(coord_vals %in% c(40, 50, 60, 61, 62, 71, 72, 80, 81, 82, 90))) / length(coord_vals)
}



###############
# NDVI data from MODIS

library(MODISTools)

PTS$NDVI.mean <- rep(NA< NPoints)
PTS$NDVI.SD <- rep(NA< NPoints)


for(x in 1:(NPoints*10000)){
	i <- which(is.na(PTS$NDVI.mean))[1]
	if(!is.na(i)){
	tryCatch(point_ndvi <- mt_subset(product = "MOD13Q1",
                    lat = PTS$lat[i],
                    lon =  PTS$long[i],
                    band = "250m_16_days_NDVI",
                    start = "2010-01-01",
                    end = "2010-12-30",
                    km_lr = 10,
                    km_ab = 10,
                    site_name = "i",
                    internal = TRUE,
                    progress = FALSE), finally = print("Hello"), error = function(e) { skip_to_next <<- TRUE})
	if(exists("point_ndvi")){
		PTS$NDVI.mean[i] <- round(mean(point_ndvi$value))
		PTS$NDVI.SD[i] <- round(sd(point_ndvi$value))
		print(paste("datpoint:", i, "   loop:", x))
	}else{}
	}else{}
	rm(point_ndvi)
	if(i==NPoints){break}else{}
}

###############
# Zoogeographic region from Holt et al.


## Use sf
pts.multi_all <- sf::st_multipoint(as.matrix(data.frame(PTS$long,PTS$lat)))
pts <- st_cast(st_sfc(pts.multi_all), "POINT")
st_crs(pts) <- st_crs(poly_Realms)

PTS$Realm <- rep(NA, NPoints)
for(r in 1:nrealm){
	polyrealm <- poly_Realms %>% filter (Realm==realm_name[r])
	PTS$Realm[which(st_intersects(pts, polyrealm, sparse = FALSE))] <- realm_name[r]
}


ggplot(poly_Realms) + geom_sf(aes(fill=Realm)) + geom_point(data=PTS, aes(x=long, y=lat))



## Safe location-specific covariate file
setwd(dir.data)
write.csv(PTS, "Data.AVM.Global_Location.Covariates_200717.csv")



###############
# NDVI data

##GIMMS NDVI dataset (https://climatedataguide.ucar.edu/climate-data/ndvi-normalized-difference-vegetation-index-3rd-generation-nasagfsc-gimms) for the year 2010

library(gimms)

library(ncdf4)
library(rgdal) 

# Download data
NDVI <- gimms::downloadGimms(x = as.Date("2010-01-01"), y = as.Date("2010-12-31"))

# Check for available files
gimms_file <- rearrangeFiles(dsn = paste0(getwd()), 
                              pattern = "ndvi3g_geo_v1_2010", full.names = TRUE)
gimms_file

## rasterize files
gimms_raster <- rasterizeGimms(gimms_file, 
                               filename = paste0(gimms_file, ".tif"))







#####################
# Host trait data

# EltonTraits v1.0 (Wilman et al. 2014), https://figshare.com/articles/Data_Paper_Data_Paper/3559887

temp <- tempfile()
download.file("https://ndownloader.figshare.com/files/5631081", temp)

Elton.dat <- read.table(temp, header = TRUE,
fill  = TRUE, quote = "\"", stringsAsFactors = FALSE,sep = "\t")
unlink(temp)


HostTraits <- data.frame(Species = sort(unique(Dat.AVM.All$Host)))
NHost <- length(HostTraits$Species)

sel.row <- match(HostTraits$Species, Elton.dat$Scientific)

HostTraits$PassNonPass <- Elton.dat$PassNonPass [sel.row]
HostTraits$IOCOrder <- Elton.dat$IOCOrder [sel.row]
HostTraits$BLFamilyLatin <- Elton.dat$BLFamilyLatin [sel.row]
HostTraits$BLFamilyEnglish <- Elton.dat$BLFamilyEnglish [sel.row]

HostTraits$Bodymass <- Elton.dat$BodyMass.Value [sel.row]
HostTraits$ForStrat.watbelowsurf <- Elton.dat$ForStrat.watbelowsurf [sel.row]
HostTraits$ForStrat.wataroundsurf <- Elton.dat$ForStrat.wataroundsurf [sel.row] 
HostTraits$ForStrat.ground <- Elton.dat$ForStrat.ground [sel.row]
HostTraits$ForStrat.understory <- Elton.dat$ForStrat.understory [sel.row]
HostTraits$ForStrat.midhigh <- Elton.dat$ForStrat.midhigh [sel.row]
HostTraits$ForStrat.canopy <- Elton.dat$ForStrat.canopy [sel.row]
HostTraits$ForStrat.aerial <- Elton.dat$ForStrat.aerial [sel.row]
HostTraits$PelagicSpecialist <- Elton.dat$PelagicSpecialist [sel.row]

HostTraits$Diet.5Cat  <- Elton.dat$Diet.5Cat [sel.row]
HostTraits$Diet.Inv  <- Elton.dat$Diet.Inv [sel.row]
HostTraits$Diet.Vend  <- Elton.dat$Diet.Vend [sel.row]
HostTraits$Diet.Vect <- Elton.dat$Diet.Vect [sel.row]
HostTraits$Diet.Vfish <- Elton.dat$Diet.Vfish [sel.row]
HostTraits$Diet.Vunk <- Elton.dat$Diet.Vunk [sel.row]
HostTraits$Diet.Scav <- Elton.dat$Diet.Scav [sel.row]
HostTraits$Diet.Fruit <- Elton.dat$Diet.Fruit [sel.row]
HostTraits$Diet.Nect <- Elton.dat$Diet.Nect [sel.row]
HostTraits$Diet.Seed <- Elton.dat$Diet.Seed [sel.row]
HostTraits$Diet.PlantO <- Elton.dat$Diet.PlantO [sel.row]


setwd(dir.data)
write.csv(HostTraits, file="Data.AVM.Global_Host.Traits_200717.csv")





