

library(tidyverse)
library(ggplot2)
library(raster)
library(sp)
library(rgeos)
library(sf)



# Avian Malaria records
DatAVM <- read.csv("AvMal_Infection.csv", header=T)

PTS <-  as_tibble(DatAVM) %>% select (c("Latitude", "Longitude", "Country")) 
PTS <-  distinct(PTS)
PTS <-  PTS %>% rename(lat=Latitude, long =Longitude, country =Country)

# ID for each point/location
PTS$pts_id <- paste0("pt.", 1: dim(PTS)[1])

# Number of points
NPoints <- dim(PTS)[1]


######
# Worldclim data 

rast_worldclim <- raster::getData("worldclim",var="bio",res=10)
rast_worldclim <- rast_worldclim[[c("bio1", "bio12", "bio14", "bio15")]] 

sp.pts <- SpatialPoints(data.frame(PTS$long, PTS$lat), proj4string = rast_worldclim@crs)
worldclim.values <- raster::extract(rast_worldclim, sp.pts)   

PTS$bio1 <- worldclim.values[, "bio1"]
PTS$bio12 <- worldclim.values[, "bio12"]
PTS$bio14 <- worldclim.values[, "bio14"] 
PTS$bio15 <- worldclim.values[, "bio15"] 


######
# Elevation data

PTS$elevation <- rep(NA, NPoints)
elevation_world <- getData('worldclim', var='alt', res=2.5)
PTS$elevation <- raster::extract(elevation_world, sp.pts)
PTS$elevation <- ifelse(PTS$elevation<1,1,PTS$elevation)

######
# Copernicus landcover data

## LCCS data and codes can be downloaded here:  
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-land-cover?tab=form  
# https://maps.elie.ucl.ac.be/CCI/viewer/download/CCI-LC_Maps_Legend.pdf

Copernicus <- raster("ESACCI-LC-L4-LCCS-Map-300m-P1Y-2010-v2.0.7cds.nc")

coord <- PTS
coordinates(coord) <- ~long + lat
proj4string(coord) <- crs(Copernicus)

coord_vals <- raster::extract(Copernicus, coord, buffer = 10000)

PTS$prop_cropland <- rep(NA, NPoints)
PTS$prop_water <- rep(NA, NPoints)
PTS$prop_urban <- rep(NA, NPoints)
PTS$prop_grassland <- rep(NA, NPoints)
PTS$prop_treecover <- rep(NA, NPoints)
PTS$prop_shrubland <- rep(NA, NPoints)
for(p in 1:NPoints){
	coord_vals <- raster::extract(Copernicus, coord[p,], buffer = 10000)
	coord_vals <- unlist(coord_vals) 
	PTS$prop_cropland[p] <- length(which(coord_vals %in% c(10, 20, 30))) / length(coord_vals)
	PTS$prop_water[p] <- length(which(coord_vals == 210)) / length(coord_vals)
	PTS$prop_urban[p] <- length(which(coord_vals == 190)) / length(coord_vals)
	PTS$prop_grassland[p] <- length(which(coord_vals == 130)) / length(coord_vals)
	PTS$prop_shrubland[p] <- length(which(coord_vals %in% c(120, 121, 122))) / length(coord_vals)
	PTS$prop_treecover[p] <- length(which(coord_vals %in% c(40, 50, 60, 61, 62, 71, 72, 80, 81, 82, 90))) / length(coord_vals)
}

###############
# NDVI data from MODIS

library(MODISTools)

PTS$NDVI.mean <- rep(NA< NPoints)
PTS$NDVI.SD <- rep(NA< NPoints)

# Run multiple loops if download interrupts
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
# Zoogeographic region from Holt et al. 2013

#Map of zoogeographical regions/realms can be downloaded here: https://macroecology.ku.dk/resources/wallace

poly_Realms <- st_read("newRealms.shp")

pts.multi_all <- sf::st_multipoint(as.matrix(data.frame(PTS$long,PTS$lat)))
pts <- st_cast(st_sfc(pts.multi_all), "POINT")
st_crs(pts) <- st_crs(poly_Realms)

# Correct wrong spelling in realm name "Oceania"
poly_Realms$Realm[which(poly_Realms$Realm=="Oceanina")] <- "Oceania"
realm_name <- unique(poly_Realms$Realm)
nrealm <- length(realm_name)

PTS$realm <- rep(NA, NPoints)
for(r in 1:nrealm){
	polyrealm <- poly_Realms %>% filter (Realm==realm_name[r])
	PTS$realm[which(st_intersects(pts, polyrealm, sparse = FALSE))] <- realm_name[r]
}

# Complete missing values
PTS$realm[which(is.na(PTS$realm) & PTS$country=="UK")] <- "Palearctic"
PTS$realm[which(is.na(PTS$realm) & PTS$country=="Greece")] <- "Palearctic"
PTS$realm[which(is.na(PTS$realm) & PTS$country=="Spain")] <- "Palearctic"
PTS$realm[which(is.na(PTS$realm) & PTS$country=="Italy")] <- "Palearctic"
PTS$realm[which(is.na(PTS$realm) & PTS$country=="Portugal")] <- "Palearctic"
PTS$realm[which(is.na(PTS$realm) & PTS$country=="Malta")] <- "Palearctic"
PTS$realm[which(is.na(PTS$realm) & PTS$country=="Japan")] <- "Palearctic"
PTS$realm[which(is.na(PTS$realm) & PTS$country=="Chile")] <- "Neotropical"
PTS$realm[which(is.na(PTS$realm) & PTS$country=="Brazil")] <- "Neotropical"
PTS$realm[which(is.na(PTS$realm) & PTS$country=="New Guinea")] <- "Oriental"
PTS$realm[which(is.na(PTS$realm) & PTS$country=="Antarctic ")] <- "Palearctic"
PTS$realm[which(is.na(PTS$realm) & PTS$country=="Falkland Islands")] <- "Neotropical"
PTS$realm[which(is.na(PTS$realm) & PTS$long<(-61.))] <- "Neotropical"
PTS$realm[which(is.na(PTS$realm) & PTS$long>110 & PTS$lat>0)] <- "Oriental"
PTS$realm[which(is.na(PTS$realm) & PTS$long>110 & PTS$lat>(-12))] <- "Oceanina"
PTS$realm[which(is.na(PTS$realm) & PTS$long>110)] <- "Australian"

ggplot(poly_Realms) + geom_sf(aes(fill=Realm)) + geom_point(data=PTS[which(is.na(PTS$realm)),], aes(x=long, y=lat))


write.csv(PTS, file = "AvMal_Locations.csv")


#####################
# Host trait data

# EltonTraits v1.0 (Wilman et al. (2014), https://figshare.com/articles/Data_Paper_Data_Paper/3559887)
temp <- tempfile()
download.file("https://ndownloader.figshare.com/files/5631081", temp)

Elton.dat <- read.table(temp, header = TRUE,
fill  = TRUE, quote = "\"", stringsAsFactors = FALSE,sep = "\t")
unlink(temp)
Elton.dat <- as.tibble(Elton.dat) 
Elton.dat <- rename(Elton.dat, Species =Scientific)

# Migration distance data from Dufour et al. 2020, https://onlinelibrary.wiley.com/doi/full/10.1111/jbi.13700
DatMigr  <- as.tibble(read.csv("Dufor_migration.data_210315.csv", header=T))
DatMigr$Species <- gsub("_", " ", DatMigr$Sp.Scien.jetz)
DatMigr  <- DatMigr  %>% select(Species, strategy_3, distance_4, distance_quanti_RES0)
DatMigr  <- rename(DatMigr, migrate.strategy.3=strategy_3, migrate.strategy.4=distance_4, migrate.distance=distance_quanti_RES0)

HostTraits <- as.tibble(data.frame(Species = sort(unique(DatAVM$Avian.host.species))))
HostTraits <- HostTraits %>% left_join(Elton.dat) %>% left_join(DatMigr)


write.csv(HostTraits, file = "AvMal_HostTrait.csv")






