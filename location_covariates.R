
library(ggplot2)
library(raster)
library(sp)
library(rgeos)
library(sf)


# WGS84 definition
WGS84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# Read location data
PTS <- read.csv("Data.AVM.Global_Location.Covariates.csv", header=T)

# Number of points
NPoints <- dim(PTS)[1]
# Spatial points
sp.pts <- sf::st_multipoint(matrix(c(PTS$long,PTS$lat), length(PTS$long), 2))
# Buffer around sample locations
buff10.pts <- st_buffer(sp.pts, d=0.1)


######
# Worldclim data 

rast_worldclim <- raster::getData("worldclim",var="bio",res=10)
rast_worldclim <- rast_worldclim[[c("bio1", "bio12", "bio14", "bio15")]] 

sp.pts <- SpatialPoints(data.frame(pts_long, pts_lat), proj4string = rast_worldclim@crs)
worldclim.values <- extract(rast_worldclim,sp.pts)   

PTS$bio1 <- worldclim.values[, "bio1"]
PTS$bio12 <- worldclim.values[, "bio12"]
PTS$bio14 <- worldclim.values[, "bio14"] 
PTS$bio15 <- worldclim.values[, "bio15"] 

######
# Elevation data

PTS$elevation <- rep(NA, NPoints)

elevation_world <- getData('worldclim', var='alt', res=2.5)

PTS$elevation <- extract(elevation_world,sp.pts)
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

# Run multiple loops if donload interrupts
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

#Map of zoogeographical regions/realms can be downloaded here: https://macroecology.ku.dk/resources/wallace

poly_Realms <- st_read("./zoo.regions_Holt2013/newRealms.shp")

pts.multi_all <- sf::st_multipoint(as.matrix(data.frame(PTS$long,PTS$lat)))
pts <- st_cast(st_sfc(pts.multi_all), "POINT")
st_crs(pts) <- st_crs(poly_Realms)

PTS$Realm <- rep(NA, NPoints)
for(r in 1:nrealm){
	polyrealm <- poly_Realms %>% filter (Realm==realm_name[r])
	PTS$Realm[which(st_intersects(pts, polyrealm, sparse = FALSE))] <- realm_name[r]
}

ggplot(poly_Realms) + geom_sf(aes(fill=Realm)) + geom_point(data=PTS, aes(x=long, y=lat))




#####################
# Host trait data

# EltonTraits v1.0 (Wilman et al. (2014), https://figshare.com/articles/Data_Paper_Data_Paper/3559887)

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






