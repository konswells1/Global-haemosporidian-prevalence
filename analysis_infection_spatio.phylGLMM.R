


##########
# Load R packages

library(Hmisc)
library(raster)
library(ape)
library(tidyverse)
library(dplyr)
library(ape)
library(sp)
library(sf)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggtree)
library(RColorBrewer)
library(viridis)
library(wesanderson)
library(INLA)

## Inverse logit function
invlogit <- function(x){ exp(x) / (1 + exp(x))}
memory.limit(size=10^13)


#############
# Load and format data

# World map
library(rnaturalearth)
worldmap <- ne_countries(scale = "medium", type = "countries", returnclass = "sf")[1]

setwd(dir.data)
# Avian Malaria records
DatAVM <- as_tibble(read.csv("AvMal_Infection.csv", header=T))
DatAVM <- rename(DatAVM, Species= Avian.host.species)
DatAVM <- rename(DatAVM, Haem_pa= Haemoproteus, Plas_pa=Plasmodium, ParaHaem_pa=Parahaemoproteus, Leuc_pa=Leucocytozoon)
DatAVM <- DatAVM %>% select(Species, Species, Latitude, Longitude, Year, Isolation.source, Protocol, Data.source, Haem_pa, Plas_pa, ParaHaem_pa, Leuc_pa)            

# Location-specific covariates
DatEnv <- as.tibble(read.csv("AvMal_Locations.csv", header=T))

# Host trait data
DatHost <- as.tibble(read.csv("AvMal_HostTrait.csv", header=T))

# Phylogenetic tree birds
setwd(dir.data0)
BirdTree100 <- read.tree("BirdTree100.tre")


# Matching of different datasets 
match_host2avm <- match(DatAVM$Species, DatHost$Species)
match_env2avm <- match(paste(DatAVM$Latitude, DatAVM$Longitude), paste(DatEnv$lat, DatEnv$long))

# Merge all data into single file
DatMerge <- data.frame(
	DatAVM,
	DatEnv[match_env2avm,],
	DatHost[match_host2avm,])


#Scaled covariates
DatMerge$abslat_sc <- as.numeric(scale(abs(DatMerge$Latitude))) 
DatMerge$bio1_sc <- as.numeric(scale(DatMerge$bio1)) 
DatMerge$bio12_sc <- as.numeric(scale(DatMerge$bio12)) 
DatMerge$bio14_sc <- as.numeric(scale(DatMerge$bio14)) 
DatMerge$bio15_sc <- as.numeric(scale(DatMerge$bio15)) 
DatMerge$elevation_sc <- as.numeric(scale(DatMerge$elevation)) 
DatMerge$prop_cropland_sc <- as.numeric(scale(DatMerge$prop_cropland)) 
DatMerge$prop_water_sc <- as.numeric(scale(DatMerge$prop_water)) 
DatMerge$prop_urban_sc <- as.numeric(scale(DatMerge$prop_urban))  
DatMerge$prop_grassland_sc <-  as.numeric(scale(DatMerge$prop_grassland)) 
DatMerge$prop_treecover_sc <-  as.numeric(scale(DatMerge$prop_treecover))
DatMerge$prop_shrubland_sc <-  as.numeric(scale(DatMerge$prop_shrubland))
DatMerge$NDVI.mean_sc <-  as.numeric(scale(DatMerge$NDVI.mean))
DatMerge$NDVI.SD_sc <-  as.numeric(scale(DatMerge$NDVI.SD))
DatMerge$Bodymass_log10.sc <-  as.numeric(scale(log10(DatMerge$BodyMass.Value)))
DatMerge$ForStrat.ground_sc <-  as.numeric(scale(DatMerge$ForStrat.ground))
DatMerge$ForStrat.understory_sc <-  as.numeric(scale(DatMerge$ForStrat.understory))
DatMerge$ForStrat.midhigh_sc <-  as.numeric(scale(DatMerge$ForStrat.midhigh))
DatMerge$ForStrat.canopy_sc <-  as.numeric(scale(DatMerge$ForStrat.canopy))
DatMerge$ForStrat.aerial_sc <-  as.numeric(scale(DatMerge$ForStrat.aerial))
DatMerge$Diet.Inv_sc <-  as.numeric(scale(DatMerge$Diet.Inv))
DatMerge$Diet.Vend_sc <-  as.numeric(scale(DatMerge$Diet.Vend))
DatMerge$Diet.Vect_sc <-  as.numeric(scale(DatMerge$Diet.Vect))
DatMerge$Diet.Vfish_sc <-  as.numeric(scale(DatMerge$Diet.Vfish))
DatMerge$Diet.Vunk_sc <-  as.numeric(scale(DatMerge$Diet.Vunk))
DatMerge$Diet.Scav_sc <-  as.numeric(scale(DatMerge$Diet.Scav))
DatMerge$Diet.Fruit_sc <-  as.numeric(scale(DatMerge$Diet.Fruit))
DatMerge$Diet.Nect_sc <-  as.numeric(scale(DatMerge$Diet.Nect))
DatMerge$bird.richness_sc <-  as.numeric(scale(DatMerge$Diet.Nect))
DatMerge$migrate.distance_sc <-  as.numeric(scale(DatMerge$migrate.distance))
  
host_spname <- sort(unique(DatMerge$Species))
host_family <- DatMerge$BLFamilyLatin[match(host_spname, DatMerge$Species)]
nhost <- length(host_spname)
host_id <- 1:nhost

family_name <- sort(unique(host_family))
nfamily <- length(family_name)
family_id <- 1:nfamily

DatMerge$family_id <- family_id[match(DatMerge$BLFamilyLatin, family_name)]

DatMerge$realm[which(DatMerge$realm=="Oceanina")] <- "Oceania"
DatMerge$realm_id <- as.numeric(as.factor(DatMerge$realm))

realm_name <- sort(unique((DatMerge$realm)))
nrealm <- length(realm_name)


# Subset of data for Leucocytozoon : select only locations/host screened
DatMergeLeuc <- DatMerge[which(!is.na(DatMerge$Leuc_pa)),]

# Subset of data for Haemoproteus: select only families known to be infected
tble_family.Haem <- table(DatMerge$BLFamilyLatin, DatMerge$Haem_pa)
sel.DatMerge_Haem <- which(!is.na(match(DatMerge$BLFamilyLatin, names(which(tble_family.Haem[,2]>1)))))
DatMergeHaem <- DatMerge[sel.DatMerge_Haem,]
tble_family.Haem[which(tble_family.Haem[,2]>1),]


# Vectors of point coordinates
pts_id <- DatEnv$pts_id
pts_lat <- DatEnv$lat
pts_long <- DatEnv$long
NPoints <- dim(DatEnv)[1]
pts_realm <- DatEnv$realm

pts_nsample.Plas <- as.numeric(table(DatMerge$pts_id)[match(pts_id, names(table(DatMerge$pts_id)))])
pts_nsample.Haem <- as.numeric( table(DatMergeHaem$pts_id)[match(pts_id, names(table(DatMergeHaem$pts_id)))])
pts_nsample.Leuc <- as.numeric(table(DatMergeLeuc$pts_id)[match(pts_id, names(table(DatMergeLeuc$pts_id)))])

rm(DatAVM); rm(DatEnv); rm(DatHost)

ggplot(worldmap) + geom_sf() +
geom_point(data = DatMerge, aes(x = long, y = lat), size = 2, shape = 19, fill = "darkred")


########
# Compute family-level phylogenetic distance matrices

# Array to store 100 species-level phylogenetic precion matrices
ntree <- 10  ##ntree <- 100 Many replicates do not change family-level phylognetic tree!
phyl.vcv <- phyl.vcv.corr <-  array(NA, dim=c(ntree, nhost, nhost))

for(x in 1:ntree){
	birdtree.sel <- BirdTree100[[x]]
	birdtree.sel$tip.label <- sub("_"," ", birdtree.sel$tip.label)

	drop.tip_host <- which(is.na(match(birdtree.sel$tip.label, host_spname)))
	phyl.tree <- ape::drop.tip(birdtree.sel, drop.tip_host)
	phyl.vcv.corr[x,,] <- ape::vcv.phylo(phyl.tree, model = "Brownian", cor=TRUE)
	phyl.vcv[x,,] <- ape::vcv(phyl.tree, model = "Brownian")
}

phyl.dist_family.ntree <- phyl.vcv_family.ntree <-
phyl.prec_family.ntree <- array(NA, dim=c(ntree, nfamily, nfamily))

for(x in 1:ntree){
	phyl.vcv.corr_family <- matrix(0, nrow=nfamily, ncol=nfamily)
	for(i in 1:nfamily){
		for(j in 1:nfamily){
			sel_fam.i <- which(host_family==family_name[i]) 
			sel_fam.j <- which(host_family==family_name[j]) 
			phyl.vcv_family.ntree[x,i,j] <- mean(phyl.vcv[x, sel_fam.i, sel_fam.j])
			phyl.vcv.corr_family[i,j] <- mean(phyl.vcv.corr[x, sel_fam.i, sel_fam.j]) 
		}
	}
	diag(phyl.vcv.corr_family) <- 1
	phyl.dist_family.ntree[x,,] <- (1 - phyl.vcv.corr_family)
}

phyl.dist_family <- phyl.vcv_family <- matrix(0, nrow=nfamily, ncol=nfamily) 
for(i in 1:nfamily){
	for(j in 1:nfamily){
		phyl.dist_family[[i,j]] <- mean(phyl.dist_family.ntree[,i,j])
		phyl.vcv_family[[i,j]] <- mean(phyl.vcv_family.ntree[,i,j])
	}
}

# Compute precision matrix
phyl.covar_family <- phyl.vcv_family/max(phyl.vcv_family)
phyl.covar_family <- phyl.covar_family / exp(determinant(phyl.covar_family)$modulus[1] /nrow(phyl.covar_family))
phyl.prec_family <- solve(phyl.covar_family)

# Family-level tree
library(phylogram)
F_as.dendrogram <- as.dendrogram(hclust(dist(phyl.dist_family)))
tree_host.family <- as.phylo.dendrogram(F_as.dendrogram)
tree_host.family$tip.label <- family_name 

plot(tree_host.family)

rm(BirdTree100)

##########################3
#
#  R-INLA models
#
############


#######
# Overall infection probabilities

model_inf0_Plas <- inla(Plas_pa ~ 1,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE))
model_inf0_ParaHaem <- inla(ParaHaem_pa ~ 1,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE))
model_inf0_Haem <- inla(Haem_pa ~ 1,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE))
model_inf0_Leuc <- inla(Leuc_pa ~ 1,family="binomial",data=DatMergeLeuc, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE))

#######
# Point-level infection probabilities
f_pts_Plas <- as.formula("Plas_pa ~ -1 + f(pts_id, model='iid')")
f_pts_ParaHaem <- as.formula("ParaHaem_pa ~ -1 + f(pts_id, model='iid')")
f_pts_Haem <- as.formula("Haem_pa ~ -1 + f(pts_id, model='iid')")
f_pts_Leuc <- as.formula("Leuc_pa ~ -1 + f(pts_id, model='iid')")

mod_pts_Plas <- inla(f_pts_Plas,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE))
mod_pts_ParaHaem <- inla(f_pts_ParaHaem,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE))
mod_pts_Haem <- inla(f_pts_Haem,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE))
mod_pts_Leuc <- inla(f_pts_Leuc,family="binomial",data=DatMergeLeuc, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE))

#######
# Realm-specific infection probabilities
f_realm_Plas <- as.formula("Plas_pa ~ -1 + f(realm_id, model='iid')")
f_realm_ParaHaem <- as.formula("ParaHaem_pa ~ -1 + f(realm_id, model='iid')")
f_realm_Haem <- as.formula("Haem_pa ~ -1 + f(realm_id, model='iid')")
f_realm_Leuc <- as.formula("Leuc_pa ~ -1 + f(realm_id, model='iid')")

mod_realm_Plas <- inla(f_realm_Plas,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE))
mod_realm_ParaHaem <- inla(f_realm_ParaHaem,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE))
mod_realm_Haem <- inla(f_realm_Haem,family="binomial",data=DatMergeHaem, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE))
mod_realm_Leuc <- inla(f_realm_Leuc,family="binomial",data=DatMergeLeuc, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE))


#######
# GLM with fixed effects only

pred_fixeff <- paste("1 + abslat_sc + elevation_sc + ",
			"bio1_sc + bio12_sc + bio14_sc + bio15_sc + ",
			"NDVI.mean_sc + NDVI.SD_sc + prop_treecover_sc + prop_water_sc +",
			"ForStrat.canopy_sc + Bodymass_log10.sc + bird.richness_sc + migrate.distance_sc")
f_fixeff_Plas <- as.formula(paste("Plas_pa ~ ", pred_fixeff))
f_fixeff_ParaHaem <- as.formula(paste("ParaHaem_pa ~ ", pred_fixeff))
f_fixeff_Haem <- as.formula(paste("Haem_pa ~ ", pred_fixeff))
f_fixeff_Leuc <- as.formula(paste("Leuc_pa ~ ", pred_fixeff))

mod_fixeff_Plas <- inla(f_fixeff_Plas,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE, hyper=pc.prec), control.compute=list(dic=TRUE,cpo=TRUE))
mod_fixeff_ParaHaem <- inla(f_fixeff_ParaHaem,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE, hyper=pc.prec), control.compute=list(dic=TRUE,cpo=TRUE))
mod_fixeff_Haem <- inla(f_fixeff_Haem,family="binomial",data=DatMergeHaem, control.predictor=list(compute=TRUE, hyper=pc.prec), control.compute=list(dic=TRUE,cpo=TRUE))
mod_fixeff_Leuc <- inla(f_fixeff_Leuc,family="binomial",data=DatMergeLeuc, control.predictor=list(compute=TRUE, hyper=pc.prec), control.compute=list(dic=TRUE,cpo=TRUE))


#######
# phylo-0-GLMM: phylogeny only as random effects

pred_phyl.0 <- "-1 + f(family_id, model='generic0', Cmatrix=phyl.prec_family)"

f_phyl.0_Plas <- as.formula(paste("Plas_pa ~ ", pred_phyl.0))
f_phyl.0_ParaHaem <- as.formula(paste("ParaHaem_pa ~ ", pred_phyl.0))
f_phyl.0_Haem <- as.formula(paste("Haem_pa ~ ", pred_phyl.0))
f_phyl.0_Leuc <- as.formula(paste("Leuc_pa ~ ", pred_phyl.0))

mod_phyl.0_Plas <- inla(f_phyl.0_Plas,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE), control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"))
mod_phyl.0_ParaHaem <- inla(f_phyl.0_ParaHaem,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE), control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"))
mod_phyl.0_Haem <- inla(f_phyl.0_Haem,family="binomial",data=DatMergeHaem, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE), control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"))
mod_phyl.0_Leuc <- inla(f_phyl.0_Leuc,family="binomial",data=DatMergeLeuc, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE), control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"))


#######
# phylo-cov-GLMM: fixed effects, with phylogeny and year, and isolation.source as random effects but NOT zooregion

pred_phyl.Cov <- paste("-1 + abslat_sc +  elevation_sc + ",
			"bio1_sc + bio12_sc + bio14_sc + bio15_sc + ",
			"NDVI.mean_sc + NDVI.SD_sc + prop_treecover_sc + prop_water_sc +",
			"ForStrat.canopy_sc + Bodymass_log10.sc + bird.richness_sc + migrate.distance_sc + f(family_id, model='generic0', Cmatrix=phyl.prec_family) +",
			"f(Year, model='iid') + f(Isolation.source, model='iid')")

f_phyl.Cov_Plas <- as.formula(paste("Plas_pa ~ ", pred_phyl.Cov))
f_phyl.Cov_ParaHaem <- as.formula(paste("ParaHaem_pa ~ ", pred_phyl.Cov))
f_phyl.Cov_Haem <- as.formula(paste("Haem_pa ~ ", pred_phyl.Cov))
f_phyl.Cov_Leuc <- as.formula(paste("Leuc_pa ~ ",pred_phyl.Cov))

mod_phyl.Cov_Plas <- inla(f_phyl.Cov_Plas,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE, hyper=pc.prec), control.compute=list(dic=TRUE,cpo=TRUE), control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"))

mod_phyl.Cov_ParaHaem <- inla(f_phyl.Cov_ParaHaem,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE, hyper=pc.prec), control.compute=list(dic=TRUE,cpo=TRUE), control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"))
mod_phyl.Cov_Haem <- inla(f_phyl.Cov_Haem,family="binomial",data=DatMergeHaem, control.predictor=list(compute=TRUE, hyper=pc.prec), control.compute=list(dic=TRUE,cpo=TRUE), control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"))
mod_phyl.Cov_Leuc <- inla(f_phyl.Cov_Leuc,family="binomial",data=DatMergeLeuc, control.predictor=list(compute=TRUE, hyper=pc.prec), control.compute=list(dic=TRUE,cpo=TRUE), control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"))


#######
# GLM with phylogeny, year and zooregion as random effects

pred_reg.phyl <- paste("-1 + abslat_sc + elevation_sc + ",
			"bio1_sc + bio12_sc + bio14_sc + bio15_sc + ",
			"NDVI.mean_sc + NDVI.SD_sc + prop_treecover_sc + prop_water_sc +",
			"ForStrat.canopy_sc + Bodymass_log10.sc + bird.richness_sc + migrate.distance_sc + f(family_id, model='generic0', Cmatrix=phyl.prec_family) +",
			"f(realm_id, model='iid') + ",
			"f(Year, model='iid') + f(Isolation.source, model='iid')")

f_reffphyl_Plas <- as.formula(paste("Plas_pa ~ ", pred_reg.phyl))
f_reffphyl_ParaHaem <- as.formula(paste("ParaHaem_pa ~ ", pred_reg.phyl))
f_reffphyl_Haem <- as.formula(paste("Haem_pa ~ ", pred_reg.phyl))
f_reffphyl_Leuc <- as.formula(paste("Leuc_pa ~ ", pred_reg.phyl))

mod_reg.phyl_Plas <- inla(f_reffphyl_Plas,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE, hyper=pc.prec), control.compute=list(dic=TRUE,cpo=TRUE), control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"))
mod_reg.phyl_ParaHaem <- inla(f_reffphyl_ParaHaem,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE, hyper=pc.prec), control.compute=list(dic=TRUE,cpo=TRUE), control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"))
mod_reg.phyl_Haem <- inla(f_reffphyl_Haem,family="binomial",data=DatMerge, control.predictor=list(compute=TRUE, hyper=pc.prec), control.compute=list(dic=TRUE,cpo=TRUE), control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"))
mod_reg.phyl_Leuc <- inla(f_reffphyl_Leuc,family="binomial",data=DatMergeLeuc, control.predictor=list(compute=TRUE, hyper=pc.prec), control.compute=list(dic=TRUE,cpo=TRUE), control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"))


#######
# Model with spatial (SPDE) and phylogenetic random effects

######
# Build mesh for spatial SPDE model

coords <- cbind(DatMerge$long, DatMerge$lat)
# Boundaries
#world.bdry <- inla.sp2segment(worldmap[1])
pts.bdry <- inla.nonconvex.hull(coords, 5, 5, resolution=100)
# Make global mesh
mesh.world <- inla.mesh.2d(coords, boundary = pts.bdry, max.edge = c(1, 3), offset = c(1e-5, 1.5), cutoff = 0.1)

dim(mesh.world$loc)
 
plot(worldmap[1], col = 'grey'); plot(mesh.world, add=T)
plot(mesh.world)
	
plot_mesh <- function(){
	plot(mesh.world, main="")
	#points(DatMerge$long, DatMerge$lat, pch=19, col = brewer.pal(n = nrealm, name = "Set3")[DatMerge$realm_id])
	points(pts_long, pts_lat, pch=19, cex= log(pts_nsample.Plas)/5, col = brewer.pal(n = nrealm, name = "Set3")[as.numeric(as.factor(pts_realm))])
	legend("bottom", legend=realm_name, pch = 19, col = brewer.pal(n = nrealm, name = "Set3"), xpd = TRUE, horiz = TRUE, inset = c(0, -0.05), cex= 1 )
	}
plot_mesh ()

tiff("plot_mesh.tiff", width = 3.307*3.5, height = 3.307*2, units = 'in', res = 600); plot_mesh()
dev.off()

# Compute projector matrix A
A <- inla.spde.make.A(mesh.world, loc = coords)
dim(A)  # dimension of number of locations and number of spatial mesh points

# Matern object
world.spde <- inla.spde2.matern(mesh = mesh.world, alpha = 2)

field.indices <- inla.spde.make.index(name = "field", n.spde = world.spde$n.spde)
field.indices <- inla.spde.make.index(name = "field", n.spde = mesh.world$n)

spde <- inla.spde2.pcmatern(mesh = mesh.world,
  prior.range = c(0.05, 0.01), # P(practic.range < 0.05) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

# Mesh for Haemoproteus
coords_Haem <- cbind(DatMergeHaem$long, DatMergeHaem$lat)
pts.bdry_Haem <- inla.nonconvex.hull(coords_Haem, 5, 5, resolution=100)
# Make global mesh
mesh.world_Haem <- inla.mesh.2d(coords_Haem, boundary = pts.bdry_Haem, max.edge = c(1, 3), offset = c(1e-5, 1.5), cutoff = 0.1)
# Compute projector matrix A
A_Haem <- inla.spde.make.A(mesh.world_Haem, loc = coords_Haem)
# Matern object
world.spde_Haem <- inla.spde2.matern(mesh = mesh.world_Haem, alpha = 2)
spde_Haem <- inla.spde2.pcmatern(mesh = mesh.world_Haem,
  prior.range = c(0.05, 0.01), # P(practic.range < 0.05) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

# Mesh for Leucocytozoon
coords_Leuc <- cbind(DatMergeLeuc$long, DatMergeLeuc$lat)
pts.bdry_Leuc <- inla.nonconvex.hull(coords_Leuc, 5, 5, resolution=100)
# Make global mesh
mesh.world_Leuc <- inla.mesh.2d(coords_Leuc, boundary = pts.bdry_Leuc, max.edge = c(1, 3), offset = c(1e-5, 1.5), cutoff = 0.1)
# Compute projector matrix A
A_Leuc <- inla.spde.make.A(mesh.world_Leuc, loc = coords_Leuc)
# Matern object
world.spde_Leuc <- inla.spde2.matern(mesh = mesh.world_Leuc, alpha = 2)
spde_Leuc <- inla.spde2.pcmatern(mesh = mesh.world_Leuc,
  prior.range = c(0.05, 0.01), # P(practic.range < 0.05) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01


# Create stack for SPDE models
stack_Plas <- inla.stack(
  data = list(y = DatMerge$Plas_pa), 
  A = list(A,1,1,1,1,1,1),   # One A for each effect
  effects = list(
	sp_id = 1:spde$n.spde,
	family_id = DatMerge$family_id,
	year = DatMerge$Year, 
	isolation.source = DatMerge$Isolation.source,
	zooreg = DatMerge$realm_id,
      intercept = rep(1, nrow(DatMerge)), 
      list(abslat = DatMerge$abslat_sc,  
      elevation = DatMerge$elevation_sc,   
	bio1 = DatMerge$bio1_sc,
      bio12 = DatMerge$bio12_sc,
      bio14 = DatMerge$bio14_sc,
      bio15 = DatMerge$bio15_sc,
      ndvi.mean = DatMerge$NDVI.mean_sc,
      ndvi.sd = DatMerge$NDVI.SD_sc,
      LC.forest = DatMerge$prop_treecover_sc,
      LC.wetland = DatMerge$prop_water_sc,
	forstrat.canopy = DatMerge$ForStrat.canopy_sc,
	log.bm = DatMerge$Bodymass_log10.sc,
	bird.richness_sc = DatMerge$bird.richness_sc,
	migrate.distance_sc = DatMerge$migrate.distance_sc,
	r.1 = as.factor(DatMerge$realm),
	r.2 = as.factor(DatMerge$realm),
	r.3 = as.factor(DatMerge$realm),
	r.4 = as.factor(DatMerge$realm),
	r.5 = as.factor(DatMerge$realm),
	r.6 = as.factor(DatMerge$realm),
	r.7 = as.factor(DatMerge$realm),
	r.8 = as.factor(DatMerge$realm),
	r.9 = as.factor(DatMerge$realm),
	r.10 = as.factor(DatMerge$realm),
	r.11 = as.factor(DatMerge$realm),
	r.12 = as.factor(DatMerge$realm),
	r.13 = as.factor(DatMerge$realm))),
  tag = 'Plas_pa') 

stack_ParaHaem <- inla.stack(
  data = list(y = DatMerge$ParaHaem_pa), 
  A = list(A,1,1,1,1,1,1),   # One A for each effect
  effects = list(sp_id = 1:spde$n.spde, 
	family_id = DatMerge$family_id,
	year = DatMerge$Year, 
	isolation.source = DatMerge$Isolation.source,
	zooreg = DatMerge$realm_id,
      intercept = rep(1, nrow(DatMerge)), 
      list(abslat = DatMerge$abslat_sc,  
      elevation = DatMerge$elevation_sc,      
	bio1 = DatMerge$bio1_sc,
      bio12 = DatMerge$bio12_sc,
      bio14 = DatMerge$bio14_sc,
      bio15 = DatMerge$bio15_sc,
      ndvi.mean = DatMerge$NDVI.mean_sc,
      ndvi.sd = DatMerge$NDVI.SD_sc,
      LC.forest = DatMerge$prop_treecover_sc,
      LC.wetland = DatMerge$prop_water_sc,
	forstrat.canopy = DatMerge$ForStrat.canopy_sc,
	log.bm = DatMerge$Bodymass_log10.sc,
	bird.richness_sc = DatMerge$bird.richness_sc,
	migrate.distance_sc = DatMerge$migrate.distance_sc,
	r.1 = as.factor(DatMerge$realm),
	r.2 = as.factor(DatMerge$realm),
	r.3 = as.factor(DatMerge$realm),
	r.4 = as.factor(DatMerge$realm),
	r.5 = as.factor(DatMerge$realm),
	r.6 = as.factor(DatMerge$realm),
	r.7 = as.factor(DatMerge$realm),
	r.8 = as.factor(DatMerge$realm),
	r.9 = as.factor(DatMerge$realm),
	r.10 = as.factor(DatMerge$realm),
	r.11 = as.factor(DatMerge$realm),
	r.12 = as.factor(DatMerge$realm),
	r.13 = as.factor(DatMerge$realm))),
  tag = 'ParaHaem_pa') 

stack_Haem <- inla.stack(
  data = list(y = DatMergeHaem$Haem_pa), 
  A = list(A_Haem,1,1,1,1,1,1),   # One A for each effect
  effects = list(sp_id = 1:spde_Haem$n.spde, 
	family_id = DatMergeHaem$family_id,
	year = DatMergeHaem$Year, 
	isolation.source = DatMerge$Isolation.source,
	zooreg = DatMergeHaem$realm_id,
      intercept = rep(1, nrow(DatMergeHaem)), 
      list(abslat = DatMergeHaem$abslat_sc,  
      elevation = DatMergeHaem$elevation_sc,       
	bio1 = DatMergeHaem$bio1_sc,
      bio12 = DatMergeHaem$bio12_sc,
      bio14 = DatMergeHaem$bio14_sc,
      bio15 = DatMergeHaem$bio15_sc,
      ndvi.mean = DatMergeHaem$NDVI.mean_sc,
      ndvi.sd = DatMergeHaem$NDVI.SD_sc,
      LC.forest = DatMergeHaem$prop_treecover_sc,
      LC.wetland = DatMergeHaem$prop_water_sc,
	forstrat.canopy = DatMergeHaem$ForStrat.canopy_sc,
	log.bm = DatMergeHaem$Bodymass_log10.sc,
	bird.richness_sc = DatMergeHaem$bird.richness_sc,
	migrate.distance_sc = DatMergeHaem$migrate.distance_sc,
	r.1 = as.factor(DatMergeHaem$realm),
	r.2 = as.factor(DatMergeHaem$realm),
	r.3 = as.factor(DatMergeHaem$realm),
	r.4 = as.factor(DatMergeHaem$realm),
	r.5 = as.factor(DatMergeHaem$realm),
	r.6 = as.factor(DatMergeHaem$realm),
	r.7 = as.factor(DatMergeHaem$realm),
	r.8 = as.factor(DatMergeHaem$realm),
	r.9 = as.factor(DatMergeHaem$realm),
	r.10 = as.factor(DatMergeHaem$realm),
	r.11 = as.factor(DatMergeHaem$realm),
	r.12 = as.factor(DatMergeHaem$realm),
	r.13 = as.factor(DatMergeHaem$realm))),
  tag = 'Haem_pa') 
 
stack_Leuc <- inla.stack(
  data = list(y = DatMergeLeuc$Leuc_pa), 
  A = list(A_Leuc,1,1,1,1,1,1),   # One A for each effect
  effects = list(sp_id = 1:spde_Leuc$n.spde, 
	family_id = DatMergeLeuc$family_id,
	year = DatMergeLeuc$Year,
	isolation.source = DatMerge$Isolation.source, 
	zooreg = DatMergeLeuc$realm_id,
      intercept = rep(1, nrow(DatMergeLeuc)), 
      list(abslat = DatMergeLeuc$abslat_sc,  
      elevation = DatMergeLeuc$elevation_sc,        
	bio1 = DatMergeLeuc$bio1_sc,
      bio12 = DatMergeLeuc$bio12_sc,
      bio14 = DatMergeLeuc$bio14_sc,
      bio15 = DatMergeLeuc$bio15_sc,
      ndvi.mean = DatMergeLeuc$NDVI.mean_sc,
      ndvi.sd = DatMergeLeuc$NDVI.SD_sc,
      LC.forest = DatMergeLeuc$prop_treecover_sc,
      LC.wetland = DatMergeLeuc$prop_water_sc,
	forstrat.canopy = DatMergeLeuc$ForStrat.canopy_sc,
	log.bm = DatMergeLeuc$Bodymass_log10.sc,
	bird.richness_sc = DatMergeLeuc$bird.richness_sc,
	migrate.distance_sc = DatMergeLeuc$migrate.distance_sc,
	r.1 = as.factor(DatMergeLeuc$realm),
	r.2 = as.factor(DatMergeLeuc$realm),
	r.3 = as.factor(DatMergeLeuc$realm),
	r.4 = as.factor(DatMergeLeuc$realm),
	r.5 = as.factor(DatMergeLeuc$realm),
	r.6 = as.factor(DatMergeLeuc$realm),
	r.7 = as.factor(DatMergeLeuc$realm),
	r.8 = as.factor(DatMergeLeuc$realm),
	r.9 = as.factor(DatMergeLeuc$realm),
	r.10 = as.factor(DatMergeLeuc$realm),
	r.11 = as.factor(DatMergeLeuc$realm),
	r.12 = as.factor(DatMergeLeuc$realm),
	r.13 = as.factor(DatMergeLeuc$realm))),
  tag = 'Leuc_pa')


f_spde <- y ~ -1 + abslat + elevation + bio1 + bio12 + bio14 + bio15 + 
			ndvi.mean + ndvi.sd + LC.forest + LC.wetland + forstrat.canopy + log.bm + bird.richness_sc + migrate.distance_sc + 
			f(family_id, model='generic0', Cmatrix=phyl.prec_family) + 
			f(sp_id, model = spde) +
			f(zooreg, model='iid') + 
			f(year, model='iid') + 
			f(isolation.source, model='iid')

mod_spde_Plas <- inla(f_spde, family = 'Binomial',
  data = inla.stack.data(stack_Plas, spde = spde), 
  control.predictor = list(A = inla.stack.A(stack_Plas), compute = TRUE, link = 1),
  control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"),
  control.compute = list(dic=TRUE,cpo=TRUE), num.threads = 3, verbose = TRUE)

mod_spde_ParaHaem <- inla(f_spde, family = 'Binomial',
  data = inla.stack.data(stack_ParaHaem, spde = spde), 
  control.predictor = list(A = inla.stack.A(stack_ParaHaem), compute = TRUE, link = 1), 
  control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"),
  control.compute = list(dic=TRUE,cpo=TRUE), num.threads = 3, verbose = TRUE)

mod_spde_Haem <- inla(f_spde, family = 'Binomial',
  data = inla.stack.data(stack_Haem, spde = spde), 
  control.predictor = list(A = inla.stack.A(stack_Haem), compute = TRUE, link = 1), 
  control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"),
  control.compute = list(dic=TRUE,cpo=TRUE), num.threads = 3, verbose = TRUE)

mod_spde_Leuc <- inla(f_spde, family = 'Binomial',
  data = inla.stack.data(stack_Leuc, spde = spde), 
  control.predictor = list(A = inla.stack.A(stack_Leuc), compute = TRUE, link = 1), 
  control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"),
  control.compute = list(dic=TRUE,cpo=TRUE), num.threads = 3, verbose = TRUE)


#######
# Model with spatial and phylogenetic random effects and varying coefficient

f_spde.varcoeff.rw <- y ~ -1 + abslat + elevation + bio1 + bio12 + bio14 + bio15 +
			ndvi.mean + ndvi.sd + LC.forest + LC.wetland + forstrat.canopy + 
			log.bm + bird.richness_sc + migrate.distance_sc +
			f(r.1, elevation, model='rw1') +
			f(r.2, bio1, model='rw1') + 
			f(r.3, bio12, model='rw1') + 
			f(r.4, bio14, model='rw1') + 
			f(r.5, bio15, model='rw1') + 
			f(r.6, ndvi.mean, model='rw1') + 
			f(r.7, ndvi.sd, model='rw1') + 
			f(r.8, LC.forest, model='rw1') + 
			f(r.9, LC.wetland, model='rw1') +
			f(r.10, forstrat.canopy, model='rw1') + 
			f(r.11, log.bm, model='rw1') +
			f(r.12, bird.richness_sc, model='rw1') +
			f(r.13, migrate.distance_sc, model='rw1') +
			f(family_id, model='generic0', Cmatrix=phyl.prec_family) + 
			f(sp_id, model = spde) +
			f(zooreg, model='iid') +
			f(year, model='rw1') +
			f(isolation.source, model='iid')


mod_spde.varcoeff_Plas <- inla(f_spde.varcoeff.rw, family = 'binomial',
  data = inla.stack.data(stack_Plas, spde = spde), 
  control.predictor = list(A = inla.stack.A(stack_Plas), hyper=pc.prec,compute = TRUE, link = 1),
  control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"),
  control.compute = list(dic=TRUE,cpo=TRUE), num.threads = 3, verbose = TRUE)

mod_spde.varcoeff_ParaHaem <- inla(f_spde.varcoeff.rw, family = 'binomial',
  data = inla.stack.data(stack_ParaHaem, spde = spde), 
  control.predictor = list(A = inla.stack.A(stack_ParaHaem), hyper=pc.prec,compute = TRUE, link = 1),
  control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"),
  control.compute = list(dic=TRUE,cpo=TRUE), num.threads = 3, verbose = TRUE)

mod_spde.varcoeff_Haem <- inla(f_spde.varcoeff.rw, family = 'binomial',
  data = inla.stack.data(stack_Haem, spde = spde), 
  control.predictor = list(A = inla.stack.A(stack_Haem), hyper=pc.prec,compute = TRUE, link = 1),
  control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"),
  control.compute = list(dic=TRUE,cpo=TRUE), num.threads = 3, verbose = TRUE)

mod_spde.varcoeff_Leuc <- inla(f_spde.varcoeff.rw, family = 'binomial',
  data = inla.stack.data(stack_Leuc, spde = spde), 
  control.predictor = list(A = inla.stack.A(stack_Leuc), hyper=pc.prec, compute = TRUE, link = 1),
  control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"),
  control.compute = list(dic=TRUE,cpo=TRUE), num.threads = 3, verbose = TRUE)

Penalized complexity prior
 pc.prec = list(prec = list(prior = "pc.prec", param = c(5, 0.01)))

mod_spde.varcoeff.pc.prior_Plas <- inla(f_spde.varcoeff.rw, family = 'binomial',
  data = inla.stack.data(stack_Plas, spde = spde), 
  control.predictor = list(A = inla.stack.A(stack_Plas), hyper=pc.prec, compute = TRUE, link = 1),
  control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"),
  control.compute = list(dic=TRUE,cpo=TRUE), num.threads = 3, verbose = TRUE)

mod_spde.varcoeff.pc.prior_ParaHaem <- inla(f_spde.varcoeff.rw, family = 'binomial',
  data = inla.stack.data(stack_ParaHaem, spde = spde), 
  control.predictor = list(A = inla.stack.A(stack_ParaHaem), hyper=pc.prec, compute = TRUE, link = 1),
  control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"),
  control.compute = list(dic=TRUE,cpo=TRUE), num.threads = 3, verbose = TRUE)

mod_spde.varcoeff.pc.prior_Haem <- inla(f_spde.varcoeff.rw, family = 'binomial',
  data = inla.stack.data(stack_Haem, spde = spde), 
  control.predictor = list(A = inla.stack.A(stack_Haem), hyper=pc.prec, compute = TRUE, link = 1),
  control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"),
  control.compute = list(dic=TRUE,cpo=TRUE), num.threads = 3, verbose = TRUE)

mod_spde.varcoeff.pc.prior_Leuc <- inla(f_spde.varcoeff.rw, family = 'binomial',
  data = inla.stack.data(stack_Leuc, spde = spde), 
  control.predictor = list(A = inla.stack.A(stack_Leuc), hyper=pc.prec, compute = TRUE, link = 1),
  control.inla=list(cmin=0, int.strategy="eb", strategy="adaptive"),
  control.compute = list(dic=TRUE,cpo=TRUE), num.threads = 3, verbose = TRUE)


#############
# Functions to extract selected INLA model output

fu_extract.fixeff <- function(mod){
	mod <- mod
	fixeff <- mod$summary.fixed[, c("mean","0.025quant","0.975quant")]
	colnames(fixeff) <- c("mean", "CI.1", "CI.2")
	fixeff$significant <- 0
	fixeff$significant[which(fixeff$CI.1>0 & fixeff$CI.2>0)] <- 1
	fixeff$significant[which(fixeff$CI.1<0 & fixeff$CI.2<0)] <- 1 
	rm(mod)
	fixeff
}

fu_extract.randeff.family <- function(mod){
	mod <- mod
	reffno <- which(names(mod$summary.random)=="family_id" | names(mod$summary.random)=="family_id")
	randeff.family <- cbind(family= family_name, 
	round(invlogit(mod$summary.random[[reffno]][, c("mean","0.025quant","0.975quant")]),3))
	colnames(randeff.family) <- c("family", "mean", "CI.1", "CI.2") 
	rm(mod)
	randeff.family
}

fu_extract.randeff.realm <- function(mod){
	mod <- mod
	reffno <- which(names(mod$summary.random)=="realm_id" | names(mod$summary.random)=="realm")
	randeff.realm <-  data.frame(realm= realm_name, mean=rep(NA, nrealm), CI.1=rep(NA, nrealm), CI.2=rep(NA, nrealm))
	randeff <- round(invlogit(mod$summary.random[[reffno]][, c("mean","0.025quant","0.975quant")]),3)
	realm.IDs <- mod$summary.random$realm_id$ID
	randeff.realm[realm.IDs,c("mean", "CI.1", "CI.2")] <- randeff
	rm(mod)
	randeff.realm
}

fu_extract.randeff.pts <- function(mod){
	mod <- mod
	npts <- length(pts_id)
	reffno <- which(names(mod$summary.random)=="ID.pts" | names(mod$summary.random)=="pts_id")
	randeff <- round(invlogit(mod$summary.random[[reffno]][, c("mean","0.025quant","0.975quant")]),3)
	randeff.pts <- data.frame(pts_id= pts_id, mean=rep(NA,npts), CI.1=rep(NA, npts), CI.2=rep(NA, npts))
   	pts.match <- match(mod$summary.random$pts_id$ID, pts_id)
	randeff.pts[pts.match,c("mean", "CI.1", "CI.2")] <- randeff
	rm(mod)
	randeff.pts
}

fu_extract.var.randeff <- function(mod, name_reff){
	mod <- mod
	reff <- which(names(mod$marginals.hyperpar)==name_reff)
	var.randeff <- round(inla.tmarginal(function(x) 1 / x, mod$marginals.hyperpa[[reff]], method = "linear") %>% inla.qmarginal(c(0.5, 0.025, 0.975), .), 3)
	as.vector(var.randeff)
	names(var.randeff) <- c("mean", "CI.1", "CI.2")
	var.randeff
}

fu_plot_randeff.family.1 <- function(species){
  sel <- which(posterior_randeff.family$Species==species)
  ggplot(posterior_randeff.family[sel,], aes(x=family, y=mean)) + 	
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_point(position=position_dodge(0.4), size=2)+ xlab("Covariate") + ylab("Posterior estimate (infection risk)") +
  geom_errorbar(aes(ymin = CI1, ymax = CI2), size=0.8, width=0.3, position=position_dodge(0.4)) + scale_color_brewer(palette="Dark2") +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1)) + ylim(0,1) + ggtitle(paste(species))
}

fu_plot_randeff.family.2 <- function(mod, species){
	mod <- mod
	species = species
	randeff.family <- fu_extract.randeff.family(mod)
  ggplot(randeff.family, aes(x=family, y=mean)) + 	
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_point(position=position_dodge(0.4), size=2)+ xlab("Covariate") + ylab("Posterior estimate (infection risk)") +
  geom_errorbar(aes(ymin = CI.1, ymax = CI.2), size=0.8, width=0.3, position=position_dodge(0.4)) + scale_color_brewer(palette="Dark2") +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1)) + ylim(0,1) + ggtitle(paste(species))
}


#################
# Output summary

# Vector of covariate (fixed effect names)
covar_str <- c("abslat", "elevation", "bio1", "bio12", "bio14", "bio15", "ndvi.mean", "ndvi.sd", "LC.forest", "LC.wetland", "forstrat.canopy", "log.bm", "bird.richness_sc", "migrate.distance_sc")
covar_name <- c("Equator dist", "Elevation", "Temp mean", "Rain total", "Rain min", "Rain var", "NDVI mean", "NDVI sd", "Forest cover", "Wetland cover", "Canopy foraging", "Body mass", "Bird species richness", "Migration distance")
ncovar <- length(covar_str)

# Vector of model names 
models_str <- c("mod_realm_", "mod_pts_", "mod_fixeff_", "mod_phyl.0_", "mod_phyl.Cov_", "mod_reg.phyl_", "mod_spde_", "mod_spde.varcoeff.pc.prior_")
models_name <- c("regional GLM", "location GLM", "covariate GLM", "phyl-0 GLMM", "phyl-cov GLMM", "regional phyl GLMM", "spatio-phyl GLMM", "spatio-phyl-varcoeff GLMM")
nmodels <- length(models_str)

par_str  <- c("Haem", "Leuc", "ParaHaem", "Plas")
par_name <- c("Haemoproteus", "Leucocytozoon", "Parahaemoproteus", "Plasmodium")
npar <- length(par_str )

# Tables of model DICs and cpo 
tble_DIC <- tble_cpo.mean <- tble_cpo.sd <-matrix(NA, nrow=nmodels, ncol=npar)
colnames(tble_DIC) <- colnames(tble_cpo.mean) <- colnames(tble_cpo.sd) <- par_name
rownames(tble_DIC) <- rownames(tble_cpo.mean) <- rownames(tble_cpo.sd) <- models_name

for(i in 1:nmodels){
	for(j in 1:npar){
		mod <- get(paste0(models_str[i], par_str [j]))
		tble_DIC[i,j] <- mod$dic$dic	
		tble_cpo.mean[i,j] <- round(mean(mod$cpo$cpo, na.rm=T),2)
		tble_cpo.sd[i,j] <- round(sd(mod$cpo$cpo, na.rm=T),2)
	}
}


#### 
#  Overall average infection probabilities for different parasite genera

tble_inf0 <- rbind(
round(100*invlogit(fu_extract.fixeff(get(paste0("model_inf0_", par_str[1])))[,1:3]),1),
round(100*invlogit(fu_extract.fixeff(get(paste0("model_inf0_", par_str[2])))[,1:3]),1),
round(100*invlogit(fu_extract.fixeff(get(paste0("model_inf0_", par_str[3])))[,1:3]),1),
round(100*invlogit(fu_extract.fixeff(get(paste0("model_inf0_", par_str[4])))[,1:3]),1))
rownames(tble_inf0) <- par_name


#### 
#  Infection probs in different regions

post_infprop.realm <- rbind(
	cbind(fu_extract.randeff.realm(get(paste0("mod_realm_", par_str [1]))), genus=par_name[1]),
	cbind(fu_extract.randeff.realm(get(paste0("mod_realm_", par_str [2]))), genus=par_name[2]),
	cbind(fu_extract.randeff.realm(get(paste0("mod_realm_", par_str [3]))), genus=par_name[3]),
	cbind(fu_extract.randeff.realm(get(paste0("mod_realm_", par_str [4]))), genus=par_name[4]))

# Highest regional infection probabilties 
arrange(post_infprop.realm %>% filter(genus==par_name[1]), CI.1)
arrange(post_infprop.realm %>% filter(genus==par_name[2]), CI.1)
arrange(post_infprop.realm %>% filter(genus==par_name[3]), CI.1)
arrange(post_infprop.realm %>% filter(genus==par_name[4]), CI.1)

post_infprop.realm %>% filter(realm=="Saharo-Arabian"  & genus!="Haemoproteus")
# Infection < 10% for common genera
post_infprop.realm %>% filter(genus!="Haemoproteus" & CI.2<0.1)
# Highest/lowest regional infection with Haemoproteus
arrange(post_infprop.realm %>% filter(genus=="Haemoproteus"), CI.1)
arrange(post_infprop.realm %>% filter(genus=="Haemoproteus"), CI.2)

plot_inf.realm <- ggplot(post_infprop.realm, aes(x=realm, y=mean, group=genus, color=genus)) + 
  geom_hline(yintercept=0, linetype=1, color = "lightgrey") +
  geom_point(position=position_dodge(0.4), size=2, pch=sort(rep(c(19,15,17,8),nrealm)))+ xlab("Region") + ylab("Infection probability") +
  geom_errorbar(aes(ymin = CI.1, ymax = CI.2), size=0.8, width=0.3, position=position_dodge(0.4)) + scale_color_brewer(palette="Dark2") +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1)) +
  theme(text = element_text(size=12),
          title = element_text(size=12),
          strip.text.x = element_text(size=12),
          strip.text.y = element_text(size=15),
	    axis.text.y.left = element_text(size=12))
plot_inf.realm


#### 
# Infection probs in different locations
post_infprop.pts <- rbind(
	cbind(fu_extract.randeff.pts(get(paste0("mod_pts_", par_str[1])))[,2:4], genus=par_name[1], pointID=pts_id, long = pts_long, lat = pts_lat, realm=pts_realm, nsamp = pts_nsample.Haem),
	cbind(fu_extract.randeff.pts(get(paste0("mod_pts_", par_str[2])))[,2:4], genus=par_name[2], pointID=pts_id, long = pts_long, lat = pts_lat, realm=pts_realm, nsamp = pts_nsample.Leuc),
	cbind(fu_extract.randeff.pts(get(paste0("mod_pts_", par_str[3])))[,2:4], genus=par_name[3], pointID=pts_id, long = pts_long, lat = pts_lat, realm=pts_realm, nsamp = pts_nsample.Plas),
	cbind(fu_extract.randeff.pts(get(paste0("mod_pts_", par_str[4])))[,2:4], genus=par_name[4], pointID=pts_id, long = pts_long, lat = pts_lat, realm=pts_realm, nsamp = pts_nsample.Plas))
post_infprop.pts$CIsize <- abs(post_infprop.pts$CI.1 - post_infprop.pts$CI.2)

# Number of locations with < 10% uncertainty in estimates
post_infprop.pts %>% group_by(genus) %>% summarise(n = length(which(abs(CIsize)<0.10)))


# Which threshold in CI size to use for identifying hotspots
threshold_CI.hotspots <- 0.1

# Locations with highest infection probabilities
hotspots_local <- rbind(
top_n(post_infprop.pts %>% filter(genus==par_name[1] & CIsize<=threshold_CI.hotspots), 3, CI.1),
top_n(post_infprop.pts %>% filter(genus==par_name[2] & CIsize<=threshold_CI.hotspots), 3, CI.1),
top_n(post_infprop.pts %>% filter(genus==par_name[3] & CIsize<=threshold_CI.hotspots), 3, CI.1),
top_n(post_infprop.pts %>% filter(genus==par_name[4] & CIsize<=threshold_CI.hotspots), 3, CI.1))
# infection rates at hotspots
range(c(hotspots_local$CI.1, hotspots_local$CI.2)); range(hotspots_local$CI.2)

plotfu_infprop.pts <- function(genus, threshold) {
	sel_infprop <- which(post_infprop.pts$genus==genus & post_infprop.pts$CIsize<=threshold)
	datpts <- post_infprop.pts[sel_infprop,]
	datpts2 <- post_infprop.pts[-sel_infprop,]
	ggplot(data = worldmap) + geom_sf() + theme_classic() + ylim(-55, 80) +
	geom_point(data=datpts, aes(x=long, y=lat, color =mean)) + scale_color_gradient(low = "blue", high = "red", name="Posterior mean") +
	ggtitle(paste(genus)) + ylab("") + xlab("")
} 

plot_infprop.pts <- ggarrange(plotfu_infprop.pts(par_name[1], threshold_CI.hotspots ), plotfu_infprop.pts(par_name[2], threshold_CI.hotspots ),
	plotfu_infprop.pts(par_name[3], threshold_CI.hotspots), plotfu_infprop.pts(par_name[4], threshold_CI.hotspots ), ncol=2, nrow=2, common.legend = F)
plot_infprop.pts

#### 
#  Infection probs of different families
post_phyl.0_family <- rbind(
	cbind(fu_extract.randeff.family(get(paste0("mod_phyl.0_", par_str [1])))[,1:4], genus=par_name[1]),
	cbind(fu_extract.randeff.family(get(paste0("mod_phyl.0_", par_str [2])))[,1:4], genus=par_name[2]),
	cbind(fu_extract.randeff.family(get(paste0("mod_phyl.0_", par_str [3])))[,1:4], genus=par_name[3]),
	cbind(fu_extract.randeff.family(get(paste0("mod_phyl.0_", par_str [4])))[,1:4], genus=par_name[4]))
post_phyl.0_family$CIsize <- abs(post_phyl.0_family$CI.1 - post_phyl.0_family$CI.2)

post_spde.varcoeff_family <- rbind(
	cbind(fu_extract.randeff.family(get(paste0(tail(models_str,1), par_str [1])))[,1:4], genus=par_name[1]),
	cbind(fu_extract.randeff.family(get(paste0(tail(models_str,1), par_str [2])))[,1:4], genus=par_name[2]),
	cbind(fu_extract.randeff.family(get(paste0(tail(models_str,1), par_str [3])))[,1:4], genus=par_name[3]),
	cbind(fu_extract.randeff.family(get(paste0(tail(models_str,1), par_str [4])))[,1:4], genus=par_name[4]))
post_spde.varcoeff_family$CIsize <- abs(post_spde.varcoeff_family$CI.1 - post_spde.varcoeff_family$CI.2)

plot_comp.phyleffect <- plot(post_phyl.0_family$mean, post_spde.varcoeff_family$mean, xlim=c(0,1), ylim=c(0,1)); abline(0,1)
plot_comp.phyleffect

## Variance of phylogenetic random effect
post_spde.varcoeff_var.family <- rbind(
	round(c(fu_extract.var.randeff( get(paste0(tail(models_str,1), par_str [1])), "Precision for family_id")),1),
	round(c(fu_extract.var.randeff( get(paste0(tail(models_str,1), par_str [2])), "Precision for family_id")),1),
	round(c(fu_extract.var.randeff( get(paste0(tail(models_str,1), par_str [3])), "Precision for family_id")),1),
	round(c(fu_extract.var.randeff( get(paste0(tail(models_str,1), par_str [4])), "Precision for family_id")),1))
rownames(post_spde.varcoeff_var.family) <- par_name
post_spde.varcoeff_var.family

# Number of phylogenetic random effects with < 10% uncertainty in estimates
post_spde.varcoeff_family %>% group_by(species) %>% summarise(n = length(which(abs(CIsize)<0.10)))

# Highest/lowest infection probabilities
arrange(post_spde.varcoeff_family %>% filter(genus==par_name[1]), CI.1)
arrange(post_spde.varcoeff_family %>% filter(genus==par_name[1]), CI.2)

arrange(post_spde.varcoeff_family %>% filter(genus==par_name[2]), CI.1)
arrange(post_spde.varcoeff_family %>% filter(genus==par_name[2]), CI.2)

arrange(post_spde.varcoeff_family %>% filter(genus==par_name[3]), CI.1)
arrange(post_spde.varcoeff_family %>% filter(genus==par_name[3]), CI.2)

arrange(post_spde.varcoeff_family %>% filter(genus==par_name[4]), CI.1)
arrange(post_spde.varcoeff_family %>% filter(genus==par_name[4]), CI.2)

# Infection probabilites < 1%
arrange(post_spde.varcoeff_family %>% filter(genus==par_name[1] & CIsize<0.1 & CI.2<0.01), CI.2)
arrange(post_spde.varcoeff_family %>% filter(genus==par_name[2] & CIsize<0.1 & CI.2<0.01), CI.2)
arrange(post_spde.varcoeff_family %>% filter(genus==par_name[3] & CIsize<0.1 & CI.2<0.01), CI.2)
arrange(post_spde.varcoeff_family %>% filter(genus==par_name[4] & CIsize<0.1 & CI.2<0.01), CI.2)

# Infection probability above global average
arrange(post_spde.varcoeff_family %>% filter(genus==par_name[1] & CIsize<0.1 & CI.1>tble_inf0$CI.2[1]), CI.2)
arrange(post_spde.varcoeff_family %>% filter(genus==par_name[2] & CIsize<0.1 & CI.1>tble_inf0$CI.2[2]), CI.2)
arrange(post_spde.varcoeff_family %>% filter(genus==par_name[3] & CIsize<0.1 & CI.1>tble_inf0$CI.2[3]), CI.2)
arrange(post_spde.varcoeff_family %>% filter(genus==par_name[4] & CIsize<0.1 & CI.1>tble_inf0$CI.2[4]), CI.2)



plot_phyleffect <- 
ggarrange(
fu_plot_randeff.family.2(get(paste0( models_str[8], par_str [1])), paste(par_name[1])),
fu_plot_randeff.family.2(get(paste0( models_str[8], par_str [2])), paste(par_name[2])),
fu_plot_randeff.family.2(get(paste0(models_str[8], par_str [3])), paste(par_name[3])),
fu_plot_randeff.family.2(get(paste0(models_str[8], par_str [4])), paste(par_name[4])))


#### 
#  Covariate/fixed effects
post_fixeff <- rbind(
	cbind(fu_extract.fixeff(get(paste0(tail(models_str,1), par_str [1]))), genus=par_name[1], covar = covar_name),
	cbind(fu_extract.fixeff(get(paste0(tail(models_str,1), par_str [2]))), genus=par_name[2], covar = covar_name),
	cbind(fu_extract.fixeff(get(paste0(tail(models_str,1), par_str [3]))), genus=par_name[3], covar = covar_name),
	cbind(fu_extract.fixeff(get(paste0(tail(models_str,1), par_str [4]))), genus=par_name[4], covar = covar_name))

# Odds ratios for 'significant' effects
post_fixeff_sign <- post_fixeff %>% filter(significant==1)
post_fixeff_sign.odds <- post_fixeff_sign[, c("covar", "genus", "mean", "CI.1", "CI.2")] 
post_fixeff_sign.odds[, c("mean", "CI.1", "CI.2")] <- round((exp(post_fixeff_sign.odds[, c("mean", "CI.1", "CI.2")])),2)


pf=post_fixeff_sign.odds
post_fixeff_sign.odds_print <- data.frame(Variable=pf$covar, Species=pf$genus, Odds_ratio=paste0(pf$mean, " (", pf$CI.1, " - ", pf$CI.2, ")"))

save(post_fixeff_sign.odds_print, file ="post_fixeff_sign.odds_print.rds")
write.csv(post_fixeff_sign.odds_print, file ="post_fixeff_sign.odds_print.csv")


plot_fixeff <- ggplot(post_fixeff_sign, aes(x=covar, y=mean, group=genus, color=genus)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_point(position=position_dodge(0.4), size=2, alpha=(0.2+post_fixeff_sign$significant*3))+ xlab("Covariate") + ylab("Posterior estimate") +
  geom_errorbar(aes(ymin = CI.1, ymax = CI.2), size=0.8, width=0.3, position=position_dodge(0.4), alpha=(0.2+post_fixeff_sign$significant*3)) + scale_color_brewer(palette="Dark2") +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1)) + ylim(range(post_fixeff_sign$CI.1, post_fixeff_sign$CI.2))
plot_fixeff


######
## Variation in random effect for fixed effect coefficients
post_coefficient_reffvar <- data.frame(par_str= rep(par_str, each=13), genus= rep(par_name, each=13), r_str = rep(paste0("r.", 1:13), 4) , covar = rep(covar_name[-1], 4), variance_mean=rep(NA, 4*13), variance_CI.1=rep(NA, 4*13), variance_CI.2=rep(NA, 4*13))
for(i in 1:dim(post_coefficient_reffvar)[1]){
	mod <- get(paste0(tail(models_str,1), post_coefficient_reffvar$par_str[i]))
	effname <- paste0("Precision for ", post_coefficient_reffvar$r_str[i])
	post_coefficient_reffvar$variance_mean[i] <- fu_extract.var.randeff(mod, effname)["mean"]
	post_coefficient_reffvar$variance_CI.1[i] <- fu_extract.var.randeff(mod, effname)["CI.1"]
	post_coefficient_reffvar$variance_CI.2[i] <- fu_extract.var.randeff(mod, effname)["CI.2"]
}


table_covariate.variance <- as_tibble(post_coefficient_reffvar) %>% filter(genus!="Haemoproteus") %>% select(genus, covar, variance_mean, variance_CI.1, variance_CI.2)


## Effect sizes for covariates with variable coefficients
post_coefficient_reffvar.sel <-  post_coefficient_reffvar %>% filter( variance_CI.1>0)

nvarcoeff <- dim(post_coefficient_reffvar.sel)[1]

# Generate dataframe with posterior estimates of variable coefficient estiamates
for(i in 1:nvarcoeff){
	gen <- post_coefficient_reffvar.sel$genus[i]
	pstr <- post_coefficient_reffvar.sel$par_str[i]
	cov <- post_coefficient_reffvar.sel$covar[i]
	rstr <- post_coefficient_reffvar.sel$r_str[i]
	effect <- paste0(gen, " ~ ", cov)
	mod <- get(paste0(tail(models_str,1), pstr))
	fixeff.mean <- post_fixeff %>% filter(genus==gen & covar==cov)
	fixeff.mean <- fixeff.mean$mean
	reff <- which(names(mod$summary.random)==rstr)
	coef.adj_mean <- fixeff.mean + mod$summary.random[[reff]][, "mean"]
	coef.adj_CI.1 <- fixeff.mean + mod$summary.random[[reff]][, "0.025quant"]
	coef.adj_CI.2  <- fixeff.mean + mod$summary.random[[reff]][, "0.975quant"]
	if(pstr =="Leuc"){
		realms <- sort(unique(DatMergeLeuc$realm))
	}else{
		realms <- realm_name
	}
	df_varcoeff <- data.frame(
			effect.no = i,
			genus = gen,
			par_str = pstr,
			covar = cov,
			r_str= rstr,
			realm = realms,
			effect = effect,	
			coef.adj_mean = coef.adj_mean,
			coef.adj_CI.1 = coef.adj_CI.1,
			coef.adj_CI.2 = coef.adj_CI.2)
	if(i==1){
		post_varcoeff <- df_varcoeff
	}else{
		post_varcoeff <- rbind(post_varcoeff, df_varcoeff)	
	}	
}


# Select variable coefficient with clearly distinct CIs
varcoeff_sel <- rep(NA, nvarcoeff)
varcoeff.opposite_sel <- rep(NA, nvarcoeff)
for(i in 1:nvarcoeff){
  	eff <- post_varcoeff %>% filter(effect.no==i)
	# distinct lower CIs > any other upper CIs OR upper CIS < any other lower CIs
	varcoeff_sel[i] <- ifelse(max(eff$coef.adj_CI.1) > min(eff$coef.adj_CI.2), 1,0)	 
	varcoeff.opposite_sel[i] <- ifelse((max(eff$coef.adj_CI.1)>0 & min(eff$coef.adj_CI.2)<0), 1,0) 
}


post_varcoeff <- post_varcoeff %>% filter(!is.na(match(effect.no, which(varcoeff_sel==1))))
post_varcoeff <- post_varcoeff %>% filter(par_str!="Haem")

unique(post_varcoeff$effect)
write.csv(post_varcoeff, file="post_varcoeff.csv")

unique(post_varcoeff$covar [which(post_varcoeff$genus==par_name[1])])
unique(post_varcoeff$covar [which(post_varcoeff$genus==par_name[2])])
unique(post_varcoeff$covar [which(post_varcoeff$genus==par_name[3])])
unique(post_varcoeff$covar [which(post_varcoeff$genus==par_name[4])])


plot_varcoeff <- ggplot(post_varcoeff , aes(x=realm, y=coef.adj_mean, color=genus)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_point(position=position_dodge(0.4), size=2)+ xlab("Realm") + ylab("Posterior estimate") +
  geom_errorbar(aes(ymin = coef.adj_CI.1, ymax = coef.adj_CI.2), size=0.8, width=0.3, position=position_dodge(0.4)) + scale_color_brewer(palette="Dark2") +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1)) +
	facet_wrap(~ effect)
plot_varcoeff


## Variables with opposing variable effects
post_varcoeff.oppose <- post_varcoeff %>% filter(!is.na(match(effect.no, which(varcoeff.opposite_sel==1))))

plot_varcoeff.oppose <- ggplot(post_varcoeff.oppose , aes(x=realm, y=coef.adj_mean, color=genus)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_point(position=position_dodge(0.4), size=2)+ xlab("Realm") + ylab("Posterior estimate") +
  geom_errorbar(aes(ymin = coef.adj_CI.1, ymax = coef.adj_CI.2), size=0.8, width=0.3, position=position_dodge(0.4)) + scale_color_brewer(palette="Dark2") +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1)) +
	facet_wrap(~ effect)
plot_varcoeff.oppose




