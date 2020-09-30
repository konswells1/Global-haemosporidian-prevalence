# Global-haemosporidian-prevalence

Some code to analyse to compile environmental data and analyse the infection probability of birds with avian haemosporidian parasites from multiple field survey studies conducted at > 500 geo-referenced locations.

Relevant climatic variables at sample locations were obtained from the WorldClim database of gridded climate data at 0.01 degree resolution (http://world clim.org/version2) with an interest in the follwoing variables: annual mean temperature (bio1), annual rainfall (bio12), rainfall of driest month (bio14), and rainfall seasonality (coefficient of variation in rainfall over year, bio15) to
Elevation for all locations were quantified using Shuttle Radar Topography Mission (SRTM) data as accessible through the raster package in R. 
The proportion of cover with forest and wetland in buffers of 10km radius around sample locations were estimated based on Copernicus landcover data from 2010 (map version 2.07; https://cds.climate.copernicus.eu).
The normalized difference vegetation index (NDVI) for the year 2010 in buffers of 10km radius around all sampling locations was computed from the Terra Moderate Resolution Imaging Spectroradiometer (MODIS, MOD13Q1 version 6, https://lpdaac.usgs.gov/products/mod13q1v006/) and calculated the mean and 1 standard deviation of the NDVI as measures of the vegetation density and its variation within a year.

Bird host-specific species traits from the EltonTraits v1.0 database (https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/13-1917.1). In particular, we considered host body mass and the proportion of time species spend foraging in the upper canopy as relevant covariates inour analysis.

In order to identify key drivers of bird infection with haemosporidian parasites while accounting for possible spatio-temporal and phylogenetic pattern underpinning the global dataset, we used a Bayesian statistical model to jointly estimate the posterior distribution of fixed and random effect model parameters. 
