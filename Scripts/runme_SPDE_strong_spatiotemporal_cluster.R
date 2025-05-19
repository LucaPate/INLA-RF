# Title: spatio_temporal_clustering
# Author: Mario Figueira
# Date: 2025-04-24
# Description: Create a spatio-temporal clustering base on an quasi equally-distribution of the data along the partitions

# Last update: 2025-04-24

remove(list=ls())

# Loading libraries ----

library(Matrix)
library(parallel)

library(INLA)
library(fmesher)
library(inlabru)

library(ggplot2)
library(gridExtra)
library(ggtext)
library(dplyr)
library(sf)

seed <- 1234 # set a seed for reproducibility
set.seed(seed = seed)

# Custom functions ----

split_poly <- function(sf_poly, sf_obs_points, n_areas, seed){
  #' Split polygon
  #' 
  #' @sf_poly Geometry with the polygon to be splitted
  #' @param n_areas The number of resulting areas
  #' @param n_points The number of points randomly simulated for the k-means algorithm
  
  if(!missing(seed)){set.seed(seed)}
  
  # k-means clustering using sf_obs_points
  points <- do.call(rbind, st_geometry(sf_obs_points))
  k_means <- kmeans(points, centers = n_areas)
  # Create voronoi polygons
  voronoi_polys <- st_cast(x = st_multipolygon(x = lapply(st_voronoi(x = st_multipoint(x = k_means$centers)), FUN = function(x){x})) %>% st_geometry(obj = .), to = "POLYGON")
  
  # Intersect to set the right boundaries and compute the area for each polygon
  st_crs(voronoi_polys) <- st_crs(sf_poly)
  clust_areas <- st_intersection(voronoi_polys, sf_poly) # sf_poly$geometry
  clust_areas <- st_sf(data.frame(id = as.factor(1:n_areas), area = st_area(clust_areas)), geometry = clust_areas)
  return(clust_areas)
}

# Loading data ----

DF_strong_sim <- readRDS(file = "./Code_completed/Data_strong.RDS")

## Building the spatial mesh for SPDE-FEM inference ----
data(PRborder)
sf_PRborder <- st_sfc(st_polygon(x = list(PRborder))) # boundary of the sr (study region)
PR_nchull <- fm_nonconvex_hull_inla(x =  sf_PRborder, convex = -0.05)
st_PR_int_nchull <-  st_polygon(x = list(PR_nchull$loc[c(PR_nchull$idx[1:which(PR_nchull$idx[1,1]==PR_nchull$idx[,2]),1],1),1:2]))

mesh_inf <- fm_mesh_2d_inla(boundary = list(st_PR_int_nchull), max.edge = c(0.3,0.6), offset = c(-0.01,-0.1))
sf_poly_int <- st_polygon(x = list(mesh_inf$loc[c(mesh_inf$segm$int$idx[,1],1),1:2]))

ggplot() + geom_sf(data = st_sf(geometry = st_geometry(sf_poly_int)), mapping = aes()) + theme_bw()

table(DF_strong_sim$id_time)
t_each_group <- 2
id_tgroup <- rep(1:4, each = t_each_group*300)
table(id_tgroup)

list_clust_polygons <- list()
id_spt_group <- c()
n_areas <- 4
for(i in unique(id_tgroup)){
  df_sf_points <- DF_strong_sim$geometry[DF_strong_sim$id_time %in% ((i-1)*t_each_group+1):(i*t_each_group),]
  clust_polygons <- split_poly(sf_poly = sf_poly_int, sf_obs_points = df_sf_points, n_areas = n_areas, seed = seed)
  list_clust_polygons[[i]] <- clust_polygons
  id_sp <- st_intersects(x = df_sf_points, y = clust_polygons) %>% as.numeric(.)
  id_spt_group <- c(id_spt_group,(i-1)*n_areas+id_sp)
}

DF_strong_sim$id_spt_group <- id_spt_group

ggplot() + 
  geom_sf(data = list_clust_polygons[[1]], mapping = aes(fill = id)) +
  geom_sf(data = DF_strong_sim[DF_strong_sim$id_spt_group == 1,], mapping = aes(), colour = "red") +
  geom_sf(data = DF_strong_sim[DF_strong_sim$id_spt_group == 2,], mapping = aes(), colour = "blue") +
  geom_sf(data = DF_strong_sim[DF_strong_sim$id_spt_group == 3,], mapping = aes(), colour = "green") +
  geom_sf(data = DF_strong_sim[DF_strong_sim$id_spt_group == 4,], mapping = aes(), colour = "orange") +
  theme_bw()
