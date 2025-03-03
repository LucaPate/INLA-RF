# Title: runme_SPDE_strong_datasim.R
# Author: Mario Figueira
# Date: 2025-02-18
# Description: Code to simulate a strong non-linearity effect for the covariates
#              with a spatio-temporal structure, using: Q_st = Q_s(SPDE) \otimes Q_t(AR1)

# Last update: 2025-02-18

remove(list=ls())

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

library(ranger)

seed <- 1234 # set a seed for reproducibility
set.seed(seed = seed)

# Custom functions ----

## Function to fix the precision matrix from inla.call.object$misc$configs$config[[k]]$Q (or Qprior or Qinv)
fix.Q <- function(Q) {
  # Q: the precision matrix from the inla.call.object; inla.call.object <- inla(...) 
  d <- diag(Q)
  Q <- Q + t(Q)
  diag(Q) <- d
  return (Q)
}

## Function to simulate a SPDE effect in R2
simulation_SPDE <- function(mesh, sigma = 1, range = 1, constr = TRUE, n_rep = 1, n_cores = 1, seed){
  # mesh: mesh to create the precision matrix of the SPDE
  # sigma: the marginal precision of the CAR prior (the variability between first-order neighbors)
  # range: the spatial range; it means the range at which the correlation falls to 0.1 (in the correl)
  # constr: an argument that, if TRUE, forces the GMRF to have a global mean of zero
  # seed: seed to reproduce the result
  if(!missing(seed)){set.seed(seed)}
  Q = fmesher::fm_matern_precision(x = mesh, alpha = 2, rho = range, sigma = sigma)
  C0 = diag(fmesher::fm_fem(mesh = mesh, Q = Q)$c0)
  L = Matrix::chol(Q)
  u_sp <- parallel::mclapply(X = seq_len(n_rep), mc.cores = n_cores, FUN = function(i){
    w = rnorm(nrow(Q))
    u_sp = backsolve(L, w)
    if(constr){
      u_sp = u_sp - mean(C0 * u_sp) / C0
    }
    return(u_sp)
  }) %>% do.call(what = cbind, .)
  return(list(u=u_sp, Q=Q))
}

# Function to repeat rows rep(A, times = n) = (A, A, A) by rows

rep_row <- function(x, times){
  lapply(X = seq_len(times), FUN = function(i){x}) %>% do.call(what = rbind, .)
}

## Defining the study region ----
data(PRborder)
sf_PRborder <- st_sfc(st_polygon(x = list(PRborder))) # boundary of the sr (study region)
PR_nchull <- fm_nonconvex_hull_inla(x =  sf_PRborder, convex = -0.05)
coord_PR_nchull <- PR_nchull$loc[c(PR_nchull$idx[1:which(PR_nchull$idx[1,1]==PR_nchull$idx[,2]),1],1),1:2]
st_PR_int_nchull <-  st_polygon(x = list(coord_PR_nchull))

# Plot of the simplified boundary and the original one 
gg_PR_int <- ggplot() + 
  geom_sf(st_PR_int_nchull, mapping = aes(color = "Internal")) + 
  geom_sf(sf_PRborder, mapping = aes(color = "Original")) + 
  theme_bw() + labs(title = "Original and internal boundary", color = "Boundaries") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

## Mesh creation and 
mesh_sim <- fm_mesh_2d_inla(boundary = list(st_PR_int_nchull), max.edge = c(0.2,0.5), offset = c(-0.01,-0.1))

gg_mesh <- ggplot() + 
  gg(mesh_sim) + labs(title = "Mesh for the SPDE-FEM simulation") +
  coord_sf() + theme_bw() + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank())

# grid.arrange(arrangeGrob(grobs = list(gg_PR_int, gg_mesh), ncol = 2))
gg_PR_int; gg_mesh

## Spatio-temporal effect simulation
st_bbox_sr <- st_bbox(obj = st_PR_int_nchull)
loc_grid <- expand.grid(x = seq(from = st_bbox_sr[1], to = st_bbox_sr[3], length.out = 1E2), y = seq(from = st_bbox_sr[2], to = st_bbox_sr[4], length.out = 1E2))
st_loc_grid_t <- st_multipoint(x = as.matrix(loc_grid)) %>% st_geometry(.) %>% st_cast("POINT")
idx <- st_intersects(x = st_loc_grid_t, st_PR_int_nchull, sparse = FALSE)
st_loc_grid <- st_loc_grid_t[idx,] # Grid points that are inside the study region

size_spatial <- 1E3 # Spatial sample size per temporal node (time points)
size_time <- 8 # Number of time points
sample_size  <- size_spatial * size_time # Total sample size

st_sample_points <- st_sample(x = st_PR_int_nchull, size = sample_size) # Sample points

# simulation_SPDE(mesh = mesh_sim, sigma = 1, range = 0.3, constr = TRUE, seed = seed) # Test of the simulation_SPDE function for one temporal node
sp_sim <- simulation_SPDE(mesh = mesh_sim, n_rep = size_time, sigma = 1, range = coord_PR_nchull %>% dist(.) %>% max(.)/2, constr = TRUE)$u # %>% as.vector(.)

spt_sim <- sp_sim
rho_t <- 0.7 # Autoregressive parameter for the temporal structure
for(i in 2:size_time){
  # spt_sim[((i-1)*size_spatial+1):(i*size_spatial)] <- rho_t * spt_sim[((i-2)*size_spatial+1):((i-1)*size_spatial)] + sqrt(1 - rho_t^2) * sp_sim[((i-1)*size_spatial+1):(i*size_spatial)]
  spt_sim[,i] <- rho_t*spt_sim[,i-1] + sqrt(1 - rho_t**2) * sp_sim[,i]
}

spt_sim <- spt_sim %>% as.vector(.)

# df_sp <- st_sf(data.frame(sp = drop(fm_basis(x = mesh_sim, loc = rep_row(x = st_coordinates(st_loc_grid)[,1:2], times = size_time)) %*% sp_sim$u)), geometry = rep(st_loc_grid, times = size_time))
A_spt_grid <- INLA::inla.spde.make.A(mesh  = mesh_sim, loc = rep_row(x = st_coordinates(st_loc_grid)[,1:2], times = size_time), group = rep(seq_len(size_time), each = length(st_loc_grid)), n.group = size_time)
df_sp_grid <- st_sf(data.frame(sp = drop(A_spt_grid %*% spt_sim), id_time = rep(seq_len(size_time), each = length(st_loc_grid))), geometry = rep(st_loc_grid, times = size_time))

A_spt_sample <- INLA::inla.spde.make.A(mesh  = mesh_sim, loc = st_sample_points, group = rep(seq_len(size_time), each = size_spatial), n.group = size_time)
df_sp_sample <- st_sf(data.frame(sp = drop(A_spt_sample %*% spt_sim), id_time = rep(seq_len(size_time), each = size_spatial)), geometry = st_sample_points)

## Plot of the spatial effect (grid)
fwrap_names <- paste("Temporal node:", unique(df_sp_sample$id_time)); names(fwrap_names) <- unique(df_sp_sample$id_time)
ggsp_grid <- ggplot() + 
  geom_sf(data = df_sp_grid, mapping = aes(color = sp), pch = 15) +
  scale_color_viridis_c(option = "turbo") + 
  coord_sf() + theme_bw() + labs(color = "Values", title = "Spatial effect (SPDE)") +
  facet_wrap(facets = ~ id_time, labeller = labeller(id_time = fwrap_names), ncol = 4) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), strip.text = element_text(face = "bold", hjust = 0.5))

## Plot of the spatial effect (sample locations)
ggsp_sample <- ggplot() +
  geom_sf(data = df_sp_sample, mapping = aes(color = sp)) +
  scale_color_viridis_c(option = "turbo") +
  coord_sf() + theme_bw() + labs(color = "Values", title = "Spatial effect (SPDE)") +
  facet_wrap(facets = ~ id_time, labeller = labeller(id_time = fwrap_names), ncol = 4) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), strip.text = element_text(face = "bold", hjust = 0.5))

## Simulation of the covariates ----

X1 <- factor(sample(LETTERS[1:3], sample_size, replace = TRUE))
# table(X1) # distribution of the different levels of the categorical variable

# Coefficients for the qualitative regressor (categorical variables) with a sum zero constraint
beta_cat <- cbind(rnorm(n = 3, mean = 0, sd = 2)) %>% apply(X = ., MARGIN = 2, FUN = function(x){x - mean(x)})

ggcat <- ggplot() + 
  geom_sf(data = st_sf(data.frame(X1 = unclass(X1)), geometry = st_sample_points), mapping = aes(color = as.factor(X1))) + 
  coord_sf() + theme_bw() + labs(color = "Categories", title = "Categories in the sample locations") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# grid.arrange(arrangeGrob(grobs = list(ggsp_grid, ggsp_sample, ggcat), ncol = 3))
ggsp_grid; ggsp_sample; ggcat

## Simulation of the response variable under the strong non-linearity case ----
X2 <- rnorm(sample_size)
fX2 <- sin(2*X2)*(2*X2)
ggcov2_strong <- ggplot(data.frame(X2, fX2)) + 
  geom_line(mapping = aes(x = X2, y = fX2)) + 
  labs(title = "Variable X2") +
  theme_bw() + theme(plot.title = element_text(face = "bold", hjust = 0.5))

X3 <- runif(sample_size)
fX3 <- sin(X3^4) + cos(2.5*pi*X3)  
ggcov3_strong <- ggplot(data.frame(X3, fX3)) + 
  geom_line(mapping = aes(x = X3, y = fX3)) + 
  labs(title = "Variable X3") +
  theme_bw() + theme(plot.title = element_text(face = "bold", hjust = 0.5))

sigma_gauss <- (50)**(-1/2) # Defining the stdev. through the precision
y <- rnorm(n = sample_size, mean = df_sp_sample$sp + beta_cat[unclass(X1)] + fX2 + fX3, sd = sigma_gauss) # Simulation of the response variable
dfsf_sample_strong <- st_sf(data.frame(y = y, sp = df_sp_sample$sp, id_time = rep(seq_len(size_time), each = size_spatial),
                                       X1 = X1, fX1 = beta_cat[unclass(X1)], 
                                       X2 = X2, fX2 = fX2, X3 = X3, fX3, fX3), 
                            geometry = st_geometry(df_sp_sample)) # Data frame for the simulation with strong non-linearity  

gg_y_strong <- ggplot() + # Response variable under weak linearity in the sample locations
  geom_sf(data = dfsf_sample_strong, mapping = aes(color = y)) + 
  scale_color_viridis_c(option = "turbo") +
  coord_sf() + theme_bw() + labs(color = "Values", title = "Response variable in the sample locations") + 
  facet_wrap(facets = ~ id_time, labeller = labeller(id_time = fwrap_names), ncol = 4) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), strip.text = element_text(face = "bold", hjust = 0.5))

ggcov2_strong; ggcov3_strong; gg_y_strong

list_spt_sim <- list(spt_sim = spt_sim, mesh = mesh_sim)
saveRDS(object = dfsf_sample_strong, file = "./Scripts/Data_strong.RDS")
saveRDS(object = list_spt_sim, file = "./Scripts/list_spt_sim.RDS")
