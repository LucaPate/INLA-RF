# Title: runme_SPDE_strong_Offset_and_OffsetUnvertainty_CV
# Author: Mario Figueira
# Date: 2025-05-23
# Description: Code with simulated data implementing both algorithm to combine INLA and RF and to perform a spatio-temporal Cross-Validation analysis: 
#              A. Correcting exclusively the offset 
#              B. and correcting the offset along with tranferring the uncertainty

# Last update: 2025-05-26

remove(list=ls())

# Loading libraries ----

library(Matrix)
library(Rcpp)
library(RcppEigen)
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

## Loading a cpp custom function to compute the diagonal terms of the product of two  sparse matrices ("dgCMatrix" class)
Rcpp::sourceCpp("./Code_completed/diagonal_product.cpp")

## Function to fix the precision matrix from inla.call.object$misc$configs$config[[k]]$Q (or Qprior or Qinv)
fix.Q <- function(Q) {
  # Q: the precision matrix from the inla.call.object; inla.call.object <- inla(...) 
  d <- diag(Q)
  Q <- Q + t(Q)
  diag(Q) <- d
  return (Q)
}

diag_Mprod <- function(A,B, num_cores = 1){
  # A: a dense matrix
  # B: a sparse matrix in Compressed Sparse Column (CSC) Format
  if(class(B)=="dgCMatrix"){
    res <- parallel::mclapply(X = seq_len(B@Dim[1]), mc.cores = num_cores, FUN = function(i){
      idx_x <- ((B@p[i]+1):B@p[i+1])
      idx_r <- B@i[idx_x]+1
      return(sum(A[i,idx_r]*B@x[idx_x]))
    }) %>% do.call(what = c, .)
  } else{
    stop("B is not 'dgCMatrix' class.")
  }
  return(res)
}

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

spt_partition <- function(DF_input, mesh_inf, t_each_group = 2, n_areas = 4){
  sf_poly_int <- st_polygon(x = list(mesh_inf$loc[c(mesh_inf$segm$int$idx[,1],1),1:2]))
  id_tgroup <- rep(seq_len(length(unique(DF_input$id_time))/t_each_group), each = t_each_group*unique(table(DF_input$id_time)))
  
  if(length(unique(table(id_tgroup)))!=1 | !(unique(table(id_tgroup))[1] %>% is.integer(.))){
    stop("Partitions with different number of point per temporal group or non integer number of observations for each partition.")
  }
  
  list_clust_polygons <- list()
  id_spt_group <- c()
  
  for(i in unique(id_tgroup)){
    df_sf_points <- DF_input$geometry[DF_input$id_time %in% ((i-1)*t_each_group+1):(i*t_each_group),]
    clust_polygons <- split_poly(sf_poly = sf_poly_int, sf_obs_points = df_sf_points, n_areas = n_areas, seed = seed)
    list_clust_polygons[[i]] <- clust_polygons
    id_sp <- st_intersects(x = df_sf_points, y = clust_polygons) %>% as.numeric(.)
    id_spt_group <- c(id_spt_group,(i-1)*n_areas+id_sp)
  }
  
  DF_input$id_spt_group <- id_spt_group
  return(list(DF_input = DF_input, list_clust_polygons = list_clust_polygons))
}

# Loading data ----

DF_strong_sim <- readRDS(file = "./Code_completed/Data_strong.RDS")

# Analysing strong data ----

## Building the spatial mesh for SPDE-FEM inference ----
data(PRborder)
sf_PRborder <- st_sfc(st_polygon(x = list(PRborder))) # boundary of the sr (study region)
PR_nchull <- fm_nonconvex_hull_inla(x =  sf_PRborder, convex = -0.05)
st_PR_int_nchull <-  st_polygon(x = list(PR_nchull$loc[c(PR_nchull$idx[1:which(PR_nchull$idx[1,1]==PR_nchull$idx[,2]),1],1),1:2]))

mesh_inf <- fm_mesh_2d_inla(boundary = list(st_PR_int_nchull), max.edge = c(0.3,0.6), offset = c(-0.01,-0.1))

## Spliting the data into a train and test set

DF_polys_spt_partitions <- spt_partition(DF_input = DF_strong_sim, mesh_inf = mesh_inf, t_each_group = 2, n_areas = 4)
DF_spt_part <- DF_polys_spt_partitions$DF_input

ggplot() + 
  geom_sf(data = st_sf(cbind(IID = "A. Temporal nodes from 1 to 4", DF_polys_spt_partitions$list_clust_polygons[[1]])), mapping = aes(fill = id), alpha = 0.25) +
  geom_sf(data = st_sf(cbind(IID = "B. Temporal nodes from 5 to 8", DF_polys_spt_partitions$list_clust_polygons[[2]] %>% mutate(id = factor(as.numeric(id) + 3)))), mapping = aes(fill = id), alpha = 0.25) +
  geom_sf(data = st_sf(cbind(IID = "A. Temporal nodes from 1 to 4", DF_polys_spt_partitions$DF_input[DF_polys_spt_partitions$DF_input$id_spt_group %in% 1:3,])), mapping = aes(color = factor(id_spt_group))) +
  geom_sf(data = st_sf(cbind(IID = "B. Temporal nodes from 5 to 8", DF_polys_spt_partitions$DF_input[DF_polys_spt_partitions$DF_input$id_spt_group %in% 4:6,])), mapping = aes(color = factor(id_spt_group))) +
  labs(colour = "ID cluster", fill = "ID cluster") +
  facet_wrap(facets = ~ IID) +
  theme_bw() + theme(strip.text.x = element_text(size = 20, face = "bold", h = 0.5), text = element_text(size = 16))

ggplot() + 
  geom_sf(data = st_sf(cbind(IID = "A. Temporal nodes from 1 to 4", DF_polys_spt_partitions$list_clust_polygons[[1]])), mapping = aes(fill = id), alpha = 0.25) +
  geom_sf(data = st_sf(cbind(IID = "B. Temporal nodes from 5 to 8", DF_polys_spt_partitions$list_clust_polygons[[2]] %>% mutate(id = factor(as.numeric(id) + 4)))), mapping = aes(fill = id), alpha = 0.25) +
  geom_sf(data = st_sf(cbind(IID = "A. Temporal nodes from 1 to 4", DF_polys_spt_partitions$DF_input[DF_polys_spt_partitions$DF_input$id_spt_group %in% 1:4,])), mapping = aes(color = factor(id_spt_group))) +
  geom_sf(data = st_sf(cbind(IID = "B. Temporal nodes from 5 to 8", DF_polys_spt_partitions$DF_input[DF_polys_spt_partitions$DF_input$id_spt_group %in% 5:8,])), mapping = aes(color = factor(id_spt_group))) +
  labs(colour = "ID cluster", fill = "ID cluster") +
  facet_wrap(facets = ~ IID) +
  theme_bw() + theme(strip.text.x = element_text(size = 20, face = "bold", h = 0.5), text = element_text(size = 16))

ggplot() + 
  geom_sf(data = st_sf(cbind(IID = "A. Temporal nodes (1 & 2)", DF_polys_spt_partitions$list_clust_polygons[[1]])), mapping = aes(fill = id), alpha = 0.25) +
  geom_sf(data = st_sf(cbind(IID = "B. Temporal nodes (3 & 4)", DF_polys_spt_partitions$list_clust_polygons[[2]] %>% mutate(id = factor(as.numeric(id) + 4)))), mapping = aes(fill = id), alpha = 0.25) +
  geom_sf(data = st_sf(cbind(IID = "C. Temporal nodes (5 & 6)", DF_polys_spt_partitions$list_clust_polygons[[3]] %>% mutate(id = factor(as.numeric(id) + 8)))), mapping = aes(fill = id), alpha = 0.25) +
  geom_sf(data = st_sf(cbind(IID = "D. Temporal nodes (7 & 8)", DF_polys_spt_partitions$list_clust_polygons[[4]] %>% mutate(id = factor(as.numeric(id) + 12)))), mapping = aes(fill = id), alpha = 0.25) +
  geom_sf(data = st_sf(cbind(IID = "A. Temporal nodes (1 & 2)", DF_polys_spt_partitions$DF_input[DF_polys_spt_partitions$DF_input$id_spt_group %in% 1:4,])), mapping = aes(color = factor(id_spt_group))) +
  geom_sf(data = st_sf(cbind(IID = "B. Temporal nodes (3 & 4)", DF_polys_spt_partitions$DF_input[DF_polys_spt_partitions$DF_input$id_spt_group %in% 5:8,])), mapping = aes(color = factor(id_spt_group))) +
  geom_sf(data = st_sf(cbind(IID = "C. Temporal nodes (5 & 6)", DF_polys_spt_partitions$DF_input[DF_polys_spt_partitions$DF_input$id_spt_group %in% 9:12,])), mapping = aes(color = factor(id_spt_group))) +
  geom_sf(data = st_sf(cbind(IID = "D. Temporal nodes (7 & 8)", DF_polys_spt_partitions$DF_input[DF_polys_spt_partitions$DF_input$id_spt_group %in% 13:16,])), mapping = aes(color = factor(id_spt_group))) +
  labs(colour = "ID cluster", fill = "ID cluster") +
  facet_wrap(facets = ~ IID) +
  theme_bw() + theme(strip.text.x = element_text(size = 20, face = "bold", h = 0.5), text = element_text(size = 16))


list_DF_cv <- list()

for(id_cv in seq_along(unique(DF_spt_part$id_spt_group))){
  
  idx_train <- which(DF_spt_part$id_spt_group != id_cv)
  idx_test <- which(DF_spt_part$id_spt_group == id_cv)
  DF_train <- DF_strong_sim; DF_train[idx_test,"y"] <- NA
  DF_test <- DF_strong_sim; DF_test[idx_train,"y"] <- NA
  
  ## Construction of the SPDE-FEM effect----
  
  max_dist <- lapply(fm_bbox(mesh_inf), "[") %>% unlist(.) %>% matrix(data = ., ncol = 2, byrow = FALSE) %>% dist(.)
  spde_spt <- inla.spde2.pcmatern(mesh = mesh_inf, alpha = 2, prior.range = c(max_dist/5, 0.5), prior.sigma = c(1,0.5), constr = TRUE)
  spde_spt_idx <- inla.spde.make.index(name = "spt", n.spde = spde_spt$n.spde, n.group = length(unique(DF_strong_sim$id_time)))
  A_spt_inf <- inla.spde.make.A(mesh = mesh_inf, loc = st_coordinates(DF_strong_sim), group = DF_strong_sim$id_time, n.group = length(unique(DF_strong_sim$id_time)))
  
  x2_group <- inla.group(x = DF_strong_sim$X2, n = 50)
  x3_group <- inla.group(x = DF_strong_sim$X3, n = 50)
  
  # We can analyse the covariates as fixed effects or as random effects, using a "non-linear" structure as a Random Walk of second order.
  inf_stk <- inla.stack(data = list(y = DF_train$y),
                        A = list(1, A_spt_inf),
                        effects = list(
                          list(beta0 = rep(1, nrow(DF_train)),
                               f_x1 = as.numeric(DF_train$X1),
                               # f_x2 = x2_group,
                               # f_x3 = x3_group,
                               f_x2 = DF_strong_sim$X2,
                               f_x3 = DF_strong_sim$X3
                          ),
                          spde_spt_idx
                        ),
                        remove.unused = FALSE,
                        compress = FALSE,
                        tag = "inf_stk")
  
  spt_formula <- y ~ -1 + beta0 + f(f_x1, model = "iid", constr = TRUE) + 
    f_x2 + f_x3 +
    # f(f_x2, model = "rw2", constr = TRUE) + f(f_x3, model = "rw2", constr = TRUE) +
    f(spt, model = spde_spt, group = spt.group, control.group = list(model = "ar1"))
  
  rinla <- inla(data = inla.stack.data(stack = inf_stk),
                family = "gaussian",
                formula = spt_formula,
                control.predictor = list(A = inla.stack.A(inf_stk)),
                control.compute = list(config = TRUE),
                num.threads = 24,
                verbose = FALSE)
  
  rinla_orig <- rinla # saving the first analysis
  mu_new <- rinla$misc$configs$config[[1]]$improved.mean
  Q_new <- fix.Q(rinla$misc$configs$config[[1]]$Q)
  
  # Implementing Algorithm 1 ----
  
  ## Defining the measures for predictive performance ----
  
  DFrmse_train <- data.frame(INLA = NA, INLA.RFnoun = NA, INLA.RFun = NA)
  DFrmse_test <- data.frame(INLA = NA, INLA.RFnoun = NA, INLA.RFun = NA)
  
  DFmae_train <- data.frame(INLA = NA, INLA.RFnoun = NA, INLA.RFun = NA)
  DFmae_test <- data.frame(INLA = NA, INLA.RFnoun = NA, INLA.RFun = NA)
  
  DFcp_train <- data.frame(INLA = NA, INLA.RFnoun = NA, INLA.RFun = NA)
  DFcp_test <- data.frame(INLA = NA, INLA.RFnoun = NA, INLA.RFun = NA)
  
  DFaiw_train <- data.frame(INLA = NA, INLA.RFnoun = NA, INLA.RFun = NA)
  DFaiw_test <- data.frame(INLA = NA, INLA.RFnoun = NA, INLA.RFun = NA)
  
  ## Computing the measure for the train/test sets of the standard INLA approach ----
  
  DFrmse_train$INLA <- sqrt(mean((DF_strong_sim$y[idx_train] - rinla_orig$summary.fitted.values[idx_train,"mean"])**2))
  DFrmse_test$INLA <- sqrt(mean((DF_strong_sim$y[idx_test] - rinla_orig$summary.fitted.values[idx_test,"mean"])**2))
  
  DFmae_train$INLA <- mean(abs(DF_strong_sim$y[idx_train] - rinla_orig$summary.fitted.values[idx_train,"mean"]))
  DFmae_test$INLA <- mean(abs(DF_strong_sim$y[idx_test] - rinla_orig$summary.fitted.values[idx_test,"mean"]))
  
  DFcp_train$INLA <- mean(DF_strong_sim$y[idx_train] >= rinla_orig$summary.fitted.values[idx_train,"0.025quant"] & DF_strong_sim$y[idx_train] <= rinla_orig$summary.fitted.values[idx_train,"0.975quant"])
  DFcp_test$INLA <- mean(DF_strong_sim$y[idx_test] >= rinla_orig$summary.fitted.values[idx_test,"0.025quant"] & DF_strong_sim$y[idx_test] <= rinla_orig$summary.fitted.values[idx_test,"0.975quant"])
  
  DFaiw_train$INLA <- mean(rinla_orig$summary.fitted.values[idx_train,"0.975quant"] - rinla_orig$summary.fitted.values[idx_train,"0.025quant"])
  DFaiw_test$INLA <- mean(rinla_orig$summary.fitted.values[idx_test,"0.975quant"] - rinla_orig$summary.fitted.values[idx_test,"0.025quant"])
  
  ## Setting the configuration for the INLA-RF Algorithm ----
  
  ysim <- DF_strong_sim$y
  e <- y.e_hat <- rep(0, nrow(DF_strong_sim)) # Offset of the residual estimations from the RF. Initialized with y.e_hat = 0.
  offx <- rep(0, times = nrow(DF_strong_sim))
  KLD_GMRF_den <- 10 # Initial value for the KLD
  i <- 0 # Number of loops = i
  
  KLD_comp <- "fast_full"
  verbose <- TRUE
  
  t1 <- Sys.time()
  while(KLD_GMRF_den > 1E-2){ # Using the KLD as condition
    if(i>0){
      mu_old <- mu_new
      Q_old <- Q_new
      list_control_mode <- list(x = c(rep(0, rinla$misc$configs$mnpred), mu_old), theta = rinla$mode$theta, restart = TRUE, fixed = FALSE)
      
      inf_stk_loop <- inf_stk
      inf_stk_loop$data$data[idx_train,] <- inf_stk$data$data[idx_train,] - y.e_hat[idx_train]
      
      # offx[idx_train] <- offx[idx_train] + y.e_hat[idx_train]
      rinla <- inla(data = inla.stack.data(stack = inf_stk_loop), 
                    family = "gaussian",
                    formula = spt_formula,
                    # offset = offx,
                    control.predictor = list(A = inla.stack.A(inf_stk_loop)),
                    control.compute = list(config = TRUE),
                    control.mode = list_control_mode,
                    num.threads = 24,
                    verbose = FALSE)
      
      mu_new <- rinla$misc$configs$config[[1]]$improved.mean
      Q_new <- fix.Q(rinla$misc$configs$config[[1]]$Q)
      Q_inv_new <- inla.qinv(Q_new) %>% as(., Class = "CsparseMatrix") # The original Compress Sparse format is COO. Therefore, we need to convert it to CSC, aka "CsparseMatrix", for faster computational performance.
      
      if(KLD_comp == "fast_full"){ # Computing the KLD suing the whole matrix
        # t1 <- Sys.time()
        idx_reord <- inla.qreordering(graph = Q_new)
        L_ir <- chol(Q_new[idx_reord$ireordering,idx_reord$ireordering])
        L_ir_old <- chol(Q_old[idx_reord$ireordering,idx_reord$ireordering])
        KLD_GMRF_den <- 1/2*(sum(diagonal_producto_csc(A = Q_inv_new, B = Q_old)) - nrow(Q_new) + drop((mu_new-mu_old) %*% Q_new %*% cbind(mu_new-mu_old)) + 2*sum(log(diag(L_ir))) - 2*sum(log(diag(L_ir_old))))/nrow(Q_new)
        # t2 <- Sys.time()
      } else if(KLD_comp == "diag"){
        # t1 <- Sys.time()
        idx_reord <- inla.qreordering(graph = Q_new)
        L_ir <- chol(Q_new[idx_reord$ireordering,idx_reord$ireordering])
        L_ir_old <- chol(Q_old[idx_reord$ireordering,idx_reord$ireordering])
        KLD_GMRF_den <- 1/2*(sum(diag(Q_inv_new) * diag(Q_old)) - nrow(Q_new) + drop(t(mu_new-mu_old) %*% Diagonal(x=diag(Q_new)) %*% cbind(mu_new-mu_old)) + 2*sum(log(diag(L_ir))) - 2*sum(log(diag(L_ir_old))))/nrow(Q_new)
        # t2 <- Sys.time()
      } else{ # Slow, pretty slow
        KLD_GMRF_den <- 1/2*(sum(diag(solve(Q_new)%*%Q_old)) - nrow(Q_new) + drop((mu_new-mu_old) %*% Q_new %*% cbind(mu_new-mu_old)) + determinant(Q_new, logarithm = TRUE)$modulus - determinant(Q_old, logarithm = TRUE)$modulus)/nrow(Q_new)
      }
    }
    
    e[idx_train] <- ysim[idx_train] - rinla$summary.fitted.values[idx_train,"mean"]
    sp_coord <- st_coordinates(DF_strong_sim$geometry)
    df_e <- data.frame(e = e, id_time = DF_strong_sim$id_time, X1 = DF_strong_sim$X1, X2 = DF_strong_sim$X2, X3 = DF_strong_sim$X3, sp_x = sp_coord[,1], sp_y = sp_coord[,2])
    
    fit_rf_e <- 
      ranger(formula = e ~ id_time + X1 + X2 + X3 + sp_x + sp_y,
             data = df_e[idx_train,],
             importance = "none",
             replace = FALSE,
             seed = seed,
             oob.error = TRUE)
    
    y.e_hat[idx_train] <- fit_rf_e$predictions
    y.e_hat[idx_test] <- predict(fit_rf_e, data = data.frame(df_e[idx_test,]))$predictions
    
    i = i + 1 # Increasing in 1 the value of the index for the loop
    # DFrmse_A1_train[i,] <- c(sqrt(mean((ysim[idx_train] - rinla$summary.fitted.values[idx_train,"mean"])**2)), sqrt(mean((ysim[idx_train] - rinla$summary.fitted.values[idx_train,"mean"] - y.e_hat[idx_train])**2)))
    # DFrmse_A1_test[i,] <- c(sqrt(mean((ysim[idx_test] - rinla$summary.fitted.values[idx_test,"mean"])**2)), sqrt(mean((ysim[idx_test] - rinla$summary.fitted.values[idx_test,"mean"] - y.e_hat[idx_test])**2)))
    
    if(verbose && i > 0){
      sprintf("KLD: %.4f. \n", KLD_GMRF_den) %>% cat(.)
      sprintf("Iteration: %i. \n", i-1) %>% cat(.)
      cat("---------------------------------------------- \n")
    }
  }
  t2 <- Sys.time()
  difftime(t2, t1)
  
  ## Computing the measure for the train/test sets of the SPDE-RF (with no uncertainty transferring) approach ----
  
  DFrmse_train$INLA.RFnoun <- sqrt(mean((DF_strong_sim$y[idx_train] - rinla$summary.fitted.values[idx_train,"mean"] - y.e_hat[idx_train])**2))
  DFrmse_test$INLA.RFnoun <- sqrt(mean((DF_strong_sim$y[idx_test] - rinla$summary.fitted.values[idx_test,"mean"] - y.e_hat[idx_test])**2))
  
  DFmae_train$INLA.RFnoun <- mean(abs(DF_strong_sim$y[idx_train] - rinla$summary.fitted.values[idx_train,"mean"] - y.e_hat[idx_train]))
  DFmae_test$INLA.RFnoun <- mean(abs(DF_strong_sim$y[idx_test] - rinla$summary.fitted.values[idx_test,"mean"] - y.e_hat[idx_test]))
  
  DFcp_train$INLA.RFnoun <- mean(DF_strong_sim$y[idx_train] >= (rinla$summary.fitted.values[idx_train,"0.025quant"] + y.e_hat[idx_train]) & DF_strong_sim$y[idx_train] <= (rinla$summary.fitted.values[idx_train,"0.975quant"] + y.e_hat[idx_train]))
  DFcp_test$INLA.RFnoun <- mean(DF_strong_sim$y[idx_test] >= (rinla$summary.fitted.values[idx_test,"0.025quant"] + y.e_hat[idx_test]) & DF_strong_sim$y[idx_test] <= (rinla$summary.fitted.values[idx_test,"0.975quant"] + y.e_hat[idx_test]))
  
  DFaiw_train$INLA.RFnoun <- mean(rinla$summary.fitted.values[idx_train,"0.975quant"] - rinla$summary.fitted.values[idx_train,"0.025quant"]) # Neutralized the effect of the offset for the interval calculation
  DFaiw_test$INLA.RFnoun <- mean(rinla$summary.fitted.values[idx_test,"0.975quant"] - rinla$summary.fitted.values[idx_test,"0.025quant"]) # Neutralized the effect of the offset for the interval calculation
  
  # Implementing Algorithm 2 ----
  
  ## Setting the configuration for the INLA-RF Algorithm ----
  rinla <- rinla_orig
  ysim <- DF_strong_sim$y
  e <- y.e_hat <- rep(0, nrow(DF_strong_sim)) # Offset of the residual estimations from the RF. Initialized with y.e_hat = 0.
  offx <- rep(0, times = nrow(DF_strong_sim))
  off_variance <- 0
  KLD_GMRF_den <- 10 # Initial value for the KLD
  i <- 0 # Number of loops = i
  
  KLD_comp <- "fast_full"
  verbose <- TRUE
  
  t1 <- Sys.time()
  while(KLD_GMRF_den > 1E-2){ # Using the KLD as condition
    if(i>0){
      mu_old <- mu_new
      Q_old <- Q_new
      list_control_mode <- list(x = c(rep(0, rinla$misc$configs$mnpred), mu_old), theta = rinla$mode$theta, restart = TRUE, fixed = FALSE)
      
      off_variance <- off_variance + fit_rf_e$prediction.error
      
      inf_stk_loop <- inf_stk
      inf_stk_loop$data$data[idx_train,] <- inf_stk$data$data[idx_train,] - y.e_hat[idx_train]
      
      # offx[idx_train] <- offx[idx_train] + y.e_hat[idx_train]
      rinla <- inla(data = inla.stack.data(stack = inf_stk_loop), 
                    family = "gaussian",
                    formula = spt_formula,
                    # offset = offx,
                    control.predictor = list(A = inla.stack.A(inf_stk_loop)),
                    control.compute = list(config = TRUE),
                    control.mode = list_control_mode,
                    control.family = list(hyper = list(precoffset = list(initial = -log(off_variance), fixed = TRUE))),
                    num.threads = 24,
                    verbose = FALSE)
      
      mu_new <- rinla$misc$configs$config[[1]]$improved.mean
      Q_new <- fix.Q(rinla$misc$configs$config[[1]]$Q)
      Q_inv_new <- inla.qinv(Q_new) %>% as(., Class = "CsparseMatrix") # The original Compress Sparse format is COO. Therefore, we need to convert it to CSC, aka "CsparseMatrix", for faster computational performance.
      
      if(KLD_comp == "fast_full"){ # Computing the KLD suing the whole matrix
        # t1 <- Sys.time()
        idx_reord <- inla.qreordering(graph = Q_new)
        L_ir <- chol(Q_new[idx_reord$ireordering,idx_reord$ireordering])
        L_ir_old <- chol(Q_old[idx_reord$ireordering,idx_reord$ireordering])
        KLD_GMRF_den <- 1/2*(sum(diagonal_producto_csc(A = Q_inv_new, B = Q_old)) - nrow(Q_new) + drop((mu_new-mu_old) %*% Q_new %*% cbind(mu_new-mu_old)) + 2*sum(log(diag(L_ir))) - 2*sum(log(diag(L_ir_old))))/nrow(Q_new)
        # t2 <- Sys.time()
      } else if(KLD_comp == "diag"){
        # t1 <- Sys.time()
        idx_reord <- inla.qreordering(graph = Q_new)
        L_ir <- chol(Q_new[idx_reord$ireordering,idx_reord$ireordering])
        L_ir_old <- chol(Q_old[idx_reord$ireordering,idx_reord$ireordering])
        KLD_GMRF_den <- 1/2*(sum(diag(Q_inv_new) * diag(Q_old)) - nrow(Q_new) + drop(t(mu_new-mu_old) %*% Diagonal(x=diag(Q_new)) %*% cbind(mu_new-mu_old)) + 2*sum(log(diag(L_ir))) - 2*sum(log(diag(L_ir_old))))/nrow(Q_new)
        # t2 <- Sys.time()
      } else{ # Slow, pretty slow, ...
        KLD_GMRF_den <- 1/2*(sum(diag(solve(Q_new)%*%Q_old)) - nrow(Q_new) + drop((mu_new-mu_old) %*% Q_new %*% cbind(mu_new-mu_old)) + determinant(Q_new, logarithm = TRUE)$modulus - determinant(Q_old, logarithm = TRUE)$modulus)/nrow(Q_new)
      }
    }
    
    e[idx_train] <- ysim[idx_train] - rinla$summary.fitted.values[idx_train,"mean"]
    sp_coord <- st_coordinates(DF_strong_sim$geometry)
    df_e <- data.frame(e = e, id_time = DF_strong_sim$id_time, X1 = DF_strong_sim$X1, X2 = DF_strong_sim$X2, X3 = DF_strong_sim$X3, sp_x = sp_coord[,1], sp_y = sp_coord[,2])
    
    fit_rf_e <- 
      ranger(formula = e ~ id_time + X1 + X2 + X3 + sp_x + sp_y,
             data = df_e[idx_train,],
             importance = "none",
             replace = FALSE,
             seed = seed,
             oob.error = TRUE)
    
    y.e_hat[idx_train] <- fit_rf_e$predictions
    y.e_hat[idx_test] <- predict(fit_rf_e, data = data.frame(df_e[idx_test,]))$predictions
    
    i = i + 1 # Increasing in 1 the value of the index for the loop
    # DFrmse_A2_train[i,] <- c(sqrt(mean((ysim[idx_train] - rinla$summary.fitted.values[idx_train,"mean"])**2)), sqrt(mean((ysim[idx_train] - rinla$summary.fitted.values[idx_train,"mean"] - y.e_hat[idx_train])**2)))
    # DFrmse_A2_test[i,] <- c(sqrt(mean((ysim[idx_test] - rinla$summary.fitted.values[idx_test,"mean"])**2)), sqrt(mean((ysim[idx_test] - rinla$summary.fitted.values[idx_test,"mean"] - y.e_hat[idx_test])**2)))
    
    if(verbose && i > 1){
      sprintf("KLD: %.4f. \n", KLD_GMRF_den) %>% cat(.)
      sprintf("Iteration: %i. \n", i-1) %>% cat(.)
      cat("---------------------------------------------- \n")
    }
  }
  t2 <- Sys.time()
  difftime(t2, t1)
  
  ## Computing the measure for the train/test sets of the SPDE-RF (with no uncertainty transferring) approach ----
  
  DFrmse_train$INLA.RFun <- sqrt(mean((DF_strong_sim$y[idx_train] - rinla$summary.fitted.values[idx_train,"mean"] - y.e_hat[idx_train])**2))
  DFrmse_test$INLA.RFun <- sqrt(mean((DF_strong_sim$y[idx_test] - rinla$summary.fitted.values[idx_test,"mean"] - y.e_hat[idx_test])**2))
  
  DFmae_train$INLA.RFun <- mean(abs(DF_strong_sim$y[idx_train] - rinla$summary.fitted.values[idx_train,"mean"] - y.e_hat[idx_train]))
  DFmae_test$INLA.RFun <- mean(abs(DF_strong_sim$y[idx_test] - rinla$summary.fitted.values[idx_test,"mean"] - y.e_hat[idx_test]))
  
  inlarf_train_q13 <- mclapply(X = idx_train, mc.cores = 1, FUN = function(idx){
    DF_mean_sd <- data.frame(mean = rinla$summary.fitted.values[idx,"mean"] + y.e_hat[idx], sd = sqrt(rinla$summary.fitted.values[idx,"sd"]**2 + fit_rf_e$prediction.error))
    return(data.frame(idx = idx, q1 = qnorm(p = 0.025, mean = DF_mean_sd$mean, sd = DF_mean_sd$sd), q3 = qnorm(p = 0.975, mean = DF_mean_sd$mean, sd = DF_mean_sd$sd)))
  }) %>% do.call(what = rbind, .)
  
  inlarf_test_q13 <- mclapply(X = idx_test, mc.cores = 1, FUN = function(idx){
    DF_mean_sd <- data.frame(mean = rinla$summary.fitted.values[idx,"mean"] + y.e_hat[idx], sd = sqrt(rinla$summary.fitted.values[idx,"sd"]**2 + fit_rf_e$prediction.error))
    return(data.frame(idx = idx, q1 = qnorm(p = 0.025, mean = DF_mean_sd$mean, sd = DF_mean_sd$sd), q3 = qnorm(p = 0.975, mean = DF_mean_sd$mean, sd = DF_mean_sd$sd)))
  }) %>% do.call(what = rbind, .)
  
  DFcp_train$INLA.RFun <- mean(DF_strong_sim$y[idx_train] >= (inlarf_train_q13$q1) & DF_strong_sim$y[idx_train] <= (inlarf_train_q13$q3))
  DFcp_test$INLA.RFun <- mean(DF_strong_sim$y[idx_test] >= (inlarf_test_q13$q1) & DF_strong_sim$y[idx_test] <= (inlarf_test_q13$q3))
  
  DFaiw_train$INLA.RFun <- mean(inlarf_train_q13$q3 - inlarf_train_q13$q1) # Neutralized the effect of the offset for the interval calculation
  DFaiw_test$INLA.RFun <- mean(inlarf_test_q13$q3 - inlarf_test_q13$q1) # Neutralized the effect of the offset for the interval calculation
  
  list_DF_cv[[id_cv]] <- list(DFrmse = list(DFrmse_train, DFrmse_test), DFmae = list(DFmae_train, DFmae_test), DFcp = list(DFcp_train, DFcp_test), DFaiw = list(DFaiw_train, DFaiw_test))
}

idx <- 6
cat(
  paste0(
    "RMSE & ",
    "($", round(list_DF_cv[[idx]]$DFrmse[[1]][[1]], digits = 2), "$, $", round(list_DF_cv[[idx]]$DFrmse[[2]][[1]], digits = 2), "$) ", "& ",
    "($", round(list_DF_cv[[idx]]$DFrmse[[1]][[2]], digits = 2), "$, $", round(list_DF_cv[[idx]]$DFrmse[[2]][[2]], digits = 2), "$) ", "& ",
    "($", round(list_DF_cv[[idx]]$DFrmse[[1]][[3]], digits = 2), "$, $", round(list_DF_cv[[idx]]$DFrmse[[2]][[3]], digits = 2), "$) \\\\ \n",
    "MAE & ",
    "($", round(list_DF_cv[[idx]]$DFmae[[1]][[1]], digits = 2), "$, $", round(list_DF_cv[[idx]]$DFmae[[2]][[1]], digits = 2), "$) ", "& ",
    "($", round(list_DF_cv[[idx]]$DFmae[[1]][[2]], digits = 2), "$, $", round(list_DF_cv[[idx]]$DFmae[[2]][[2]], digits = 2), "$) ", "& ",
    "($", round(list_DF_cv[[idx]]$DFmae[[1]][[3]], digits = 2), "$, $", round(list_DF_cv[[idx]]$DFmae[[2]][[3]], digits = 2), "$) \\\\ \n",
    "CP & ",
    "($", round(list_DF_cv[[idx]]$DFcp[[1]][[1]], digits = 2), "$, $", round(list_DF_cv[[idx]]$DFcp[[2]][[1]], digits = 2), "$) ", "& ",
    "($", round(list_DF_cv[[idx]]$DFcp[[1]][[2]], digits = 2), "$, $", round(list_DF_cv[[idx]]$DFcp[[2]][[2]], digits = 2), "$) ", "& ",
    "($", round(list_DF_cv[[idx]]$DFcp[[1]][[3]], digits = 2), "$, $", round(list_DF_cv[[idx]]$DFcp[[2]][[3]], digits = 2), "$) \\\\ \n",
    "AIW & ",
    "($", round(list_DF_cv[[idx]]$DFaiw[[1]][[1]], digits = 2), "$, $", round(list_DF_cv[[idx]]$DFaiw[[2]][[1]], digits = 2), "$) ", "& ",
    "($", round(list_DF_cv[[idx]]$DFaiw[[1]][[2]], digits = 2), "$, $", round(list_DF_cv[[idx]]$DFaiw[[2]][[2]], digits = 2), "$) ", "& ",
    "($", round(list_DF_cv[[idx]]$DFaiw[[1]][[3]], digits = 2), "$, $", round(list_DF_cv[[idx]]$DFaiw[[2]][[3]], digits = 2), "$) \\\\"
  ))


mat_mean <- matrix(data = 0, ncol = 6, nrow = 4)
ks <- rep(1:2, times = 3)
km <- rep(1:3, each = 2) 
for(ii in 1:length(list_DF_cv)){
  for(k in 1:ncol(mat_mean)){
    for(l in 1:4){
      mat_mean[l,k] <- mat_mean[l,k] + list_DF_cv[[ii]][[l]][[ks[k]]][[km[k]]]
    }
  }
}
mat_mean/6
