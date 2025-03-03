# Title: runme_SPDE_strong_lowRankNodeCorrection.R
# Author: Mario Figueira
# Date: 2025-02-22
# Description: Code combining INLA and RF to perform a low-rank correction of some nodes for a chosen structure of the model.
#              In this case some nodes of the spatio-temporal structure.

# Last update: 2025-02-26

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
library(ggridges)
library(ggtext)
library(dplyr)
library(sf)

library(ranger)

seed <- 1234 # set a seed for reproducibility
set.seed(seed = seed)

# Custom functions ----

## Loading a cpp custom function to compute the diagonal terms of the product of two  sparse matrices ("dgCMatrix" class)
Rcpp::sourceCpp("./Scripts/diagonal_product.cpp")

## Function to fix the precision matrix from inla.call.object$misc$configs$config[[k]]$Q (or Qprior or Qinv)
fix.Q <- function(Q) {
  # Q: the precision matrix from the inla.call.object; inla.call.object <- inla(...) 
  d <- diag(Q)
  Q <- Q + t(Q)
  diag(Q) <- d
  return (Q)
}

# A function to compute only the diagonal elements of a product of two matrices (A and B) using parallelized code in R (it is better to use the C++ function)
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

# Loading data ----

DF_strong_sim <- readRDS(file = "./Scripts/Data_strong.RDS")
list_spt_sim <- readRDS(file = "./Scripts/list_spt_sim.RDS")

## Spliting the data into a train and test set

idx_train <- sample(x = seq_len(nrow(DF_strong_sim)), size = nrow(DF_strong_sim)*0.8) %>% sort(.)
idx_test <- setdiff(seq_len(nrow(DF_strong_sim)), idx_train)
DF_train <- DF_strong_sim; DF_train[idx_test,"y"] <- NA
DF_test <- DF_strong_sim; DF_test[idx_train,"y"] <- NA

# Analysing strong data ----

## Building the spatial mesh for SPDE-FEM inference ----
data(PRborder)
sf_PRborder <- st_sfc(st_polygon(x = list(PRborder))) # boundary of the sr (study region)
PR_nchull <- fm_nonconvex_hull_inla(x =  sf_PRborder, convex = -0.05)
st_PR_int_nchull <-  st_polygon(x = list(PR_nchull$loc[c(PR_nchull$idx[1:which(PR_nchull$idx[1,1]==PR_nchull$idx[,2]),1],1),1:2]))

mesh_inf <- fm_mesh_2d_inla(boundary = list(st_PR_int_nchull), max.edge = c(0.3,0.6), offset = c(-0.01,-0.1))

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
              num.threads = 48,
              verbose = FALSE)

rinla_orig <- rinla # saving the first analysis
mu_new <- rinla$misc$configs$config[[1]]$improved.mean
Q_new <- fix.Q(rinla$misc$configs$config[[1]]$Q)

# Nodes to be corrected and true values ----

nstress <- 100 # The number of stress nodes to correct. Here we will corect only nodes related to the spatio-temporal effect
idx_linked <- which(apply(X = A_spt_inf[idx_train,], MARGIN = 2, FUN = sum) != 0)

## Choosing the nodes to correct ----
# Choosing the nodes to correct by higher variance or related to higher RMSE its different when we have a projection from the latent field to the linear predictor.
# Given that the RMSE is computed for the observations, which implies that it is related to (\eta) and not directly linked to the latent field (x). 
# Meanwhile, variance is linked directly to the latent field nodes.

idx_linked_spt <- order(rinla$summary.random$spt$sd[idx_linked], decreasing = TRUE)[1:nstress]
selected_spt_nodes <- rinla$summary.random$spt[idx_linked[idx_linked_spt],]
idx_obs.linked <- which(apply(X = A_spt_inf[idx_train,selected_spt_nodes$ID+1], MARGIN = 1, FUN = sum) != 0)


selected_mesh_nodes <- (selected_spt_nodes$ID + 1) %% mesh_inf$n
group_mesh <- ceiling((selected_spt_nodes$ID + 1) / mesh_inf$n)
if(any(selected_mesh_nodes == 0)){selected_mesh_nodes[selected_mesh_nodes==0] <- mesh_inf$n}

sim_node_values <- mclapply(X = unique(group_mesh), mc.cores = 10, FUN = function(i){
  res <- fm_basis(x = list_spt_sim$mesh, loc = mesh_inf$loc[selected_mesh_nodes[group_mesh == i],1:2]) %*% list_spt_sim$spt_sim[((i-1)*list_spt_sim$mesh$n+1):(i*list_spt_sim$mesh$n)] %>% drop(.)
  return(res)
}) %>% do.call(what = c, .)

n_group <- 10
nodes_sel_marg <- mclapply(X = 1:nstress, mc.cores = 1, FUN = function(i){cbind(data.frame(rinla$marginals.random$spt[[selected_spt_nodes$ID[i]+1]]), group = i%%n_group + if(i%%n_group==0){n_group}else{0}, height = 0, ID_group = ceiling(i/n_group))}) %>% 
  do.call(what = rbind, .)

ngroup_idx_to_plot <- 1:10
height_error_bar <- 1:n_group
scale_ridge <- 1
ggplot_post.correction <- ggplot() +
  geom_ridgeline_gradient(data = nodes_sel_marg[nodes_sel_marg$ID_group %in% ngroup_idx_to_plot,], mapping = aes(x = x, y = height + group, group = group, height = y, scale = scale_ridge, fill = y), linewidth = 1, colour = "blue", alpha = 0.5) +
  geom_linerange(data = data.frame(x = sim_node_values,
                                   y_max = height_error_bar+0.85, y_min = height_error_bar,
                                   ID_group = ceiling(seq_len(last(ngroup_idx_to_plot)*n_group)/n_group)),
                 mapping = aes(x = x, ymin = y_min, ymax = y_max),
                 linewidth = 1, colour = "black") +
  scale_y_continuous(name = "Posterior distribution", breaks = 1:10, 
                     labels = 1:10) +
  scale_fill_viridis_c(name = "Density", option = "mako") +
  facet_wrap(facets = ~ ID_group, ncol = 5) +
  theme_bw()

## Setting the configuration for the INLA-RF Algorithm ----

### Setting the configuration of the inla.stack and the new model formula ----

# selected_spt_nodes: Nodes selected from the spatio-temporal structure
A_u.iid <- matrix(data = 0, nrow = nrow(DF_strong_sim), ncol = nstress) %>% as(., "TsparseMatrix")
for(k in (selected_spt_nodes$ID+1)){
  A_idx_train <- idx_train[which(A_spt_inf[idx_train,k] != 0)]
  i <- A_idx_train - 1 # 'sparseMatrix-class' uses the python indexx notation 0:(length-1) 
  j <- rep(which(k == (selected_spt_nodes$ID + 1)) - 1, length(i)) # 'sparseMatrix-class' uses the python indexx notation 0:(length-1)
  A_u.iid@i <- c(A_u.iid@i, as.integer(i))
  A_u.iid@j <- c(A_u.iid@j, as.integer(j))
  A_u.iid@x <- c(A_u.iid@x, A_spt_inf[A_idx_train, k])
}

spt_formula_loop <- y ~ -1 + offset(offx) + beta0 + f(f_x1, model = "iid", constr = TRUE) + 
  f_x2 + f_x3 +
  # f(f_x2, model = "rw2", constr = TRUE) + f(f_x3, model = "rw2", constr = TRUE) +
  f(spt, model = spde_spt, group = spt.group, control.group = list(model = "ar1")) +
  f(u.iid, model = "iid", hyper = list(prec = list(initial = log(tau.iid), fixed = TRUE)))

inf_stk_loop <- inla.stack(data = list(y = DF_train$y),
                           A = list(1, A_spt_inf, A_u.iid),
                           effects = list(
                             list(beta0 = rep(1, nrow(DF_train)),
                                  f_x1 = as.numeric(DF_train$X1),
                                  # f_x2 = x2_group,
                                  # f_x3 = x3_group,
                                  f_x2 = DF_strong_sim$X2,
                                  f_x3 = DF_strong_sim$X3
                             ),
                             spde_spt_idx,
                             list(u.iid = 1:nstress)
                           ),
                           remove.unused = FALSE,
                           compress = FALSE,
                           tag = "inf_stk")

### Setting the initial configuration of the algorithm and offsets ----

ysim <- DF_strong_sim$y
e <- y.e_hat <- rep(0, nrow(DF_strong_sim)) # Offset of the residual estimations from the RF. Initialized with y.e_hat = 0.
# offx <- rep(0, times = nrow(DF_strong_sim))
offx <- rep(0, times = rinla_orig$misc$configs$npred + nstress)
# off_variance <- 0
KLD_GMRF_den <- 10 # Initial value for the KLD
i <- 0 # Number of loops = i

DFrmse_train <- data.frame(INLA = NA, INLA.RF = NA)
DFrmse_test <- data.frame(INLA = NA, INLA.RF = NA)
KLD_comp <- "fast_full"
verbose <- TRUE

t1 <- Sys.time()
while(KLD_GMRF_den > 1E-2){ # Using the KLD as condition
  if(i>0){
    
    mu_old <- mu_new
    Q_old <- Q_new
    
    if(i == 1){
      tau.iid <- as.numeric(fit_rf_e$prediction.error)**(-1)
      # mnpred = sum(dim(inf_stk$A)) = size(y_obs) + size(letent_field), without any reduction.
      list_control_mode <- list(x = c(rep(0, rinla$misc$configs$mpred + dim(inf_stk_loop$A)[2]), rinla$misc$configs$config[[1]]$improved.mean[1:(length(rinla$misc$configs$config[[1]]$improved.mean)-nrow(rinla$summary.fixed))], rep(0, times = nstress), rev(length(rinla$misc$configs$config[[1]]$improved.mean) - 0:(nrow(rinla$summary.fixed)-1))), 
                                theta = rinla$mode$theta, restart = TRUE, fixed = FALSE)
    } else{
      tau.iid <- as.numeric(fit_rf_e$prediction.error)**(-1)
      # mnpred = sum(dim(inf_stk$A)) = size(y_obs) + size(letent_field)
      list_control_mode <- list(x = c(rep(0, sum(dim(inf_stk_loop$A))), rinla$misc$configs$config[[1]]$improved.mean),
                                theta = rinla$mode$theta, restart = TRUE, fixed = FALSE)
    }
    
    offx[(length(offx) - nrow(rinla$summary.fixed) - nstress) + 1:nstress] <- y.e_hat_spt
    
    # offx[idx_train] <- offx[idx_train] + y.e_hat[idx_train]
    rinla <- inla(data = inla.stack.data(stack = inf_stk_loop), 
                  family = "gaussian",
                  formula = spt_formula_loop,
                  # offset = offx,
                  control.predictor = list(A = inla.stack.A(inf_stk_loop)),
                  control.compute = list(config = TRUE),
                  control.mode = list_control_mode,
                  num.threads = 48,
                  verbose = FALSE)
    
    size_imprmean <- length(rinla$misc$configs$config[[1]]$improved.mean)
    size_fixed <- nrow(rinla$summary.fixed)
    idx_latent_orig <- setdiff(1:size_imprmean, (size_imprmean-nstress-size_fixed+1):(size_imprmean-size_fixed))
    mu_new <- rinla$misc$configs$config[[1]]$improved.mean[idx_latent_orig]
    Q_new <- fix.Q(rinla$misc$configs$config[[1]]$Q[idx_latent_orig, idx_latent_orig])
    Q_inv_new <- inla.qinv(Q_new) %>% as(., Class = "CsparseMatrix") # The original Compress Sparse format is COO. Therefore, we need to convert it to CSC, aka "CsparseMatrix", for faster computational performance.
    
    if(KLD_comp == "fast_full"){# Computing the KLD using the whole Q precision matrix (using the whole information from both GMRFs)
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
    } else{ # Slow, pretty slow, ... (we don't need to compute the KLD this way, as the fast_full is the computational right manner of compute the KLD for huge GMRF)
      KLD_GMRF_den <- 1/2*(sum(diag(solve(Q_new)%*%Q_old)) - nrow(Q_new) + drop((mu_new-mu_old) %*% Q_new %*% cbind(mu_new-mu_old)) + determinant(Q_new, logarithm = TRUE)$modulus - determinant(Q_old, logarithm = TRUE)$modulus)/nrow(Q_new)
    }
  }
  
  e[idx_train] <- ysim[idx_train] - rinla$summary.fitted.values[idx_train,"mean"]
  sp_coord <- st_coordinates(DF_strong_sim$geometry)
  df_e <- data.frame(e = e, id_time = DF_strong_sim$id_time, X1 = DF_strong_sim$X1, X2 = DF_strong_sim$X2, X3 = DF_strong_sim$X3, sp_x = sp_coord[,1], sp_y = sp_coord[,2])
  
  fit_rf_e <- 
    ranger(formula = e ~ id_time + sp_x + sp_y,
           data = df_e[idx_train,],
           importance = "none",
           replace = FALSE,
           seed = seed,
           oob.error = TRUE)
  
  y_hat_pred <- predict(fit_rf_e, data = data.frame(df_e[idx_test,]))$predictions
  
  idx_nnodes.obs <- mclapply(X = seq_len(nrow(A_u.iid)), mc.cores = 10, FUN = function(i){sum(A_u.iid[i,] != 0)}) %>% do.call(what = c, .)
  
  y.e_hat[idx_train] <- fit_rf_e$predictions
  # ei_hat.u <- y.e_hat[A_u.iid@i+1]*A_u.iid@x/idx_nnodes.obs[A_u.iid@i+1]
  ei_hat.u <- y.e_hat[A_u.iid@i+1]/idx_nnodes.obs[A_u.iid@i+1]
  
  y.e_hat_spt <- mclapply(X = seq_len(nstress), mc.cores = 10, FUN = function(i){mean(ei_hat.u[(A_u.iid@j+1)==i])}) %>% do.call(what = c, .)
  
  i = i + 1 # Increasing in 1 the value of the index for the loop
  DFrmse_train[i,] <- c(sqrt(mean((ysim[idx_train] - rinla$summary.fitted.values[idx_train,"mean"])**2)), sqrt(mean((ysim[idx_train] - rinla$summary.fitted.values[idx_train,"mean"] - y.e_hat[idx_train])**2)))

  if(verbose && i > 1){
    sprintf("KLD: %.4f. \n", KLD_GMRF_den) %>% cat(.)
    sprintf("Iteration: %i. \n", i-1) %>% cat(.)
    cat("---------------------------------------------- \n")
  }
}
t2 <- Sys.time()
difftime(t2, t1)

# Graphical results for the node correction ----

selected_spt_nodes$ID+1
corr_mean <- offx[(length(offx) - nrow(rinla$summary.fixed) - nstress) + 1:nstress] + rinla$summary.random$u.iid$mean

marginals_spt_rf <- mclapply(X = selected_spt_nodes$ID+1, mc.cores = 20, FUN = function(i){
  inla.tmarginal(marginal = rinla$marginals.random$spt[[i]], fun = function(x){(x - rinla$summary.random$spt$mean[i])*tau.iid**(-1/2) + corr_mean[which((selected_spt_nodes$ID+1)==i)]})
})

n_group <- 10
nodes_sel_marg_rf <- mclapply(X = 1:nstress, mc.cores = 1, FUN = function(i){cbind(data.frame(marginals_spt_rf[[i]]), group = i%%n_group + if(i%%n_group==0){n_group}else{0}, height = 0, ID_group = ceiling(i/n_group))}) %>% 
  do.call(what = rbind, .)

ggplot() +
  geom_ridgeline_gradient(data = nodes_sel_marg[nodes_sel_marg$ID_group %in% ngroup_idx_to_plot,], mapping = aes(x = x, y = height + group, group = group, height = y, scale = scale_ridge, fill = y), linewidth = 1, colour = "blue", alpha = 0.5) +
  geom_ridgeline_gradient(data = nodes_sel_marg_rf[nodes_sel_marg$ID_group %in% ngroup_idx_to_plot,], mapping = aes(x = x, y = height + group, group = group, height = y, scale = scale_ridge, fill = y), linewidth = 1, colour = "red", alpha = 0.5) +
  geom_linerange(data = data.frame(x = sim_node_values,
                                   y_max = height_error_bar+0.85, y_min = height_error_bar,
                                   ID_group = ceiling(seq_len(last(ngroup_idx_to_plot)*n_group)/n_group)),
                 mapping = aes(x = x, ymin = y_min, ymax = y_max),
                 linewidth = 1, colour = "black") +
  scale_y_continuous(name = "Posterior distribution", breaks = 1:10, 
                     labels = 1:10) +
  scale_fill_viridis_c(name = "Density", option = "mako") +
  facet_wrap(facets = ~ ID_group, ncol = 5) +
  theme_bw()
