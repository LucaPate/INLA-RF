# Title: runme_SPDE_strong_offset
# Author: Mario Figueira
# Date: 2025-02-18
# Description: Code combining INLA and RF, correcting the offset

# Last update: 2025-02-20

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

# Loading data ----

DF_strong_sim <- readRDS(file = "./Code_completed/Data_strong.RDS")

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
              num.threads = 24,
              verbose = FALSE)

rinla_orig <- rinla # saving the first analysis
mu_new <- rinla$misc$configs$config[[1]]$improved.mean
Q_new <- fix.Q(rinla$misc$configs$config[[1]]$Q)

# ggplot() + 
#   geom_line(data = DF_strong_sim, mapping = aes(x = X2, y = fX2), color = "red") + 
#   geom_line(data = rinla_orig$summary.random$f_x2, mapping = aes(x = ID, y = mean), color = "blue") +
#   theme_bw()
# 
# ggplot() + 
#   geom_line(data = DF_strong_sim, mapping = aes(x = X3, y = fX3), color = "red") + 
#   geom_line(data = rinla_orig$summary.random$f_x3, mapping = aes(x = ID, y = mean), color = "blue") +
#   theme_bw()

## Setting the configuration for the INLA-RF Algorithm ----
ysim <- DF_strong_sim$y
e <- y.e_hat <- rep(0, nrow(DF_strong_sim)) # Offset of the residual estimations from the RF. Initialized with y.e_hat = 0.
offx <- rep(0, times = nrow(DF_strong_sim))
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
  DFrmse_train[i,] <- c(sqrt(mean((ysim[idx_train] - rinla$summary.fitted.values[idx_train,"mean"])**2)), sqrt(mean((ysim[idx_train] - rinla$summary.fitted.values[idx_train,"mean"] - y.e_hat[idx_train])**2)))
  DFrmse_test[i,] <- c(sqrt(mean((ysim[idx_test] - rinla$summary.fitted.values[idx_test,"mean"])**2)), sqrt(mean((ysim[idx_test] - rinla$summary.fitted.values[idx_test,"mean"] - y.e_hat[idx_test])**2)))
  
  if(verbose && i > 1){
    sprintf("KLD: %.4f. \n", KLD_GMRF_den) %>% cat(.)
    sprintf("Iteration: %i. \n", i-1) %>% cat(.)
    cat("---------------------------------------------- \n")
  }
}
t2 <- Sys.time()
difftime(t2, t1)

