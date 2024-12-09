# Algorithm for the INLA-RF approach. Applied to Simulated data with Strong non-linearity
# Description: Combining INLA and RF for correcting the predictions. The correction from the random forest is 
# integrated as an offset, without transferring the uncertainty from the RF to the correction in INLA.

# Last update: 15/11/2024

remove(list=ls())

library(Matrix)
library(parallel)

library(INLA)
library(inlabru)
library(fmesher)

library(ggplot2)
library(gridExtra)
library(ggtext)
library(dplyr)
library(sf)

library(ranger)

seed <- 1234 # set a seed for reproducibility
set.seed(seed = seed)

# Custom functions 

## Function to fix the precision matrix from inla.call.object$misc$configs$config[[k]]$Q (or Qprior or Qinv)
fix.Q <- function(Q) {
  # Q: the precision matrix from the inla.call.object; inla.call.object <- inla(...) 
  d <- diag(Q)
  Q <- Q + t(Q)
  diag(Q) <- d
  return (Q)
}

## Function to simulate a SPDE effect in R2
simulation_SPDE <- function(mesh, sigma = 1, range = 1, constr=TRUE, seed){
  # mesh: mesh to create the precision matrix of the SPDE
  # sigma: the marginal precision of the CAR prior (the variability between first-order neighbors)
  # range: the spatial range; it means the range at which the correlation falls to 0.1 (in the correl)
  # constr: an argument that, if TRUE, forces the GMRF to have a global mean of zero
  # seed: seed to reproduce the result
  if(!missing(seed)){set.seed(seed)}
  Q = fmesher::fm_matern_precision(x = mesh, alpha = 2, rho = range, sigma = sigma)
  L = Matrix::chol(Q)
  w = rnorm(nrow(Q))
  u_sp = backsolve(L, w)
  if(constr){
    u_sp = u_sp - mean(u_sp)
  }
  return(list(u=u_sp, Q=Q))
}

## Defining the study region ----
data(PRborder)
sf_PRborder <- st_sfc(st_polygon(x = list(PRborder))) # boundary of the sr (study region)
PR_nchull <- fm_nonconvex_hull_inla(x =  sf_PRborder, convex = -0.05)
st_PR_int_nchull <-  st_polygon(x = list(PR_nchull$loc[c(PR_nchull$idx[1:which(PR_nchull$idx[1,1]==PR_nchull$idx[,2]),1],1),1:2]))

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

## Spatial effect simulation
sp_sim <- simulation_SPDE(mesh = mesh_sim, sigma = 1, range = 0.3, constr = TRUE, seed = seed)

st_bbox_sr <- st_bbox(obj = st_PR_int_nchull)
loc_grid <- expand.grid(x = seq(from = st_bbox_sr[1], to = st_bbox_sr[3], length.out = 1E2), y = seq(from = st_bbox_sr[2], to = st_bbox_sr[4], length.out = 1E2))
st_loc_grid_t <- st_multipoint(x = as.matrix(loc_grid)) %>% st_geometry(.) %>% st_cast("POINT")
idx <- st_intersects(x = st_loc_grid_t, st_PR_int_nchull, sparse = FALSE)

sample_size <- 1E3
st_loc_grid <- st_loc_grid_t[idx,] # Grid points that are inside the study region
st_sample_points <- st_sample(x = st_PR_int_nchull, size = sample_size) # Sample points

df_sp <- st_sf(data.frame(sp = drop(fm_basis(x = mesh_sim, loc = st_coordinates(st_loc_grid)[,1:2]) %*% sp_sim$u)), geometry = st_loc_grid)
df_sp_sample <- st_sf(data.frame(sp = drop(fm_basis(x = mesh_sim, loc = st_coordinates(st_sample_points)[,1:2]) %*% sp_sim$u)), geometry = st_sample_points)

## Plot of the spatial effect (grid)
ggsp_grid <- ggplot() + 
  geom_sf(data = df_sp, mapping = aes(color = sp), pch = 15) + 
  scale_color_viridis_c(option = "turbo") + 
  coord_sf() + theme_bw() + labs(color = "Values", title = "Spatial effect (SPDE)") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

## Plot of the spatial effect (sample locations)
ggsp_sample <- ggplot() + 
  geom_sf(data = df_sp_sample, mapping = aes(color = sp)) + 
  scale_color_viridis_c(option = "turbo") + 
  coord_sf() + theme_bw() + labs(color = "Values", title = "Spatial effect (SPDE)") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

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
dfsf_sample_strong <- st_sf(data.frame(y = y, sp = df_sp_sample$sp,  
                                       X1 = X1, fX1 = beta_cat[unclass(X1)], 
                                       X2 = X2, fX2 = fX2, X3 = X3, fX3, fX3), 
                            geometry = st_geometry(df_sp_sample)) # Data frame for the simulation with strong non-linearity  

gg_y_strong <- ggplot() + # Response variable under weak linearity in the sample locations
  geom_sf(data = dfsf_sample_strong, mapping = aes(color = y)) + 
  scale_color_viridis_c(option = "turbo") +
  coord_sf() + theme_bw() + labs(color = "Values", title = "Response variable in the sample locations") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggcov2_strong; ggcov3_strong; gg_y_strong

## Creating the train/test groups ----
idx.test <- sample(x = 1:nrow(dfsf_sample_strong), size = round(nrow(dfsf_sample_strong)*0.2)) 
simdata_train = dfsf_sample_strong[-idx.test,]
simdata_test = dfsf_sample_strong[idx.test,]

## INLA settings ----

coord_strong_df <- st_coordinates(dfsf_sample_strong)
coord_strong_train_df <- st_coordinates(simdata_train)
coord_strong_test_df <- st_coordinates(simdata_test)

idx_chull <- chull(x = coord_strong_df)
st_chull <- st_polygon(x = list(coord_strong_df[c(idx_chull, idx_chull[1]),]))

mesh_inf <- fm_mesh_2d_inla(boundary = list(st_chull), max.edge = c(0.15, 0.5), offset = c(0.01, -0.1))

rho_0 <- st_coordinates(st_chull)[,1:2] %>% dist(.) %>% max(.)/5
spde_inf <- inla.spde2.pcmatern(mesh = mesh_inf, alpha = 2, prior.range = c(rho_0, 0.5), prior.sigma = c(1, 0.5), constr = TRUE)
A_sp <- fm_basis(x = mesh_inf, loc = coord_strong_train_df)

y.e_hat <- rep(0, nrow(simdata_train)) # Offset of the residual estimations from the RF. Initialized with y.e_hat = 0.
KLD_GMRF_den <- 10 # Initial value for the KLD
i <- 1 # Number of loops = i-1

DFrmse <- data.frame(idx = NA, INLA = NA, INLA.RF = NA)

t1 <- Sys.time()
while(KLD_GMRF_den>1E-2){ # Using the KLD as condition
    
  if(i > 1){
    mu_old <- mu_new
    Q_old <- Q_new
    list_control_mode <- list(x = inla_model_1$mode$x, theta = inla_model_1$mode$theta, restart = TRUE, fixed = FALSE) 
  } else{
    list_control_mode <- inla.set.control.mode.default()
  }
  
  inf_stk_1 <- inla.stack(data = list(y = simdata_train$y - y.e_hat),
                          A = list(A_sp, 1), 
                          effects = list(
                            list(sp = 1:mesh_inf$n),
                            list(
                              beta_0 = rep(1, nrow(simdata_train)),
                              idx_cat = as.numeric(simdata_train$X1),
                              beta_2 = simdata_train$X2,
                              beta_3 = simdata_train$X3
                            )
                          ),
                          tag = "inf_stk_1")
  
  formula_inla_1 <- y ~ -1 + beta_0 + beta_2 + beta_3 + f(idx_cat, model = "iid", constr = TRUE) + f(sp, model = spde_inf)
  
  inla_model_1 <- inla(data = inla.stack.data(inf_stk_1), formula = formula_inla_1, family = "gaussian",
                       control.predictor = list(A = inla.stack.A(inf_stk_1)),
                       control.compute = list(config = TRUE),
                       control.mode = list_control_mode,
                       verbose = FALSE)
  
  mu_new <- inla_model_1$misc$configs$config[[1]]$improved.mean
  Q_new <- inla_model_1$misc$configs$config[[1]]$Q
  Q_inv_new <- inla_model_1$misc$configs$config[[1]]$Qinv
  
  if(i == 1){
    mu_old <- rep(0, times = length(mu_new))
    Q_old <- inla_model_1$misc$configs$config[[1]]$Qprior
  }
  
  KLD_GMRF_den <- 1/2*(sum(diag(solve(Q_new)%*%Q_old)) - nrow(Q_new) + drop((mu_new-mu_old) %*% Q_new %*% cbind(mu_new-mu_old)) + determinant(Q_new, logarithm = TRUE)$modulus - determinant(Q_old, logarithm = TRUE)$modulus)/nrow(Q_new)
  
  idx.fitted <- inla_model_1$summary.fitted.values %>% rownames(.) %>% grepl(pattern = "fitted.APredictor.", x = .)
  e.train <- simdata_train$y - inla_model_1$summary.fitted.values$mean[idx.fitted]
  e.test <- simdata_test$y - (inla_model_1$summary.fixed$mean[1] + 
                                inla_model_1$summary.fixed$mean[2]*simdata_test$X2 +
                                inla_model_1$summary.fixed$mean[3]*simdata_test$X3 +
                                inla_model_1$summary.random$idx_cat$mean[as.numeric(simdata_test$X1)] +
                                drop(fm_basis(x = mesh_inf, loc = coord_strong_test_df) %*% inla_model_1$summary.random$sp$mean))
  
  # Prediction on test data from INLA model
  
  df_e.train <- data.frame(e = e.train, xcoord = coord_strong_train_df[,1], ycoord = coord_strong_train_df[,2], 
                           X1 = dfsf_sample_strong$X1[-idx.test], 
                           X2 = dfsf_sample_strong$X2[-idx.test], 
                           X3 = dfsf_sample_strong$X3[-idx.test])
  
  # Fitting the residual information with a RF
  
  fit_rf_e.train <- 
    ranger(formula = e ~ X1 + X2 + X3 + xcoord + ycoord,
           data = df_e.train,
           importance = "none",
           replace = FALSE,
           seed = seed,
           oob.error = TRUE) 
  
  # Compute predictions on the left out station (this need to be changed, the code is wrong as it is now)
  y_hat <- predict(fit_rf_e.train, data = data.frame(simdata_test, xcoord = coord_strong_test_df[,1], ycoord = coord_strong_test_df[,2]))$predictions # (this need to be changed, the code is wrong as it is now)
  
  # Compute prediction intervals
  oob_error <- df_e.train$e - fit_rf_e.train$predictions #oob predictions
  
  alpha = 0.05 #for prediction intervals
  lowerPred <- y_hat + quantile(oob_error, alpha/2)
  upperPred <- y_hat + quantile(oob_error, 1 - alpha/2)
  
  dens_oob_error <- density(oob_error)
  dens_oob_error$x <- dens_oob_error$x[dens_oob_error$y!=0]
  dens_oob_error$y <- dens_oob_error$y[dens_oob_error$y!=0]
  
  # Optimazing only the standard deviation of the optimal Gaussian
  # optim_sd <- optim(par = sd(oob_error), method = "Brent", fn = function(x){pracma::trapz(x = dens_oob_error$x, y = abs(dens_oob_error$y * (log(dens_oob_error$y) - log(dnorm(x = dens_oob_error$x,mean = mean(oob_error), sd = x)))))}, lower = 0, upper = 5*sd(oob_error))
  
  # Optimazing the mean and the standard deviation of the optimal Gaussian
  optim_mu_sd <- optim(par = c(mean(oob_error), sd(oob_error)), method = "BFGS", fn = function(x){pracma::trapz(x = dens_oob_error$x, y = abs(dens_oob_error$y * (log(dens_oob_error$y) - log(dnorm(x = dens_oob_error$x, mean = x[1], sd = x[2])))))})
  
  KLD_emp_the <-  pracma::trapz(x = dens_oob_error$x, y = abs(dens_oob_error$y * (log(dens_oob_error$y) - log(dnorm(x = dens_oob_error$x,mean = optim_mu_sd$par[1], sd = optim_mu_sd$par[2])))))
  
  ggplot() + 
    geom_line(data = data.frame(dens_oob_error[c("x", "y")]), mapping = aes(x = x, y = y)) +
    geom_line(data = data.frame(x = seq(min(dens_oob_error$x),max(dens_oob_error$x),length.out=1E3), y = dnorm(x = seq(min(dens_oob_error$x), max(dens_oob_error$x), length.out = 1E3), mean = optim_mu_sd$par[1], sd = optim_mu_sd$par[2])),
              mapping = aes(x = x, y = y), colour = "red") + 
    labs(title = "Real and optimal Guassian PDF", subtitle = sprintf("KLD = %s", formatC(KLD_emp_the, format = "e", digits = 3))) +
    theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"),
                       plot.subtitle = element_markdown(hjust=0.5, face="bold", size=12, color = "black"))
  
  ggplot() + geom_line(data = data.frame(x = dens_oob_error$x, y = dens_oob_error$y * (log(dens_oob_error$y) - log(dnorm(x = dens_oob_error$x,mean = optim_mu_sd$par[1], sd = optim_mu_sd$par[2])))), 
                       mapping = aes(x = x, y = y), colour = "blue") + 
    labs(title = "Derivative of KLD(p||q)") + 
    theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"))
  
  e.test_irf <- simdata_test$y - (inla_model_1$summary.fixed$mean[1] +
                                    inla_model_1$summary.fixed$mean[2] * dfsf_sample_strong$X2 +
                                    inla_model_1$summary.fixed$mean[3] * dfsf_sample_strong$X3 +
                                    inla_model_1$summary.random$idx_cat$mean[as.numeric(simdata_test$X1)] + 
                                    drop(fm_basis(x = mesh_inf, loc = coord_strong_test_df) %*% inla_model_1$summary.random$sp$mean) + 
                                    y_hat)  # Prediction on test data (INLA+RF) 
  
  df_res <- data.frame(idx = 1, matrix(data = c(sqrt(mean(e.train**2)), sqrt(mean((oob_error)**2)), sqrt(mean(e.test**2)), sqrt(mean(e.test_irf**2))), byrow = TRUE, ncol = 2))
  colnames(df_res) <- colnames(DFrmse); rownames(df_res) <- c(paste("Train", i), paste("Test", i))
  
  if(i == 1){DFrmse <- df_res} else{DFrmse <- rbind(DFrmse, df_res)}
  
  y.e_hat <- fit_rf_e.train$predictions
  i = i + 1
  
  sprintf("KLD: %.4f. \n", KLD_GMRF_den) %>% cat(.)
  sprintf("Iteration: %i. \n", i-1) %>% cat(.)
  cat("---------------------------------------------- \n")
}
t2 <- Sys.time()

difftime(t2, t1)
DFrmse
