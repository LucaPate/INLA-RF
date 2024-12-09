# Algorithm for the INLA-RF approach. Applied to <<meuse>> data
# Last update: 9/12/2024

remove(list=ls())

library(Matrix)
library(INLA)
library(fmesher)
library(ggplot2)
library(dplyr)
library(ranger)
library(sf)

seed <- 1234
set.seed(seed = seed)

# Loading the Meuse data (from the 'sp' package) 

data("meuse")
dfsf_meuse <- st_as_sf(x = meuse, coords = c("x", "y"))

# Splitting data in train/test sets

idx.test <- sample(x = 1:nrow(dfsf_meuse), size = round(nrow(dfsf_meuse)*0.2))
df_meuse_train = dfsf_meuse[-idx.test,]
df_meuse_test = dfsf_meuse[idx.test,]

## INLA settings ----

coords_meuse_sf <- st_coordinates(dfsf_meuse)
coord_meuse_train_df <- st_coordinates(df_meuse_train)
coord_meuse_test_df <- st_coordinates(df_meuse_test)

idx_chull <- chull(x = coords_meuse_sf)
st_chull <- st_polygon(x = list(coords_meuse_sf[c(idx_chull, idx_chull[1]),]))

min_dist_chull <- min(st_bbox(st_chull) %>% matrix(data = ., ncol = 2) %>% apply(X = ., MARGIN = 1, FUN = diff))
mesh_inf <- fm_mesh_2d_inla(boundary = list(st_chull), max.edge = c(min_dist_chull/25, min_dist_chull/10), offset = c(0.01, 0.2*min_dist_chull))

rho_0 <- st_coordinates(st_chull)[,1:2] %>% dist(.) %>% max(.)/5
spde_inf <- inla.spde2.pcmatern(mesh = mesh_inf, alpha = 2, prior.range = c(rho_0, 0.5), prior.sigma = c(1, 0.5), constr = TRUE)
A_sp <- fm_basis(x = mesh_inf, loc = coord_meuse_train_df)

y.e_hat <- rep(0, nrow(df_meuse_train)) # Offset of the residual estimations from the RF. Initialized with y.e_hat = 0.
KLD_GMRF_den <- 10 # Initial value for the KLD
i <- 1 # Number of loops = i-1

DF_RMSE <- data.frame(Train = NA, Test = NA)
verbose <- TRUE

# Defining the INLA-RF algorithm
while(KLD_GMRF_den>1E-3){ # Using the KLD as condition
  ### Fitting the model with INLA ----
  
  inf_stk_2 <- inla.stack(data = list(y = df_meuse_train$cadmium - y.e_hat), # list(y = rep(NA, times = length(dfsf_meuse$y))),
                          A = list(A_sp, 1), 
                          effects = list(
                            list(sp = 1:mesh_inf$n),
                            list(
                              beta_0 = rep(1, nrow(df_meuse_train)),
                              idx_soil = as.numeric(df_meuse_train$soil),
                              beta_1 = df_meuse_train$elev,
                              beta_2 = df_meuse_train$dist
                            )
                          ),
                          remove.unused = FALSE,
                          tag = "inf_stk_2")
  
  formula_inla_2 <- y ~ -1 + 
    beta_0 + beta_1 + beta_2 +
    f(idx_soil, model = "iid", constr = TRUE) + f(sp, model = spde_inf)
  
  if(i>1){
    mu_old <- inla_model_2$misc$configs$config[[1]]$improved.mean
    Q_old <- inla_model_2$misc$configs$config[[1]]$Q
  }
  
  inla_model_2 <- inla(data = inla.stack.data(inf_stk_2),
                       formula = formula_inla_2, family = "gaussian",
                       control.predictor = list(A = inla.stack.A(inf_stk_2)),
                       control.mode = if(i==1){list()}else{list(theta = inla_model_2$mode$theta, x = inla_model_2$mode$x)},
                       control.compute = list(config = TRUE),
                       verbose = FALSE)
  
  mu_new <- inla_model_2$misc$configs$config[[1]]$improved.mean
  Q_new <- inla_model_2$misc$configs$config[[1]]$Q
  Q_inv_new <- inla_model_2$misc$configs$config[[1]]$Qinv
  
  if(i>1){
    KLD_GMRF_den <- 1/2*(sum(diag(solve(Q_new)%*%Q_old)) - nrow(Q_new) + drop((mu_new-mu_old) %*% Q_new %*% cbind(mu_new-mu_old)) + determinant(Q_new, logarithm = TRUE)$modulus - determinant(Q_old, logarithm = TRUE)$modulus)/nrow(Q_new)
  }
  
  idx.fitted <- inla_model_2$summary.fitted.values %>% rownames(.) %>% grepl(pattern = "fitted.APredictor.", x = .)
  e.train <- df_meuse_train$cadmium - inla_model_2$summary.fitted.values$mean[idx.fitted]
  e.test <- df_meuse_test$cadmium - (inla_model_2$summary.fixed$mean[1] +
                                       df_meuse_test$elev*inla_model_2$summary.fixed$mean[2] + df_meuse_test$dist*inla_model_2$summary.fixed$mean[3] +
                                       inla_model_2$summary.random$idx_soil$mean[as.numeric(df_meuse_test$soil)] +
                                       drop(fm_basis(x = mesh_inf, loc = coord_meuse_test_df) %*% inla_model_2$summary.random$sp$mean)) # Prediction on test data 
  
  df_e.train <- data.frame(e = e.train, sp = drop(A_sp %*% inla_model_2$summary.random$sp$mean), 
                           elev = dfsf_meuse$elev[-idx.test], 
                           dist = dfsf_meuse$dist[-idx.test], 
                           soil = dfsf_meuse$soil[-idx.test])
  
  ### 2. Residual analysis by Random Forest  ----
  ### leveraging the posterior information geometry from the INLA object
  
  fit_rf_e.train <- 
    ranger(formula = e ~ elev + dist + soil + sp,
           data = df_e.train,
           importance = "none",
           replace = FALSE,
           seed = seed,
           oob.error = TRUE) 
  
  # Compute predictions on the left out station (this need to be changed, the code is wrong as it is now)
  y_hat <- predict(fit_rf_e.train, data = data.frame(df_meuse_test, sp = drop(fm_basis(x = mesh_inf, loc = coord_meuse_test_df) %*% inla_model_2$summary.random$sp$mean)))$predictions # (this need to be changed, the code is wrong as it is now)
  
  # Compute prediction intervals
  oob_error <- df_e.train$e - fit_rf_e.train$predictions #oob predictions
  
  alpha = 0.05 #for prediction intervals
  lowerPred <- y_hat + quantile(oob_error, alpha/2)
  upperPred <- y_hat + quantile(oob_error, 1 - alpha/2)
  
  dens_oob_error <- density(oob_error)
  optim_sd <- optim(par = sd(oob_error), method = "Brent", fn = function(x){pracma::trapz(x = dens_oob_error$x, y = abs(dens_oob_error$y * (log(dens_oob_error$y) - log(dnorm(x = dens_oob_error$x,mean = mean(oob_error), sd = x)))))}, lower = 0, upper = 5*sd(oob_error))
  
  ggplot() + 
    geom_line(data = data.frame(dens_oob_error[c("x", "y")]), mapping = aes(x = x, y = y)) +
    geom_line(data = data.frame(x = seq(min(dens_oob_error$x),max(dens_oob_error$x),length.out=1E3), y = dnorm(x = seq(min(dens_oob_error$x), max(dens_oob_error$x), length.out = 1E3), mean = mean(oob_error), sd = optim_sd$par)),
              mapping = aes(x = x, y = y), colour = "red") + 
    theme_bw()
  
  ggplot() + geom_line(data = data.frame(x = dens_oob_error$x, y = dens_oob_error$y * (log(dens_oob_error$y) - log(dnorm(x = dens_oob_error$x,mean = mean(oob_error), sd = optim_sd$par)))), 
                       mapping = aes(x = x, y = y), colour = "blue") + labs(title = "Derivative of KLD(p||q)") + theme_bw()
  KLD_emp_the <-  pracma::trapz(x = dens_oob_error$x, y = abs(dens_oob_error$y * (log(dens_oob_error$y) - log(dnorm(x = dens_oob_error$x,mean = mean(oob_error), sd = optim_sd$par)))))
  
  e.test_irf <- df_meuse_test$cadmium - (inla_model_2$summary.fixed$mean[1] +
                                           df_meuse_test$elev*inla_model_2$summary.fixed$mean[2] + df_meuse_test$dist*inla_model_2$summary.fixed$mean[3] +
                                           inla_model_2$summary.random$idx_soil$mean[as.numeric(df_meuse_test$soil)] +
                                           drop(fm_basis(x = mesh_inf, loc = coord_meuse_test_df) %*% inla_model_2$summary.random$sp$mean) + 
                                           y_hat)  # Prediction on test data (INLA+RF) 
  
  DF_RMSE[i,1] <- sqrt(mean((oob_error)**2))/sqrt(mean(e.train**2)) # Train RMSE
  DF_RMSE[i,2] <- sqrt(mean(e.test**2))/sqrt(mean(e.test_irf**2)) # Test RMSE
  
  y.e_hat <- fit_rf_e.train$predictions
  i <- i+1 
  if(verbose){
    sprintf("KLD: %.4f. \n", KLD_GMRF_den) %>% cat(.)
    sprintf("Iteration: %i. \n", i-1) %>% cat(.)
    cat("---------------------------------------------- \n")
  }
}

DF_RMSE
