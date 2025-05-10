# Algorithm for the INLA-RF approach. Applied to Simulated data with jumps or stress points 
# Description: Combining INLA and RF for correcting the predictions. The correction from the random forest is 
# integrated as an offset, transferring the uncertainty from the RF to the correction in INLA. 
# In this case, there are only some points that are corrected by the RF in the INLA latent field, those points 
# are called stress points. This avoids an overfitting from correcting the whole data.

# Date (Creation): 9/12/2024
# Last update: 10/05/2025

remove(list=ls())

library(Matrix)
library(INLA)
library(ranger)

library(ggplot2)
library(ggtext)
library(gridExtra)
library(ggmagnify)

library(dplyr)

seed <- 123
set.seed(seed)

nsize <- 2E3
prec.rw1 <- prec.gauss <- 20

n.strsp <- 10 # number of stress points
loc.stress <- 1:n.strsp * floor(nsize/(n.strsp + 1)) # location of the stress points 
vec.strsp <- rep(0, times = nsize)
vec.strsp[loc.stress] <- sign(rbinom(n = n.strsp, size = 1, prob = 0.5) - 0.5) * rnorm(n = n.strsp, mean = 5, sd = 0.1)

u.rw1 <- cumsum(rnorm(nsize, mean = 0, sd = (prec.rw1)**(-1/2)) + vec.strsp); u.rw1 <- u.rw1 - mean(u.rw1)
u.rw1 %>% plot(., type = "l")

ysim <- rnorm(n = nsize, mean = u.rw1 + 2, sd = (prec.gauss)**(-1/2))
ysim %>% plot(., type = "l")

ggplot() + 
  geom_line(data = data.frame(Time = seq_along(ysim), y = ysim), mapping = aes(x = Time, y = y)) +
  theme_bw() + theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), 
                     axis.text = element_text(size = 14))

rinla <- inla(data = list(y = ysim, beta0 = rep(1, nsize), u = 1:nsize), 
              family = "gaussian",
              formula = y ~ -1 + beta0 + f(u, model = "rw2", constr = TRUE),
              control.compute = list(config = TRUE),
              verbose = FALSE)
rinla_orig <- rinla # saving the first analysis

ggrw1 <- ggplot() +
  geom_ribbon(data = data.frame(rinla_orig$summary.random$u, id = "Temporal Effect") %>%
                rename(., all_of(c(q1 = 'X0.025quant', q3 = 'X0.975quant'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "blue", alpha = 0.4) +
  geom_line(data = data.frame(rinla_orig$summary.random$u, id = "Temporal Effect"),
            mapping = aes(x = ID, y = mean), color = "blue") +
  geom_line(data = data.frame(ID = rinla_orig$summary.random$u$ID, mean = u.rw1, id = "Temporal Effect"),
            mapping = aes(x = ID, y = mean), color = "black") +
  geom_vline(xintercept = loc.stress, color = "red") +
  theme_bw() + labs(title = "A. Temporal effect (marginals)") +
  theme(plot.title = element_text(size = 20, h = 0, face = "bold"),
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 14))

ggvar <- ggplot() +
  geom_point(data = rinla_orig$summary.random$u, mapping = aes(x = ID, y = sd**2), color = "blue") +
  geom_vline(xintercept = loc.stress, color = "red") + 
  labs(title = "B. Latent-field (variance)") +
  theme_bw() + ylab(label = expression(sigma^2~(u[i]))) +
  theme(plot.title = element_text(size = 20, h = 0, face = "bold"),
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 14))

gglp_var <- ggplot() +
  geom_point(data = data.frame(ID = 1:nrow(rinla_orig$summary.fitted.values), rinla_orig$summary.fitted.values), mapping = aes(x = ID, y = sd**2), color = "blue") +
  geom_vline(xintercept = loc.stress, color = "red") + 
  labs(title = "C. Linear predictor (variance)") +
  theme_bw() + ylab(label = expression(sigma^2~(eta[i]))) +
  theme(plot.title = element_text(size = 20, h = 0, face = "bold"),
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 14))

grid.arrange(arrangeGrob(grobs = list(ggrw1, ggvar, gglp_var), ncol = 1))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors_base <- gg_color_hue(2)

nstress <- 1E2
sel_stressp <- (rinla_orig$summary.random$u$sd**2) %>% order(., decreasing = TRUE) %>% .[1:nstress]
stress_col <- rep("blue", times = nsize); stress_col[sel_stressp] <- "red"
ggplot() +
  geom_ribbon(data = data.frame(ID = 1:nsize, rinla_orig$summary.random$u, id = "Temporal Effect") %>%
                rename(., all_of(c(q1 = 'X0.025quant', q3 = 'X0.975quant'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "blue", alpha = 0.4) +
  geom_line(data = data.frame(ID = 1:nsize, rinla_orig$summary.random$u, id = "Temporal Effect"),
            mapping = aes(x = ID, y = mean), color = "blue") +
  geom_line(data = data.frame(ID = 1:nsize, mean = u.rw1, id = "Temporal Effect"),
            mapping = aes(x = ID, y = mean), color = "black") +
  geom_point(data = data.frame(ID = 1:nsize, mean = u.rw1), mapping = aes(x = ID, y = mean, colour = stress_col), size = 2) +
  theme_bw() + labs(title = "Temporal effect nodes") +
  labs(colour = "Stress points") +
  scale_color_manual(labels = c("Not stressed", "Stressed"),
                     values = colors_base) +
  theme(plot.title = element_text(size = 20, h = 0.5, face = "bold"),
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16))

# Setting the configuration for the Algorithm ----

y.e_hat <- rep(0, length(ysim)) # Offset of the residual estimations from the RF. Initialized with y.e_hat = 0.
KLD_GMRF_den <- 10 # Initial value for the KLD
i <- 1 # Number of loops = i-1

DFrmse <- data.frame(INLA = NA, INLA.RF = NA)
verbose <- TRUE

rinla <- rinla
sel <- sel_stressp

A_u.iid <- matrix(data = 0, nrow = nsize, ncol = length(sel)) %>% as(., "sparseMatrix")
for(k in sel){
  A_u.iid[k,which(k == sel)] <- 1
}

inf_stack <- inla.stack(data = list(y = ysim),
                        A = list(1, A_u.iid),
                        effects = list(
                          list(beta0 = rep(1, nsize),
                               u.rw1 = 1:nsize),
                          list(u.iid = 1:nstress)
                        ),
                        compress = FALSE,
                        remove.unused = FALSE,
                        tag = "inf.stk")

t1 <- Sys.time()
while(KLD_GMRF_den>1E-2){ # Using the KLD as condition
  if(i == 2){
    mu_old <- mu_new
    Q_old <- Q_new
    tau.iid <- as.numeric(optim_mu_sd[2])**(-2)
    list_control_mode <- list(theta = rinla$mode$theta, restart = TRUE, fixed = FALSE)
  } else if(i > 2){
    mu_old <- mu_new
    Q_old <- Q_new
    tau.iid <- as.numeric(optim_mu_sd[2])**(-2)
    list_control_mode <- list(x = c(rep(0, rinla$misc$configs$mnpred), rinla$misc$configs$config[[1]]$improved.mean), 
                              theta = rinla$mode$theta, restart = TRUE, fixed = FALSE)
  } else{
    tau.iid <- 0.001
    list_control_mode <- inla.set.control.mode.default()
  }
  
  if(!("rinla" %in% ls()) | i != 1){
    if(length(rinla$misc$configs$offsets) != (dim(inf_stack$A)[2] + length(ysim))){
      offx <- rep(0, times = dim(inf_stack$A)[2])
      offx[(length(offx) - nrow(rinla$summary.fixed) - length(sel) + 1) + 1:length(sel)] <- y.e_hat[sel]
    } else{
      offx[(length(offx) - nrow(rinla$summary.fixed) - length(sel) + 1) + 1:length(sel)] <- offx[(length(offx) - nrow(rinla$summary.fixed) - length(sel) + 1) + 1:length(sel)] + rinla$summary.random$u.iid$mean + y.e_hat[sel]
    }
    
    rinla <- inla(data = inla.stack.data(inf_stack), 
                  family = "gaussian",
                  formula = y ~ -1 + offset(offx) + beta0 + f(u.rw1, model = "rw2", constr = TRUE) + f(u.iid, model = "iid", hyper = list(prec = list(initial = log(tau.iid), fixed = TRUE))),
                  control.compute = list(config = TRUE),
                  control.predictor = list(A = inla.stack.A(inf_stack)),
                  control.mode = list_control_mode,
                  verbose = FALSE)
    
    idx.uiid <- rinla$misc$configs$contents$length[3:(which(rinla$misc$configs$contents$tag == "u.iid")-1)] + (1:rinla$misc$configs$contents$length[which(rinla$misc$configs$contents$tag == "u.iid")])
    mu_new <- rinla$misc$configs$config[[1]]$improved.mean[-c(idx.uiid)]
    Q_new <- rinla$misc$configs$config[[1]]$Q[-c(idx.uiid), -c(idx.uiid)]
    Q_inv_new <- rinla$misc$configs$config[[1]]$Qinv[-c(idx.uiid), -c(idx.uiid)]
    if(i == 1){
      mu_old <- rep(0, times = length(mu_new))
      Q_old <- rinla$misc$configs$config[[1]]$Qprior[-c(idx.uiid), -c(idx.uiid)]
    }
  } else{
    mu_new <- rinla$misc$configs$config[[1]]$improved.mean
    Q_new <- rinla$misc$configs$config[[1]]$Q
    Q_inv_new <- rinla$misc$configs$config[[1]]$Qinv
    
    mu_old <- rep(0, times = length(mu_new))
    Q_old <- rinla$misc$configs$config[[1]]$Qprior
  }
  
  KLD_GMRF_den <- 1/2*(sum(diag(solve(Q_new)%*%Q_old)) - nrow(Q_new) + drop((mu_new-mu_old) %*% Q_new %*% cbind(mu_new-mu_old)) + determinant(Q_new, logarithm = TRUE)$modulus - determinant(Q_old, logarithm = TRUE)$modulus)/nrow(Q_new)
  
  e <- ysim - rinla$summary.fitted.values[1:length(ysim),"mean"]
  df_e <- data.frame(e = e, rw1 = 1:nsize)
  
  fit_rf_e <- 
    ranger(formula = e ~ rw1,
           data = df_e,
           importance = "none",
           replace = FALSE,
           seed = seed,
           oob.error = TRUE)
  
  oob_error <- df_e$e - fit_rf_e$predictions #oob predictions
  
  dens_oob_error <- density(oob_error)
  dens_oob_error$x <- dens_oob_error$x[dens_oob_error$y!=0]
  dens_oob_error$y <- dens_oob_error$y[dens_oob_error$y!=0]
  
  # Optimazing the mean and the standard deviation of the optimal Gaussian
  # optim_mu_sd <- optim(par = c(mean(oob_error), sd(oob_error)), method = "BFGS", fn = function(x){pracma::trapz(x = dens_oob_error$x, y = abs(dens_oob_error$y * (log(dens_oob_error$y) - log(dnorm(x = dens_oob_error$x, mean = x[1], sd = x[2])))))})
  # KLD_emp_the <-  pracma::trapz(x = dens_oob_error$x, y = abs(dens_oob_error$y * (log(dens_oob_error$y) - log(dnorm(x = dens_oob_error$x,mean = optim_mu_sd$par[1], sd = optim_mu_sd$par[2])))))
  optim_mu_sd <- c(mean(oob_error), sd(oob_error)) # alternative to avoid the optimization step (just assume that the empirical mean and variance are used
  
  if("sel" %in% ls()){
    y.e_hat[sel] <- fit_rf_e$predictions[sel]
  } else{
    y.e_hat <- fit_rf_e$predictions
  }
  
  DFrmse[i,] <- c(sqrt(mean((ysim - rinla$summary.fitted.values[1:length(ysim),"mean"])**2)), sqrt(mean((ysim - rinla$summary.fitted.values[1:length(ysim),"mean"] - y.e_hat)**2)))
  
  i = i + 1
  if(verbose){
    sprintf("KLD: %.4f. \n", KLD_GMRF_den) %>% cat(.)
    sprintf("Iteration: %i. \n", i-1) %>% cat(.)
    cat("---------------------------------------------- \n")
  }
}
t2 <- Sys.time()
difftime(t2, t1)

# RMSE for the different steps: INLA means the inla output after incorporating the RF corrections (for i > 1), and INLA.RF means the rmse for the RF output
DFrmse

sqrt(mean((ysim[sel_stressp] - rinla_orig$summary.fitted.values[1:length(ysim),"mean"][sel_stressp])**2)) # RMSE for the stress points using the first INLA without RF corrections
sqrt(mean((ysim[sel_stressp] - rinla$summary.fitted.values[1:length(ysim),"mean"][sel_stressp])**2)) # RMSE for the stress points after the RF corrections

gg_lpred_INLA <- ggplot() + # ILA without RF corrections 
  geom_ribbon(data = data.frame(ID = 1:nsize, rinla_orig$summary.fitted.values[1:length(ysim),], id = "Temporal Effect") %>%
                rename(., all_of(c(q1 = 'X0.025quant', q3 = 'X0.975quant'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "blue", alpha = 0.4) +
  geom_line(data = data.frame(ID = 1:nsize, rinla_orig$summary.fitted.values[1:length(ysim),], id = "Temporal Effect"),
            mapping = aes(x = ID, y = mean), color = "blue") +
  geom_line(data = data.frame(ID = 1:nsize, mean = ysim, id = "Temporal Effect"),
            mapping = aes(x = ID, y = mean), color = "black") +
  geom_point(data = data.frame(x = sel_stressp, y = ysim[sel_stressp]), mapping = aes(x = x, y = y), color = "red", alpha = 0.25) +
  labs(title = "Linear predictor (INLA)") +
  theme_bw() + theme(plot.title = element_text(face = "bold", h = 0.5))


gg_lpred_INLA.RF <- ggplot() + # INLA without RF corrections 
  geom_ribbon(data = data.frame(ID = 1:nsize, rinla$summary.fitted.values[1:length(ysim),], id = "Temporal Effect") %>%
                rename(., all_of(c(q1 = 'X0.025quant', q3 = 'X0.975quant'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "blue", alpha = 0.4) +
  geom_line(data = data.frame(ID = 1:nsize, rinla$summary.fitted.values[1:length(ysim),], id = "Temporal Effect"),
            mapping = aes(x = ID, y = mean), color = "blue") +
  geom_line(data = data.frame(ID = 1:nsize, mean = ysim, id = "Temporal Effect"),
            mapping = aes(x = ID, y = mean), color = "black") +
  geom_point(data = data.frame(x = sel_stressp, y = ysim[sel_stressp]), mapping = aes(x = x, y = y), color = "red", alpha = 0.25) +
  labs(title = "Linear predictor (INLA-RF)") +
  theme_bw() + theme(plot.title = element_text(face = "bold", h = 0.5))

grid.arrange(arrangeGrob(grobs = list(gg_lpred_INLA, gg_lpred_INLA.RF), ncol = 1))

gg_origcor_INLA <- ggplot() +
  geom_point(data = rinla_orig$summary.fitted.values[1:length(ysim),][sel_stressp[order(sel_stressp)],], mapping = aes(x = 1:length(sel_stressp)-0.25, y = mean)) +
  geom_errorbar(data = rinla_orig$summary.fitted.values[1:length(ysim),][sel_stressp[order(sel_stressp)],], aes(x = 1:length(sel_stressp)-0.25, ymin=mean-sd, ymax=mean+sd), width=.2, position = position_dodge(0.05)) +
  geom_point(data = rinla$summary.fitted.values[1:length(ysim),][sel_stressp[order(sel_stressp)],], mapping = aes(x = 1:length(sel_stressp) + 0.25, y = mean), colour = "red") +
  geom_errorbar(data = rinla$summary.fitted.values[1:length(ysim),][sel_stressp[order(sel_stressp)],], aes(x = 1:length(sel_stressp) + 0.25, ymin=mean-sd, ymax=mean+sd), width=.2, position = position_dodge(0.05), colour = "red") +
  geom_segment(data.frame(xi = 1:length(sel_stressp)-0.5, xe = 1:length(sel_stressp)+0.5, y = ysim[sel_stressp[order(sel_stressp)]]), mapping = aes(x = xi, xend = xe, y = y), lty = "solid", colour = "blue", alpha = 0.75) +
  labs(title = "Stress points (original vs corrected by RF)") +
  theme_bw() + theme(plot.title = element_text(face = "bold", h = 0.5))

# In black the original estimations from INLA, 
# in red the estimations from INLA corrected by the RF, 
# and in blue the true value for each stress point that is plotted
gg_origcor_INLA

grid.arrange(arrangeGrob(grobs = list(gg_lpred_INLA, gg_lpred_INLA.RF, gg_origcor_INLA), layout_matrix = matrix(data = c(1,2,3,3), ncol = 2)))

rinla_u.rw1_sel <- rinla$summary.random$u.rw1[sel_stressp[order(sel_stressp)],]
rinla_u.rf_sel <- rinla$summary.random$u.iid[order(sel_stressp),]

rinla_u.rw1_sel$mean <- rinla_u.rw1_sel$mean + rinla_u.rf_sel$mean + + offx[(length(offx) - nrow(rinla$summary.fixed) - length(sel) + 1) + 1:length(sel)][order(sel_stressp)]
rinla_u.rw1_sel$sd <- sqrt(rinla_u.rw1_sel$sd**2 + rinla_u.rf_sel$sd**2)

rw1_corr <- rinla$summary.random$u.rw1
rw1_corr$mean[sel_stressp[order(sel_stressp)]] <- rinla_u.rw1_sel$mean
rw1_corr$sd[sel_stressp[order(sel_stressp)]] <- rinla_u.rw1_sel$sd
rw1_corr[,"0.025quant"] <- qnorm(p = 0.025, mean = rw1_corr$mean, sd = rw1_corr$sd)
rw1_corr[,"0.975quant"] <- qnorm(p = 0.975, mean = rw1_corr$mean, sd = rw1_corr$sd)

gg_rw_INLA <- ggplot(rw1_corr, mapping = aes(ID, mean)) +
  geom_ribbon(data = rinla_orig$summary.random$u %>% data.frame(., id = "Original") %>%
                rename(., all_of(c(q1 = 'X0.025quant', q3 = 'X0.975quant'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "blue", alpha = 0.4) +
  geom_line(data = data.frame(rinla_orig$summary.random$u, id = "Original"),
            mapping = aes(x = ID, y = mean), color = "blue") +
  geom_line(data = data.frame(ID = 1:nsize, mean = u.rw1, id = "Original"),
            mapping = aes(x = ID, y = mean), color = "black") +
  geom_ribbon(data = rw1_corr %>% data.frame(., id = "Corrected") %>%
                rename(., all_of(c(q1 = 'X0.025quant', q3 = 'X0.975quant'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "blue", alpha = 0.4) +
  geom_line(data = data.frame(rw1_corr, id = "Corrected"),
            mapping = aes(x = ID, y = mean), color = "blue") +
  geom_line(data = data.frame(ID = 1:nsize, mean = u.rw1, id = "Corrected"),
            mapping = aes(x = ID, y = mean), color = "black") +
  geom_vline(xintercept = sel_stressp[order(sel_stressp)], colour = "salmon", alpha = 0.1) +
  facet_wrap(facets = ~ id, ncol = 1) +
  labs(title = "Latent field (orginal and corrected by RF)") +
  theme_bw() + # theme(plot.title = element_text(face = "bold", h = 0.5)) +
  theme(plot.title = element_text(size = 20, h = 0.5, face = "bold"),
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16))

gg_rw_INLA_mag <- gg_rw_INLA + 
  geom_magnify(aes(from = ID > 352 & ID < 373), to = c(150, 310, -1, 9)) +
  geom_magnify(aes(from = ID > 534 & ID < 554), to = c(390, 540, 2, 13)) +
  geom_magnify(aes(from = ID > 895 & ID < 920), to = c(730, 870, -2, 9)) +
  geom_magnify(aes(from = ID > 1073 & ID < 1094), to = c(940, 1070, 5, 16)) +
  geom_magnify(aes(from = ID > 1436 & ID < 1456), to = c(1300, 1410, 5, 16)) +
  geom_magnify(aes(from = ID > 1620 & ID < 1637), to = c(1650, 1800, -2, 9))
  

gg_rw_origcor_INLA <- ggplot() +
  geom_point(data = rinla_orig$summary.random$u[1:length(ysim),][sel_stressp[order(sel_stressp)],], mapping = aes(x = 1:length(sel_stressp)-0.25, y = mean)) +
  geom_errorbar(data = rinla_orig$summary.random$u[1:length(ysim),][sel_stressp[order(sel_stressp)],], aes(x = 1:length(sel_stressp)-0.25, ymin=mean-sd, ymax=mean+sd), width=.2, position = position_dodge(0.05)) +
  geom_point(data = rinla_u.rw1_sel, mapping = aes(x = 1:length(sel_stressp) + 0.25, y = mean), colour = "red") +
  geom_errorbar(data = rinla_u.rw1_sel, aes(x = 1:length(sel_stressp) + 0.25, ymin=mean-sd, ymax=mean+sd), width=.2, position = position_dodge(0.05), colour = "red") +
  geom_segment(data.frame(xi = 1:length(sel_stressp)-0.5, xe = 1:length(sel_stressp)+0.5, y = u.rw1[sel_stressp[order(sel_stressp)]]), mapping = aes(x = xi, xend = xe, y = y), lty = "solid", linewidth = 1, colour = "blue", alpha = 0.75) +
  labs(title = "Latent field stress points (original vs corrected by RF)") + xlab(label = "ID") +
  theme_bw() + # theme(plot.title = element_text(face = "bold", h = 0.5))
  theme(plot.title = element_text(size = 20, h = 0.5, face = "bold"),
        axis.title.x = element_blank(), # element_text(size = 18), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16))

grid.arrange(arrangeGrob(grobs = list(gg_rw_INLA, gg_rw_origcor_INLA), ncol = 2))
