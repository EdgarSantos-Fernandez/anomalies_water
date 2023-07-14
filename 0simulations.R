#Author: Edgar Santos-Fernandez
# edgar.santosfdez@gmail.com

# R Code:
# This set of R codes performs the following tasks:
# 1. Simulates network data and anomalies.
# 2. Fits the models using Markov Chain Monte Carlo (MCMC) methods.
# 3. Predicts anomalies based on the fitted models.
# 4. Computes performance measures to assess the accuracy of the anomaly detection.

# Types of anomalies:
# - Persistent anomalies: These anomalies persist across time at the same spatial location and are caused by measurement errors.

# Note: To avoid re-fitting the MCMC models, please download the necessary files from the following GitHub repository: https://github.com/EdgarSantos-Fernandez/anomalies_water
# After downloading the files, place them in your working directory.

# Important Note: The files "fit_mix1_2022-06-01 120505_.rds" and "fit_mix2_2022-06-09 153355_.rds" are too large to be hosted on GitHub (>100MB).

# To obtain the final outcomes from all the models, you can use the file "ypred_2022-06-10 094358_.rds".

# Please note that running the Stan models may take a few minutes. Your patience is appreciated.






library('Rcpp')
library('dplyr')
library('ggplot2')
library('RColorBrewer')
library('rstan')
library('bayesplot')
library('nlme')
library('SSN')
library('abind')
library('viridis')
library('SSNbayes')
library('qcc')

source('functions.R') # functions to fit the models and utils

RNGkind(sample.kind = "Rounding")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')



# simulated dataset. Reproducible example.
obs_data <- readRDS('obs_data2022-02-22 130519_.rds')


seed <- 20210614 

set.seed(seed)

path <- paste0("./b_", seed  ,".ssn")

if(file.exists(path)){
  n <- importSSN(path, "preds")
  
}  else{ n <- createSSN(n = c(120), 
                        obsDesign = binomialDesign(30), 
                        predDesign = binomialDesign(50),
                        importToR = TRUE, path = path,
                        treeFunction = iterativeTreeLayout)}

x11(); plot(n, lwdLineCol = "addfunccol",  lwdLineEx = 8,
            lineCol = 4,  col = 1,  pch = 16,  xlab = "x-coordinate",  ylab = "y-coordinate")

createDistMat(n, "preds", o.write=TRUE, amongpred = TRUE)



# set to False to generate a new dataset
use_dataset <- T


if(use_dataset == F) { # this will generate new space-time data
 #======================================================
# generating a spatial stream network process with 30 locations, 

rawDFobs <- getSSNdata.frame(n, "Obs")

# generating 3 continous covariates
rawDFobs[,"X1"] <- rnorm(length(rawDFobs[,1]))
rawDFobs[,"X2"] <- rnorm(length(rawDFobs[,1]))
rawDFobs[,"X3"] <- rnorm(length(rawDFobs[,1]))


n <- putSSNdata.frame(rawDFobs,n, Name = 'Obs')

set.seed(seed)

# spatial parameters
#var_tu <- 3
#alpha_tu <- 20
#var_nug <- 0.1

#covariates
#betas <- c(10,-1,0,1) 


sim.out <- SimulateOnSSN( n, ObsSimDF = rawDFobs,
                          formula = ~ X1 + X2 + X3,  coefficients = c(10, 1, 0, -1),
                          CorModels = c("Exponential.taildown"), use.nugget = TRUE,
                          CorParms = c(3, 10, .1), addfunccol = "addfunccol")

sim.ssn <- sim.out$ssn.object
simDFobs <- getSSNdata.frame(sim.ssn, "Obs")



# temporal part
#phi <- 0.8 and t <- 30

df <- getSSNdata.frame(sim.ssn, "Obs")
t <- 30 # days
df <- do.call("rbind", replicate(t, df, simplify = FALSE))# replicating the df
df$date <- rep(1:t, each = (nrow(df)/t)) 
obs_data <- df


set.seed(seed + 1000)
phi <- 0.8 
ar1 <- corAR1(form = ~ unique(obs_data$date), value = phi) # can use corExp for 
AR1 <- Initialize(ar1, data = data.frame(unique(obs_data$date)))

epsilon <- t(chol(corMatrix(AR1))) %*% rnorm(length(unique(obs_data$date)), 0, 1)   #NB AR1 error

#x11(); plot(epsilon[,1] , type = 'l')

epsilon <- rep(epsilon, each = length(unique(obs_data$locID)) ) + 
  rnorm(length(epsilon)*length(unique(obs_data$locID)), 0, 0.5) # for all the dates

obs_data$epsilon <- epsilon

obs_data$y <- obs_data$Sim_Values + obs_data$epsilon
obs_data$pid <- rep(1:nrow(obs_data))



# simulating anomalies
mean_anom <- 5
sd_anom <- 1

obs_data$locID <- as.numeric(as.character(obs_data$locID))

set.seed(seed + 1)

obs_data$ind <-
  sensor_anom( #HERE
    obs = obs_data,
    prop = 0.05,
    lambda= .8
  )

table(obs_data$ind)
prop.table(table(obs_data$ind))

locs <- unique(obs_data$locID)

obs_data$ind_type <- NA
obs_data[obs_data$ind  == 1,]$ind_type <- 'large_spike'
table(obs_data$ind_type)

set.seed(seed+10) 
for(i in 1:length(locs)){
  print(i)
  oo <- obs_data[obs_data$locID == locs[i],]
  for(j in 2:nrow(oo)){
    print(j)
    oo[j,]$ind_type <- ifelse( (oo[j,]$ind == 1 & (oo[j + 1,]$ind == 1 | oo[j-1,]$ind == 1) ), 
                               ifelse( oo[j - 1,]$ind == 0, 
           sample(c( 'shift', 'high_var', 'drift'), 1, replace = T, prob = c(0.333, 0.333, 0.334)), #c(0.2, 0.35, 0.45)
           oo[j - 1,]$ind_type ),oo[j,]$ind_type)
  }
  obs_data[obs_data$locID == i,'ind_type'] <- oo$ind_type 
  
}

obs_data$ind_sd <- 0.0001
obs_data$ind_sd <- ifelse(obs_data$ind_type  == 'high_var', 3 ,  
                          ifelse(obs_data$ind_type  == 'large_spike', 1, obs_data$ind_sd) )

obs_data[is.na(obs_data$ind_sd),]$ind_sd <- 0.0001


set.seed(seed)
obs_data$yobs <- obs_data$ind * rnorm(nrow(obs_data), mean = mean_anom, sd = sd_anom * obs_data$ind_sd) +
  obs_data$y
obs_data$anomaly <- ifelse(obs_data$ind == 1, 'anomalous', 'no_anomalous')


# drifts 
drift_redo <- 0
if(drift_redo == 1){
  obs_data[obs_data$ind_type  == 'drift' & !is.na(obs_data$ind_type), ]$yobs <- 
    obs_data[obs_data$ind_type  == 'drift' & !is.na(obs_data$ind_type), ]$y + 
    obs_data[obs_data$ind_type  == 'drift' & !is.na(obs_data$ind_type), ]$ind * 
    rnorm(nrow(obs_data[obs_data$ind_type  == 'drift' & !is.na(obs_data$ind_type),]), 
          mean = 2, sd = sd_anom * obs_data[obs_data$ind_type  == 'drift' & !is.na(obs_data$ind_type), ]$ind_sd) 
}

# making drifts
set.seed(seed)
for(i in 1:length(locs)){
  print(i)
  oo <- obs_data[obs_data$locID == locs[i],]
  
  if(any(oo$ind_type == 'drift' & !is.na(oo$ind_type))){
    
      oo[oo$ind_type == 'drift' & !is.na(oo$ind_type),]$yobs <- sort(oo[oo$ind_type == 'drift' & !is.na(oo$ind_type),]$yobs)

    obs_data[obs_data$locID == i & obs_data$ind_type == 'drift' & !is.na(obs_data$ind_type) ,'yobs'] <- sort(obs_data[obs_data$locID == i & obs_data$ind_type == 'drift' & !is.na(obs_data$ind_type) ,'yobs']) 
    
    if( sum( oo$ind_type == 'drift' & !is.na(oo$ind_type) ) < 4) {
      print(i)
      obs_data[obs_data$locID == i & obs_data$ind_type == 'drift' & !is.na(obs_data$ind_type) ,'ind_type'] <- 'shift'
    }
    
  }
}


table(obs_data$ind_type)





obs_data_coord <- data.frame(n@obspoints@SSNPoints[[1]]@point.coords)
obs_data_coord$locID <- as.numeric(as.character(factor(1:nrow(obs_data_coord)) ))
obs_data$locID <- as.numeric(as.character(obs_data$locID))
obs_data <- obs_data %>% left_join(obs_data_coord, by = c('locID'))


obs_data$anom <- ifelse(obs_data$anomaly == 'no_anomalous',
                           'normal',
                           obs_data$anomaly) 



}



colo <- brewer.pal(4, 'Set1')
x11(width = 12, height = 9); ggplot(obs_data) + 
  geom_line(aes(x = date, y = yobs, group = locID)) +
  geom_point(aes(x = date, y = yobs, col = ind_type), size = 1.35)+
  facet_wrap(~locID, scales = 'free')+ theme_bw() +
  scale_color_manual(values = colo)+
  ylab('')+ 
  xlab('time point')+
  theme(
    axis.text=element_text(size=12), 
    axis.title=element_text(size=13),
    legend.text=element_text(size=13),
    legend.title=element_text(size=13),
    #axis.text.x = element_text(angle = 45, hjust=1),
    strip.background =element_rect(fill='white'),
    legend.position = "bottom")


nets <- SSNbayes::collapse(n, par = 'addfunccol')
nets$afv_cat <- cut(nets$computed_afv, 
                    breaks = seq(min(nets$computed_afv),
                                 max(nets$computed_afv),
                                 length.out=6),
                    labels = 1:5,
                    include.lowest = T)

library(ggrepel)



x11(width = 12, height = 11); ggplot(nets) + 
  geom_path(aes(X1, X2, group = slot, size = afv_cat), 
            lineend = 'round', linejoin = 'round', col = 'gray')+
  geom_point(data = dplyr::filter(obs_data, date %in% 1:4) ,
             aes(x = coords.x1, y = coords.x2, col = anom),
             size = 3) + # factor(ind)
  geom_text_repel(data = dplyr::filter(obs_data, date %in% 1:4) ,
                  aes(x = coords.x1, y = coords.x2 , label = locID), size = 2.5) + # factor(ind)
  scale_size_manual(values = seq(0.2,2,length.out = 5))+
  scale_color_manual(values = colo)+
  facet_wrap(~date, nrow = 2)+
  xlab("x-coordinate") +
  ylab("y-coordinate")+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    axis.text=element_text(size=12), 
    axis.title=element_text(size=13),
    legend.text=element_text(size=13),
    strip.background =element_rect(fill='white'),
    legend.position = "bottom")+
  guides(size = F)+
  labs(size="", colour="anomaly")



obs_data$date_num <- as.numeric(as.factor(obs_data$date))
obs_data <- obs_data %>% dplyr::select(locID, pid, date, date_num, y, everything())


#saveRDS(obs_data, paste0('obs_data', gsub(":", "", Sys.time()),'_','.rds'))



                


# fitting the model using SSNbayes


# set to False to generate fit a new model
use_saved_model <- T

if(use_saved_model == F) { # 
  
fit <- ssnbayes_ppd(formula = yobs ~ X1 + X2 + X3, 
                   data = obs_data,
                   path = path,
                   time_method = list("ar", "date"),
                   space_method = list('use_ssn', c("Exponential.tailup")),
                   iter = 4000,
                   warmup = 2000,
                   chains = 3,
                   refresh = 50,
                   addfunccol='addfunccol',
                   loglik = T) 
#saveRDS(fit, paste0('fit', gsub(":", "", Sys.time()),'_','.rds'))
} else{  fit <- readRDS('fit_2022-02-22 131437_.rds') }

      

# Summary stats
stats <- summary(fit)$summary %>%
  data.frame()
stats <- data.frame(Variable = rownames(stats), stats)
stats[1:15,]

data.frame(stats[grep("beta\\[", row.names(stats)),'mean'])
data.frame(stats[grep("phi", row.names(stats)),'mean'])


plot_post <- F # plot posterior distr?

if(plot_post == T) {

    ## Create plots of chain results for seven betas
    array <- as.array(fit)
    x11(width = 10, height = 2); mcmc_dens_overlay(
      array,
      pars = c(
        "beta[1]",
        "beta[2]",
        "beta[3]" ,
        "beta[4]" ),
      facet_args = list(nrow = 1)
    ) +
      theme(axis.text = element_text(size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.x = element_text(size = 14))
    
    
    x11(width = 8, height = 2); mcmc_dens_overlay(
      array,
      pars = c(
        "var_tu",
        "var_nug",
        "alpha_tu" ),
      facet_args = list(nrow = 1)
    ) + 
      expand_limits(x = 0)+
      theme(axis.text = element_text(size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.x = element_text(size = 14))
    
    
    ## Plot the distribution of phi
    x11(width = 6, height = 1);mcmc_intervals(
      array,
      pars = paste0("phi"),
      point_size = .1,
      prob_outer = 0.95
    ) +
      theme(
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        strip.text.x = element_text(size = 12)
      )
    
    
    
    x11(width = 6, height = 3); mcmc_intervals(
      array,
      pars = paste0("phi[", 1:30, "]"),
      point_size = .1,
      prob_outer = 0.95
    ) +
      theme(
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        strip.text.x = element_text(size = 12)
      )
    
    
    x11(width = 10, height = 3); mcmc_dens_overlay(
      array,
      pars = c(
        "y_pred[1,1]",
        "y_pred[1,2]",
        "y_pred[1,3]" ,
        "y_pred[1,4]" ),
      facet_args = list(nrow = 1)
    ) +
      theme(axis.text = element_text(size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.x = element_text(size = 14))

}




ypred <- data.frame(stats[grep("y_pred\\[", row.names(stats)),])

ypred$yobs <- obs_data$yobs
ypred$ypred <- ypred$mean
ypred$ind <- obs_data$ind
ypred$date <- obs_data$date
ypred$locID <- obs_data$locID
ypred$anom_ind <- ifelse(ypred$yobs < ypred$X2.5. | 
                       ypred$yobs > ypred$X97.5., 1, 0 )

ypred$anomaly <- ifelse(ypred$ind == 1, 'anomalous', 'no_anomalous')

ypred$y <- obs_data$y


ypred$ind_type <- obs_data$ind_type


table(ypred$anomaly, ypred$anom_ind)

t1 <- table(ypred$ind, ypred$anom_ind) #ypred$ind_type,   
#acc
(acc <- (t1[1,1] + t1[2,2]) / sum(t1) )
#se
(se <- t1[2,2] / (t1[2,2] + t1[2,1]) )
#sp
(sp <- t1[1,1] / (t1[1,1] + t1[1,2]))

(acc_adj <- mean(c(t1[2,2] / (t1[2,2] + t1[2,1]), t1[1,1] / (t1[1,1] + t1[1,2]))) ) # adj Acc

(mcc <-  ((t1[2,2] * t1[1,1] ) - (t1[1,2] * t1[2,1])) / sqrt( (t1[2,2]+t1[1,2])*as.numeric((t1[2,2]+t1[2,1])*(t1[1,1]+t1[1,2])*(t1[1,1]+t1[2,1])  ) ) )






#-ITERATIVE ------------------------------------
# will refit the model after removing the anomalies

# setting the anomalies to NA
obs_data$anom_ind <- NULL
obs_data <- obs_data %>% 
  left_join(ypred[,c('locID', 'date' , 'anom_ind')],
                       by = c('locID', 'date') )

obs_data$yobs2 <- ifelse(obs_data$anom_ind == 1, NA, obs_data$yobs )



# set to False to generate fit a new model
use_saved_model <- T

if(use_saved_model == F) { #
  
  fit2 <- ssnbayes_ppd(formula = yobs2 ~ X1 + X2 + X3,
                       data = obs_data,
                       path = path,
                       time_method = list("ar", "date"), # var
                       space_method = list('use_ssn', c("Exponential.tailup")), # NB:  Exponential.taildown
                       iter = 4000,
                       warmup = 2000,
                       chains = 3,
                       refresh = 50,
                       net = 1, # second network on the ssn object
                       addfunccol='addfunccol',
                       loglik = T) 
  #saveRDS(fit2, paste0('fit2_', gsub(":", "", Sys.time()),'_','.rds'))


} else{  fit2 <- readRDS('fit2_2022-02-22 134939_.rds')}


stats2 <- summary(fit2)$summary %>%
  data.frame()
stats2 <- data.frame(Variable = rownames(stats2), stats2)


data.frame(stats2[grep("beta\\[", row.names(stats2)),'mean'])
data.frame(stats2[grep("phi", row.names(stats2)),'mean'])
data.frame(stats2[grep("alpha", row.names(stats2)),'mean'])
data.frame(stats2[grep("var_tu", row.names(stats2)),'mean'])
data.frame(stats2[grep("var_nug", row.names(stats2)),'mean'])



yp2 <- data.frame(stats2[grep("y_pred\\[", row.names(stats2)),])

ypred$X2.5_2 <- yp2$X2.5.
ypred$X97.5_2 <- yp2$X97.5.
ypred$ypred2 <- yp2$mean


ypred$anom_ind2 <- ifelse(ypred$yobs < ypred$X2.5_2 | 
                           ypred$yobs > ypred$X97.5_2, 1, 0 )

t2 <- table(ypred$ind, ypred$anom_ind2)

(acc2 <- (t2[1,1] + t2[2,2]) / sum(t2) )
#se
(se2 <- t2[2,2] / (t2[2,2] + t2[2,1]) )
#sp
(sp2 <- t2[1,1] / (t2[1,1] + t2[1,2]))

(acc_adj2 <- mean(c(t2[2,2] / (t2[2,2] + t2[2,1]), t2[1,1] / (t2[1,1] + t2[1,2]))) ) # adj Acc

(mcc2 <-  ((t2[2,2] * t2[1,1] ) - (t2[1,2] * t2[2,1])) / sqrt( (t2[2,2]+t2[1,2])*as.numeric((t2[2,2]+t2[2,1])*(t2[1,1]+t2[1,2])*(t2[1,1]+t2[2,1])  ) ) )



# imputed by the model
ys2 <- data.frame(stats2[grep("y\\[", row.names(stats2)),])
ys2$ypred_mod <- ifelse( ys2$sd == 0, NA,  ys2$mean)
ypred$ypred_mod <- ys2$ypred_mod 

#rmse
dplyr::filter(ypred, !is.na(ypred_mod)) %>%
  summarize( rmse = sqrt(mean((yobs - ypred2)^2)),
             rmse2 = sqrt(mean((yobs - ypred_mod)^2)))


plot_post <- F # plot posterior distr?

if(plot_post == T) {
      
    array2 <- as.array(fit2)
    x11(width = 10, height = 3); mcmc_dens_overlay(
      array2,
      pars = c(
        "beta[1]",
        "beta[2]",
        "beta[3]" ,
        "beta[4]" ),
      facet_args = list(nrow = 1)
    ) +
      theme(axis.text = element_text(size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.x = element_text(size = 14))
    
    
    x11(width = 6, height = 1);mcmc_intervals(
      array2,
      pars = paste0("phi"),
      point_size = .1,
      prob_outer = 0.95
    ) +
      theme(
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        strip.text.x = element_text(size = 12)
      )
    
    
    
    x11(width = 8, height = 2); mcmc_dens_overlay(
      array2,
      pars = c(
        "var_tu",
        "var_nug",
        "alpha_tu" ),
      facet_args = list(nrow = 1)
    ) +
      theme(axis.text = element_text(size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.x = element_text(size = 14))
    
}


obs_data <- obs_data %>% left_join(ypred[,c('locID', 'date', 'ypred', 'ypred2')],
                       by = c('locID', 'date') )

cols <- brewer.pal(4, 'Set1')









# Error distribution
# empirical residuals

ypred$er <- ypred$yobs - ypred$ypred 
ypred$er2 <- ypred$yobs - ypred$ypred2 




# Mixtures models

# On iteration 1. 

# set to False to generate fit a new model
use_saved_model <- T

if(use_saved_model == F) {  
    
    clus <- 'data {
      int<lower=0> N;  // number of data points
      int<lower=1> D;  // number of dimensions
      int<lower=1> K;  // number of clusters
      vector[D] y[N];  // observations
    }
    transformed data {
      real<upper=0> neg_log_K;
      neg_log_K = -log(K);
    }
    parameters {
      vector[D] mu[K]; // cluster means
    }
    transformed parameters {
      //real<upper=0> soft_z[N, K]; // log unnormalized clusters
      matrix[N, K] soft_z ;
      matrix[N, K] z; // ESF: normalized clusters
    
      for (n in 1:N)
        for (k in 1:K)
          soft_z[n, k] = neg_log_K  - 0.5 * dot_self(mu[k] - y[n]);
    
      for (n in 1:N){
        z[n, 1] = exp(soft_z[n, 1]) / (exp(soft_z[n, 1]) + exp(soft_z[n, 2]));
        z[n, 2] = exp(soft_z[n, 2]) / (exp(soft_z[n, 1]) + exp(soft_z[n, 2]));
      }
    }
    model {
      // prior
        mu[1] ~ normal(0,1) ; 
        mu[2] ~ normal(5,0.5) ; 
      
      // likelihood
      for (n in 1:N)
        target += log_sum_exp(soft_z[n]);
    }
    '
    
    data_mix2 <- list(N = nrow(ypred),
                     D = 1,
                      K = 2,
                     y = data.frame(y = ypred$er)) 
    

    fit_mix <- rstan::stan(model_code = clus,
                           model_name = "clus",
                           data = data_mix2,
                           iter = 8000,
                           warmup = 6000,
                           chains = 3)
    
   #saveRDS(fit_mix, paste0('fit_mix1', gsub(":", "", Sys.time()),'_','.rds'))

} else{  fit_mix <- readRDS('fit_mix1_2022-06-01 120505_.rds') } #fit_mix2_2022-06-01 114348_

# 



# fit_mix <- readRDS('fit_mix2021-12-01 140931_.RDS')

stats_mix1 <- summary(fit_mix)$summary %>%
  data.frame()
stats_mix1 <- data.frame(Variable = rownames(stats_mix1), stats_mix1)


mix_prop1 <- data.frame(stats_mix1[grep("z\\[", row.names(stats_mix1)),])
#mix_prop2 <- mix_prop2[((nrow(mix_prop2)/2)+1): nrow(mix_prop2),]

mix_prop1p <- data.frame(p1 = mix_prop1[seq(1,nrow(mix_prop1),2),2], 
           p2 = mix_prop1[seq(2,nrow(mix_prop1),2),2])

mix_prop1p$clus <- ifelse(mix_prop1p$p1 < mix_prop1p$p2, 0,1)

table(mix_prop1p$clus)


ypred$mix_ind <- mix_prop1p$clus 


te <- table(ypred$ind, ypred$mix_ind)

(acce <- (te[1,1] + te[2,2]) / sum(te) )
#se
(see <- te[2,2] / (te[2,2] + te[2,1]) )
#sp
(spe <- te[1,1] / (te[1,1] + te[1,2]))

(acc_adje <- mean(c(te[2,2] / (te[2,2] + te[2,1]), te[1,1] / (te[1,1] + te[1,2]))) ) # adj Acc

(mcce <-  ((te[2,2] * te[1,1] ) - (te[1,2] * te[2,1])) / sqrt( (te[2,2]+te[1,2])*as.numeric((te[2,2]+te[2,1])*(te[1,1]+te[1,2])*(te[1,1]+te[2,1])  ) ) )




# mixture
# On iteration 2. 


use_saved_model <- T

if(use_saved_model == F) {  
  
    clus <- 'data {
      int<lower=0> N;  // number of data points
      int<lower=1> D;  // number of dimensions
      int<lower=1> K;  // number of clusters
      vector[D] y[N];  // observations
    }
    transformed data {
      real<upper=0> neg_log_K;
      neg_log_K = -log(K);
    }
    parameters {
      vector[D] mu[K]; // vector[D] mu[K]  cluster means
    }
    transformed parameters {
      //real<upper=0> soft[N, K]; // log unnormalized clusters
      matrix[N, K] soft ;
      matrix[N, K] z; // ESF: normalized clusters
    
      for (n in 1:N)
        for (k in 1:K)
          soft[n, k] = neg_log_K  - 0.5 * dot_self(mu[k] - y[n]);
    
      for (n in 1:N){
        z[n, 1] = exp(soft[n, 1]) / (exp(soft[n, 1]) + exp(soft[n, 2]));
        z[n, 2] = exp(soft[n, 2]) / (exp(soft[n, 1]) + exp(soft[n, 2]));
      }
    }
    model {
      // prior
        mu[1,1] ~ normal(-1,0.5) ; // std_normal();
        mu[2,1] ~ normal(10,0.25) ; // std_normal();
      
      // likelihood
      for (n in 1:N)
        target += log_sum_exp(soft[n]);
    }
    '
    
    
    data_mix2 <- list(N = nrow(ypred),
                      D = 1,
                      K = 2,
                      y = data.frame(y = ypred$er2)) # NB er_abs
    
    
    fit_mix2 <- rstan::stan(model_code = clus,
                            model_name = "clus",
                            data = data_mix2,
                            pars = c('z', 'mu'),
                            iter = 8000,
                            warmup = 6000,
                            chains = 2,
                            seed = 10)
    #saveRDS(fit_mix2, paste0('fit_mix2_', gsub(":", "", Sys.time()),'_','.rds'))
    
    
    } else{ 
      fit_mix2 <- readRDS('fit_mix2_2022-06-09 153355_.rds') } #fit_mix2_2022-06-01 114348_
    


#save.image('image.RData') 
# load('image.RData')

stats_mix2 <- summary(fit_mix2)$summary %>%
  data.frame()
stats_mix2 <- data.frame(Variable = rownames(stats_mix2), stats_mix2)


mix_prop2 <- data.frame(stats_mix2[grep("z\\[", row.names(stats_mix2)),])
#mix_prop3 <- mix_prop3[((nrow(mix_prop3)/2)+1): nrow(mix_prop3),]

mix_prop2p <- data.frame(p1 = mix_prop2[seq(1,nrow(mix_prop2),2),2], 
                         p2 = mix_prop2[seq(2,nrow(mix_prop2),2),2])

mix_prop2p$clus <- ifelse(mix_prop2p$p1 < mix_prop2p$p2, 0,1)

table(mix_prop2p$clus)

ypred$mix_ind2 <- mix_prop2p$clus 





te <- table(ypred$ind, ypred$mix_ind2)


(acce1 <- (te[1,1] + te[2,2]) / sum(te) )
#se
(see1 <- te[2,2] / (te[2,2] + te[2,1]) )
#sp
(spe1 <- te[1,1] / (te[1,1] + te[1,2]))

(acc_adje1 <- mean(c(te[2,2] / (te[2,2] + te[2,1]), te[1,1] / (te[1,1] + te[1,2]))) ) # adj Acc

(mcce1 <-  ((te[2,2] * te[1,1] ) - (te[1,2] * te[2,1])) / sqrt( (te[2,2]+te[1,2])*as.numeric((te[2,2]+te[2,1])*(te[1,1]+te[1,2])*(te[1,1]+te[2,1])  ) ) )






# Hidden Markov model 


library(depmixS4)
set.seed(1)
mod <- depmix(response = er ~ 1, data = ypred, nstates = 2) # turb  NB: cond  #,trstart = runif(4)

fm <- fit(mod, verbose = FALSE) 
est.states <- posterior(fm)

ypred$anom <- ifelse(est.states$state == 1, 'anom', 'normal')
table(ypred$anom )


ypred$hmm_ind <- ifelse(ypred$anom == 'anom', 1, 0)

tab_hmm <- table(ypred$ind, ypred$hmm_ind)


(acce2 <- (tab_hmm[1,1] + tab_hmm[2,2]) / sum(tab_hmm) )
#se
(see2 <- tab_hmm[2,2] / (tab_hmm[2,2] + tab_hmm[2,1]) )
#sp
(spe2 <- tab_hmm[1,1] / (tab_hmm[1,1] + tab_hmm[1,2]))

(acc_adje2 <- mean(c(tab_hmm[2,2] / (tab_hmm[2,2] + tab_hmm[2,1]), tab_hmm[1,1] / (tab_hmm[1,1] + tab_hmm[1,2]))) ) # adj Acc

(mcce2 <-  ((tab_hmm[2,2] * tab_hmm[1,1] ) - (tab_hmm[1,2] * tab_hmm[2,1])) / sqrt( (tab_hmm[2,2]+tab_hmm[1,2])*as.numeric((tab_hmm[2,2]+tab_hmm[2,1])*(tab_hmm[1,1]+tab_hmm[1,2])*(tab_hmm[1,1]+tab_hmm[2,1])  ) ) )



# HMM Bayesian with sigma specific  
# credit: sections of this code come from the Stan guide

use_saved_model <- T

if(use_saved_model == F) {  
  
    hmm <- 'data {
      int<lower=0> N;
      int<lower=0> K;
      real y[N];
    }
    
    parameters {
      simplex[K] theta[K];
      positive_ordered[K] mu;
      vector<lower=0> [K] sigma; // SD of the state
    real<lower=0> v;
      }
    
    transformed parameters {
    }
    model {
      // priors
      sigma[1] ~ uniform(0,2);
      sigma[2] ~ uniform(2,10);
      v ~ normal(5,1); // ESF
      target+= normal_lpdf(mu[1] | 0, .5); 
      mu[2] ~ normal(mu[1] + v, 0.25); // ESF
      // forward algorithm
      {
        real acc[K];
        real gamma[N, K];
        for (k in 1:K)
          gamma[1, k] = normal_lpdf(y[1] | mu[k], sigma[k]); //
        for (t in 2:N) {
          for (k in 1:K) {
            for (j in 1:K)
              acc[j] = gamma[t-1, j] + log(theta[j, k]) + normal_lpdf(y[t] | mu[k], sigma[k]);
            gamma[t, k] = log_sum_exp(acc);
          }
        }
        target += log_sum_exp(gamma[N]);
      }
    }
    
    generated quantities {
      int<lower=1,upper=K> z_star[N];
      real log_p_z_star;
      {
        int back_ptr[N, K];
        real best_logp[N, K];
        for (k in 1:K)
          best_logp[1, k] = normal_lpdf(y[1] | mu[k], sigma[k]);
        for (t in 2:N) {
          for (k in 1:K) {
            best_logp[t, k] = negative_infinity();
            for (j in 1:K) {
              real logp;
              logp = best_logp[t-1, j] + log(theta[j, k]) + normal_lpdf(y[t] | mu[k], sigma[k]);
              if (logp > best_logp[t, k]) {
                back_ptr[t, k] = j ;
                best_logp[t, k] = logp;
              }
            }
          }
        }
        log_p_z_star = max(best_logp[N]);
        for (k in 1:K)
          if (best_logp[N, k] == log_p_z_star)
            z_star[N] = k;
        for (t in 1:(N - 1))
          z_star[N - t] = back_ptr[N - t + 1, z_star[N - t + 1]];
      }
    }'
    
    
    
    rstan_options(auto_write = TRUE)
    
    
    stan_data <- list(N = length(ypred$er2), # NB: er_abs2  was er_abs
                      K = 2,
                      y = ypred$er2)
    
    hmm_fit <- stan(model_code = hmm,
                    data = stan_data, iter = 4000, 
                    warmup = 2000, 
                    chains = 3,
                    refresh = 100,
                    seed = 100)
    
    
    saveRDS(hmm_fit, "hmm_sd_20220601.RDS")
} else{ 
  hmm_fit <- readRDS("hmm_sd_20220601.RDS") } 



stats_hmm <- summary(hmm_fit)$summary %>%
  data.frame()
stats_hmm <- data.frame(Variable = rownames(stats_hmm), stats_hmm)


stats_hmm2 <- data.frame(stats_hmm[grep("z_star\\[", row.names(stats_hmm)),])

stats_hmm2$hmm_state <- ifelse(stats_hmm2$mean < 1.5, 0, 1)  

ypred$hmm_state <- stats_hmm2$hmm_state 



th <- table(ypred$ind, ypred$hmm_state)


(acce6 <- (th[1,1] + th[2,2]) / sum(th) )
#se
(see6 <- th[2,2] / (th[2,2] + th[2,1]) )
#sp
(spe6 <- th[1,1] / (th[1,1] + th[1,2]))

(acc_adje6 <- mean(c(th[2,2] / (th[2,2] + th[2,1]), th[1,1] / (th[1,1] + th[1,2]))) ) # adj Acc

(mcc6 <-  ((th[2,2] * th[1,1] ) - (th[1,2] * th[2,1])) / sqrt( (th[2,2]+th[1,2])*as.numeric((th[2,2]+th[2,1])*(th[1,1]+th[1,2])*(th[1,1]+th[2,1])  ) ) )




# HMM Bayesian END







# ARIMA START
library(forecast)

locs <- unique(obs_data$locID)

obs_data$anom_arima <- NA

for(j in 1:length(locs)) {
  print(j)
  s1 <- data.frame(yobs = obs_data[obs_data$locID == locs[j],]$yobs) 
  print(auto.arima(s1, stepwise = F, approx = F, D=0, seasonal = TRUE)) # best parameters in the ARIMA model
  for(i in 5:nrow(s1)){ 
    print(i)
    s1_ <- s1[1:(i - 1),'yobs']  
    pred <- Arima(s1_, order = c(0,1,0), lambda = NULL) %>%
      forecast(h=1) 
    obs_data[obs_data$locID == locs[j],]$anom_arima[i] <- ifelse(s1[i,'yobs'] > pred$lower[2] & s1[i,'yobs'] < pred$upper[2] , 0, 1)
  }
}


ypred$anom_arima <- obs_data$anom_arima

ta <- table(obs_data$ind, obs_data$anom_arima)

(acce5 <- (ta[1,1] + ta[2,2]) / sum(ta) )
#se
(see5 <- ta[2,2] / (ta[2,2] + ta[2,1]) )
#sp
(spe5 <- ta[1,1] / (ta[1,1] + ta[1,2]))

(acc_adje5 <- mean(c(ta[2,2] / (ta[2,2] + ta[2,1]), ta[1,1] / (ta[1,1] + ta[1,2]))) ) # adj Acc

(mcc5 <-  ((ta[2,2] * ta[1,1] ) - (ta[1,2] * ta[2,1])) / sqrt( (ta[2,2]+ta[1,2])*as.numeric((ta[2,2]+ta[2,1])*(ta[1,1]+ta[1,2])*(ta[1,1]+ta[2,1])  ) ) )



# ARIMA method END





# univariate control charts 




locs <- unique(ypred$locID)

ypred$cusum <- 0
df_all <- NULL 


ypred$pid <- obs_data$pid

for(i in 1:length(locs)){
  print(i)
  
  qccsub <- filter(ypred, locID == locs[i], !is.na(er2)) %>% dplyr::select(pid, er2, yobs)

  cu <- cusum(qccsub$er2, chart.all=FALSE, sizes = 1, center=0, std.dev = "MR", k = 2)

  vcu <- cu$violations
  vcu <- unlist(vcu)
  
  qccsub$cusum <- 0
  qccsub[vcu,'cusum'] <- 1
  
  qccsub$locID <- i
  df_all <- rbind(df_all, qccsub)
  }

df_all %>% group_by(locID) %>% summarise(mean(cusum))


ypred$cusum <- NULL
ypred_back <- ypred
ypred <- ypred %>% left_join(df_all, b = c('pid', "locID", 'er2', 'yobs'))

ypred %>% group_by(locID) %>% summarise(mean(cusum, na.rm = T))
table(ypred$cusum)

ypred$cusum2 <- ifelse(ypred$cusum == 1, 'anom', 'non_anom' )

# combining PPD2 with cusum
ypred$anom_ind_qcc <- ifelse(ypred$anom_ind2 == 1 | ypred$cusum == 1, 1, 0)


t_qcc <- table(ypred$ind, ypred$anom_ind_qcc) # table(ypred$ind, ypred$xbar_mr2) # 


#se
(se2qcc <- t_qcc[2,2] / (t_qcc[2,2] + t_qcc[2,1]) )
#sp
(sp2qcc <- t_qcc[1,1] / (t_qcc[1,1] + t_qcc[1,2]))

(acc2qcc <- (t_qcc[1,1] + t_qcc[2,2]) / sum(t_qcc) )

(acc_adj2qcc <- mean(c(t_qcc[2,2] / (t_qcc[2,2] + t_qcc[2,1]), t_qcc[1,1] / (t_qcc[1,1] + t_qcc[1,2]))) ) # adj Acc

(mcc2qcc <-  ((t_qcc[2,2] * t_qcc[1,1] ) - (t_qcc[1,2] * t_qcc[2,1])) / sqrt( (t_qcc[2,2]+t_qcc[1,2])*as.numeric((t_qcc[2,2]+t_qcc[2,1])*(t_qcc[1,1]+t_qcc[1,2])*(t_qcc[1,1]+t_qcc[2,1])  ) ) )



# performance measures
(pm <- data.frame(iter = c('arima','iter_1', 'iter_2', 'mix', 'mix_2','hmm', 'hmm_2', 'ppd_i2_qcc'),
                   se = c(see5, se,se2,see, see1, see2,see6, se2qcc),
                   sp = c(spe5, sp, sp2,spe,spe1, spe2,spe6,sp2qcc),
                   acc = c(acce5, acc, acc2,acce,acce1, acce2,acce6,acc2qcc),
                   acc_adj = c(acc_adje5, acc_adj, acc_adj2,acc_adje,acc_adje1, acc_adje2,acc_adje6, acc_adj2qcc),
                   mcc = c(mcc5, mcc, mcc2, mcce, mcce1, mcce2, mcc6, mcc2qcc)) )
library(xtable)
print(xtable(
  pm, type = "latex", 
  digits=c(4,4,4,4,4,4,4),
), include.rownames=F,  file = "perf_measures_ALL_20220609.tex") 


saveRDS(ypred, paste0('ypred_', gsub(":", "", Sys.time()),'_','.rds'))





# performance by anomaly group START 

table(ypred$ind_type)



ind_types <- c('large_spike',  'high_var', 'shift',  'drift')

ALL <- NULL

for(i in 1:length(ind_types)) {
  print(i)
  ypred2 <- dplyr::filter(ypred, ind_type == ind_types[i]| is.na(ind_type))
  table(ypred2$ind_type, useNA = 'always')
  
  #PPD 1
  
  t1 <- table(ypred2$ind, ypred2$anom_ind)   
  #acc
  (acc <- (t1[1,1] + t1[2,2]) / sum(t1) )
  #se
  (se <- t1[2,2] / (t1[2,2] + t1[2,1]) )
  #sp
  (sp <- t1[1,1] / (t1[1,1] + t1[1,2]))
  
  (acc_adj <- mean(c(t1[2,2] / (t1[2,2] + t1[2,1]), t1[1,1] / (t1[1,1] + t1[1,2]))) ) # adj Acc
  
  (mcc <-  ((t1[2,2] * t1[1,1] ) - (t1[1,2] * t1[2,1])) / sqrt( (t1[2,2]+t1[1,2])*as.numeric((t1[2,2]+t1[2,1])*(t1[1,1]+t1[1,2])*(t1[1,1]+t1[2,1])  ) ) )
  
  
  
  
  #PPD 2
  
  t2 <- table(ypred2$ind, ypred2$anom_ind2)
  #acc
  (t2[1,1] + t2[2,2]) / sum(t2) 
  #se
  t2[2,2] / (t2[2,2] + t2[2,1])
  #sp
  t2[1,1] / (t2[1,1] + t2[1,2])
  
  
  (acc2 <- (t2[1,1] + t2[2,2]) / sum(t2) )
  #se
  (se2 <- t2[2,2] / (t2[2,2] + t2[2,1]) )
  #sp
  (sp2 <- t2[1,1] / (t2[1,1] + t2[1,2]))
  
  (acc_adj2 <- mean(c(t2[2,2] / (t2[2,2] + t2[2,1]), t2[1,1] / (t2[1,1] + t2[1,2]))) ) # adj Acc
  
  (mcc2 <-  ((t2[2,2] * t2[1,1] ) - (t2[1,2] * t2[2,1])) / sqrt( (t2[2,2]+t2[1,2])*as.numeric((t2[2,2]+t2[2,1])*(t2[1,1]+t2[1,2])*(t2[1,1]+t2[2,1])  ) ) )
  
  
  
  
  # mixtures 1
  
  te <- table(ypred2$ind, ypred2$mix_ind)

  
  (acce <- (te[1,1] + te[2,2]) / sum(te) )
  #se
  (see <- te[2,2] / (te[2,2] + te[2,1]) )
  #sp
  (spe <- te[1,1] / (te[1,1] + te[1,2]))
  
  (acc_adje <- mean(c(te[2,2] / (te[2,2] + te[2,1]), te[1,1] / (te[1,1] + te[1,2]))) ) # adj Acc
  
  (mcce <-  ((te[2,2] * te[1,1] ) - (te[1,2] * te[2,1])) / sqrt( (te[2,2]+te[1,2])*as.numeric((te[2,2]+te[2,1])*(te[1,1]+te[1,2])*(te[1,1]+te[2,1])  ) ) )
  
  
  
  
  
  # mixtures 2
  te <- table(ypred2$ind, ypred2$mix_ind2)
  
  (acce1 <- (te[1,1] + te[2,2]) / sum(te) )
  #se
  (see1 <- te[2,2] / (te[2,2] + te[2,1]) )
  #sp
  (spe1 <- te[1,1] / (te[1,1] + te[1,2]))
  
  (acc_adje1 <- mean(c(te[2,2] / (te[2,2] + te[2,1]), te[1,1] / (te[1,1] + te[1,2]))) ) # adj Acc
  
  (mcce1 <-  ((te[2,2] * te[1,1] ) - (te[1,2] * te[2,1])) / sqrt( (te[2,2]+te[1,2])*as.numeric((te[2,2]+te[2,1])*(te[1,1]+te[1,2])*(te[1,1]+te[2,1])  ) ) )
  
  
  
  
  # hmm
  tab_hmm <- table(ypred2$ind, ypred2$hmm_ind)
  
  (acce2 <- (tab_hmm[1,1] + tab_hmm[2,2]) / sum(tab_hmm) )
  #se
  (see2 <- tab_hmm[2,2] / (tab_hmm[2,2] + tab_hmm[2,1]) )
  #sp
  (spe2 <- tab_hmm[1,1] / (tab_hmm[1,1] + tab_hmm[1,2]))
  
  (acc_adje2 <- mean(c(tab_hmm[2,2] / (tab_hmm[2,2] + tab_hmm[2,1]), tab_hmm[1,1] / (tab_hmm[1,1] + tab_hmm[1,2]))) ) # adj Acc
  
  (mcce2 <-  ((tab_hmm[2,2] * tab_hmm[1,1] ) - (tab_hmm[1,2] * tab_hmm[2,1])) / sqrt( (tab_hmm[2,2]+tab_hmm[1,2])*as.numeric((tab_hmm[2,2]+tab_hmm[2,1])*(tab_hmm[1,1]+tab_hmm[1,2])*(tab_hmm[1,1]+tab_hmm[2,1])  ) ) )
  
  
  
  # arima
  ta <- table(ypred2$ind, ypred2$anom_arima)
  
  (acce5 <- (ta[1,1] + ta[2,2]) / sum(ta) )
  #se
  (see5 <- ta[2,2] / (ta[2,2] + ta[2,1]) )
  #sp
  (spe5 <- ta[1,1] / (ta[1,1] + ta[1,2]))
  
  (acc_adje5 <- mean(c(ta[2,2] / (ta[2,2] + ta[2,1]), ta[1,1] / (ta[1,1] + ta[1,2]))) ) # adj Acc
  
  (mcc5 <-  ((ta[2,2] * ta[1,1] ) - (ta[1,2] * ta[2,1])) / sqrt( (ta[2,2]+ta[1,2])*as.numeric((ta[2,2]+ta[2,1])*(ta[1,1]+ta[1,2])*(ta[1,1]+ta[2,1])  ) ) )
  
  
  
  #
  
  # HMM bayes
  
  th <- table(ypred2$ind, ypred2$hmm_state)
  
  (acce6 <- (th[1,1] + th[2,2]) / sum(th) )
  #se
  (see6 <- th[2,2] / (th[2,2] + th[2,1]) )
  #sp
  (spe6 <- th[1,1] / (th[1,1] + th[1,2]))
  
  (acc_adje6 <- mean(c(th[2,2] / (th[2,2] + th[2,1]), th[1,1] / (th[1,1] + th[1,2]))) ) # adj Acc
  
  (mcc6 <-  ((th[2,2] * th[1,1] ) - (th[1,2] * th[2,1])) / sqrt( (th[2,2]+th[1,2])*as.numeric((th[2,2]+th[2,1])*(th[1,1]+th[1,2])*(th[1,1]+th[2,1])  ) ) )
  
  
  # Cusum
  t_qcc <- table(ypred2$ind, ypred2$anom_ind_qcc)

  #se
  (se2qcc <- t_qcc[2,2] / (t_qcc[2,2] + t_qcc[2,1]) )
  #sp
  (sp2qcc <- t_qcc[1,1] / (t_qcc[1,1] + t_qcc[1,2]))
  
  (acc2qcc <- (t_qcc[1,1] + t_qcc[2,2]) / sum(t_qcc) )
  
  (acc_adj2qcc <- mean(c(t_qcc[2,2] / (t_qcc[2,2] + t_qcc[2,1]), t_qcc[1,1] / (t_qcc[1,1] + t_qcc[1,2]))) ) # adj Acc
  
  (mcc2qcc <-  ((t_qcc[2,2] * t_qcc[1,1] ) - (t_qcc[1,2] * t_qcc[2,1])) / sqrt( (t_qcc[2,2]+t_qcc[1,2])*as.numeric((t_qcc[2,2]+t_qcc[2,1])*(t_qcc[1,1]+t_qcc[1,2])*(t_qcc[1,1]+t_qcc[2,1])  ) ) )
  
  
  (pm6 <- data.frame(iter = c('arima','iter_1', 'iter_2', 'mix', 'mix_2','hmm', 'hmm_2', 'iter_2_qcc'),
                     se = c(see5, se,se2,see, see1, see2,see6, se2qcc),
                     sp = c(spe5, sp, sp2,spe,spe1, spe2,spe6,sp2qcc),
                     acc = c(acce5, acc, acc2,acce,acce1, acce2,acce6,acc2qcc),
                     acc_adj = c(acc_adje5, acc_adj, acc_adj2,acc_adje,acc_adje1, acc_adje2,acc_adje6,acc_adj2qcc),
                     mcc = c(mcc5, mcc, mcc2, mcce, mcce1, mcce2, mcc6,mcc2qcc)) )
  
  pm6$ind_type <- ind_types[i]
  
  pm6 <- dplyr::select(pm6, ind_type, everything())
  
  ALL <- rbind(ALL, pm6)
  
}



print(xtable(
  ALL, type = "latex", 
  digits=c(4, 4,4,4,4,4,4,4),
), include.rownames=F,  file = "by_type_anom_ALL_perf_measures_20220609.tex") #digits = 4


# performance by anomaly group END *******************


#saveRDS(ypred, paste0('ypred_', gsub(":", "", Sys.time()),'_','.rds'))


# EOF ****************************************************************



