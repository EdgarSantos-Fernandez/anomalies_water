#Author: Edgar Santos-Fernandez
# edgar.santosfdez@gmail.com


#persistent anomalies
sensor_anom <- function(obs = obs, prop = 0.05, lambda = 3, seed = 1){
  
  obs <- arrange(obs,locID)
  
  set.seed(seed+1)
  obs$IND <- sample(0:1, nrow(obs), replace = T, prob = c(1-prop,prop))
  
  
  for(i in 1:nrow(obs)){
    print(i)
    if(obs$IND[i] == 1) {
      
      loc <- obs$locID[i]
      dat <- obs$date[i]
      
      obs[obs$locID %in% loc & 
          obs$date %in% (dat:(min(dat + rpois(1,lambda), max(obs$date))) )
          ,'IND'] <- 1
      i = i + max(obs$locID)
    }
  
  }
  obs <- arrange(obs,pid)
  
  obs$IND
}


# generates the posterior predictive distributions from an SSNbayes model

ssnbayes_ppd <- function(formula = formula,
                         data = data,
                         path = path,
                         time_method = time_method, # list("ar", "date")
                         space_method = space_method, #list('use_ssn', 'Exponential.tailup'),
                         iter = 3000,
                         warmup = 1500,
                         chains = 3,
                         refresh = max(iter/100, 1),
                         net = 1,
                         addfunccol = addfunccol,
                         loglik = F,
                         seed = seed
){
  
  # checks
  if(missing(time_method)){
    stop("Need to define the method (ar or var) and the column associated with time")
  }
  
  if(length(time_method) == 1){
    stop("Need to specify the column in the the data with the time variable")
  }
  
  time_points <- time_method[[2]]
  
  #if('date' %in% names(data) == F) stop("There is no column date on the data. Please, set a column called date with the time")
  if('locID' %in% names(data) == F) stop("There is no column locID on the data. Please, set a column called locID with the observation locations")
  
  if(missing(seed)) seed <- sample(1:1E6,1,replace=T)
  
  if(!missing(space_method)){
    print('using SSN object')
    if(space_method[[1]] == 'use_ssn'){
      ssn_object <- T
      
      
      if(length(space_method) > 1){
        if(space_method[[2]] %in% c("Exponential.tailup", "LinearSill.tailup" , "Spherical.tailup" ,
                                    "Exponential.taildown" ,"LinearSill.taildown" ,"Spherical.taildown",
                                    "Exponential.Euclid") == F) {stop("Need to specify one or more of the following covariance matrices: Exponential.tailup, LinearSill.tailup , Spherical.tailup ,
		Exponential.taildown, LinearSill.taildown, Spherical.taildown or Exponential.Euclid")}
        CorModels <- space_method[[2]]
      }
      if(length(space_method) == 1){
        CorModels <- "Exponential.tailup"
        print('using an Exponential.tailup model')
      }
      
    }
    if(space_method[[1]] == 'no_ssn'){
      print('no SSN object defined')
      ssn_object <- F
      if(space_method[[2]] %in% c("Exponential.Euclid") == F) {stop("Need to specify Exponential.Euclid")}
      # when using Euclidean distance, need to specify the columns with long and lat.
      if(length(space_method) < 3){ stop("Please, specify the columns in the data frame with the latitude and longitude (c('lon', 'lat'))") }
      
      data$lon <- data[,names(data) == space_method[[3]][1]]
      data$lat <- data[,names(data) == space_method[[3]][2]]
      CorModels <- space_method[[2]]
    }
    
    
  }
  
  if(missing(space_method)) {space_method <- 'no_ssn'; ssn_object <- F; CorModels <- "Exponential.Euclid" }# if missing use Euclidean distance
  
  
  
  
  # Cov
  cor_tu <- case_when(CorModels == "Exponential.tailup" ~ 1,
                      CorModels == "LinearSill.tailup" ~ 2,
                      CorModels == "Spherical.tailup" ~ 3,
                      TRUE ~ 5)
  cor_tu <- sort(cor_tu)[1]
  
  cor_td <- case_when(CorModels == "Exponential.taildown" ~ 1,
                      CorModels == "LinearSill.taildown" ~ 2,
                      CorModels == "Spherical.taildown" ~ 3,
                      TRUE ~ 5)
  cor_td <- sort(cor_td)[1]
  
  
  cor_ed <- case_when(CorModels == "Exponential.Euclid" ~ 1,
                      TRUE ~ 5)
  #CorModels == "Spherical.Euclid" ~ 2,
  #CorModels == "Gaussian.Euclid" ~ 3)
  cor_ed <- sort(cor_ed)[1]
  
  cor_re <- case_when(CorModels == "RE1" ~ 1,
                      TRUE ~ 5)
  cor_re <- sort(cor_re)[1]
  
  
  
  data_com <-  'data {
    int<lower=1> N;
    int<lower=1> K;
    int<lower=1> T;
    matrix[N,K] X[T] ; // real X[N,K,T]; //

    int<lower = 0> N_y_obs; // number observed values
    int<lower = 0> N_y_mis; // number missing values

    int<lower = 1> i_y_obs[N_y_obs] ;  //[N_y_obs,T]
    int<lower = 1> i_y_mis[N_y_mis] ;  // N_y_mis,T]

    vector[N_y_obs] y_obs;    //matrix[N_y_obs,1] y_obs[T];

    matrix[N, N] W ; // spatial weights
    matrix[N, N] h ; // total hydrological dist
    matrix[N, N] I ; // diag matrix

    matrix[N, N] D ; // downstream hydrological dist matrix
    matrix[N, N] flow_con_mat; // flow conected matrix

    matrix[N, N] e ; // Euclidean dist mat
    real<lower=0.01> alpha_max ;

  }'
  
  
  param_com <- '
  parameters {
    vector[K] beta;
    real<lower=0> sigma_nug;

    vector[N_y_mis] y_mis;//declaring the missing y
  '
  
  
  param_phi_ar <- '
    real <lower=-1, upper = 1> phi; // NBNB
  '
  param_phi_var <- '
    vector<lower=-1, upper = 1> [N] phi  ; // vector of autoregresion pars
  '
  
  
  
  param_tu <- '
    real<lower=0> sigma_tu;
    real<lower=0> alpha_tu;
'
  
  param_td <- '
    real<lower=0> sigma_td; // sd of tail-down
    real<lower=0> alpha_td; // range of the tail-down model
'
  
  param_ed <- '
    real<lower=0> sigma_ed;
    real<lower=0> alpha_ed; // range of the Euclidean dist model
'
  
  param_re <- '
    real<lower=0> sigma_RE1;
'
  
  tparam_com <- '
  transformed parameters {
    vector[N * T] y;
    vector[N] Y[T];

    vector[N] epsilon[T]; // error term
    vector[N] mu[T]; // mean

   real<lower=0> var_nug; // nugget

   matrix[N, N] C_tu; //tail-up cov
   matrix[N, N] C1; //tail-up cov
   matrix[N, N] Ind; //tail-up indicator

   matrix[N, N] C_td; //tail-down cov
   matrix[N, N] Ind2; //tail-down indicator
   matrix[2,1] iji;

   matrix[N, N] C_ed ;// Euclidean cov

   matrix[N, N] C_re ;// random effect cov
   matrix[N, N] RE1; // random effect 1

   '
  
  tparam_tu <- '
   // tail up exponential
   real<lower=0> var_tu; // parsil tail-down
  '
  
  tparam_td <- '
    real<lower=0> var_td; // parsil tail-down
 '
  
  tparam_ed <- '
    real<lower=0> var_ed; //  Euclidean dist var
'
  
  tparam_re <- '
    real<lower=0> var_RE1; // Random effect 1
'
  
  
  tparam_com2 <- '
    y[i_y_obs] = y_obs;
    y[i_y_mis] = y_mis;
    for (t in 1:T){
      //y[ i_y_obs[,t] ] = y_obs[((t - 1) * N_y_obs + 1):(t * N_y_obs)];
      //y[ i_y_mis[,t] ] = y_mis[((t - 1) * N_y_mis + 1):(t * N_y_mis)];
      Y[t] = y[((t - 1) * N + 1):(t * N)];
    }

    var_nug = sigma_nug ^ 2; // variance nugget
        mu[1] = X[1] * beta;
    epsilon[1] = Y[1] - mu[1];

'
  
  tparam_com_ar <- '
        for (t in 2:T){
        mu[t] = X[t] * beta;
        epsilon[t] = Y[t] - mu[t];
        mu[t] = mu[t] + phi * epsilon[t-1]; //
    }

  '
  
  tparam_com_var <- '
        for (t in 2:T){
        mu[t] = X[t] * beta;
        epsilon[t] = Y[t] - mu[t];
        mu[t] = mu[t] + phi .* epsilon[t-1]; // element wise mult two vectors
    }
  '
  
  
  tparam_tu2_exp <- '
    // tail up exponential
    var_tu = sigma_tu ^ 2; // variance tail-up
      C1 = var_tu * exp(- 3 * h / alpha_tu); // tail up exponential model
    C_tu = C1 .* W; // Hadamard (element-wise) product
'
  
  tparam_tu2_lin <- '
     //Tail-up linear-with-sill model
     var_tu = sigma_tu ^ 2; // variance tail-up
      for (i in 1:N) {
        for (j in 1:N) {
        Ind[i,j] = (h[i,j] / alpha_tu) <= 1 ? 1 : 0; // indicator
        }
      }
      C1 = var_tu * (1 - (h / alpha_tu)) .* Ind ; //Tail-up linear-with-sill model
    C_tu = C1 .* W; // Hadamard (element-wise) product
'
  
  tparam_tu2_sph <- '
      // Tail-up spherical model
      var_tu = sigma_tu ^ 2; // variance tail-up
      for (i in 1:N) {// Tail-up spherical model
        for (j in 1:N) {
          Ind[i,j] = (h[i,j] / alpha_tu) <= 1 ? 1 : 0; // indicator
        }
      }
    C1 = var_tu * (1 - (1.5 * h / alpha_tu) + (h .* h .* h / (2 * alpha_tu ^ 3))) .* Ind ; // Tail-up spherical model
    C_tu = C1 .* W; // Hadamard (element-wise) product
'
  #tail-up models end
  
  
  #tail-down models start
  
  # tail-down exponential
  tparam_td2_exp <- '
    var_td= sigma_td ^ 2; // variance tail-down
     	for (i in 1:N) {// Tail-down exponential model
          for (j in 1:N) {
            if(flow_con_mat[i,j] == 1){ // if points are flow connected
               C_td[i,j] = var_td * exp(- 3 * h[i,j] / alpha_td);
            }
            else{// if points are flow unconnected
              C_td[i,j] = var_td * exp(- 3 * (D[i,j] + D[j,i]) / alpha_td);
            }
          }
        }

'
  
  #Tail-down linear-with-sill model
  tparam_td2_lin <- '
    var_td= sigma_td ^ 2; // variance tail-down
        for (i in 1:N) {// Tail-down linear-with-sill model
          for (j in 1:N) {
            if(flow_con_mat[i,j] == 1){ // if points are flow connected
              Ind2[i,j] = (h[i,j] / alpha_td) <= 1 ? 1 : 0; // indicator
              C_td[i,j] = var_td * (1 - (h[i,j] / alpha_td)) .* Ind2[i,j] ; //Tail-up linear-with-sill model
            }
            else{// if points are flow unconnected
              iji[1,1] = D[i,j];
              iji [2,1] = D[j,i];
              Ind2[i,j] = (max(iji) / alpha_td) <= 1 ? 1 : 0; // indicator
              C_td[i,j] = var_td * (1 - (max(iji) / alpha_td)) * Ind2[i,j] ;
            }
          }
        }

  '
  
  #tail-down spherical model
  tparam_td2_sph <- '
    var_td= sigma_td ^ 2; // variance tail-down
        for (i in 1:N) {// tail-down spherical model
          for (j in 1:N) {
            if(flow_con_mat[i,j] == 1){ // if points are flow connected
              Ind2[i,j] = (h[i,j] / alpha_td) <= 1 ? 1 : 0; // indicator
              C_td[i,j] = var_td * (1 - (1.5 * h[i,j] / alpha_td) + ( (h[i,j] ^ 3) / (2 * alpha_td ^ 3))) * Ind2[i,j];
            }
            else{// if points are flow unconnected
              iji[1,1] = D[i,j];
              iji [2,1] = D[j,i];
              Ind2[i,j] = (max(iji) / alpha_td) <= 1 ? 1 : 0; // indicator
              C_td[i,j] = var_td * (1 - (1.5 * min(iji) / alpha_td) + ( max(iji)/(2 * alpha_td)   )) * (1 - (max(iji) / alpha_td) ) ^ 2 * Ind2[i,j];
            }
          }
        }

'
  
  #tail-down models end
  
  
  tparam_ed2 <- '
	  //Euclidean distance models start
    var_ed = sigma_ed ^ 2; // var Euclidean dist
      C_ed = var_ed * exp(- 3 * e / alpha_ed); // exponential model
    //Euclidean distance models end
    //print("C_ed: ", C_ed)
'
  
  tparam_re2 <- '
    // random effect
    var_RE1 = sigma_RE1 ^ 2;
      C_re = var_RE1 * RE1;

'
  
  model_com <- '
  model {
    for (t in 1:T){
      target += multi_normal_cholesky_lpdf(Y[t] | mu[t], cholesky_decompose(C_tu + C_td + C_re + C_ed + var_nug * I + 1e-6) );
    }

    sigma_nug ~ uniform(0,50); // cauchy(0,1) prior nugget effect
    //phi ~ uniform(-1, 1); //
    phi ~ normal(0.5,0.3); //NB informative   .666
'
  
  model_tu <- '
    sigma_tu ~ uniform(0,100);  // cauchy(0,2) prior sd  tail-up model
    alpha_tu ~ uniform(0, alpha_max);
'
  
  model_td <- '
    sigma_td ~ uniform(0,100); // sd tail-down
    alpha_td ~ uniform(0, alpha_max);
'
  
  model_ed <- '
    sigma_ed ~ uniform(0,100); // sd Euclidean dist
    alpha_ed ~ uniform(0, alpha_max); // Euclidean dist range
'
  
  model_re <- '
    sigma_RE1 ~ uniform(0,5);
'
  
  
  gen_quant <- '
  generated quantities {
   //vector[T] log_lik;
   vector[N] y_pred[T];
     for (t in 1:T){
     y_pred[t] = multi_normal_cholesky_rng( mu[t], cholesky_decompose(C_tu + C_td + C_re + C_ed + var_nug * I + 1e-6) );
        }
  }
  '
  
  ssn_ar <- paste(
    data_com,
    
    param_com,
    
    if(cor_tu %in% 1:3) param_tu,
    if(cor_td %in% 1:3) param_td,
    if(cor_ed %in% 1:3) param_ed,
    if(cor_re %in% 1:3) param_re,
    
    if(time_method[[1]] == 'ar') param_phi_ar,
    if(time_method[[1]] == 'var') param_phi_var,
    
    '}',
    
    tparam_com,
    if(cor_tu %in% 1:3)tparam_tu,
    if(cor_td %in% 1:3)tparam_td,
    if(cor_ed %in% 1:3)tparam_ed,
    if(cor_re %in% 1:3) tparam_re,
    tparam_com2,
    
    
    if(time_method[[1]] == 'ar') tparam_com_ar,
    if(time_method[[1]] == 'var') tparam_com_var,
    
    
    #ifelse(cor_tu %in% 1:3, tparam_tu2, 'C_tu = rep_matrix(0, N, N);'),
    case_when(cor_tu == 1 ~ tparam_tu2_exp,
              cor_tu == 2 ~ tparam_tu2_lin,
              cor_tu == 3 ~ tparam_tu2_sph,
              cor_tu >= 4 | cor_tu <= 0 ~ 'C_tu = rep_matrix(0, N, N);'),
    
    
    case_when(cor_td == 1 ~ tparam_td2_exp,
              cor_td == 2 ~ tparam_td2_lin,
              cor_td == 3 ~ tparam_td2_sph,
              cor_td >= 4 | cor_td <= 0 ~ 'C_td = rep_matrix(0, N, N);'),
    
    case_when(cor_ed == 1 ~ tparam_ed2,
              cor_ed >= 2 | cor_ed <= 0 ~ 'C_ed = rep_matrix(0, N, N);'),
    
    case_when(cor_re == 1 ~ tparam_re2,
              cor_re >= 2 | cor_re <= 0 ~ 'C_re = rep_matrix(0, N, N);'),
    
    '}',
    model_com,
    if(cor_tu %in% 1:3)model_tu,
    if(cor_td %in% 1:3)model_td,
    if(cor_ed %in% 1:3)model_ed,
    if(cor_re %in% 1:3) model_re,
    '}',
    
    if(loglik == T) gen_quant
  )
  
  `%notin%` <- Negate(`%in%`)
  
  pars <- c(
    case_when(cor_tu %in% 1:3 ~ c('var_tu', 'alpha_tu'),
              cor_tu %notin% 1:3 ~ ""),
    
    case_when(cor_td %in% 1:3 ~ c('var_td', 'alpha_td'),
              cor_td %notin% 1:3 ~ ""),
    
    case_when(cor_ed %in% 1:3 ~ c('var_ed', 'alpha_ed'),
              cor_ed %notin% 1:3 ~ ""),
    
    case_when(cor_re %in% 1:3 ~ c('var_re', 'alpha_re'),
              cor_re %notin% 1:3 ~ ""),
    
    if(loglik == T) 'y_pred',
    
    'var_nug',
    'beta',
    'phi',
    'y'
  )
  pars <- pars[pars != '']
  
  
  # data part
  options(na.action='na.pass') # to preserve the NAs
  out_list <- mylm(formula = formula, data = data) # produces the design matrix
  
  response <- out_list$y # response variable
  design_matrix <- out_list$X # design matrix
  
  obs_data <- data
  #ndays <- length(unique(obs_data$date))
  
  ndays <- length(unique(obs_data[, names(obs_data) %in% time_points] ))
  
  
  N <- nrow(obs_data)/ndays #nobs
  
  
  nobs <- nrow(obs_data)/ndays #nobs
  
  #obs_data$date_num <- as.numeric(factor(obs_data$date))
  
  obs_data$date_num <- as.numeric(factor(obs_data[, names(obs_data) %in% time_points]	))
  
  
  resp_var_name <- gsub("[^[:alnum:]]", " ", formula[2])
  obs_data$y <- obs_data[,names(obs_data) %in% resp_var_name]
  
  #train <- nrow(obs_data[obs_data$date_num==1 & !is.na(obs_data$y),])
  #test <- nrow(obs_data[obs_data$date_num==1 & is.na(obs_data$y),])
  
  
  
  # array structure
  X <- design_matrix #cbind(1,obs_data[, c("X1", "X2", "X3")]) # design matrix
  #Xarray <- abind(matrix(X[,1],nrow = ndays, ncol = nobs, byrow = T),
  #                matrix(X[,2],nrow = ndays, ncol = nobs, byrow = T),
  #                matrix(X[,3],nrow = ndays, ncol = nobs, byrow = T),
  #                matrix(X[,4],nrow = ndays, ncol = nobs, byrow = T), along = 3)
  
  # NB: this array order is Stan specific
  Xarray <- aperm(array( c(X), dim=c(N, ndays, ncol(X)) ),c(2, 1, 3))
  
  
  y_obs <- response[!is.na(response)]#obs_data[!is.na(obs_data$y),]$y
  #Yarray <- array(matrix(ys, nrow = ndays, ncol = train, byrow = T),
  #                dim = c(ndays, train))
  
  # index for observed values
  i_y_obs <- obs_data[!is.na(obs_data$y),]$pid
  # i_y_obs <- matrix(i_y_obs, nrow = train , ncol = ndays, byrow = F)
  
  
  # index for missing values
  i_y_mis <- obs_data[is.na(obs_data$y),]$pid
  #  i_y_mis <- matrix(i_y_mis, nrow = test, ncol = ndays, byrow = F)
  
  
  
  if(ssn_object == T){ # the ssn object exist?
    mat_all <- dist_wei_mat(path = path, net = net, addfunccol = addfunccol)
  }
  
  if(ssn_object == F){ # the ssn object does not exist- purely spatial
    
    first_date <- unique(obs_data[, names(obs_data) %in% time_points])[1]
    
    di <- dist(obs_data[obs_data$date == first_date, c('lon', 'lat')], #data$date == 1
               method = "euclidean",
               diag = FALSE,
               upper = FALSE) %>% as.matrix()
    mat_all <-  list(e = di, D = di, H = di, w.matrix = di, flow.con.mat = di)
  }
  
  
  data_list <- list(N = N, # obs + preds  points
                    T = ndays, # time points
                    K = ncol(X),  # ncol of design matrix
                    y_obs = y_obs,# y values in the obs df
                    
                    N_y_obs = length(i_y_obs),  #nrow(i_y_obs) numb obs points
                    N_y_mis = length(i_y_mis), #nrow(i_y_mis) numb preds points
                    
                    i_y_obs = i_y_obs, # index of obs points
                    i_y_mis = i_y_mis, # index of preds points
                    
                    X = Xarray, # design matrix
                    mat_all = mat_all,
                    alpha_max = 4 * max(mat_all$H) ) # a list with all the distance/weights matrices
  
  
  data_list$e = data_list$mat_all$e #Euclidean dist
  #for tail-up
  data_list$h = data_list$mat_all$H # total stream distance
  data_list$W = data_list$mat_all$w.matrix # spatial weights
  
  #for tail-down
  data_list$flow_con_mat = data_list$mat_all$flow.con.mat #flow connected matrix
  data_list$D = data_list$mat_all$D #downstream hydro distance matrix
  
  #RE1 = RE1mm # random effect matrix
  
  data_list$I = diag(1, nrow(data_list$W), nrow(data_list$W))  # diagonal matrix
  
  # phi_ini <- ifelse(time_method == "ar", 0.5, rep(0.5,N))
  
  
  ini <- function(){list(var_nug =  .1#,
                         #y = rep( mean(obs_data$Y, na.rm = T),T*N)
                         #, phi = phi_ini #phi= rep(0.5,N)
                         #y = rep( mean(obs_data$temp, na.rm = T),T*N)
  )}
  
  # saveRDS(data_list, 'data_list_20210805.RDS')
  fit <- rstan::stan(model_code = ssn_ar,
                     model_name = "ssn_ar",
                     data = data_list,
                     pars = pars,
                     iter = iter,
                     warmup = warmup,
                     init = ini,
                     chains = chains,
                     verbose = F,
                     seed = seed,
                     refresh = refresh
  )
  attributes(fit)$formula <- formula
  
  fit
}



