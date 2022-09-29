library(dplyr)
library(utils)
library(readr)
library(invgamma)
library(goftest)
##### library(matlib)
library(geosphere)
library(ape)
library(zoo)
library(matrixcalc)
library(MASS)
library(mnormt)
library(matrixStats)
library(truncnorm)
library(mcmc)
#library(tmvtnorm)
library(sandwich) 
#library(gmm)
library(Matrix)
library(evd)
library(truncdist)
library(stats4)
library(parallel)
##### library(caret)
library(Rcpp)
library(RcppArmadillo)
##### library(stableGR)
library(coda)
library(diversitree)
library(SamplerCompare)
library("readxl")
library(tidyverse)

### Determine which of these are not installed on the Stats machines
### Run library() on server

### Look at Siddharth's data - is it the same? Any new variables?

# df = read.csv("D:\\Year_3\\Term2\\Research\\changepoint_analysis\\changepoint_gitrepo\\Code\\US_COVID_weekly_vaccine_NY_data_updated.csv", row.names = NULL)
#df$year_week_num = df$week

df <- read.csv("currentjoineddata.csv")

df$week = as.Date(df$week)

df = df %>% 
  arrange( week, country, state, key, lat, long, county,population) %>% 
  filter(Recip_County == "Albany County")
  

# df = df %>% mutate(first_dose_prevalence = total_dose_one/population,
#                    second_dose_prevalence = total_complete/population)

# View(df %>% group_by(county) %>% summarise(total_tp = n()))
# loc_not_all_time_points = c("New York", "St. Lawrence")
#test_loc = sample(unique(df$key),2)
start_training_wk = max(df$week) - 0
test_wk_min = max(df$week) - 0
validation_wk_min = max(df$week) - 0
training_df = df %>% filter(week <= validation_wk_min )#& !(county %in% loc_not_all_time_points)

#Scaling the covariates

training_df$scaled_prev_log_new_death = as.numeric(scale(training_df$prev_log_new_death))
training_df$scaled_first_dose_prevalence = as.numeric(scale(training_df$first_dose_prevalence))
#training_df$scaled_second_dose_prevalence = as.numeric(scale(training_df$second_dose_prevalence))

# if(is_all_loc){

#   training_df = df %>% filter(week <= validation_wk_min & week >= start_training_wk)
#   validation_df = df %>% filter(week > validation_wk_min & week <= test_wk_min)

# }else{
#   training_df = df %>% filter(! state %in% test_loc & week <= validation_wk_min)
#   validation_df = df %>% filter(! state %in% test_loc & week > validation_wk_min & week <= test_wk_min)
# }


unique_county_lat_long = unique(training_df %>% dplyr::select(key,lat,long))
distance_mat_training = as.matrix(distm(cbind(unique_county_lat_long$long,unique_county_lat_long$lat), 
                                        fun = distHaversine))/(1000)

#training values
#y_training = training_df$category
n_total_training = nrow(training_df)
n_sp_training = nrow(unique_county_lat_long)
n_tmp_training = n_total_training/n_sp_training

week_diff_mat_training = as.matrix(dist(1:n_tmp_training, diag =  TRUE,upper = TRUE))
#No of categories
#no_of_categ = length(unique(training_df$category))

#Number of covariates
p = 2

accuracy_pred_vec = vector()
pred_category_mat = vector()

# burn_period_training = 5000
sample_size_training = 1000
diff_in_random_draws_training = 10
#constant parameters values#############################################################################
training_df_clone = training_df
no_of_categ = length(unique(training_df_clone$category))
a_training = 2
lambda_training = 1
# stdnorm = 100
# alpha_vec = rep(1, no_of_categ)


gap_bw_changepoints = 5
changepoint_vec = c(0, seq(11,n_tmp_training - 11,by=1), n_tmp_training)# - gap_bw_changepoints

changepoint_vec = changepoint_vec[c(-39,-40)] #As 49 is creating problem for X'X inverse

changepoint_week_vec = c(min(training_df$week) - 7,unique(training_df$week))

#changepoint_week_vec = changepoint_week_vec[-40]
#training_df_clone = training_df
#Giving values to different parameters of our model
training_df_clone$v_vec1 = rep(0, n_total_training)#training_df_clone$v_vec_art
training_df_clone$v_vec2 = rep(0, n_total_training)
training_df_clone$v_vec3 = rep(0, n_total_training)

training_df_clone$u_vec1 = rep(0, n_total_training)#training_df_clone$v_vec_art
training_df_clone$u_vec2 = rep(0, n_total_training)
training_df_clone$u_vec3 = rep(0, n_total_training)

# training_df_clone$w_vec1 = rep(0, n_total_training)#training_df_clone$v_vec_art
# training_df_clone$w_vec2 = rep(0, n_total_training)
# training_df_clone$w_vec3 = rep(0, n_total_training)


#time taken analysis for parameters
total_tt_v = 0
total_tt_w = 0
total_tt_u = 0
total_tt_beta = 0
total_tt_sigma_u_sq = 0
#total_tt_sigma_w_sq = 0
total_tt_sigma_eps_sq = 0

total_tt_phi_u_sp = 0
total_tt_phi_u_tmp = 0
# total_tt_phi_w_sp = 0
# total_tt_phi_w_tmp = 0
total_tt_phi_v_sp = 0
total_tt_phi_v_tmp = 0

total_tt_t_0 = 0
total_tt_pi = 0
total_tt_delta = 0

total_tt_mat_inv = 0

#C++ code##########################################################################################
arma_code <-
  "arma::vec arma_mm(const arma::mat& m, const arma::vec& v, const arma::rowvec& tv) {
    return tv * m * v;
};"

arma_mm = cppFunction(code = arma_code, depends = "RcppArmadillo")
#Setting values for markov chains before convergence########################################################
p = 2

set.seed(8)
beta_mat_preconvergence_training = matrix(c(seq(1,p+1,1),
                                            seq(1,p+1,1) + rnorm(p+1, sd =0.4),
                                            seq(1,p+1,1) + rnorm(p+1, sd =0.4)), 
                                          ncol = 3, nrow = p+1)
set.seed(9)
beta_star_mat_preconvergence_training = matrix(c(seq(1,p+1,1),
                                                 seq(1,p+1,1)+ rnorm(p+1, sd =0.6),
                                                 seq(1,p+1,1) + rnorm(p+1, sd =0.6)), 
                                               ncol = 3, nrow = p+1)

#sigma_w_sq_preconvergence_training = rep(1,3)#c(1, 1.5, 2)

sigma_v_sq_preconvergence_training = rep(1,3)#c(1, 1.5, 2)

sigma_u_sq_preconvergence_training = rep(1,3)

sigma_eps_sq_preconvergence_training = rep(1,3)#c(1, 1.5, 2)

sigma_eps_sq_star_preconvergence_training = rep(1,3)#c(1, 1.5, 2)

phi_u_sp_preconvergence_training = c(0.01, 0.03, 0.05)

phi_u_tmp_preconvergence_training = c(0.1, 0.2, 0.3)

# phi_w_sp_preconvergence_training = c(0.01, 0.03, 0.05)
# 
# phi_w_tmp_preconvergence_training = c(0.1, 0.2, 0.3)

phi_v_sp_preconvergence_training = c(0.01, 0.03, 0.05)

phi_v_tmp_preconvergence_training = c(0.1, 0.2, 0.3)

phi_sp_vec_training = c(0.001, 0.002, 0.005, 0.01, 0.05, 0.1, 0.5)

phi_tmp_vec_training = c(0.05, 0.1, 0.25, 0.5, 1)


t0_preconvergence_training = c(61,60,62)#c(changepoint_vec[floor(length(changepoint_vec)/2)-6], #rep(t0_art, 3)
                             #changepoint_vec[floor(length(changepoint_vec)/2) -5],
                             #changepoint_vec[floor(length(changepoint_vec)/2) - 4])

changepoint_t0_week_preconvergence_training = changepoint_week_vec[t0_preconvergence_training + 1] #+1 because added zero

set.seed(5)
delta_temp = c(0, sort(runif(no_of_categ - 2, 1, 40)))#taking values from min and max pi vec art
delta_mat_preconvergence_training = cbind(delta_temp, sort(delta_temp + c(0, rnorm(2))),
                                          sort(delta_temp + c(0, rnorm(2))))

delta_mat_preconvergence_training =  rbind(c(-Inf, -Inf, -Inf), delta_mat_preconvergence_training,
                                           c(Inf, Inf, Inf)) #cbind(delta_art, delta_art, delta_art)



############################################################################################################
#Samples collected after convergence of markov chain
pi_vec_sample_training = vector()
delta_sample_training = vector()
beta_sample_training = vector()
beta_star_sample_training = vector()
sigma_u_sq_sample_training = vector()
#sigma_w_sq_sample_training = vector()
sigma_v_sq_sample_training = vector()
sigma_eps_sq_sample_training = vector()
sigma_eps_sq_star_sample_training = vector()

phi_u_sp_sample_training = vector()
phi_u_tmp_sample_training = vector()
phi_v_sp_sample_training = vector()
phi_v_tmp_sample_training = vector()
# phi_w_sp_sample_training = vector()
# phi_w_tmp_sample_training = vector()

v_vec_sample_training = vector()
u_vec_sample_training = vector()
#w_vec_sample_training = vector()
t_0_sample_training = vector()
#Initialtzing for different markov chain samples#########################################
beta_sample_chains_list_training = list()
beta_star_sample_chains_list_training = list()
sigma_u_sq_sample_chains_list_training = list()
sigma_w_sq_sample_chains_list_training = list()
sigma_v_sq_sample_chains_list_training = list()
sigma_eps_sq_sample_chains_list_training = list()
sigma_eps_sq_star_sample_chains_list_training = list()
phi_u_sp_sample_chains_list_training = list()
phi_u_tmp_sample_chains_list_training = list()
phi_v_sp_sample_chains_list_training = list()
phi_v_tmp_sample_chains_list_training = list()
phi_w_sp_sample_chains_list_training = list()
phi_w_tmp_sample_chains_list_training = list()


v_vec_sample_chains_list_training = list()
w_vec_sample_chains_list_training = list()
u_vec_sample_chains_list_training = list()

pi_vec_sample_chains_list_training = list()
t0_sample_chains_list_training = list()
delta_sample_chains_list_training = list()
mcmc_obj_list_training = list()

for (init in 1:3) {
  
  beta_sample_chains_list_training[[init]] = vector()
  beta_star_sample_chains_list_training[[init]] = vector()
  sigma_u_sq_sample_chains_list_training[[init]] = vector()
  sigma_w_sq_sample_chains_list_training[[init]] = vector()
  sigma_v_sq_sample_chains_list_training[[init]] = vector()
  sigma_eps_sq_sample_chains_list_training[[init]] = vector()
  sigma_eps_sq_star_sample_chains_list_training[[init]] = vector()
  
  phi_u_sp_sample_chains_list_training[[init]] = vector()
  phi_u_tmp_sample_chains_list_training[[init]] = vector()
  phi_v_sp_sample_chains_list_training[[init]] = vector()
  phi_v_tmp_sample_chains_list_training[[init]] = vector()
  phi_w_sp_sample_chains_list_training[[init]] = vector()
  phi_w_tmp_sample_chains_list_training[[init]] = vector()
  
  v_vec_sample_chains_list_training[[init]] = vector()
  #w_vec_sample_chains_list_training[[init]] = vector()
  u_vec_sample_chains_list_training[[init]] = vector()
  
  pi_vec_sample_chains_list_training[[init]] = vector()
  t0_sample_chains_list_training[[init]] = vector()
  delta_sample_chains_list_training[[init]] = vector()
  mcmc_obj_list_training[[init]] = vector()
}

#X matrix############################################################################
X_mat_training = cbind(rep(1,n_total_training), training_df$scaled_prev_log_new_death, training_df$scaled_first_dose_prevalence)
p = ncol(X_mat_training) - 1

ident_n_sp_mat_training = diag(n_sp_training)
#Function definition####################################################################################
phi_u_sp_log_density_func = function(phi_u_sp){
  
  SIGMA_u_sp = exp(- phi_u_sp * distance_mat_training)
  log_det_SIGMA_u_sp = (determinant(SIGMA_u_sp, logarithm = TRUE))$modulus[1]
  inv_SIGMA_u_sp = chol2inv(chol(SIGMA_u_sp))
  inv_SIGMA_u_tmp_kro_inv_SIGMA_u_sp = inv_SIGMA_u_tmp_training %x% inv_SIGMA_u_sp
  
  (- (as.numeric(arma_mm(inv_SIGMA_u_tmp_kro_inv_SIGMA_u_sp,
                         u_vec_training, u_vec_training)))/(2*sigma_u_sq_training) - 
      n_tmp_training*log_det_SIGMA_u_sp/2)
  
}

phi_u_tmp_log_density_func = function(phi_u_tmp){
  
  SIGMA_u_tmp = exp(- phi_u_tmp * week_diff_mat_training)
  log_det_SIGMA_u_tmp = (determinant(SIGMA_u_tmp, logarithm = TRUE))$modulus[1]
  inv_SIGMA_u_tmp = chol2inv(chol(SIGMA_u_tmp))
  inv_SIGMA_u_tmp_kro_inv_SIGMA_u_sp = inv_SIGMA_u_tmp %x% inv_SIGMA_u_sp_training
  
  (- (as.numeric(arma_mm(inv_SIGMA_u_tmp_kro_inv_SIGMA_u_sp,
                         u_vec_training,u_vec_training)))/(2*sigma_u_sq_training) - 
      n_sp_training*log_det_SIGMA_u_tmp/2)
  
}

# phi_w_sp_log_density_func = function(phi_w_sp){
#   
#   SIGMA_w_sp = exp(- phi_w_sp * distance_mat_training)
#   log_det_SIGMA_w_sp = (determinant(SIGMA_w_sp, logarithm = TRUE))$modulus[1]
#   inv_SIGMA_w_sp = chol2inv(chol(SIGMA_w_sp))
#   inv_SIGMA_w_tmp_kro_inv_SIGMA_w_sp = inv_SIGMA_w_tmp_training %x% inv_SIGMA_w_sp
#   
#   (- (as.numeric(arma_mm(inv_SIGMA_w_tmp_kro_inv_SIGMA_w_sp,
#                          w_vec_training,w_vec_training)))/(2*sigma_w_sq_training) - 
#       n_tmp_training*log_det_SIGMA_w_sp/2)
#   
# }
# 
# phi_w_tmp_log_density_func = function(phi_w_tmp){
#   
#   SIGMA_w_tmp = exp(- phi_w_tmp * week_diff_mat_training)
#   log_det_SIGMA_w_tmp = (determinant(SIGMA_w_tmp, logarithm = TRUE))$modulus[1]
#   inv_SIGMA_w_tmp = chol2inv(chol(SIGMA_w_tmp))
#   inv_SIGMA_w_tmp_kro_inv_SIGMA_w_sp = inv_SIGMA_w_tmp %x% inv_SIGMA_w_sp_training
#   
#   (- (as.numeric(arma_mm(inv_SIGMA_w_tmp_kro_inv_SIGMA_w_sp,
#                          w_vec_training,w_vec_training)))/(2*sigma_w_sq_training) - 
#       n_sp_training*log_det_SIGMA_w_tmp/2)
#   
# }

phi_v_sp_log_density_func = function(phi_v_sp){
  
  SIGMA_v_sp = exp(- phi_v_sp * distance_mat_training)
  log_det_SIGMA_v_sp = (determinant(SIGMA_v_sp, logarithm = TRUE))$modulus[1]
  inv_SIGMA_v_sp = chol2inv(chol(SIGMA_v_sp))
  inv_SIGMA_v_tmp_kro_inv_SIGMA_v_sp = inv_SIGMA_v_tmp_training %x% inv_SIGMA_v_sp
  
  (- (as.numeric(arma_mm(inv_SIGMA_v_tmp_kro_inv_SIGMA_v_sp,
                         v_vec_training,v_vec_training)))/(2*sigma_v_sq_training) - 
      n_tmp_training*log_det_SIGMA_v_sp/2)
  
}

phi_v_tmp_log_density_func = function(phi_v_tmp){
  
  SIGMA_v_tmp = exp(- phi_v_tmp * week_diff_mat_training)
  log_det_SIGMA_v_tmp = (determinant(SIGMA_v_tmp, logarithm = TRUE))$modulus[1]
  inv_SIGMA_v_tmp = chol2inv(chol(SIGMA_v_tmp))
  inv_SIGMA_v_tmp_kro_inv_SIGMA_v_sp = inv_SIGMA_v_tmp %x% inv_SIGMA_v_sp_training
  
  (- (as.numeric(arma_mm(inv_SIGMA_v_tmp_kro_inv_SIGMA_v_sp,
                         v_vec_training,v_vec_training)))/(2*sigma_v_sq_training) - 
      n_sp_training*log_det_SIGMA_v_tmp/2)
  
}

delta3_dens_func = function(delta3){
  
  delta4 = delta_training[4]
  y_vec = training_df_clone$category
  density_log_y = 0
  
  
  u_val_minus_vec = training_df_clone[[u_vec_mc_col]][1:(t_0_training* n_sp_training)]
  u_val_plus_vec = training_df_clone[[u_vec_mc_col]][(t_0_training* n_sp_training + 1):
                                                       (n_total_training)]
  
  #w_val_minus_vec = training_df_clone[[w_vec_mc_col]][1:(tp* n_sp_training)]
  v_val_plus_vec = training_df_clone[[v_vec_mc_col]][(t_0_training* n_sp_training + 1):
                                                       (n_total_training)]
  
  x_beta_minus_vec = training_df_clone$x_beta_vec[1:(t_0_training* n_sp_training)]
  x_beta_plus_vec = training_df_clone$x_beta_vec[(t_0_training* n_sp_training + 1):
                                                   (n_total_training)]
  
  y_minus_vec = y_vec[1:(t_0_training * n_sp_training)]
  y_plus_vec = y_vec[(t_0_training* n_sp_training + 1):
                       (n_total_training)]
  
  delta2_minus_vec = rep(delta_training[2], t_0_training* n_sp_training)
  delta2_plus_vec = rep(delta_training[2], (n_tmp_training - t_0_training)* n_sp_training)
  
  delta3_minus_vec = rep(delta3, t_0_training* n_sp_training)
  delta3_plus_vec = rep(delta3, (n_tmp_training - t_0_training)* n_sp_training)
  
  delta4_minus_vec = rep(delta4, t_0_training* n_sp_training)
  delta4_plus_vec = rep(delta4, (n_tmp_training - t_0_training)* n_sp_training)
  
  idx_categ1_minus = which(y_minus_vec == 1)
  idx_categ1_plus = which(y_plus_vec == 1)
  
  idx_categ2_minus = which(y_minus_vec == 2)
  idx_categ2_plus = which(y_plus_vec == 2)
  
  idx_categ3_minus = which(y_minus_vec == 3)
  idx_categ3_plus = which(y_plus_vec == 3)
  
  idx_categ4_minus = which(y_minus_vec == 4)
  idx_categ4_plus = which(y_plus_vec == 4)
  
  delta2_dash_minus_vec = ((delta2_minus_vec - x_beta_minus_vec - u_val_minus_vec)/# - w_val_minus_vec
                             sqrt(sigma_eps_sq_training))
  
  delta2_dash_plus_vec = ((delta2_plus_vec - x_beta_plus_vec - u_val_plus_vec - v_val_plus_vec)/
                            sqrt(sigma_eps_sq_star_training))
  
  
  delta3_dash_minus_vec = ((delta3_minus_vec - x_beta_minus_vec - u_val_minus_vec)/# - w_val_minus_vec
                             sqrt(sigma_eps_sq_training))
  
  delta3_dash_plus_vec = ((delta3_plus_vec - x_beta_plus_vec - u_val_plus_vec - v_val_plus_vec)/
                            sqrt(sigma_eps_sq_star_training))
  
  delta4_dash_minus_vec = ((delta4_minus_vec - x_beta_minus_vec - u_val_minus_vec)/# - w_val_minus_vec
                             sqrt(sigma_eps_sq_training))
  
  delta4_dash_plus_vec = ((delta4_plus_vec - x_beta_plus_vec - u_val_plus_vec - v_val_plus_vec)/
                            sqrt(sigma_eps_sq_star_training))
  
  
  #Categ2 minus density##############################################################################################
  categ2_minus_log_prob_vec = log(pnorm(delta3_dash_minus_vec[idx_categ2_minus], lower.tail = TRUE, log.p = FALSE) - 
                                    pnorm(delta2_dash_minus_vec[idx_categ2_minus], lower.tail = TRUE, log.p = FALSE))
  
  if(length(which(categ2_minus_log_prob_vec == -Inf)) > 0){
    categ2_minus_log_prob_vec[which(categ2_minus_log_prob_vec == -Inf)] = log(0.000000001)
  }
  
  density_log_y = density_log_y + sum(categ2_minus_log_prob_vec)
  
  #Categ2 plus density###############################################################################################
  categ2_plus_log_prob_vec = log(pnorm(delta3_dash_plus_vec[idx_categ2_plus], lower.tail = TRUE, log.p = FALSE) - 
                                   pnorm(delta2_dash_plus_vec[idx_categ2_plus], lower.tail = TRUE, log.p = FALSE))
  
  if(length(which(categ2_plus_log_prob_vec == -Inf)) > 0){
    categ2_plus_log_prob_vec[which(categ2_plus_log_prob_vec == -Inf)] = log(0.000000001)
  }
  
  density_log_y = density_log_y + sum(categ2_plus_log_prob_vec)
  
  #Categ3 minus density##############################################################################################
  categ3_minus_log_prob_vec = log(pnorm(delta4_dash_minus_vec[idx_categ3_minus], lower.tail = TRUE, log.p = FALSE) - 
                                    pnorm(delta3_dash_minus_vec[idx_categ3_minus], lower.tail = TRUE, log.p = FALSE) )
  
  if(length(which(categ3_minus_log_prob_vec == -Inf)) > 0){
    categ3_minus_log_prob_vec[which(categ3_minus_log_prob_vec == -Inf)] = log(0.000000001)
  }
  
  
  density_log_y = density_log_y + sum(categ3_minus_log_prob_vec)
  
  #Categ3 plus density###############################################################################################
  categ3_plus_log_prob_vec = log(pnorm(delta4_dash_plus_vec[idx_categ3_plus], lower.tail = TRUE, log.p = FALSE) - 
                                   pnorm(delta3_dash_plus_vec[idx_categ3_plus], lower.tail = TRUE, log.p = FALSE) )
  
  if(length(which(categ3_plus_log_prob_vec == -Inf)) > 0){
    categ3_plus_log_prob_vec[which(categ3_plus_log_prob_vec == -Inf)] = log(0.000000001)
  }
  
  density_log_y = density_log_y + sum(categ3_plus_log_prob_vec)
  
  density_log_y 
  
}


delta4_dens_func = function(delta4){
  
  delta3 = delta_training[3]
  y_vec = training_df_clone$category
  density_log_y = 0
  
  u_val_minus_vec = training_df_clone[[u_vec_mc_col]][1:(t_0_training* n_sp_training)]
  u_val_plus_vec = training_df_clone[[u_vec_mc_col]][(t_0_training* n_sp_training + 1):
                                                       (n_total_training)]
  
  #w_val_minus_vec = training_df_clone[[w_vec_mc_col]][1:(tp* n_sp_training)]
  v_val_plus_vec = training_df_clone[[v_vec_mc_col]][(t_0_training* n_sp_training + 1):
                                                       (n_total_training)]
  
  x_beta_minus_vec = training_df_clone$x_beta_vec[1:(t_0_training* n_sp_training)]
  x_beta_plus_vec = training_df_clone$x_beta_vec[(t_0_training* n_sp_training + 1):
                                                   (n_total_training)]
  
  y_minus_vec = y_vec[1:(t_0_training * n_sp_training)]
  y_plus_vec = y_vec[(t_0_training* n_sp_training + 1):
                       (n_total_training)]
  
  # delta2_minus_vec = rep(delta_training[2], t_0_training* n_sp_training)
  # delta2_plus_vec = rep(delta_training[2], (n_tmp_training - t_0_training)* n_sp_training)
  
  delta3_minus_vec = rep(delta3, t_0_training* n_sp_training)
  delta3_plus_vec = rep(delta3, (n_tmp_training - t_0_training)* n_sp_training)
  
  delta4_minus_vec = rep(delta4, t_0_training* n_sp_training)
  delta4_plus_vec = rep(delta4, (n_tmp_training - t_0_training)* n_sp_training)
  
  idx_categ1_minus = which(y_minus_vec == 1)
  idx_categ1_plus = which(y_plus_vec == 1)
  
  idx_categ2_minus = which(y_minus_vec == 2)
  idx_categ2_plus = which(y_plus_vec == 2)
  
  idx_categ3_minus = which(y_minus_vec == 3)
  idx_categ3_plus = which(y_plus_vec == 3)
  
  idx_categ4_minus = which(y_minus_vec == 4)
  idx_categ4_plus = which(y_plus_vec == 4)
  
  # delta2_dash_minus_vec = ((delta2_minus_vec - x_beta_minus_vec - u_val_minus_vec)/# - w_val_minus_vec
  #                            sqrt(sigma_eps_sq_training))
  # 
  # delta2_dash_plus_vec = ((delta2_plus_vec - x_beta_plus_vec - u_val_plus_vec - v_val_plus_vec)/
  #                           sqrt(sigma_eps_sq_star_training))
  
  
  delta3_dash_minus_vec = ((delta3_minus_vec - x_beta_minus_vec - u_val_minus_vec)/# - w_val_minus_vec
                             sqrt(sigma_eps_sq_training))
  
  delta3_dash_plus_vec = ((delta3_plus_vec - x_beta_plus_vec - u_val_plus_vec - v_val_plus_vec)/
                            sqrt(sigma_eps_sq_star_training))
  
  delta4_dash_minus_vec = ((delta4_minus_vec - x_beta_minus_vec - u_val_minus_vec)/# - w_val_minus_vec
                             sqrt(sigma_eps_sq_training))
  
  delta4_dash_plus_vec = ((delta4_plus_vec - x_beta_plus_vec - u_val_plus_vec - v_val_plus_vec)/
                            sqrt(sigma_eps_sq_star_training))
  
  
  #Categ3 minus density##############################################################################################
  categ3_minus_log_prob_vec = log(pnorm(delta4_dash_minus_vec[idx_categ3_minus], lower.tail = TRUE, log.p = FALSE) - 
                                    pnorm(delta3_dash_minus_vec[idx_categ3_minus], lower.tail = TRUE, log.p = FALSE) )
  
  if(length(which(categ3_minus_log_prob_vec == -Inf)) > 0){
    categ3_minus_log_prob_vec[which(categ3_minus_log_prob_vec == -Inf)] = log(0.000000001)
  }
  
  
  density_log_y = density_log_y + sum(categ3_minus_log_prob_vec)
  
  #Categ3 plus density###############################################################################################
  categ3_plus_log_prob_vec = log(pnorm(delta4_dash_plus_vec[idx_categ3_plus], lower.tail = TRUE, log.p = FALSE) - 
                                   pnorm(delta3_dash_plus_vec[idx_categ3_plus], lower.tail = TRUE, log.p = FALSE) )
  
  if(length(which(categ3_plus_log_prob_vec == -Inf)) > 0){
    categ3_plus_log_prob_vec[which(categ3_plus_log_prob_vec == -Inf)] = log(0.000000001)
  }
  
  density_log_y = density_log_y + sum(categ3_plus_log_prob_vec)
  
  #Categ4 minus density################################################################################################
  categ4_minus_log_prob_vec = pnorm(delta4_dash_minus_vec[idx_categ4_minus], lower.tail = FALSE, log.p = TRUE)
  
  if(length(which(categ4_minus_log_prob_vec == -Inf)) > 0){
    categ4_minus_log_prob_vec[which(categ4_minus_log_prob_vec == -Inf)] = log(0.000000001)
  }
  
  
  density_log_y = density_log_y + sum(categ4_minus_log_prob_vec)
  
  #Categ4 plus density#################################################################################################
  categ4_plus_log_prob_vec = pnorm(delta4_dash_plus_vec[idx_categ4_plus], lower.tail = FALSE, log.p = TRUE)
  
  if(length(which(categ4_plus_log_prob_vec == -Inf)) > 0){
    categ4_plus_log_prob_vec[which(categ4_plus_log_prob_vec == -Inf)] = log(0.000000001)
  }
  
  density_log_y = density_log_y + sum(categ4_plus_log_prob_vec)
  
  density_log_y 
  
}

#Iteration Starts here
converged = FALSE
converged_i = -1
i = 1

while(i > 0){
  if(!converged){
    
    for (mc in (1:3)) {
      
      v_vec_mc_col = paste("v_vec", as.character(mc), sep = "")
      #w_vec_mc_col = paste("w_vec", as.character(mc), sep = "")
      u_vec_mc_col = paste("u_vec", as.character(mc), sep = "")
      pi_vec_mc_col = paste("pi_vec", as.character(mc), sep = "")
      #training Df pre and post change point
      # training_pre_changepoint_df = training_df_clone %>% filter(week <= changepoint_t_0_week)
      # training_post_changepoint_df = training_df_clone %>% filter(week > changepoint_t_0_week)
      beta_training = beta_mat_preconvergence_training[,mc]
      beta_star_training = beta_star_mat_preconvergence_training[,mc]
      delta_training = delta_mat_preconvergence_training[,mc]
      
      sigma_u_sq_training = sigma_u_sq_preconvergence_training[mc]
      #sigma_w_sq_training = sigma_w_sq_preconvergence_training[mc]
      sigma_v_sq_training = sigma_v_sq_preconvergence_training[mc]
      sigma_eps_sq_training = sigma_eps_sq_preconvergence_training[mc]
      sigma_eps_sq_star_training = sigma_eps_sq_star_preconvergence_training[mc]
      
      phi_u_sp_training = phi_u_sp_preconvergence_training[mc]
      phi_u_tmp_training = phi_u_tmp_preconvergence_training[mc]
      # phi_w_sp_training = phi_w_sp_preconvergence_training[mc]
      # phi_w_tmp_training = phi_w_tmp_preconvergence_training[mc]
      phi_v_sp_training = phi_v_sp_preconvergence_training[mc]
      phi_v_tmp_training = phi_v_tmp_preconvergence_training[mc]
      
      t_0_training = t0_preconvergence_training[mc]
      changepoint_t_0_week = changepoint_week_vec[t_0_training + 1] #+1 because added zero
      n_tmp_pre_changepoint_training = t_0_training
      n_tmp_post_changepoint_training = n_tmp_training - n_tmp_pre_changepoint_training
      #Need to run for the first iteration#########################################################################
      if(i == 1){
        
        st_mat_inv = Sys.time()
        #Getting inv Sigma_sp_22 to use for every iteration
        
        SIGMA_v_sp_training = exp(- phi_v_sp_training * distance_mat_training)
        inv_SIGMA_v_sp_training = chol2inv(chol(SIGMA_v_sp_training))
        #log_det_SIGMA_v_sp_training = (determinant(SIGMA_v_sp_training, logarithm = TRUE))$modulus[1]
        
        SIGMA_v_tmp_training = exp(- phi_v_tmp_training * week_diff_mat_training)
        inv_SIGMA_v_tmp_training = chol2inv(chol(SIGMA_v_tmp_training))
        #log_det_SIGMA_v_tmp_training = (determinant(SIGMA_v_tmp_training, logarithm = TRUE))$modulus[1]
        
        #inv_SIGMA_v_tmp_kro_inv_SIGMA_v_sp_training = inv_SIGMA_v_tmp_training %x% inv_SIGMA_v_sp_training
        
        inv_SIGMA_v_tmp_trans_22_matlist_training = list()
        mat_mult_temp_SIGMA_v_tmp_vec = vector()
        # inv_SIGMA_cond_pre_changepoint_matlist_training = list()
        # inv_SIGMA_cond_post_changepoint_matlist_training = list()
        SIGMA_v_tmp_trans_21_mat_training = vector()
        
        for (tmp in (1:n_tmp_training)) {
          idx = 1:n_tmp_training
          
          if(tmp != 1){
            idx[1] = tmp
            idx[tmp] = 1
          }
          SIGMA_v_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% SIGMA_v_tmp_training %*% 
                                          diag(n_tmp_training)[idx,])
          
          inv_SIGMA_v_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% inv_SIGMA_v_tmp_training %*% 
                                              diag(n_tmp_training)[idx,])
          
          inv_SIGMA_v_tmp_trans_a_training = inv_SIGMA_v_tmp_trans_training[1,1]
          inv_SIGMA_v_tmp_trans_b_training = inv_SIGMA_v_tmp_trans_training[1,-1]
          inv_SIGMA_v_tmp_trans_c_training = inv_SIGMA_v_tmp_trans_training[-1,1]
          inv_SIGMA_v_tmp_trans_d_training = inv_SIGMA_v_tmp_trans_training[-1,-1]
          
          inv_SIGMA_v_tmp_trans_22_training = (inv_SIGMA_v_tmp_trans_d_training - (inv_SIGMA_v_tmp_trans_c_training %*%
                                                                                     t(inv_SIGMA_v_tmp_trans_b_training))/
                                                 inv_SIGMA_v_tmp_trans_a_training)
          
          SIGMA_v_tmp_trans_21_training = SIGMA_v_tmp_trans_training[-1,1]
          
          mat_mult_temp_SIGMA_v_tmp_vec[tmp] = as.numeric(1 - t(SIGMA_v_tmp_trans_21_training) %*%
                                                            inv_SIGMA_v_tmp_trans_22_training %*%
                                                            SIGMA_v_tmp_trans_21_training)
          
          inv_SIGMA_v_tmp_trans_22_matlist_training[[tmp]] = inv_SIGMA_v_tmp_trans_22_training
          SIGMA_v_tmp_trans_21_mat_training = cbind(SIGMA_v_tmp_trans_21_mat_training, SIGMA_v_tmp_trans_21_training)
          
        }
        #Getting the matrix division for W calculations##################################################################
        # SIGMA_w_sp_training = exp(- phi_w_sp_training * distance_mat_training)
        # inv_SIGMA_w_sp_training = chol2inv(chol(SIGMA_w_sp_training))
        # #log_det_SIGMA_w_sp_training = (determinant(SIGMA_w_sp_training, logarithm = TRUE))$modulus[1]
        # 
        # 
        # SIGMA_w_tmp_training = exp(- phi_w_tmp_training * week_diff_mat_training)
        # inv_SIGMA_w_tmp_training = chol2inv(chol(SIGMA_w_tmp_training))
        # #log_det_SIGMA_w_tmp_training = (determinant(SIGMA_w_tmp_training, logarithm = TRUE))$modulus[1]
        # 
        # #inv_SIGMA_w_tmp_kro_inv_SIGMA_w_sp_training = inv_SIGMA_w_tmp_training %x% inv_SIGMA_w_sp_training
        # 
        # inv_SIGMA_w_tmp_trans_22_matlist_training = list()
        # mat_mult_temp_SIGMA_w_tmp_vec = vector()
        # 
        # SIGMA_w_tmp_trans_21_mat_training = vector()
        # 
        # for (tmp in (1:n_tmp_training)) {
        #   idx = 1:n_tmp_training
        #   
        #   if(tmp != 1){
        #     idx[1] = tmp
        #     idx[tmp] = 1
        #   }
        #   SIGMA_w_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% SIGMA_w_tmp_training %*% 
        #                                   diag(n_tmp_training)[idx,])
        #   
        #   inv_SIGMA_w_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% inv_SIGMA_w_tmp_training %*% 
        #                                       diag(n_tmp_training)[idx,])
        #   
        #   inv_SIGMA_w_tmp_trans_a_training = inv_SIGMA_w_tmp_trans_training[1,1]
        #   inv_SIGMA_w_tmp_trans_b_training = inv_SIGMA_w_tmp_trans_training[1,-1]
        #   inv_SIGMA_w_tmp_trans_c_training = inv_SIGMA_w_tmp_trans_training[-1,1]
        #   inv_SIGMA_w_tmp_trans_d_training = inv_SIGMA_w_tmp_trans_training[-1,-1]
        #   
        #   inv_SIGMA_w_tmp_trans_22_training = (inv_SIGMA_w_tmp_trans_d_training - (inv_SIGMA_w_tmp_trans_c_training %*%
        #                                                                              t(inv_SIGMA_w_tmp_trans_b_training))/
        #                                          inv_SIGMA_w_tmp_trans_a_training)
        #   
        #   SIGMA_w_tmp_trans_21_training = SIGMA_w_tmp_trans_training[-1,1]
        #   
        #   mat_mult_temp_SIGMA_w_tmp_vec[tmp] = as.numeric(1 - t(SIGMA_w_tmp_trans_21_training) %*%
        #                                                     inv_SIGMA_w_tmp_trans_22_training %*%
        #                                                     SIGMA_w_tmp_trans_21_training)
        #   
        #   
        #   inv_SIGMA_w_tmp_trans_22_matlist_training[[tmp]] = inv_SIGMA_w_tmp_trans_22_training
        #   SIGMA_w_tmp_trans_21_mat_training = cbind(SIGMA_w_tmp_trans_21_mat_training, SIGMA_w_tmp_trans_21_training)
        #   
        # }
        
        #Getting the matrix division for U calculations##################################################################
        
        SIGMA_u_sp_training = exp(- phi_u_sp_training * distance_mat_training)
        inv_SIGMA_u_sp_training = chol2inv(chol(SIGMA_u_sp_training))
        #log_det_SIGMA_u_sp_training = (determinant(SIGMA_u_sp_training, logarithm = TRUE))$modulus[1]
        
        SIGMA_u_tmp_training = exp(- phi_u_tmp_training * week_diff_mat_training)
        inv_SIGMA_u_tmp_training = chol2inv(chol(SIGMA_u_tmp_training))
        #log_det_SIGMA_u_tmp_training = (determinant(SIGMA_u_tmp_training, logarithm = TRUE))$modulus[1]
        
        #inv_SIGMA_u_tmp_kro_inv_SIGMA_u_sp_training = inv_SIGMA_u_tmp_training %x% inv_SIGMA_u_sp_training
        
        inv_SIGMA_u_tmp_trans_22_matlist_training = list()
        mat_mult_temp_SIGMA_u_tmp_vec = vector()
        # inv_SIGMA_cond_pre_changepoint_matlist_training = list()
        # inv_SIGMA_cond_post_changepoint_matlist_training = list()
        SIGMA_u_tmp_trans_21_mat_training = vector()
        
        for (tmp in (1:n_tmp_training)) {
          idx = 1:n_tmp_training
          
          if(tmp != 1){
            idx[1] = tmp
            idx[tmp] = 1
          }
          SIGMA_u_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% SIGMA_u_tmp_training %*% 
                                          diag(n_tmp_training)[idx,])
          
          inv_SIGMA_u_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% inv_SIGMA_u_tmp_training %*% 
                                              diag(n_tmp_training)[idx,])
          
          inv_SIGMA_u_tmp_trans_a_training = inv_SIGMA_u_tmp_trans_training[1,1]
          inv_SIGMA_u_tmp_trans_b_training = inv_SIGMA_u_tmp_trans_training[1,-1]
          inv_SIGMA_u_tmp_trans_c_training = inv_SIGMA_u_tmp_trans_training[-1,1]
          inv_SIGMA_u_tmp_trans_d_training = inv_SIGMA_u_tmp_trans_training[-1,-1]
          
          inv_SIGMA_u_tmp_trans_22_training = (inv_SIGMA_u_tmp_trans_d_training - (inv_SIGMA_u_tmp_trans_c_training %*%
                                                                                     t(inv_SIGMA_u_tmp_trans_b_training))/
                                                 inv_SIGMA_u_tmp_trans_a_training)
          
          SIGMA_u_tmp_trans_21_training = SIGMA_u_tmp_trans_training[-1,1]
          
          mat_mult_temp_SIGMA_u_tmp_vec[tmp] = as.numeric(1 - t(SIGMA_u_tmp_trans_21_training) %*%
                                                            inv_SIGMA_u_tmp_trans_22_training %*%
                                                            SIGMA_u_tmp_trans_21_training)
          
          
          inv_SIGMA_u_tmp_trans_22_matlist_training[[tmp]] = inv_SIGMA_u_tmp_trans_22_training
          SIGMA_u_tmp_trans_21_mat_training = cbind(SIGMA_u_tmp_trans_21_mat_training, SIGMA_u_tmp_trans_21_training)
          
        }
        
        total_tt_mat_inv = total_tt_mat_inv + as.numeric(difftime(Sys.time(), st_mat_inv), units="secs")
      }
      
      #pi posterior################################################################################################
      st_pi = Sys.time()
      pi_vec = vector()
      for (st in 1:n_total_training) {
        
        st_week = training_df_clone$week[st]
        y_val_st = training_df_clone$category[st]
        v_val_st = training_df_clone[[v_vec_mc_col]][st]
        u_val_st = training_df_clone[[u_vec_mc_col]][st]
        #w_val_st = training_df_clone[[w_vec_mc_col]][st]
        x_vec_st = c(1, training_df_clone$scaled_prev_log_new_death[st],
                     training_df_clone$scaled_first_dose_prevalence[st])
        
        if(st_week <= changepoint_t_0_week){
          
          pi_st = truncnorm::rtruncnorm(1, a = delta_training[y_val_st], b = delta_training[y_val_st + 1],
                                        mean = t(x_vec_st) %*% beta_training + u_val_st ,#+ w_val_st
                                        sd = sqrt(sigma_eps_sq_training))
          
        }else{
          
          pi_st = truncnorm::rtruncnorm(1, a = delta_training[y_val_st], b = delta_training[y_val_st + 1],
                                        mean = t(x_vec_st) %*% beta_star_training + u_val_st + v_val_st,
                                        sd = sqrt(sigma_eps_sq_star_training))
          
        }
        
        pi_vec = c(pi_vec, pi_st)
      }
      training_df_clone[[pi_vec_mc_col]] = NULL
      training_df_clone[[pi_vec_mc_col]] = pi_vec
      
      pi_vec_training = training_df_clone[[pi_vec_mc_col]]
      
      training_pre_changepoint_df = training_df_clone %>% filter(week <= changepoint_t_0_week)
      training_post_changepoint_df = training_df_clone %>% filter(week > changepoint_t_0_week)
      
      #Getting X matrix for pre and post change point
      X_mat_pre_changepoint_training = cbind(rep(1,n_sp_training * n_tmp_pre_changepoint_training),
                                             training_pre_changepoint_df$scaled_prev_log_new_death,
                                             training_pre_changepoint_df$scaled_first_dose_prevalence)
      
      X_mat_post_changepoint_training = cbind(rep(1,n_sp_training * n_tmp_post_changepoint_training),
                                              training_post_changepoint_df$scaled_prev_log_new_death,
                                              training_post_changepoint_df$scaled_first_dose_prevalence)
      
      p = ncol(X_mat_pre_changepoint_training) - 1
      
      total_tt_pi = total_tt_pi + as.numeric(difftime(Sys.time(), st_pi), units="secs")
      ###############################################################################################################
      #v posterior#############################################################
      st_v = Sys.time()
      #Xbeta matrix
      v_mat_training = matrix((training_df_clone[[v_vec_mc_col]]),
                              nrow = n_sp_training, ncol = n_tmp_training,
                              byrow = FALSE)
      
      u_mat_training = matrix((training_df_clone[[u_vec_mc_col]]),
                              nrow = n_sp_training, ncol = n_tmp_training,
                              byrow = FALSE)
      
      # w_mat_training = matrix((training_df_clone[[w_vec_mc_col]]),
      #                         nrow = n_sp_training, ncol = n_tmp_training,
      #                         byrow = FALSE)
      
      pi_mat_training = matrix((training_df_clone[[pi_vec_mc_col]]),
                               nrow = n_sp_training,
                               ncol = n_tmp_training, byrow = FALSE)
      
      xbeta_mat_training = matrix(c(as.numeric(X_mat_pre_changepoint_training %*% beta_training),
                                    as.numeric(X_mat_post_changepoint_training %*%beta_star_training)),
                                  nrow = n_sp_training,
                                  ncol = n_tmp_training, byrow = FALSE)
      
      for (tmp in (1:n_tmp_training)){
        idx = 1:n_tmp_training
        
        if(tmp!=1){
          idx[1] = tmp
          idx[tmp] = 1
        }
        #Switching first and jth index column for V matrix and then taking values from the 2nd column
        v_cond_training = (v_mat_training %*%
                             diag(n_tmp_training)[idx,])[(1 + n_sp_training):(n_total_training)]
        
        
        
        #Can be made faster by doing this for all change points or instead of doing it for all check if for
        #the change point this is already done if so get the value from there
        mu_cond_training = ((t(SIGMA_v_tmp_trans_21_mat_training[,tmp]) %*%
                               inv_SIGMA_v_tmp_trans_22_matlist_training[[tmp]]) %x%
                              ident_n_sp_mat_training) %*% v_cond_training
        
        #Can be made faster by dynamically saving values in a vector
        
        
        
        
        if(tmp > t_0_training){
          
          inv_SIGMA_cond_training = (inv_SIGMA_v_sp_training/mat_mult_temp_SIGMA_v_tmp_vec[tmp])
          
          v_covar_training = chol2inv(chol(inv_SIGMA_cond_training/sigma_v_sq_training +
                                             ident_n_sp_mat_training/sigma_eps_sq_star_training))
          
          
          v_mean_training  = v_covar_training %*% (inv_SIGMA_cond_training %*%
                                                     mu_cond_training/sigma_v_sq_training +
                                                     (pi_mat_training[,tmp] -
                                                        xbeta_mat_training[,tmp] - u_mat_training[,tmp])/
                                                     sigma_eps_sq_star_training)
          
          
          v_mat_training[,tmp] = v_mean_training + t(chol(v_covar_training))  %*% rnorm(n_sp_training)
          
        }else{
          
          SIGMA_cond_training = (SIGMA_v_sp_training * mat_mult_temp_SIGMA_v_tmp_vec[tmp])
          
          v_mat_training[,tmp] = (mu_cond_training + t(chol(SIGMA_cond_training))  %*% 
                                    rnorm(n_sp_training, sd = sqrt(sigma_v_sq_training)))
          
        }
        
        
      }
      training_df_clone[[v_vec_mc_col]] = v_mat_training[1:n_total_training]
      v_vec_training = training_df_clone[[v_vec_mc_col]]
      
      total_tt_v = total_tt_v + as.numeric(difftime(Sys.time(), st_v), units="secs")
      ########################################################################################################################
      #w posterior###########################################################################################
      # st_w = Sys.time()
      # 
      # for (tmp in (1:n_tmp_training)){
      #   idx = 1:n_tmp_training
      #   
      #   if(tmp!=1){
      #     idx[1] = tmp
      #     idx[tmp] = 1
      #   }
      #   #Switching first and jth index column for V matrix and then taking values from the 2nd column
      #   w_cond_training = (w_mat_training %*%
      #                        diag(n_tmp_training)[idx,])[(1 + n_sp_training):(n_total_training)]
      #   
      #   
      #   
      #   #Can be made faster by doing this for all change points or instead of doing it for all check if for
      #   #the change point this is already done if so get the value from there
      #   mu_cond_training = ((t(SIGMA_w_tmp_trans_21_mat_training[,tmp]) %*%
      #                          inv_SIGMA_w_tmp_trans_22_matlist_training[[tmp]]) %x%
      #                         ident_n_sp_mat_training) %*% w_cond_training
      #   
      #   
      #   if(tmp <= t_0_training){
      #     
      #     inv_SIGMA_cond_training = (inv_SIGMA_w_sp_training/mat_mult_temp_SIGMA_w_tmp_vec[tmp])
      #     
      #     w_covar_training = chol2inv(chol(inv_SIGMA_cond_training/sigma_w_sq_training +
      #                                        ident_n_sp_mat_training/sigma_eps_sq_training))
      #     
      #     
      #     w_mean_training  = w_covar_training %*% (inv_SIGMA_cond_training %*%
      #                                                mu_cond_training/sigma_w_sq_training +
      #                                                (pi_mat_training[,tmp] -
      #                                                   xbeta_mat_training[,tmp] - u_mat_training[,tmp])/
      #                                                sigma_eps_sq_training)
      #     
      #     
      #     w_mat_training[,tmp] = w_mean_training + t(chol(w_covar_training))  %*% rnorm(n_sp_training)
      #     
      #   }else{
      #     
      #     SIGMA_cond_training = (SIGMA_w_sp_training * mat_mult_temp_SIGMA_w_tmp_vec[tmp])
      #     
      #     w_mat_training[,tmp] = (mu_cond_training + t(chol(SIGMA_cond_training))  %*% 
      #                               rnorm(n_sp_training, sd = sqrt(sigma_w_sq_training)))
      #     
      #   }
      #   
      #   
      # }
      # training_df_clone[[w_vec_mc_col]] = w_mat_training[1:n_total_training]
      # w_vec_training = training_df_clone[[w_vec_mc_col]]
      # 
      # total_tt_w = total_tt_w + as.numeric(difftime(Sys.time(), st_w), units="secs")
      #u posterior##################################################################################################
      st_u = Sys.time()
      
      for (tmp in (1:n_tmp_training)){
        idx = 1:n_tmp_training
        
        if(tmp!=1){
          idx[1] = tmp
          idx[tmp] = 1
        }
        #Switching first and jth index column for V matrix and then taking values from the 2nd column
        u_cond_training = (u_mat_training %*%
                             (diag(n_tmp_training)[idx,]))[(1 + n_sp_training):(n_total_training)]
        
        
        
        #Can be made faster by doing this for all change points or instead of doing it for all check if for
        #the change point this is already done if so get the value from there
        mu_cond_training = ((t(SIGMA_u_tmp_trans_21_mat_training[,tmp]) %*%
                               inv_SIGMA_u_tmp_trans_22_matlist_training[[tmp]]) %x%
                              ident_n_sp_mat_training) %*% u_cond_training
        
        #Can be made faster by dynamically saving values in a vector
        inv_SIGMA_cond_training = (inv_SIGMA_u_sp_training/mat_mult_temp_SIGMA_u_tmp_vec[tmp])
        
        
        
        if(tmp <= t_0_training){
          
          
          
          u_covar_training = chol2inv(chol(inv_SIGMA_cond_training/sigma_u_sq_training +
                                             ident_n_sp_mat_training/sigma_eps_sq_training))
          
          
          u_mean_training  = u_covar_training %*% (inv_SIGMA_cond_training %*%
                                                     mu_cond_training/sigma_u_sq_training +
                                                     (pi_mat_training[,tmp] -
                                                        xbeta_mat_training[,tmp])/
                                                     sigma_eps_sq_training)
          
          
          u_mat_training[,tmp] = u_mean_training + t(chol(u_covar_training))  %*% rnorm(n_sp_training)
          
        }else{
          
          u_covar_training = chol2inv(chol(inv_SIGMA_cond_training/sigma_u_sq_training +
                                             ident_n_sp_mat_training/sigma_eps_sq_star_training))
          
          u_mean_training  = u_covar_training %*% (inv_SIGMA_cond_training %*%
                                                     mu_cond_training/sigma_u_sq_training +
                                                     (pi_mat_training[,tmp] -
                                                        xbeta_mat_training[,tmp] - v_mat_training[,tmp])/
                                                     sigma_eps_sq_star_training)
          
          u_mat_training[,tmp] = u_mean_training + t(chol(u_covar_training))  %*% rnorm(n_sp_training)
          
        }
        
        
      }
      training_df_clone[[u_vec_mc_col]] = u_mat_training[1:n_total_training]
      u_vec_training = training_df_clone[[u_vec_mc_col]]
      
      total_tt_u = total_tt_u + as.numeric(difftime(Sys.time(), st_u), units="secs")
      #beta posterior##########################################################################################
      st_beta = Sys.time()
      
      pi_pre_changepoint_training = pi_mat_training[1:(n_tmp_pre_changepoint_training*n_sp_training)]
      pi_post_changepoint_training = pi_mat_training[(n_tmp_pre_changepoint_training*n_sp_training+1):n_total_training]
      
      #w_minus_training = w_mat_training[1:(n_tmp_pre_changepoint_training*n_sp_training)]
      v_plus_training = v_mat_training[(n_tmp_pre_changepoint_training*n_sp_training+1):n_total_training]
      u_minus_training = u_mat_training[1:(n_tmp_pre_changepoint_training*n_sp_training)]
      u_plus_training = u_mat_training[(n_tmp_pre_changepoint_training*n_sp_training+1):n_total_training]
      #dividing by 10^4 for large variance for beta prior
      beta_covar_training = chol2inv(chol( (t(X_mat_pre_changepoint_training) %*% 
                                              X_mat_pre_changepoint_training)/sigma_eps_sq_training)) #+ diag(p+1)/10^4
      
      beta_star_covar_training = chol2inv(chol( (t(X_mat_post_changepoint_training) %*% 
                                                   X_mat_post_changepoint_training)/sigma_eps_sq_star_training))
      
      beta_mean_training = beta_covar_training %*% ( t(X_mat_pre_changepoint_training)%*% 
                                                       (pi_pre_changepoint_training - u_minus_training))/sigma_eps_sq_training
      
      beta_star_mean_training = (beta_star_covar_training %*% (t(X_mat_post_changepoint_training)%*% 
                                                                 (pi_post_changepoint_training - u_plus_training 
                                                                  - v_plus_training))/sigma_eps_sq_star_training)
      
      beta_training = beta_mean_training + t(chol(beta_covar_training)) %*% rnorm(p+1)
      
      beta_star_training = beta_star_mean_training + t(chol(beta_star_covar_training)) %*% rnorm(p+1)
      
      beta_mat_preconvergence_training[,mc] = beta_training
      beta_star_mat_preconvergence_training[,mc] = beta_star_training
      
      training_df_clone$x_beta_vec = c(as.numeric(X_mat_pre_changepoint_training %*% beta_training),
                                       as.numeric(X_mat_post_changepoint_training %*%beta_star_training))
      
      
      total_tt_beta = total_tt_beta + as.numeric(difftime(Sys.time(), st_beta), units="secs")
      #######################################################################################################################
      #Sigma1^2 posterior######################################################
      #st_sigma_w_sq = Sys.time()
      
      # sigma_w_sq_a_training = a_training + (n_sp_training * n_tmp_pre_changepoint_training)/2
      # 
      # sigma_v_sq_a_training = a_training + (n_sp_training * n_tmp_post_changepoint_training)/2
      # 
      # sigma_w_sq_lambda_training = (t(v_training) %*%
      #                                inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[t_0_training/gap_bw_changepoints - 1]]
      #                              %*% v_training)/2 + lambda_training
      # 
      # sigma_v_sq_lambda_training = (t(v_star_training) %*%
      #                                     inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_mat_list[[t_0_training/gap_bw_changepoints - 1]]
      #                                   %*% v_star_training)/2 + lambda_training
      
      #sigma_w_sq_training = 1#invgamma::rinvgamma(1,sigma_w_sq_a_training,rate = sigma_w_sq_lambda_training)
      
      sigma_v_sq_training = 1#invgamma::rinvgamma(1, sigma_v_sq_a_training,
      #rate = sigma_v_sq_lambda_training)
      #sigma_w_sq_preconvergence_training[mc] = sigma_w_sq_training
      sigma_v_sq_preconvergence_training[mc] = sigma_v_sq_training
      
      #total_tt_sigma_w_sq = total_tt_sigma_w_sq + as.numeric(difftime(Sys.time(), st_sigma_w_sq), units="secs")
      ##########################################################################################################################
      #Sigma2^2 posterior############################################
      st_sigma_eps_sq = Sys.time()
      
      # sigma_eps_sq_a_training = a_training + (n_sp_training * n_tmp_pre_changepoint_training)/2
      # 
      # sigma_eps_sq_star_a_training = a_training + (n_sp_training * n_tmp_post_changepoint_training)/2
      # 
      # temp_vec = pi_pre_changepoint_training - X_mat_pre_changepoint_training%*%beta_training - v_training
      # 
      # temp_vec2 = (pi_post_changepoint_training - X_mat_post_changepoint_training %*% beta_star_training -
      #                v_star_training)
      # 
      # sigma_eps_sq_lambda_training = sum(temp_vec^2)/2 + lambda_training
      # 
      # sigma_eps_sq_star_lambda_training = sum(temp_vec2^2)/2 + lambda_training
      
      sigma_eps_sq_training = 1#invgamma::rinvgamma(1,sigma_eps_sq_a_training,rate = sigma_eps_sq_lambda_training)
      
      sigma_eps_sq_star_training = 1#invgamma::rinvgamma(1, sigma_eps_sq_star_a_training,
      #rate = sigma_eps_sq_star_lambda_training)
      
      sigma_eps_sq_preconvergence_training[mc] = sigma_eps_sq_training
      sigma_eps_sq_star_preconvergence_training[mc] = sigma_eps_sq_star_training
      
      total_tt_sigma_eps_sq = total_tt_sigma_eps_sq + as.numeric(difftime(Sys.time(), st_sigma_eps_sq), units="secs")
      #phi_u_sp posterior arms########################################################################################
      if(i == 1 | i %% 20 == 0){
        
        st_phi_u_sp = Sys.time()
        
        phi_u_sp_log_density_vec_training = vector()
        for(phi_u_sp_val in  phi_sp_vec_training){
          
          phi_u_sp_log_density_vec_training = c(phi_u_sp_log_density_vec_training, 
                                                phi_u_sp_log_density_func(phi_u_sp_val))
          
          
        }
        
        phi_u_sp_log_density_vec_training = (phi_u_sp_log_density_vec_training - max(phi_u_sp_log_density_vec_training))
        
        phi_u_sp_density_vec_training = exp(phi_u_sp_log_density_vec_training)
        
        const_prop_phi_u_sp = 1/sum(phi_u_sp_density_vec_training)
        
        phi_u_sp_training = sample(x = (phi_sp_vec_training), 1, replace = FALSE,
                                   prob = round(phi_u_sp_density_vec_training * const_prop_phi_u_sp, digits = 3))
        
        phi_u_sp_preconvergence_training[mc] = phi_u_sp_training
        total_tt_phi_u_sp = total_tt_phi_u_sp + as.numeric(difftime(Sys.time(), st_phi_u_sp), units="secs")
      }
      #phi_u_tmp posterior arms########################################################################################
      if(i == 1 | i %% 20 == 0){
        
        
        st_phi_u_tmp = Sys.time()
        
        
        phi_u_tmp_log_density_vec_training = vector()
        for(phi_u_tmp_val in  phi_tmp_vec_training){
          
          phi_u_tmp_log_density_vec_training = c(phi_u_tmp_log_density_vec_training, 
                                                 phi_u_tmp_log_density_func(phi_u_tmp_val))
          
          
        }
        
        phi_u_tmp_log_density_vec_training = (phi_u_tmp_log_density_vec_training - max(phi_u_tmp_log_density_vec_training))
        
        phi_u_tmp_density_vec_training = exp(phi_u_tmp_log_density_vec_training)
        
        const_prop_phi_u_tmp = 1/sum(phi_u_tmp_density_vec_training)
        
        phi_u_tmp_training = sample(x = (phi_tmp_vec_training), 1, replace = FALSE,
                                    prob = round(phi_u_tmp_density_vec_training * const_prop_phi_u_tmp, digits = 3))
        
        phi_u_tmp_preconvergence_training[mc] = phi_u_tmp_training
        
        total_tt_phi_u_tmp = total_tt_phi_u_tmp + as.numeric(difftime(Sys.time(), st_phi_u_tmp), units="secs")
      }
      #phi_w_sp posterior arms########################################################################################
      # if(i %% 20 == 0){
      #   
      #   st_phi_w_sp = Sys.time()
      #   
      #   phi_w_sp_log_density_vec_training = vector()
      #   for(phi_w_sp_val in  phi_sp_vec_training){
      #     
      #     phi_w_sp_log_density_vec_training = c(phi_w_sp_log_density_vec_training, 
      #                                           phi_w_sp_log_density_func(phi_w_sp_val))
      #     
      #     
      #   }
      #   
      #   phi_w_sp_log_density_vec_training = (phi_w_sp_log_density_vec_training - max(phi_w_sp_log_density_vec_training))
      #   
      #   phi_w_sp_density_vec_training = exp(phi_w_sp_log_density_vec_training)
      #   
      #   const_prop_phi_w_sp = 1/sum(phi_w_sp_density_vec_training)
      #   
      #   phi_w_sp_training = sample(x = (phi_sp_vec_training), 1, replace = FALSE,
      #                              prob = round(phi_w_sp_density_vec_training * const_prop_phi_w_sp, digits = 3))
      #   
      #   phi_w_sp_preconvergence_training[mc] = phi_w_sp_training
      #   total_tt_phi_w_sp = total_tt_phi_w_sp + as.numeric(difftime(Sys.time(), st_phi_w_sp), units="secs")
      # }
      # #phi_w_tmp posterior arms########################################################################################
      # if(i %% 20 == 0){
      #   
      #   st_phi_w_tmp = Sys.time()
      #   
      #   phi_w_tmp_log_density_vec_training = vector()
      #   for(phi_w_tmp_val in  phi_tmp_vec_training){
      #     
      #     phi_w_tmp_log_density_vec_training = c(phi_w_tmp_log_density_vec_training, 
      #                                            phi_w_tmp_log_density_func(phi_w_tmp_val))
      #     
      #     
      #   }
      #   
      #   phi_w_tmp_log_density_vec_training = (phi_w_tmp_log_density_vec_training - max(phi_w_tmp_log_density_vec_training))
      #   
      #   phi_w_tmp_density_vec_training = exp(phi_w_tmp_log_density_vec_training)
      #   
      #   const_prop_phi_w_tmp = 1/sum(phi_w_tmp_density_vec_training)
      #   
      #   phi_w_tmp_training = sample(x = (phi_tmp_vec_training), 1, replace = FALSE,
      #                               prob = round(phi_w_tmp_density_vec_training * const_prop_phi_w_tmp, digits = 3))
      #   
      #   phi_w_tmp_preconvergence_training[mc] = phi_w_tmp_training
      #   total_tt_phi_w_tmp = total_tt_phi_w_tmp + as.numeric(difftime(Sys.time(), st_phi_w_tmp), units="secs")
      # }
      #phi_v_sp posterior arms########################################################################################
      if(i == 1 | i %% 20 == 0){
        
        st_phi_v_sp = Sys.time()
        
        phi_v_sp_log_density_vec_training = vector()
        for(phi_v_sp_val in  phi_sp_vec_training){
          
          phi_v_sp_log_density_vec_training = c(phi_v_sp_log_density_vec_training, 
                                                phi_v_sp_log_density_func(phi_v_sp_val))
          
          
        }
        
        phi_v_sp_log_density_vec_training = (phi_v_sp_log_density_vec_training - max(phi_v_sp_log_density_vec_training))
        
        phi_v_sp_density_vec_training = exp(phi_v_sp_log_density_vec_training)
        
        const_prop_phi_v_sp = 1/sum(phi_v_sp_density_vec_training)
        
        phi_v_sp_training = sample(x = (phi_sp_vec_training), 1, replace = FALSE,
                                   prob = round(phi_v_sp_density_vec_training * const_prop_phi_v_sp, digits = 3))
        
        phi_v_sp_preconvergence_training[mc] = phi_v_sp_training
        total_tt_phi_v_sp = total_tt_phi_v_sp + as.numeric(difftime(Sys.time(), st_phi_v_sp), units="secs")
      }
      #phi_v_tmp posterior arms########################################################################################
      if(i == 1 | i %% 20 == 0){
        
        st_phi_v_tmp = Sys.time()
        
        phi_v_tmp_log_density_vec_training = vector()
        for(phi_v_tmp_val in  phi_tmp_vec_training){
          
          phi_v_tmp_log_density_vec_training = c(phi_v_tmp_log_density_vec_training, 
                                                 phi_v_tmp_log_density_func(phi_v_tmp_val))
          
          
        }
        
        phi_v_tmp_log_density_vec_training = (phi_v_tmp_log_density_vec_training - max(phi_v_tmp_log_density_vec_training))
        
        phi_v_tmp_density_vec_training = exp(phi_v_tmp_log_density_vec_training)
        
        const_prop_phi_v_tmp = 1/sum(phi_v_tmp_density_vec_training)
        
        phi_v_tmp_training = sample(x = (phi_tmp_vec_training), 1, replace = FALSE,
                                    prob = round(phi_v_tmp_density_vec_training * const_prop_phi_v_tmp, digits = 3))
        
        phi_v_tmp_preconvergence_training[mc] = phi_v_tmp_training
        
        total_tt_phi_v_tmp = total_tt_phi_v_tmp + as.numeric(difftime(Sys.time(), st_phi_v_tmp), units="secs")
      }
      
      #Getting the matrix division for V calculations##################################################################
      #Getting inv Sigma_sp_22 to use for every iteration
      st_mat_inv = Sys.time()
      
      SIGMA_v_sp_training = exp(- phi_v_sp_training * distance_mat_training)
      inv_SIGMA_v_sp_training = chol2inv(chol(SIGMA_v_sp_training))
      #log_det_SIGMA_v_sp_training = (determinant(SIGMA_v_sp_training, logarithm = TRUE))$modulus[1]
      
      SIGMA_v_tmp_training = exp(- phi_v_tmp_training * week_diff_mat_training)
      inv_SIGMA_v_tmp_training = chol2inv(chol(SIGMA_v_tmp_training))
      #log_det_SIGMA_v_tmp_training = (determinant(SIGMA_v_tmp_training, logarithm = TRUE))$modulus[1]
      
      #inv_SIGMA_v_tmp_kro_inv_SIGMA_v_sp_training = inv_SIGMA_v_tmp_training %x% inv_SIGMA_v_sp_training
      
      inv_SIGMA_v_tmp_trans_22_matlist_training = list()
      mat_mult_temp_SIGMA_v_tmp_vec = vector()
      # inv_SIGMA_cond_pre_changepoint_matlist_training = list()
      # inv_SIGMA_cond_post_changepoint_matlist_training = list()
      SIGMA_v_tmp_trans_21_mat_training = vector()
      
      for (tmp in (1:n_tmp_training)) {
        idx = 1:n_tmp_training
        
        if(tmp != 1){
          idx[1] = tmp
          idx[tmp] = 1
        }
        SIGMA_v_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% SIGMA_v_tmp_training %*% 
                                        diag(n_tmp_training)[idx,])
        
        inv_SIGMA_v_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% inv_SIGMA_v_tmp_training %*% 
                                            diag(n_tmp_training)[idx,])
        
        inv_SIGMA_v_tmp_trans_a_training = inv_SIGMA_v_tmp_trans_training[1,1]
        inv_SIGMA_v_tmp_trans_b_training = inv_SIGMA_v_tmp_trans_training[1,-1]
        inv_SIGMA_v_tmp_trans_c_training = inv_SIGMA_v_tmp_trans_training[-1,1]
        inv_SIGMA_v_tmp_trans_d_training = inv_SIGMA_v_tmp_trans_training[-1,-1]
        
        inv_SIGMA_v_tmp_trans_22_training = (inv_SIGMA_v_tmp_trans_d_training - (inv_SIGMA_v_tmp_trans_c_training %*%
                                                                                   t(inv_SIGMA_v_tmp_trans_b_training))/
                                               inv_SIGMA_v_tmp_trans_a_training)
        
        SIGMA_v_tmp_trans_21_training = SIGMA_v_tmp_trans_training[-1,1]
        
        mat_mult_temp_SIGMA_v_tmp_vec[tmp] = as.numeric(1 - t(SIGMA_v_tmp_trans_21_training) %*%
                                                          inv_SIGMA_v_tmp_trans_22_training %*%
                                                          SIGMA_v_tmp_trans_21_training)
        
        inv_SIGMA_v_tmp_trans_22_matlist_training[[tmp]] = inv_SIGMA_v_tmp_trans_22_training
        SIGMA_v_tmp_trans_21_mat_training = cbind(SIGMA_v_tmp_trans_21_mat_training, SIGMA_v_tmp_trans_21_training)
        
      }
      #Getting the matrix division for W calculations##################################################################
      # SIGMA_w_sp_training = exp(- phi_w_sp_training * distance_mat_training)
      # inv_SIGMA_w_sp_training = chol2inv(chol(SIGMA_w_sp_training))
      # #log_det_SIGMA_w_sp_training = (determinant(SIGMA_w_sp_training, logarithm = TRUE))$modulus[1]
      # 
      # 
      # SIGMA_w_tmp_training = exp(- phi_w_tmp_training * week_diff_mat_training)
      # inv_SIGMA_w_tmp_training = chol2inv(chol(SIGMA_w_tmp_training))
      # #log_det_SIGMA_w_tmp_training = (determinant(SIGMA_w_tmp_training, logarithm = TRUE))$modulus[1]
      # 
      # #inv_SIGMA_w_tmp_kro_inv_SIGMA_w_sp_training = inv_SIGMA_w_tmp_training %x% inv_SIGMA_w_sp_training
      # 
      # inv_SIGMA_w_tmp_trans_22_matlist_training = list()
      # mat_mult_temp_SIGMA_w_tmp_vec = vector()
      # 
      # SIGMA_w_tmp_trans_21_mat_training = vector()
      # 
      # for (tmp in (1:n_tmp_training)) {
      #   idx = 1:n_tmp_training
      #   
      #   if(tmp != 1){
      #     idx[1] = tmp
      #     idx[tmp] = 1
      #   }
      #   SIGMA_w_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% SIGMA_w_tmp_training %*% 
      #                                   diag(n_tmp_training)[idx,])
      #   
      #   inv_SIGMA_w_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% inv_SIGMA_w_tmp_training %*% 
      #                                       diag(n_tmp_training)[idx,])
      #   
      #   inv_SIGMA_w_tmp_trans_a_training = inv_SIGMA_w_tmp_trans_training[1,1]
      #   inv_SIGMA_w_tmp_trans_b_training = inv_SIGMA_w_tmp_trans_training[1,-1]
      #   inv_SIGMA_w_tmp_trans_c_training = inv_SIGMA_w_tmp_trans_training[-1,1]
      #   inv_SIGMA_w_tmp_trans_d_training = inv_SIGMA_w_tmp_trans_training[-1,-1]
      #   
      #   inv_SIGMA_w_tmp_trans_22_training = (inv_SIGMA_w_tmp_trans_d_training - (inv_SIGMA_w_tmp_trans_c_training %*%
      #                                                                              t(inv_SIGMA_w_tmp_trans_b_training))/
      #                                          inv_SIGMA_w_tmp_trans_a_training)
      #   
      #   SIGMA_w_tmp_trans_21_training = SIGMA_w_tmp_trans_training[-1,1]
      #   
      #   mat_mult_temp_SIGMA_w_tmp_vec[tmp] = as.numeric(1 - t(SIGMA_w_tmp_trans_21_training) %*%
      #                                                     inv_SIGMA_w_tmp_trans_22_training %*%
      #                                                     SIGMA_w_tmp_trans_21_training)
      #   
      #   inv_SIGMA_w_tmp_trans_22_matlist_training[[tmp]] = inv_SIGMA_w_tmp_trans_22_training
      #   SIGMA_w_tmp_trans_21_mat_training = cbind(SIGMA_w_tmp_trans_21_mat_training, SIGMA_w_tmp_trans_21_training)
      #   
      # }
      
      #Getting the matrix division for U calculations##################################################################
      
      SIGMA_u_sp_training = exp(- phi_u_sp_training * distance_mat_training)
      inv_SIGMA_u_sp_training = chol2inv(chol(SIGMA_u_sp_training))
      #log_det_SIGMA_u_sp_training = (determinant(SIGMA_u_sp_training, logarithm = TRUE))$modulus[1]
      
      SIGMA_u_tmp_training = exp(- phi_u_tmp_training * week_diff_mat_training)
      inv_SIGMA_u_tmp_training = chol2inv(chol(SIGMA_u_tmp_training))
      #log_det_SIGMA_u_tmp_training = (determinant(SIGMA_u_tmp_training, logarithm = TRUE))$modulus[1]
      
      #inv_SIGMA_u_tmp_kro_inv_SIGMA_u_sp_training = inv_SIGMA_u_tmp_training %x% inv_SIGMA_u_sp_training
      
      inv_SIGMA_u_tmp_trans_22_matlist_training = list()
      mat_mult_temp_SIGMA_u_tmp_vec = vector()
      # inv_SIGMA_cond_pre_changepoint_matlist_training = list()
      # inv_SIGMA_cond_post_changepoint_matlist_training = list()
      SIGMA_u_tmp_trans_21_mat_training = vector()
      
      for (tmp in (1:n_tmp_training)) {
        idx = 1:n_tmp_training
        
        if(tmp != 1){
          idx[1] = tmp
          idx[tmp] = 1
        }
        SIGMA_u_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% SIGMA_u_tmp_training %*% 
                                        diag(n_tmp_training)[idx,])
        
        inv_SIGMA_u_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% inv_SIGMA_u_tmp_training %*% 
                                            diag(n_tmp_training)[idx,])
        
        inv_SIGMA_u_tmp_trans_a_training = inv_SIGMA_u_tmp_trans_training[1,1]
        inv_SIGMA_u_tmp_trans_b_training = inv_SIGMA_u_tmp_trans_training[1,-1]
        inv_SIGMA_u_tmp_trans_c_training = inv_SIGMA_u_tmp_trans_training[-1,1]
        inv_SIGMA_u_tmp_trans_d_training = inv_SIGMA_u_tmp_trans_training[-1,-1]
        
        inv_SIGMA_u_tmp_trans_22_training = (inv_SIGMA_u_tmp_trans_d_training - (inv_SIGMA_u_tmp_trans_c_training %*%
                                                                                   t(inv_SIGMA_u_tmp_trans_b_training))/
                                               inv_SIGMA_u_tmp_trans_a_training)
        
        SIGMA_u_tmp_trans_21_training = SIGMA_u_tmp_trans_training[-1,1]
        
        mat_mult_temp_SIGMA_u_tmp_vec[tmp] = as.numeric(1 - t(SIGMA_u_tmp_trans_21_training) %*%
                                                          inv_SIGMA_u_tmp_trans_22_training %*%
                                                          SIGMA_u_tmp_trans_21_training)
        
        inv_SIGMA_u_tmp_trans_22_matlist_training[[tmp]] = inv_SIGMA_u_tmp_trans_22_training
        SIGMA_u_tmp_trans_21_mat_training = cbind(SIGMA_u_tmp_trans_21_mat_training, SIGMA_u_tmp_trans_21_training)
        
      }
      
      total_tt_mat_inv = total_tt_mat_inv + as.numeric(difftime(Sys.time(), st_mat_inv), units="secs")
      
      ##############################################################################################################
      #delta posterior###############################
      st_delta = Sys.time()
      
      
      delta3_slice = diversitree::mcmc(lik = delta3_dens_func, x.init = c(delta_training[3]), 
                                       nsteps = 1, w= 5, lower = 0, upper = delta_training[4])
      
      delta_training[3] =delta3_slice$pars
      
      delta4_slice = diversitree::mcmc(lik = delta4_dens_func, x.init = c(delta_training[4]), 
                                       nsteps = 1, w= 5, lower = delta_training[3], upper = Inf)
      
      delta_training[4] = delta4_slice$pars
      
      #alpha_training = c(log(delta_training[3]), log(delta_training[4] - delta_training[3]))
      delta_mat_preconvergence_training[,mc] = delta_training
      
      total_tt_delta = total_tt_delta + as.numeric(difftime(Sys.time(), st_delta), units="secs")
      ###########################################################################################
      #t_0 posterior####################################
      st_t_0 = Sys.time()
      log_prop_pmf_t_0_vec_training = vector()
      for (tp in (changepoint_vec)) {
        
        week_t = changepoint_week_vec[tp + 1]#+1 because added zero
        
        
        training_post_t_df = training_df_clone %>% filter(week > week_t)
        
        
        if(tp != 0){
          training_pre_t_df = training_df_clone %>% filter(week <= week_t)
          X_mat_pre_t = cbind(rep(1,n_sp_training * tp),
                              training_pre_t_df$scaled_prev_log_new_death,
                              training_pre_t_df$scaled_first_dose_prevalence)
          
          u_vec_pre_t = training_pre_t_df[[u_vec_mc_col]]
          pi_vec_pre_t = training_pre_t_df[[pi_vec_mc_col]]
          
          temp_vec = pi_vec_pre_t - X_mat_pre_t %*% beta_training - u_vec_pre_t
        }
        
        
        if(tp != n_tmp_training){
          
          X_mat_post_t = cbind(rep(1, n_sp_training * (n_tmp_training - tp)),
                               training_post_t_df$scaled_prev_log_new_death,
                               training_post_t_df$scaled_first_dose_prevalence)
          
          #w_vec_pre_t = training_pre_t_df$w_vec
          v_vec_post_t = training_post_t_df[[v_vec_mc_col]]
          
          u_vec_post_t = training_post_t_df[[u_vec_mc_col]]
          
          pi_vec_post_t = training_post_t_df[[pi_vec_mc_col]]
          
          
          #- w_vec_pre_t
          temp_vec2 = pi_vec_post_t - X_mat_post_t %*% beta_star_training - u_vec_post_t - v_vec_post_t
          n_tmp_post_t = n_tmp_training - tp
        }
        
        if(tp == 0){
          
          log_prop_pmf_t_0_training = (((-n_sp_training*(n_tmp_post_t)/2) *
                                          log(2*pi*sigma_eps_sq_star_training)) +
                                         ((-1/(2*sigma_eps_sq_star_training)) * sum(temp_vec2^2)))
          
          
        }else if(tp == n_tmp_training){
          
          log_prop_pmf_t_0_training = ((-n_sp_training * tp/2)*log(2*pi*sigma_eps_sq_training) +
                                         ((-1/(2*sigma_eps_sq_training)) * sum(temp_vec^2)) )
          
        }else{
          
          log_prop_pmf_t_0_training = ((-n_sp_training * tp/2)*log(2*pi*sigma_eps_sq_training) +
                                         ((-1/(2*sigma_eps_sq_training)) * sum(temp_vec^2)) +
                                         
                                         ((-n_sp_training*(n_tmp_post_t)/2) *
                                            log(2*pi*sigma_eps_sq_star_training)) +
                                         ((-1/(2*sigma_eps_sq_star_training)) * sum(temp_vec2^2)))
          
        }
        
        
        log_prop_pmf_t_0_vec_training = c(log_prop_pmf_t_0_vec_training, log_prop_pmf_t_0_training)
      }
      #min_prop_pmf_t_0_vec_training = min(prop_pmf_t_0_vec_training)
      max_log_prop_pmf_t_0_vec_training = max(log_prop_pmf_t_0_vec_training)
      # if(min_prop_pmf_t_0_vec_training < 0){
      #   prop_pmf_t_0_vec_training = prop_pmf_t_0_vec_training + (-min_prop_pmf_t_0_vec_training)
      # }
      log_prop_pmf_t_0_vec_training_minus_max = log_prop_pmf_t_0_vec_training - max_log_prop_pmf_t_0_vec_training
      prop_pmf_t_0_vec_training = exp(log_prop_pmf_t_0_vec_training_minus_max)
      
      print(prop_pmf_t_0_vec_training)
      
      const_prop = 1/sum(prop_pmf_t_0_vec_training)
      
      t_0_training = sample(x = (changepoint_vec),1, replace = FALSE,
                            prob = round(prop_pmf_t_0_vec_training * const_prop, digits = 3))
      
      n_tmp_pre_changepoint_training = t_0_training
      n_tmp_post_changepoint_training = n_tmp_training - n_tmp_pre_changepoint_training
      changepoint_t_0_week = changepoint_week_vec[n_tmp_pre_changepoint_training + 1]#+1 because added zero
      
      t0_preconvergence_training[mc] = t_0_training
      total_tt_t_0 = total_tt_t_0 + as.numeric(difftime(Sys.time(), st_t_0), units="secs")
      
      #Collecting samples for Gelman Rubin stat#####################################################################
      
      v_vec_sample_chains_list_training[[mc]] = cbind(v_vec_sample_chains_list_training[[mc]],
                                                      v_vec_training)
      
      u_vec_sample_chains_list_training[[mc]] = cbind(u_vec_sample_chains_list_training[[mc]],
                                                      u_vec_training)
      
      pi_vec_sample_chains_list_training[[mc]] = cbind(pi_vec_sample_chains_list_training[[mc]],
                                                       pi_vec_training)
      
      delta_sample_chains_list_training[[mc]] = cbind(delta_sample_chains_list_training[[mc]],
                                                      delta_training)
      
      beta_sample_chains_list_training[[mc]] = cbind(beta_sample_chains_list_training[[mc]],
                                                     beta_training)
      
      beta_star_sample_chains_list_training[[mc]] = cbind(beta_star_sample_chains_list_training[[mc]],
                                                          beta_star_training)
      
      # sigma_w_sq_sample_chains_list_training[[mc]] = c(sigma_w_sq_sample_chains_list_training[[mc]],
      #                                                 sigma_w_sq_training)
      # 
      # sigma_v_sq_sample_chains_list_training[[mc]] = c(sigma_v_sq_sample_chains_list_training[[mc]],
      #                                                      sigma_v_sq_training)
      # 
      # sigma_eps_sq_sample_chains_list_training[[mc]] = c(sigma_eps_sq_sample_chains_list_training[[mc]],
      #                                                 sigma_eps_sq_training)
      # 
      # sigma_eps_sq_star_sample_chains_list_training[[mc]] = c(sigma_eps_sq_star_sample_chains_list_training[[mc]],
      #                                                      sigma_eps_sq_star_training)
      
      t0_sample_chains_list_training[[mc]] = c(t0_sample_chains_list_training[[mc]],
                                               t_0_training)
      
    }
    
    mcmc_obj_list_training[[i]] = cbind(c(beta_mat_preconvergence_training[,1],
                                          beta_star_mat_preconvergence_training[,1]), 
                                        c(beta_mat_preconvergence_training[,2],
                                          beta_star_mat_preconvergence_training[,2]), 
                                        c(beta_mat_preconvergence_training[,3],
                                          beta_star_mat_preconvergence_training[,3]))
    
    
    if( i %% 50 == 0){
      
      gelmanrubinstat = vector()
      tot_gel_var = nrow(mcmc_obj_list_training[[i]])
      
      # check for gelman rubin stat
      for (gl in 1:(tot_gel_var)){
        
        plist = lapply(apply(do.call(rbind, lapply(mcmc_obj_list_training,function(x) x[gl,])),2,
                             as.data.frame), coda::mcmc)
        gelmanrubin = gelman.diag(x = plist)
        gelmanrubinstat = c(gelmanrubinstat, gelmanrubin$psrf[,1])
        
      }
      
      no_conv_gr = length(which(gelmanrubinstat <= 1.2))
      
      if(length(which(gelmanrubinstat > 1.2)) == 0){
        converged = TRUE
        converged_i = i
        #Continuing with first chain values to collect sample
        training_df_clone$v_vec = training_df_clone$v_vec1
        #training_df_clone$w_vec = training_df_clone$w_vec1
        training_df_clone$u_vec = training_df_clone$u_vec1
        training_df_clone$pi_vec = training_df_clone$pi_vec1
        beta_training = beta_mat_preconvergence_training[,1]
        beta_star_training = beta_star_mat_preconvergence_training[,1]
        sigma_u_sq_training = sigma_u_sq_preconvergence_training[1]
        #sigma_w_sq_training = sigma_w_sq_preconvergence_training[1]
        sigma_v_sq_training = sigma_v_sq_preconvergence_training[1]
        sigma_eps_sq_training = sigma_eps_sq_preconvergence_training[1]
        sigma_eps_sq_star_training = sigma_eps_sq_star_preconvergence_training[1]
        
        phi_u_sp_training = phi_u_sp_preconvergence_training[1]
        phi_u_tmp_training = phi_u_tmp_preconvergence_training[1]
        # phi_w_sp_training = phi_w_sp_preconvergence_training[1]
        # phi_w_tmp_training = phi_w_tmp_preconvergence_training[1]
        phi_v_sp_training = phi_v_sp_preconvergence_training[1]
        phi_v_tmp_training = phi_v_tmp_preconvergence_training[1]
        
        delta_training = delta_mat_preconvergence_training[,1]
        t_0_training = t0_preconvergence_training[1]
        changepoint_t_0_week = changepoint_week_vec[t_0_training + 1]#+1 because added zero
        n_tmp_pre_changepoint_training = t_0_training
        n_tmp_post_changepoint_training = n_tmp_training - n_tmp_pre_changepoint_training
        #Saving first sample values after convergence
        v_vec_sample_training = cbind(v_vec_sample_training, v_vec_training)
        #w_vec_sample_training = cbind(w_vec_sample_training, w_vec_training)
        u_vec_sample_training = cbind(u_vec_sample_training, u_vec_training)
        pi_vec_sample_training = cbind(pi_vec_sample_training, pi_vec_training)
        delta_sample_training = cbind(delta_sample_training, delta_training)
        beta_sample_training = cbind(beta_sample_training, beta_training)
        beta_star_sample_training = cbind(beta_star_sample_training, beta_star_training)
        
        sigma_u_sq_sample_training = c(sigma_u_sq_sample_training, sigma_u_sq_training)
        #sigma_w_sq_sample_training = c(sigma_w_sq_sample_training, sigma_w_sq_training)
        sigma_v_sq_sample_training = c(sigma_v_sq_sample_training, sigma_v_sq_training)
        sigma_eps_sq_sample_training = c(sigma_eps_sq_sample_training, sigma_eps_sq_training)
        sigma_eps_sq_star_sample_training = c(sigma_eps_sq_star_sample_training, sigma_eps_sq_star_training)
        
        phi_u_sp_sample_training = c(phi_u_sp_sample_training, phi_u_sp_training)
        phi_u_tmp_sample_training = c(phi_u_tmp_sample_training, phi_u_tmp_training)
        # phi_w_sp_sample_training = c(phi_w_sp_sample_training, phi_w_sp_training)
        # phi_w_tmp_sample_training = c(phi_w_tmp_sample_training, phi_w_tmp_training)
        phi_v_sp_sample_training = c(phi_v_sp_sample_training, phi_v_sp_training)
        phi_v_tmp_sample_training = c(phi_v_tmp_sample_training, phi_v_tmp_training)
        
        t_0_sample_training = c(t_0_sample_training, t_0_training)
        
      }
    }
    
    
  }
  #Post convergence################################################################################
  #training Df pre and post change point
  # training_pre_changepoint_df = training_df_clone %>% filter(week <= changepoint_t_0_week)
  # training_post_changepoint_df = training_df_clone %>% filter(week > changepoint_t_0_week)
  
  if(converged){
    
    #pi posterior
    st_pi = Sys.time()
    pi_vec = vector()
    for (st in 1:n_total_training) {
      st_week = training_df_clone$week[st]
      y_val_st = training_df_clone$category[st]
      v_val_st = training_df_clone$v_vec[st]
      u_val_st = training_df_clone$u_vec[st]
      #w_val_st = training_df_clone$w_vec[st]
      x_vec_st = c(1, training_df_clone$scaled_prev_log_new_death[st],
                   training_df_clone$scaled_first_dose_prevalence[st])
      if(st_week <= changepoint_t_0_week){
        
        pi_st = truncnorm::rtruncnorm(1, a = delta_training[y_val_st], b = delta_training[y_val_st + 1],
                                      mean = t(x_vec_st) %*% beta_training + u_val_st ,#+ w_val_st
                                      sd = sqrt(sigma_eps_sq_training))
      }else{
        
        pi_st = truncnorm::rtruncnorm(1, a = delta_training[y_val_st], b = delta_training[y_val_st + 1],
                                      mean = t(x_vec_st) %*% beta_star_training + u_val_st + v_val_st,
                                      sd = sqrt(sigma_eps_sq_star_training))
      }
      
      pi_vec = c(pi_vec, pi_st)
    }
    training_df_clone$pi_vec = NULL
    training_df_clone$pi_vec = pi_vec
    
    pi_vec_training = training_df_clone$pi_vec
    
    training_pre_changepoint_df = training_df_clone %>% filter(week <= changepoint_t_0_week)
    training_post_changepoint_df = training_df_clone %>% filter(week > changepoint_t_0_week)
    
    #Getting X matrix for pre and post change point
    X_mat_pre_changepoint_training = cbind(rep(1,n_sp_training * n_tmp_pre_changepoint_training),
                                           training_pre_changepoint_df$scaled_prev_log_new_death,
                                           training_pre_changepoint_df$scaled_first_dose_prevalence)
    
    X_mat_post_changepoint_training = cbind(rep(1,n_sp_training * n_tmp_post_changepoint_training),
                                            training_post_changepoint_df$scaled_prev_log_new_death,
                                            training_post_changepoint_df$scaled_first_dose_prevalence)
    
    
    total_tt_pi = total_tt_pi + as.numeric(difftime(Sys.time(), st_pi), units="secs")
    ###############################################################################################################
    #v posterior#############################################################
    st_v = Sys.time()
    #Xbeta matrix
    v_mat_training = matrix((training_df_clone$v_vec),
                            nrow = n_sp_training, ncol = n_tmp_training,
                            byrow = FALSE)
    
    u_mat_training = matrix((training_df_clone$u_vec),
                            nrow = n_sp_training, ncol = n_tmp_training,
                            byrow = FALSE)
    
    # w_mat_training = matrix((training_df_clone$w_vec),
    #                         nrow = n_sp_training, ncol = n_tmp_training,
    #                         byrow = FALSE)
    
    pi_mat_training = matrix((training_df_clone$pi_vec),
                             nrow = n_sp_training,
                             ncol = n_tmp_training, byrow = FALSE)
    
    xbeta_mat_training = matrix(c(as.numeric(X_mat_pre_changepoint_training %*% beta_training),
                                  as.numeric(X_mat_post_changepoint_training %*% beta_star_training)),
                                nrow = n_sp_training,
                                ncol = n_tmp_training, byrow = FALSE)
    
    
    for (tmp in (1:n_tmp_training)){
      idx = 1:n_tmp_training
      
      if(tmp!=1){
        idx[1] = tmp
        idx[tmp] = 1
      }
      #Switching first and jth index column for V matrix and then taking values from the 2nd column
      v_cond_training = (v_mat_training %*%
                           diag(n_tmp_training)[idx,])[(1 + n_sp_training):(n_total_training)]
      
      
      
      #Can be made faster by doing this for all change points or instead of doing it for all check if for
      #the change point this is already done if so get the value from there
      mu_cond_training = ((t(SIGMA_v_tmp_trans_21_mat_training[,tmp]) %*%
                             inv_SIGMA_v_tmp_trans_22_matlist_training[[tmp]]) %x%
                            ident_n_sp_mat_training) %*% v_cond_training
      
      
      if(tmp > t_0_training){
        
        inv_SIGMA_cond_training = (inv_SIGMA_v_sp_training/mat_mult_temp_SIGMA_v_tmp_vec[tmp])
        
        v_covar_training = chol2inv(chol(inv_SIGMA_cond_training/sigma_v_sq_training +
                                           ident_n_sp_mat_training/sigma_eps_sq_star_training))
        
        
        v_mean_training  = v_covar_training %*% (inv_SIGMA_cond_training %*%
                                                   mu_cond_training/sigma_v_sq_training +
                                                   (pi_mat_training[,tmp] -
                                                      xbeta_mat_training[,tmp] - u_mat_training[,tmp])/
                                                   sigma_eps_sq_star_training)
        
        
        v_mat_training[,tmp] = v_mean_training + t(chol(v_covar_training))  %*% rnorm(n_sp_training)
        
      }else{
        
        SIGMA_cond_training = (SIGMA_v_sp_training * mat_mult_temp_SIGMA_v_tmp_vec[tmp])
        
        v_mat_training[,tmp] = (mu_cond_training + t(chol(SIGMA_cond_training))  %*% 
                                  rnorm(n_sp_training, sd = sqrt(sigma_v_sq_training)))
        
      }
      
      
    }
    
    training_df_clone$v_vec = v_mat_training[1:n_total_training]
    v_vec_training = training_df_clone$v_vec
    
    total_tt_v = total_tt_v + as.numeric(difftime(Sys.time(), st_v), units="secs")
    
    ########################################################################################################################
    #w posterior###########################################################################################
    # st_w = Sys.time()
    # 
    # for (tmp in (1:n_tmp_training)){
    #   idx = 1:n_tmp_training
    #   
    #   if(tmp!=1){
    #     idx[1] = tmp
    #     idx[tmp] = 1
    #   }
    #   #Switching first and jth index column for V matrix and then taking values from the 2nd column
    #   w_cond_training = (w_mat_training %*%
    #                        diag(n_tmp_training)[idx,])[(1 + n_sp_training):(n_total_training)]
    #   
    #   
    #   
    #   #Can be made faster by doing this for all change points or instead of doing it for all check if for
    #   #the change point this is already done if so get the value from there
    #   mu_cond_training = ((t(SIGMA_w_tmp_trans_21_mat_training[,tmp]) %*%
    #                          inv_SIGMA_w_tmp_trans_22_matlist_training[[tmp]]) %x%
    #                         ident_n_sp_mat_training) %*% w_cond_training
    #   
    #   
    #   if(tmp <= t_0_training){
    #     
    #     inv_SIGMA_cond_training = (inv_SIGMA_w_sp_training/mat_mult_temp_SIGMA_w_tmp_vec[tmp])
    #     
    #     w_covar_training = chol2inv(chol(inv_SIGMA_cond_training/sigma_w_sq_training +
    #                                        ident_n_sp_mat_training/sigma_eps_sq_training))
    #     
    #     
    #     w_mean_training  = w_covar_training %*% (inv_SIGMA_cond_training %*%
    #                                                mu_cond_training/sigma_w_sq_training +
    #                                                (pi_mat_training[,tmp] -
    #                                                   xbeta_mat_training[,tmp] - u_mat_training[,tmp])/
    #                                                sigma_eps_sq_training)
    #     
    #     
    #     w_mat_training[,tmp] = w_mean_training + t(chol(w_covar_training))  %*% rnorm(n_sp_training)
    #     
    #   }else{
    #     
    #     SIGMA_cond_training = (SIGMA_w_sp_training * mat_mult_temp_SIGMA_w_tmp_vec[tmp])
    #     
    #     w_mat_training[,tmp] = (mu_cond_training + t(chol(SIGMA_cond_training))  %*% 
    #                               rnorm(n_sp_training, sd = sqrt(sigma_w_sq_training)))
    #     
    #   }
    #   
    #   
    # }
    # training_df_clone$w_vec = w_mat_training[1:n_total_training]
    # w_vec_training = training_df_clone$w_vec
    # 
    # total_tt_w = total_tt_w + as.numeric(difftime(Sys.time(), st_w), units="secs")
    #u posterior##################################################################################################
    st_u = Sys.time()
    
    for (tmp in (1:n_tmp_training)){
      idx = 1:n_tmp_training
      
      if(tmp!=1){
        idx[1] = tmp
        idx[tmp] = 1
      }
      #Switching first and jth index column for V matrix and then taking values from the 2nd column
      u_cond_training = (u_mat_training %*%
                           diag(n_tmp_training)[idx,])[(1 + n_sp_training):(n_total_training)]
      
      
      
      #Can be made faster by doing this for all change points or instead of doing it for all check if for
      #the change point this is already done if so get the value from there
      mu_cond_training = ((t(SIGMA_u_tmp_trans_21_mat_training[,tmp]) %*%
                             inv_SIGMA_u_tmp_trans_22_matlist_training[[tmp]]) %x%
                            ident_n_sp_mat_training) %*% u_cond_training
      
      #Can be made faster by dynamically saving values in a vector
      inv_SIGMA_cond_training = (inv_SIGMA_u_sp_training/mat_mult_temp_SIGMA_u_tmp_vec[tmp])
      
      
      
      if(tmp <= t_0_training){
        
        
        
        u_covar_training = chol2inv(chol(inv_SIGMA_cond_training/sigma_u_sq_training +
                                           ident_n_sp_mat_training/sigma_eps_sq_training))
        
        
        u_mean_training  = u_covar_training %*% (inv_SIGMA_cond_training %*%
                                                   mu_cond_training/sigma_u_sq_training +
                                                   (pi_mat_training[,tmp] -
                                                      xbeta_mat_training[,tmp])/
                                                   sigma_eps_sq_training)
        
        
        u_mat_training[,tmp] = u_mean_training + t(chol(u_covar_training))  %*% rnorm(n_sp_training)
        
      }else{
        
        u_covar_training = chol2inv(chol(inv_SIGMA_cond_training/sigma_u_sq_training +
                                           ident_n_sp_mat_training/sigma_eps_sq_star_training))
        
        u_mean_training  = u_covar_training %*% (inv_SIGMA_cond_training %*%
                                                   mu_cond_training/sigma_u_sq_training +
                                                   (pi_mat_training[,tmp] -
                                                      xbeta_mat_training[,tmp] - v_mat_training[,tmp])/
                                                   sigma_eps_sq_star_training)
        
        u_mat_training[,tmp] = u_mean_training + t(chol(u_covar_training))  %*% rnorm(n_sp_training)
        
      }
      
      
    }
    training_df_clone$u_vec = u_mat_training[1:n_total_training]
    u_vec_training = training_df_clone$u_vec
    
    total_tt_u = total_tt_u + as.numeric(difftime(Sys.time(), st_u), units="secs")
    
    ########################################################################################################################
    #beta posterior###################################
    st_beta = Sys.time()
    pi_pre_changepoint_training = pi_mat_training[1:(n_tmp_pre_changepoint_training*n_sp_training)]
    pi_post_changepoint_training = pi_mat_training[(n_tmp_pre_changepoint_training*n_sp_training+1):n_total_training]
    
    #w_minus_training = w_mat_training[1:(n_tmp_pre_changepoint_training*n_sp_training)]
    v_plus_training = v_mat_training[(n_tmp_pre_changepoint_training*n_sp_training+1):n_total_training]
    u_minus_training = u_mat_training[1:(n_tmp_pre_changepoint_training*n_sp_training)]
    u_plus_training = u_mat_training[(n_tmp_pre_changepoint_training*n_sp_training+1):n_total_training]
    #dividing by 10^4 for large variance for beta prior
    beta_covar_training = chol2inv(chol( (t(X_mat_pre_changepoint_training) %*% 
                                            X_mat_pre_changepoint_training)/sigma_eps_sq_training)) #+ diag(p+1)/10^4
    
    beta_star_covar_training = chol2inv(chol( (t(X_mat_post_changepoint_training) %*% 
                                                 X_mat_post_changepoint_training)/sigma_eps_sq_star_training))
    
    beta_mean_training = beta_covar_training %*% ( t(X_mat_pre_changepoint_training)%*% 
                                                     (pi_pre_changepoint_training - u_minus_training))/sigma_eps_sq_training
    
    beta_star_mean_training = (beta_star_covar_training %*% (t(X_mat_post_changepoint_training)%*% 
                                                               (pi_post_changepoint_training - u_plus_training 
                                                                - v_plus_training))/sigma_eps_sq_star_training)
    
    beta_training = beta_mean_training + t(chol(beta_covar_training)) %*% rnorm(p+1)
    
    beta_star_training = beta_star_mean_training + t(chol(beta_star_covar_training)) %*% rnorm(p+1)
    
    training_df_clone$x_beta_vec = c(as.numeric(X_mat_pre_changepoint_training %*% beta_training),
                                     as.numeric(X_mat_post_changepoint_training %*%beta_star_training))
    
    total_tt_beta = total_tt_beta + as.numeric(difftime( Sys.time(), st_beta), units="secs")
    #######################################################################################################################
    #Sigma1^2 posterior######################################################
    #st_sigma_w_sq = Sys.time()
    
    # sigma_w_sq_a_training = a_training + (n_sp_training * n_tmp_pre_changepoint_training)/2
    # 
    # sigma_v_sq_a_training = a_training + (n_sp_training * n_tmp_post_changepoint_training)/2
    # 
    # sigma_w_sq_lambda_training = (t(v_training) %*%
    #                                inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[t_0_training/gap_bw_changepoints - 1]]
    #                              %*% v_training)/2 + lambda_training
    # 
    # sigma_v_sq_lambda_training = (t(v_star_training) %*%
    #                                     inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_mat_list[[t_0_training/gap_bw_changepoints - 1]]
    #                                   %*% v_star_training)/2 + lambda_training
    
    #sigma_w_sq_training = 1#invgamma::rinvgamma(1,sigma_w_sq_a_training,rate = sigma_w_sq_lambda_training)
    
    sigma_v_sq_training = 1#invgamma::rinvgamma(1, sigma_v_sq_a_training,
    #rate = sigma_v_sq_lambda_training)
    
    #total_tt_sigma_w_sq = total_tt_sigma_w_sq + as.numeric(difftime(Sys.time(), st_sigma_w_sq), units="secs")
    ##########################################################################################################################
    #Sigma2^2 posterior############################################
    st_sigma_eps_sq = Sys.time()
    
    # sigma_eps_sq_a_training = a_training + (n_sp_training * n_tmp_pre_changepoint_training)/2
    # 
    # sigma_eps_sq_star_a_training = a_training + (n_sp_training * n_tmp_post_changepoint_training)/2
    # 
    # temp_vec = pi_pre_changepoint_training - X_mat_pre_changepoint_training%*%beta_training - v_training
    # 
    # temp_vec2 = (pi_post_changepoint_training - X_mat_post_changepoint_training %*% beta_star_training -
    #                v_star_training)
    # 
    # sigma_eps_sq_lambda_training = sum(temp_vec^2)/2 + lambda_training
    # 
    # sigma_eps_sq_star_lambda_training = sum(temp_vec2^2)/2 + lambda_training
    
    sigma_eps_sq_training = 1#invgamma::rinvgamma(1,sigma_eps_sq_a_training,rate = sigma_eps_sq_lambda_training)
    
    sigma_eps_sq_star_training = 1#invgamma::rinvgamma(1, sigma_eps_sq_star_a_training,
    #rate = sigma_eps_sq_star_lambda_training)
    
    total_tt_sigma_eps_sq = total_tt_sigma_eps_sq + as.numeric(difftime(Sys.time(), st_sigma_eps_sq), units="secs")
    
    #phi_u_sp posterior arms########################################################################################
    if(i==1 | i %% 20 == 0){
      
      st_phi_u_sp = Sys.time()
      
      phi_u_sp_log_density_vec_training = vector()
      for(phi_u_sp_val in  phi_sp_vec_training){
        
        phi_u_sp_log_density_vec_training = c(phi_u_sp_log_density_vec_training, 
                                              phi_u_sp_log_density_func(phi_u_sp_val))
        
        
      }
      
      phi_u_sp_log_density_vec_training = (phi_u_sp_log_density_vec_training - max(phi_u_sp_log_density_vec_training))
      
      phi_u_sp_density_vec_training = exp(phi_u_sp_log_density_vec_training)
      
      const_prop_phi_u_sp = 1/sum(phi_u_sp_density_vec_training)
      
      phi_u_sp_training = sample(x = (phi_sp_vec_training), 1, replace = FALSE,
                                 prob = round(phi_u_sp_density_vec_training * const_prop_phi_u_sp, digits = 3))
      
      phi_u_sp_preconvergence_training[mc] = phi_u_sp_training
      total_tt_phi_u_sp = total_tt_phi_u_sp + as.numeric(difftime(Sys.time(), st_phi_u_sp), units="secs")
    }
    #phi_u_tmp posterior arms########################################################################################
    if(i==1 | i %% 20 == 0){
      
      st_phi_u_tmp = Sys.time()
      
      phi_u_tmp_log_density_vec_training = vector()
      for(phi_u_tmp_val in  phi_tmp_vec_training){
        
        phi_u_tmp_log_density_vec_training = c(phi_u_tmp_log_density_vec_training, 
                                               phi_u_tmp_log_density_func(phi_u_tmp_val))
        
        
      }
      
      phi_u_tmp_log_density_vec_training = (phi_u_tmp_log_density_vec_training - max(phi_u_tmp_log_density_vec_training))
      
      phi_u_tmp_density_vec_training = exp(phi_u_tmp_log_density_vec_training)
      
      const_prop_phi_u_tmp = 1/sum(phi_u_tmp_density_vec_training)
      
      phi_u_tmp_training = sample(x = (phi_tmp_vec_training), 1, replace = FALSE,
                                  prob = round(phi_u_tmp_density_vec_training * const_prop_phi_u_tmp, digits = 3))
      
      phi_u_tmp_preconvergence_training[mc] = phi_u_tmp_training
      
      total_tt_phi_u_tmp = total_tt_phi_u_tmp + as.numeric(difftime(Sys.time(), st_phi_u_tmp), units="secs")
    }
    #phi_w_sp posterior arms########################################################################################
    # if(i %% 20 == 0){
    #   
    #   st_phi_w_sp = Sys.time()
    #   
    #   phi_w_sp_log_density_vec_training = vector()
    #   for(phi_w_sp_val in  phi_sp_vec_training){
    #     
    #     phi_w_sp_log_density_vec_training = c(phi_w_sp_log_density_vec_training, 
    #                                           phi_w_sp_log_density_func(phi_w_sp_val))
    #     
    #     
    #   }
    #   
    #   phi_w_sp_log_density_vec_training = (phi_w_sp_log_density_vec_training - max(phi_w_sp_log_density_vec_training))
    #   
    #   phi_w_sp_density_vec_training = exp(phi_w_sp_log_density_vec_training)
    #   
    #   const_prop_phi_w_sp = 1/sum(phi_w_sp_density_vec_training)
    #   
    #   phi_w_sp_training = sample(x = (phi_sp_vec_training), 1, replace = FALSE,
    #                              prob = round(phi_w_sp_density_vec_training * const_prop_phi_w_sp, digits = 3))
    #   
    #   phi_w_sp_preconvergence_training[mc] = phi_w_sp_training
    #   total_tt_phi_w_sp = total_tt_phi_w_sp + as.numeric(difftime(Sys.time(), st_phi_w_sp), units="secs")
    # }
    # #phi_w_tmp posterior arms########################################################################################
    # if(i %% 20 == 0){
    #   
    #   st_phi_w_tmp = Sys.time()
    #   
    #   phi_w_tmp_log_density_vec_training = vector()
    #   for(phi_w_tmp_val in  phi_tmp_vec_training){
    #     
    #     phi_w_tmp_log_density_vec_training = c(phi_w_tmp_log_density_vec_training, 
    #                                            phi_w_tmp_log_density_func(phi_w_tmp_val))
    #     
    #     
    #   }
    #   
    #   phi_w_tmp_log_density_vec_training = (phi_w_tmp_log_density_vec_training - max(phi_w_tmp_log_density_vec_training))
    #   
    #   phi_w_tmp_density_vec_training = exp(phi_w_tmp_log_density_vec_training)
    #   
    #   const_prop_phi_w_tmp = 1/sum(phi_w_tmp_density_vec_training)
    #   
    #   phi_w_tmp_training = sample(x = (phi_tmp_vec_training), 1, replace = FALSE,
    #                               prob = round(phi_w_tmp_density_vec_training * const_prop_phi_w_tmp, digits = 3))
    #   
    #   phi_w_tmp_preconvergence_training[mc] = phi_w_tmp_training
    #   total_tt_phi_w_tmp = total_tt_phi_w_tmp + as.numeric(difftime(Sys.time(), st_phi_w_tmp), units="secs")
    # }
    #phi_v_sp posterior arms########################################################################################
    if(i==1 | i %% 20 == 0){
      
      st_phi_v_sp = Sys.time()
      
      phi_v_sp_log_density_vec_training = vector()
      for(phi_v_sp_val in  phi_sp_vec_training){
        
        phi_v_sp_log_density_vec_training = c(phi_v_sp_log_density_vec_training, 
                                              phi_v_sp_log_density_func(phi_v_sp_val))
        
        
      }
      
      phi_v_sp_log_density_vec_training = (phi_v_sp_log_density_vec_training - max(phi_v_sp_log_density_vec_training))
      
      phi_v_sp_density_vec_training = exp(phi_v_sp_log_density_vec_training)
      
      const_prop_phi_v_sp = 1/sum(phi_v_sp_density_vec_training)
      
      phi_v_sp_training = sample(x = (phi_sp_vec_training), 1, replace = FALSE,
                                 prob = round(phi_v_sp_density_vec_training * const_prop_phi_v_sp, digits = 3))
      
      phi_v_sp_preconvergence_training[mc] = phi_v_sp_training
      total_tt_phi_v_sp = total_tt_phi_v_sp + as.numeric(difftime(Sys.time(), st_phi_v_sp), units="secs")
    }
    #phi_v_tmp posterior arms########################################################################################
    if(i==1 | i %% 20 == 0){
      
      st_phi_v_tmp = Sys.time()
      
      phi_v_tmp_log_density_vec_training = vector()
      for(phi_v_tmp_val in  phi_tmp_vec_training){
        
        phi_v_tmp_log_density_vec_training = c(phi_v_tmp_log_density_vec_training, 
                                               phi_v_tmp_log_density_func(phi_v_tmp_val))
        
        
      }
      
      phi_v_tmp_log_density_vec_training = (phi_v_tmp_log_density_vec_training - max(phi_v_tmp_log_density_vec_training))
      
      phi_v_tmp_density_vec_training = exp(phi_v_tmp_log_density_vec_training)
      
      const_prop_phi_v_tmp = 1/sum(phi_v_tmp_density_vec_training)
      
      phi_v_tmp_training = sample(x = (phi_tmp_vec_training), 1, replace = FALSE,
                                  prob = round(phi_v_tmp_density_vec_training * const_prop_phi_v_tmp, digits = 3))
      
      phi_v_tmp_preconvergence_training[mc] = phi_v_tmp_training
      
      total_tt_phi_v_tmp = total_tt_phi_v_tmp + as.numeric(difftime(Sys.time(), st_phi_v_tmp), units="secs")
    }
    
    #Getting the matrix division for V calculations##################################################################
    #Getting inv Sigma_sp_22 to use for every iteration
    st_mat_inv = Sys.time()
    SIGMA_v_sp_training = exp(- phi_v_sp_training * distance_mat_training)
    inv_SIGMA_v_sp_training = chol2inv(chol(SIGMA_v_sp_training))
    #log_det_SIGMA_v_sp_training = (determinant(SIGMA_v_sp_training, logarithm = TRUE))$modulus[1]
    
    SIGMA_v_tmp_training = exp(- phi_v_tmp_training * week_diff_mat_training)
    inv_SIGMA_v_tmp_training = chol2inv(chol(SIGMA_v_tmp_training))
    #log_det_SIGMA_v_tmp_training = (determinant(SIGMA_v_tmp_training, logarithm = TRUE))$modulus[1]
    
    #inv_SIGMA_v_tmp_kro_inv_SIGMA_v_sp_training = inv_SIGMA_v_tmp_training %x% inv_SIGMA_v_sp_training
    
    inv_SIGMA_v_tmp_trans_22_matlist_training = list()
    mat_mult_temp_SIGMA_v_tmp_vec = vector()
    # inv_SIGMA_cond_pre_changepoint_matlist_training = list()
    # inv_SIGMA_cond_post_changepoint_matlist_training = list()
    SIGMA_v_tmp_trans_21_mat_training = vector()
    
    for (tmp in (1:n_tmp_training)) {
      idx = 1:n_tmp_training
      
      if(tmp != 1){
        idx[1] = tmp
        idx[tmp] = 1
      }
      SIGMA_v_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% SIGMA_v_tmp_training %*% 
                                      diag(n_tmp_training)[idx,])
      
      inv_SIGMA_v_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% inv_SIGMA_v_tmp_training %*% 
                                          diag(n_tmp_training)[idx,])
      
      inv_SIGMA_v_tmp_trans_a_training = inv_SIGMA_v_tmp_trans_training[1,1]
      inv_SIGMA_v_tmp_trans_b_training = inv_SIGMA_v_tmp_trans_training[1,-1]
      inv_SIGMA_v_tmp_trans_c_training = inv_SIGMA_v_tmp_trans_training[-1,1]
      inv_SIGMA_v_tmp_trans_d_training = inv_SIGMA_v_tmp_trans_training[-1,-1]
      
      inv_SIGMA_v_tmp_trans_22_training = (inv_SIGMA_v_tmp_trans_d_training - (inv_SIGMA_v_tmp_trans_c_training %*%
                                                                                 t(inv_SIGMA_v_tmp_trans_b_training))/
                                             inv_SIGMA_v_tmp_trans_a_training)
      
      SIGMA_v_tmp_trans_21_training = SIGMA_v_tmp_trans_training[-1,1]
      
      mat_mult_temp_SIGMA_v_tmp_vec[tmp] = as.numeric(1 - t(SIGMA_v_tmp_trans_21_training) %*%
                                                        inv_SIGMA_v_tmp_trans_22_training %*%
                                                        SIGMA_v_tmp_trans_21_training)
      
      inv_SIGMA_v_tmp_trans_22_matlist_training[[tmp]] = inv_SIGMA_v_tmp_trans_22_training
      SIGMA_v_tmp_trans_21_mat_training = cbind(SIGMA_v_tmp_trans_21_mat_training, SIGMA_v_tmp_trans_21_training)
      
    }
    #Getting the matrix division for W calculations##################################################################
    # SIGMA_w_sp_training = exp(- phi_w_sp_training * distance_mat_training)
    # inv_SIGMA_w_sp_training = chol2inv(chol(SIGMA_w_sp_training))
    # #log_det_SIGMA_w_sp_training = (determinant(SIGMA_w_sp_training, logarithm = TRUE))$modulus[1]
    # 
    # 
    # SIGMA_w_tmp_training = exp(- phi_w_tmp_training * week_diff_mat_training)
    # inv_SIGMA_w_tmp_training = chol2inv(chol(SIGMA_w_tmp_training))
    # #log_det_SIGMA_w_tmp_training = (determinant(SIGMA_w_tmp_training, logarithm = TRUE))$modulus[1]
    # 
    # #inv_SIGMA_w_tmp_kro_inv_SIGMA_w_sp_training = inv_SIGMA_w_tmp_training %x% inv_SIGMA_w_sp_training
    # 
    # inv_SIGMA_w_tmp_trans_22_matlist_training = list()
    # mat_mult_temp_SIGMA_w_tmp_vec = vector()
    # 
    # SIGMA_w_tmp_trans_21_mat_training = vector()
    # 
    # for (tmp in (1:n_tmp_training)) {
    #   idx = 1:n_tmp_training
    #   
    #   if(tmp != 1){
    #     idx[1] = tmp
    #     idx[tmp] = 1
    #   }
    #   SIGMA_w_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% SIGMA_w_tmp_training %*% 
    #                                   diag(n_tmp_training)[idx,])
    #   
    #   inv_SIGMA_w_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% inv_SIGMA_w_tmp_training %*% 
    #                                       diag(n_tmp_training)[idx,])
    #   
    #   inv_SIGMA_w_tmp_trans_a_training = inv_SIGMA_w_tmp_trans_training[1,1]
    #   inv_SIGMA_w_tmp_trans_b_training = inv_SIGMA_w_tmp_trans_training[1,-1]
    #   inv_SIGMA_w_tmp_trans_c_training = inv_SIGMA_w_tmp_trans_training[-1,1]
    #   inv_SIGMA_w_tmp_trans_d_training = inv_SIGMA_w_tmp_trans_training[-1,-1]
    #   
    #   inv_SIGMA_w_tmp_trans_22_training = (inv_SIGMA_w_tmp_trans_d_training - (inv_SIGMA_w_tmp_trans_c_training %*%
    #                                                                              t(inv_SIGMA_w_tmp_trans_b_training))/
    #                                          inv_SIGMA_w_tmp_trans_a_training)
    #   
    #   SIGMA_w_tmp_trans_21_training = SIGMA_w_tmp_trans_training[-1,1]
    #   
    #   mat_mult_temp_SIGMA_w_tmp_vec[tmp] = as.numeric(1 - t(SIGMA_w_tmp_trans_21_training) %*%
    #                                                     inv_SIGMA_w_tmp_trans_22_training %*%
    #                                                     SIGMA_w_tmp_trans_21_training)
    #   
    #   inv_SIGMA_w_tmp_trans_22_matlist_training[[tmp]] = inv_SIGMA_w_tmp_trans_22_training
    #   SIGMA_w_tmp_trans_21_mat_training = cbind(SIGMA_w_tmp_trans_21_mat_training, SIGMA_w_tmp_trans_21_training)
    #   
    # }
    
    #Getting the matrix division for U calculations##################################################################
    
    SIGMA_u_sp_training = exp(- phi_u_sp_training * distance_mat_training)
    inv_SIGMA_u_sp_training = chol2inv(chol(SIGMA_u_sp_training))
    #log_det_SIGMA_u_sp_training = (determinant(SIGMA_u_sp_training, logarithm = TRUE))$modulus[1]
    
    SIGMA_u_tmp_training = exp(- phi_u_tmp_training * week_diff_mat_training)
    inv_SIGMA_u_tmp_training = chol2inv(chol(SIGMA_u_tmp_training))
    #log_det_SIGMA_u_tmp_training = (determinant(SIGMA_u_tmp_training, logarithm = TRUE))$modulus[1]
    
    #inv_SIGMA_u_tmp_kro_inv_SIGMA_u_sp_training = inv_SIGMA_u_tmp_training %x% inv_SIGMA_u_sp_training
    
    inv_SIGMA_u_tmp_trans_22_matlist_training = list()
    mat_mult_temp_SIGMA_u_tmp_vec = vector()
    # inv_SIGMA_cond_pre_changepoint_matlist_training = list()
    # inv_SIGMA_cond_post_changepoint_matlist_training = list()
    SIGMA_u_tmp_trans_21_mat_training = vector()
    
    for (tmp in (1:n_tmp_training)) {
      idx = 1:n_tmp_training
      
      if(tmp != 1){
        idx[1] = tmp
        idx[tmp] = 1
      }
      SIGMA_u_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% SIGMA_u_tmp_training %*% 
                                      diag(n_tmp_training)[idx,])
      
      inv_SIGMA_u_tmp_trans_training = (diag(n_tmp_training)[idx,] %*% inv_SIGMA_u_tmp_training %*% 
                                          diag(n_tmp_training)[idx,])
      
      inv_SIGMA_u_tmp_trans_a_training = inv_SIGMA_u_tmp_trans_training[1,1]
      inv_SIGMA_u_tmp_trans_b_training = inv_SIGMA_u_tmp_trans_training[1,-1]
      inv_SIGMA_u_tmp_trans_c_training = inv_SIGMA_u_tmp_trans_training[-1,1]
      inv_SIGMA_u_tmp_trans_d_training = inv_SIGMA_u_tmp_trans_training[-1,-1]
      
      inv_SIGMA_u_tmp_trans_22_training = (inv_SIGMA_u_tmp_trans_d_training - (inv_SIGMA_u_tmp_trans_c_training %*%
                                                                                 t(inv_SIGMA_u_tmp_trans_b_training))/
                                             inv_SIGMA_u_tmp_trans_a_training)
      
      SIGMA_u_tmp_trans_21_training = SIGMA_u_tmp_trans_training[-1,1]
      
      mat_mult_temp_SIGMA_u_tmp_vec[tmp] = as.numeric(1 - t(SIGMA_u_tmp_trans_21_training) %*%
                                                        inv_SIGMA_u_tmp_trans_22_training %*%
                                                        SIGMA_u_tmp_trans_21_training)
      
      inv_SIGMA_u_tmp_trans_22_matlist_training[[tmp]] = inv_SIGMA_u_tmp_trans_22_training
      SIGMA_u_tmp_trans_21_mat_training = cbind(SIGMA_u_tmp_trans_21_mat_training, SIGMA_u_tmp_trans_21_training)
      
    }
    total_tt_mat_inv = total_tt_mat_inv + as.numeric(difftime(Sys.time(), st_mat_inv), units="secs")
    ###################################################################################################################
    #delta posterior###############################
    st_delta = Sys.time()
    
    
    delta3_slice = diversitree::mcmc(lik = delta3_dens_func, x.init = c(delta_training[3]), nsteps = 1,w= 5, lower = 0, 
                                     upper = delta_training[4])
    
    delta_training[3] =delta3_slice$pars
    
    delta4_slice = diversitree::mcmc(lik = delta4_dens_func, x.init = c(delta_training[4]), nsteps = 1,w= 5, lower = delta_training[3], 
                                     upper = Inf)
    
    delta_training[4] = delta4_slice$pars
    
    total_tt_delta = total_tt_delta + as.numeric(difftime(Sys.time(), st_delta), units="secs")
    ################################################################################################################
    #t_0 posterior####################################
    st_t_0 = Sys.time()
    log_prop_pmf_t_0_vec_training = vector()
    for (tp in (changepoint_vec)) {
      
      week_t = changepoint_week_vec[tp + 1]#+1 because added zero
      
      
      training_post_t_df = training_df_clone %>% filter(week > week_t)
      
      
      if(tp != 0){
        training_pre_t_df = training_df_clone %>% filter(week <= week_t)
        X_mat_pre_t = cbind(rep(1,n_sp_training * tp),
                            training_pre_t_df$scaled_prev_log_new_death,
                            training_pre_t_df$scaled_first_dose_prevalence)
        
        u_vec_pre_t = training_pre_t_df$u_vec
        pi_vec_pre_t = training_pre_t_df$pi_vec
        
        temp_vec = pi_vec_pre_t - X_mat_pre_t %*% beta_training - u_vec_pre_t
      }
      
      
      if(tp != n_tmp_training){
        
        X_mat_post_t = cbind(rep(1, n_sp_training * (n_tmp_training - tp)),
                             training_post_t_df$scaled_prev_log_new_death,
                             training_post_t_df$scaled_first_dose_prevalence)
        
        #w_vec_pre_t = training_pre_t_df$w_vec
        v_vec_post_t = training_post_t_df$v_vec
        
        u_vec_post_t = training_post_t_df$u_vec
        
        pi_vec_post_t = training_post_t_df$pi_vec
        
        #- w_vec_pre_t
        temp_vec2 = pi_vec_post_t - X_mat_post_t %*% beta_star_training - u_vec_post_t - v_vec_post_t
        n_tmp_post_t = n_tmp_training - tp
        
      }
      
      if(tp == 0){
        
        log_prop_pmf_t_0_training = (((-n_sp_training*(n_tmp_post_t)/2) *
                                        log(2*pi*sigma_eps_sq_star_training)) +
                                       ((-1/(2*sigma_eps_sq_star_training)) * sum(temp_vec2^2)))
        
        
      }else if(tp == n_tmp_training){
        
        log_prop_pmf_t_0_training = ((-n_sp_training * tp/2)*log(2*pi*sigma_eps_sq_training) +
                                       ((-1/(2*sigma_eps_sq_training)) * sum(temp_vec^2)))
        
      }else{
        
        log_prop_pmf_t_0_training = ((-n_sp_training * tp/2)*log(2*pi*sigma_eps_sq_training) +
                                       ((-1/(2*sigma_eps_sq_training)) * sum(temp_vec^2)) +
                                       
                                       ((-n_sp_training*(n_tmp_post_t)/2) *
                                          log(2*pi*sigma_eps_sq_star_training)) +
                                       ((-1/(2*sigma_eps_sq_star_training)) * sum(temp_vec2^2)))
        
      }
      
      
      log_prop_pmf_t_0_vec_training = c(log_prop_pmf_t_0_vec_training, log_prop_pmf_t_0_training)
    }
    #min_prop_pmf_t_0_vec_training = min(prop_pmf_t_0_vec_training)
    max_log_prop_pmf_t_0_vec_training = max(log_prop_pmf_t_0_vec_training)
    # if(min_prop_pmf_t_0_vec_training < 0){
    #   prop_pmf_t_0_vec_training = prop_pmf_t_0_vec_training + (-min_prop_pmf_t_0_vec_training)
    # }
    log_prop_pmf_t_0_vec_training_minus_max = log_prop_pmf_t_0_vec_training - max_log_prop_pmf_t_0_vec_training
    prop_pmf_t_0_vec_training = exp(log_prop_pmf_t_0_vec_training_minus_max)
    
    const_prop = 1/sum(prop_pmf_t_0_vec_training)
    
    t_0_training = sample(x = (changepoint_vec),1, replace = FALSE,
                          prob = round(prop_pmf_t_0_vec_training * const_prop, digits = 3))
    
    n_tmp_pre_changepoint_training = t_0_training
    n_tmp_post_changepoint_training = n_tmp_training - n_tmp_pre_changepoint_training
    changepoint_t_0_week = changepoint_week_vec[n_tmp_pre_changepoint_training + 1]#+1 because added zero
    
    total_tt_t_0 = total_tt_t_0 + as.numeric(difftime(Sys.time(), st_t_0), units="secs")
    
    #Collecting samples for prediction
    if(i %% diff_in_random_draws_training == 0){
      
      v_vec_sample_training = cbind(v_vec_sample_training, v_vec_training)
      #w_vec_sample_training = cbind(w_vec_sample_training, w_vec_training)
      u_vec_sample_training = cbind(u_vec_sample_training, u_vec_training)
      pi_vec_sample_training = cbind(pi_vec_sample_training, pi_vec_training)
      delta_sample_training = cbind(delta_sample_training, delta_training)
      beta_sample_training = cbind(beta_sample_training, beta_training)
      beta_star_sample_training = cbind(beta_star_sample_training, beta_star_training)
      
      sigma_u_sq_sample_training = c(sigma_u_sq_sample_training, sigma_u_sq_training)
      #sigma_w_sq_sample_training = c(sigma_w_sq_sample_training, sigma_w_sq_training)
      sigma_v_sq_sample_training = c(sigma_v_sq_sample_training, sigma_v_sq_training)
      sigma_eps_sq_sample_training = c(sigma_eps_sq_sample_training, sigma_eps_sq_training)
      sigma_eps_sq_star_sample_training = c(sigma_eps_sq_star_sample_training, sigma_eps_sq_star_training)
      
      phi_u_sp_sample_training = c(phi_u_sp_sample_training, phi_u_sp_training)
      phi_u_tmp_sample_training = c(phi_u_tmp_sample_training, phi_u_tmp_training)
      # phi_w_sp_sample_training = c(phi_w_sp_sample_training, phi_w_sp_training)
      # phi_w_tmp_sample_training = c(phi_w_tmp_sample_training, phi_w_tmp_training)
      phi_v_sp_sample_training = c(phi_v_sp_sample_training, phi_v_sp_training)
      phi_v_tmp_sample_training = c(phi_v_tmp_sample_training, phi_v_tmp_training)
      t_0_sample_training = c(t_0_sample_training, t_0_training)
      
    }
    
    if(length(phi_v_sp_sample_training) == sample_size_training){
      save.image(paste0("main_cp_sptmp_UV_model_gelman_delta_slice_vectorized_discrete_phi_vaccine_data_10-90_cont_tp_closeTo60_incl_0_T_fewer_covar_1st_vaccine_updated",i,"iter_completed.RData"))
      break
    }
    
  }
  
  if(i%% 1000 == 0){
    save.image(paste0(
      "main_cp_sptmp_UV_model_gelman_delta_slice_vectorized_discrete_phi_vaccine_data_10-90_cont_tp_closeTo60_incl_0_T_fewer_covar_1st_vaccine_updated",i,"iter.RData"))
  }
  i = i+1
}
#time analysis
avg_tt_v = total_tt_v/(i)
avg_tt_beta = total_tt_beta/(i)
#avg_tt_sigma_w_sq = total_tt_sigma_w_sq/(i)
avg_tt_sigma_eps_sq = total_tt_sigma_eps_sq/(i)
avg_tt_delta = total_tt_delta/(i)
avg_tt_pi = total_tt_pi/(i)
avg_tt_t_0 = total_tt_t_0/(i)

sum(avg_tt_v,avg_tt_beta, avg_tt_sigma_w_sq, avg_tt_sigma_eps_sq, avg_tt_delta, avg_tt_pi, avg_tt_t_0)
#estimation of y with the posterior parameters
# estimated_pi_training = X_mat_pre_changepoint_training %*% rowMeans(beta_sample_training) + 
#   rowMeans(v_sample_training) + rnorm(n_total_training,mean = 0, sd= sqrt(mean(sigma_eps_sq_sample_training)))
# 
# estimated_pi_star_training = X_mat_post_changepoint_training %*% rowMeans(beta_star_sample_training) + 
#   rowMeans(v_star_sample_training) + rnorm(n_total_training,mean = 0, sd= sqrt(mean(sigma_eps_sq_star_sample_training)))

#####################################################################################################

training_df_clone$pi_vec_est = rowMeans(pi_vec_sample_training)

mean_delta_training = rowMeans(delta_sample_training)

training_df_clone = training_df_clone %>% 
  mutate(pred_category = ifelse(pi_vec_est > mean_delta_training[4],4,
                                ifelse(pi_vec_est > mean_delta_training[3], 3,
                                       ifelse(pi_vec_est > mean_delta_training[2],2,1))))
pred_category = training_df_clone$pred_category

pred_category_mat = cbind(pred_category_mat, pred_category)

accuracy_pred = sum(training_df_clone$category == training_df_clone$pred_category)/n_total_training

accuracy_pred_vec = c(accuracy_pred_vec, accuracy_pred)
# mape_training = (sum(abs( (y_training - estimated_pi_training)/y_training )))/n_total_training
# mape_vec_training_training = c(mape_vec_training_training, mape_training)



#Testing for significance
beta_significance_vec = vector()
for (ts in (1:(p+1))) {
  beta_significance_vec = c(beta_significance_vec, 
                            unname(!(0 >= quantile(beta_sample_training[ts,], 0.025) & 
                                       0 <= quantile(beta_sample_training[ts,], 0.975))))
}
beta_significance_mat_training = cbind(beta_significance_mat_training, beta_significance_vec)

beta_star_significance_vec = vector()
for (ts in (1:(p+1))) {
  beta_star_significance_vec = c(beta_star_significance_vec, 
                                 unname(!(0 >= quantile(beta_star_sample_training[ts,], 0.025) & 
                                            0 <= quantile(beta_star_sample_training[ts,], 0.975))))
}
beta_star_significance_mat_training = cbind(beta_star_significance_mat_training, 
                                            beta_star_significance_vec)



###########################################################################################

min_idx = which(mape_vec_validation_training == min(mape_vec_validation_training))
if(min_idx %% length(rho_vec_sp_training) == 0){
  phi_s = length(rho_vec_sp_training)
  phi_t = min_idx/length(rho_vec_sp_training)
}else{
  phi_s = rho_vec_sp_training[min_idx %% length(rho_vec_sp_training)]
  phi_t = rho_vec_tmp_training[floor(min_idx/length(rho_vec_sp_training)) + 1]
}
return(c(phi_s, phi_t))

