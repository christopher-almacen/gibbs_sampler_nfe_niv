################################################################################
# Heading
################################################################################

# Author(s):
#   Chris Almacen

################################################################################
# Goal
################################################################################

# Estimate unassigned, assigned, and decision outcomes via Gibbs Sampler
# on synthetic data

################################################################################
# Personal Settings
################################################################################

# Set (personal) working directory
setwd("~/gibbs_sampler_nfe_niv")

# Clear variables from environment
rm(list = ls())

################################################################################
# Libraries
################################################################################

library(MCMCpack) # For inverse Wishart distribution and inverse Gamma
library(truncnorm) # For truncated normal distribution
library(Matrix) # For sparse matrices
library(Rcpp)

################################################################################
# Functions
################################################################################

# ------------------------------------------------------------------------------
# C++ functions
# ------------------------------------------------------------------------------

sourceCpp("functions/crossprod_armadillo.cpp")
sourceCpp("functions/draw_coef.cpp")
sourceCpp("functions/posterior_mean_one_outcome.cpp")

# ------------------------------------------------------------------------------
# apply row-wise the standard deviation operation
# ------------------------------------------------------------------------------

rowSd <- function(x, na.rm = FALSE) {
  apply(x, 1, sd, na.rm = na.rm)
}

################################################################################
# Load Data
################################################################################

load("input/indices_decisions_patients.RData")
load("input/indices_decisions_organs.RData")
load("input/survival_outcomes.RData")
load("input/binary_decisions.RData")
load("input/indicator_assigned.RData")
load("input/indicator_unassigned.RData")
load("input/indicator_censor.RData")
load("input/X_obs.RData")
load("input/Q_obs.RData")
load("input/interactions_obs.RData")

################################################################################
# Gibbs Sampler Set-Up
################################################################################

# ------------------------------------------------------------------------------
# Data
# ------------------------------------------------------------------------------

Q_obs <- cbind(Q_obs, interactions_obs)

# ------------------------------------------------------------------------------
# Dimensions of data
# ------------------------------------------------------------------------------

n_obs <- length(survival_outcomes) # Number of patient-organ pair observations
n_patients <- length(unique(indices_decisions_patients)) # Number of unique patients
n_organs <- length(unique(indices_decisions_organs)) # Number of unique organs (donors)
k_patients <- dim(X_obs)[2] - 1 # Number of patient characteristics
k_organs <- dim(Q_obs)[2] # Number of organ characteristics

# ------------------------------------------------------------------------------
# Initial Guess
# ------------------------------------------------------------------------------

set.seed(0)
beta_x_c <- as.numeric(rep(0, k_patients + 1))
beta_f_1_c <- 1
beta_f_2_c <- 1
alpha_x_c <- as.numeric(rep(0, k_patients + 1))
alpha_q_c <- as.numeric(rep(0, k_organs))
alpha_f_1_c <- 1
alpha_f_2_c <- 1
alpha_f_3_c <- 1
alpha_f_4_c <- 1
gamma_x_c <- as.numeric(rep(0, k_patients + 1))
gamma_q_c <- as.numeric(rep(0, k_organs))
# gamma_z_c <- 1
sigma_squared_1_c <- 1
sigma_squared_2_c <- 1
sigma_squared_4_c <- 1
sigma_squared_tilde_0_c <- 1
sigma_squared_tilde_1_c <- 1
f_1_c <- rnorm(n_patients, 0, 1)
f_2_c <- rnorm(n_patients, 0, 1)
f_3_c <- rnorm(n_obs, 0, 1)
f_4_c <- rnorm(n_organs, 0, 1)
y_unassigned_c <- survival_outcomes[indicator_unassigned == 1]
y_unassigned_full_c <- as.numeric(survival_outcomes * indicator_unassigned)
y_assigned_c <- survival_outcomes[indicator_assigned == 1]
y_assigned_full_c <- as.numeric(survival_outcomes * indicator_assigned)
y_decision_c <- rnorm(n_obs, 0, 1)

# ------------------------------------------------------------------------------
# Diffuse priors
# ------------------------------------------------------------------------------

prior_variance <- 1000 # Prior variance for coefficients
prior_shape <- 3 # Prior shape parameter for inverse-Wishart for variances
prior_scale <- 3 # Prior scale parameter for inverse-Wishart for variances

# ------------------------------------------------------------------------------
# Preliminary Computations
# ------------------------------------------------------------------------------

# Dimensions of data
n_unassigned <- sum(indicator_unassigned) # Number of patients that did not receive a "transplanted kidney"
n_assigned <- sum(indicator_assigned) # Number of patients that received a "transplanted kidney"
k_patients <- dim(X_obs)[2] - 1 # Number of patient characteristics
k_organs <- dim(Q_obs)[2] # Number of organ characteristics
k_instruments <- 1 # Number of instrumental variables

# Compute some preliminary rows and columns
col_ones_n_obs <- rep(1, n_obs) # Column of ones with length = number of patient-organ pairs

# Indicators for step 1
y_decision_accept <- binary_decisions == 1 # True if accept
y_decision_reject <- !y_decision_accept # True if reject
n_accept <- sum(y_decision_accept) # Number of observed accepts
n_reject <- sum(y_decision_reject) # Number of observed rejects
indicator_censor_unassigned <- indicator_censor * indicator_unassigned == 1
n_censor_unassigned <- sum(indicator_censor_unassigned) # Number of censored survival outcomes for unassigned
indicator_censor_assigned <- indicator_censor * indicator_assigned == 1
n_censor_assigned <- sum(indicator_censor_assigned) # Number of censored survival outcomes for assigned

# Extend RHS variables
# Z_obs <- matrix(Z_obs, ncol = 1)
X_unassigned <- X_obs[indicator_unassigned == 1, ] # Patient characteristics for observed unassigned outcomes
X_assigned <- X_obs[indicator_assigned == 1, ] # Patient characteristics given assignment
Q_assigned <- Q_obs[indicator_assigned == 1, ] # Organ characteristics given assignment

# Coefficients
beta_c <- as.numeric(c(beta_x_c, beta_f_1_c, beta_f_2_c))
alpha_c <-as.numeric(c(alpha_x_c, alpha_q_c, alpha_f_1_c, 1, alpha_f_3_c, alpha_f_4_c))

# Degrees of freedom
df_unassigned <- length(c(beta_x_c, beta_f_1_c, beta_f_2_c)) # Degrees of freedom for unassigned outcomes
df_assigned <- length(c(alpha_x_c, alpha_q_c, alpha_f_1_c, 1, alpha_f_3_c, alpha_f_4_c)) # Degrees of freedom for assigned outcomes
df_decision <- length(c(gamma_x_c, gamma_q_c, 1, 1, 1)) # Degrees of freedom for decision outcomes
# df_decision <- length(c(gamma_x_c, gamma_q_c, gamma_z_c, 1, 1, 1)) # Degrees of freedom for decision outcomes

# Compute diffuse prior matrices (for step 2)
prior_mat_beta <- diag(df_unassigned) * (1 / prior_variance)
prior_mat_alpha <- diag(df_assigned) * (1 / prior_variance)
prior_diag_alpha <- as.numeric(rep((1 / prior_variance), df_assigned))
prior_mat_gamma <- diag(df_decision) * (1 / prior_variance)

# Compute posterior shapes (for steps 3 and 5)
posterior_shape_unassigned <- prior_shape + (n_unassigned / 2)
posterior_shape_assigned <- prior_shape + (n_assigned / 2)
posterior_shape_f_1 <- prior_shape + n_patients / 2
posterior_shape_f_2 <- prior_shape + n_patients / 2
posterior_shape_f_4 <- prior_shape + n_organs / 2

# Sparse matrices of observed decisions
obs_decisions_patients <- sparseMatrix(i = seq_len(n_obs), j = indices_decisions_patients, dims = c(n_obs, n_patients), x = 1)
obs_decisions_organs <- sparseMatrix(i = seq_len(n_obs), j = indices_decisions_organs, dims = c(n_obs, n_organs), x = 1)

# Weights with respect to patients
weight_unassigned_patients <- drop(crossprod(obs_decisions_patients, indicator_unassigned))
weight_assigned_patients <- drop(crossprod(obs_decisions_patients, indicator_assigned))
weight_decision_patients <- drop(crossprod(obs_decisions_patients, col_ones_n_obs))

# Weights with respect to organs
weight_assigned_organs <- drop(crossprod(obs_decisions_organs, indicator_assigned))
weight_decision_organs <- drop(crossprod(obs_decisions_organs, col_ones_n_obs))

# Preliminary computations
X_beta_c <- X_obs %*% beta_x_c
X_alpha_c <- X_obs %*% alpha_x_c
Q_alpha_c <- Q_obs %*% alpha_q_c
X_gamma_c <- X_obs %*% gamma_x_c
Q_gamma_c <- Q_obs %*% gamma_q_c
# Z_gamma_c <- Z_obs %*% gamma_z_c
f_1_ext_c <- f_1_c[indices_decisions_patients]
f_2_ext_c <- f_2_c[indices_decisions_patients]
f_4_ext_c <- f_4_c[indices_decisions_organs]

# Construct RHS matrices
RHS_unassigned_c <- cbind(X_unassigned,
                          f_1_ext_c[indicator_unassigned == 1],
                          f_2_ext_c[indicator_unassigned == 1])
RHS_assigned_c <- cbind(X_assigned,
                        Q_assigned,
                        f_1_ext_c[indicator_assigned == 1],
                        f_2_ext_c[indicator_assigned == 1],
                        f_3_c[indicator_assigned == 1],
                        f_4_ext_c[indicator_assigned == 1])
# RHS_decision_c <- cbind(X_obs, Q_obs, Z_obs, f_1_ext_c, f_3_c, f_4_ext_c)
RHS_decision_c <- cbind(X_obs, Q_obs, f_1_ext_c, f_3_c, f_4_ext_c)

# Fix column indices for unassigned matrix
RHS_unassigned_col_f1 <- ncol(RHS_unassigned_c) - 1
RHS_unassigned_col_f2 <- RHS_unassigned_col_f1 + 1

# Fix column indices for assigned matrix
RHS_assigned_col_f1 <- ncol(RHS_assigned_c) - 3
RHS_assigned_col_f2 <- RHS_assigned_col_f1 + 1
RHS_assigned_col_f3 <- RHS_assigned_col_f1 + 2
RHS_assigned_col_f4 <- RHS_assigned_col_f1 + 3

# Fix column indices for decision matrix
RHS_decision_col_f1 <- ncol(RHS_decision_c) - 2
RHS_decision_col_f3 <- RHS_decision_col_f1 + 1
RHS_decision_col_f4 <- RHS_decision_col_f1 + 2

# ------------------------------------------------------------------------------
# Chain and Burn Settings
# ------------------------------------------------------------------------------

chain_total <- 100000 # Total chain length
chain_partition <- 10 # Number of loops before keeping it
n_chain <- chain_total / chain_partition # Chain length of loops kept
burn <- 5000 # Number of iterations to be burned
n_parameters <- length(c(beta_x_c,
                         beta_f_1_c,
                         beta_f_2_c,
                         alpha_x_c,
                         alpha_q_c,
                         alpha_f_1_c,
                         alpha_f_2_c,
                         alpha_f_3_c,
                         alpha_f_4_c,
                         gamma_x_c,
                         gamma_q_c,
                         sigma_squared_1_c,
                         sigma_squared_2_c,
                         sigma_squared_4_c,
                         sigma_squared_tilde_0_c,
                         sigma_squared_tilde_1_c))
chain <- matrix(0, nrow = n_chain, ncol = n_parameters) # Define chain

# Collect chain settings
chain_strings <- c("chain_total", "chain_partition", "chain")

################################################################################
# Gibbs Sampler: Algorithm
################################################################################

gibbs_sampler <- function() {
  
  # Gibbs-Sampler
  for (n_c in 1:chain_total) {
    
    # --------------------------------------------------------------------------
    # Step 1 Part (a): Augment latent decision outcomes
    # --------------------------------------------------------------------------
    
    # Compute mean of decision outcomes
    mean_decision <- (X_gamma_c + Q_gamma_c + f_1_ext_c + f_3_c + f_4_ext_c)
    
    # Draw latent decision outcomes conditional on observed decisions
    y_decision_c[y_decision_accept] <- rtruncnorm(n_accept, a = 0, b = Inf,
                                                  mean = mean_decision[y_decision_accept], sd = 1)
    y_decision_c[y_decision_reject] <- rtruncnorm(n_reject, a = -Inf, b = 0,
                                                  mean = mean_decision[y_decision_reject], sd = 1)
    # Print that step 1 is done
    print(sprintf("On loop %i, step 1 for latent decision outcomes is done.", n_c))
    
    # --------------------------------------------------------------------------
    # Step 1 Part (b): Augment censored unassigned outcomes
    # --------------------------------------------------------------------------
    
    # Compute mean of unassigned outcomes
    mean_unassigned <- (X_beta_c + beta_f_1_c * f_1_ext_c + beta_f_2_c * f_2_ext_c)
    
    # Draw latent assigned outcomes over censored observations
    y_unassigned_full_c[indicator_censor_unassigned] <- rtruncnorm(n_censor_unassigned, a = survival_outcomes[indicator_censor_unassigned], b = Inf,
                                                                   mean = mean_unassigned[indicator_censor_unassigned], sd = sqrt(sigma_squared_tilde_0_c))
    y_unassigned_c <- y_unassigned_full_c[indicator_unassigned == 1]
    
    # Print that step 1 is done
    print(sprintf("On loop %i, step 1 for censored unassigned outcomes is done.", n_c))
    
    # --------------------------------------------------------------------------
    # Step 1 Part (c): Augment censored assigned outcomes
    # --------------------------------------------------------------------------
    
    # Compute mean of unassigned outcomes
    mean_assigned <- (X_alpha_c + Q_alpha_c + f_1_ext_c * alpha_f_1_c
                      + f_2_ext_c + f_3_c * alpha_f_3_c + f_4_ext_c * alpha_f_4_c)
    
    # Draw latent assigned outcomes over censored observations
    y_assigned_full_c[indicator_censor_assigned] <- rtruncnorm(n_censor_assigned, a = survival_outcomes[indicator_censor_assigned], b = Inf,
                                                               mean = mean_assigned[indicator_censor_assigned], sd = sqrt(sigma_squared_tilde_1_c))
    y_assigned_c <- y_assigned_full_c[indicator_assigned == 1]
    
    # Print that step 1 is done
    print(sprintf("On loop %i, step 1 for censored assigned outcomes is done.", n_c))
    
    # --------------------------------------------------------------------------
    # Step 2: Draw posterior coefficients for unassigned outcomes
    # --------------------------------------------------------------------------
    
    # Define iteration of RHS variables for unassigned outcomes
    RHS_unassigned_c[, RHS_unassigned_col_f1] <- f_1_ext_c[indicator_unassigned == 1]
    RHS_unassigned_c[, RHS_unassigned_col_f2] <- f_2_ext_c[indicator_unassigned == 1]
    
    # Draw posterior coefficients
    beta_c <- draw_coef_rc(RHS_unassigned_c,
                           y_unassigned_c,
                           prior_mat_beta,
                           sigma_squared_tilde_0_c)
    
    # Define iteration of posterior coefficients
    beta_x_c <- beta_c[1:(k_patients + 1)]
    beta_f_1_c <- beta_c[k_patients + 2]
    beta_f_2_c <- beta_c[k_patients + 3]
    
    # Compute matrix multiplication
    X_beta_c <- X_obs %*% beta_x_c
    
    # Print that step 2 is done for beta coefficients
    print(sprintf("On loop %i, step 2 for unassigned coefficients is done.", n_c))
    
    # --------------------------------------------------------------------------
    # Step 2: Draw posterior coefficients for assigned outcomes
    # --------------------------------------------------------------------------
    
    # Define iteration of RHS variables for unassigned outcomes
    RHS_assigned_c[, RHS_assigned_col_f1] <- f_1_ext_c[indicator_assigned == 1]
    RHS_assigned_c[, RHS_assigned_col_f3] <- f_3_c[indicator_assigned == 1]
    RHS_assigned_c[, RHS_assigned_col_f4] <- f_4_ext_c[indicator_assigned == 1]
    
    # Draw posterior coefficients
    alpha_c <- draw_coef_rc(RHS_assigned_c,
                            y_assigned_c,
                            prior_mat_alpha,
                            sigma_squared_tilde_1_c)
    
    # Define iteration of posterior coefficients
    alpha_x_c <- alpha_c[1:(k_patients + 1)]
    alpha_q_c <- alpha_c[(k_patients + 2):(k_patients + 1 + k_organs)]
    alpha_f_1_c <- alpha_c[RHS_assigned_col_f1]
    alpha_c[RHS_assigned_col_f2] <- 1
    alpha_f_2_c <- 1
    alpha_f_3_c <- alpha_c[RHS_assigned_col_f3]
    alpha_f_4_c <- alpha_c[RHS_assigned_col_f4]
    
    # Compute matrix multiplication
    X_alpha_c <- X_obs %*% alpha_x_c
    Q_alpha_c <- Q_obs %*% alpha_q_c
    
    # Print that step 2 is done for alpha coefficients
    print(sprintf("On loop %i, step 2 for assigned coefficients is done.", n_c))
    
    # --------------------------------------------------------------------------
    # Step 2: Draw posterior coefficients for decision outcomes
    # --------------------------------------------------------------------------
    
    # Define iteration of RHS variables for decision outcomes
    RHS_decision_c[, RHS_decision_col_f1] <- f_1_ext_c
    RHS_decision_c[, RHS_decision_col_f3] <- f_3_c
    RHS_decision_c[, RHS_decision_col_f4] <- f_4_ext_c
    
    # Draw posterior coefficients
    gamma_c <- draw_coef_rc(RHS_decision_c, y_decision_c, prior_mat_gamma, 1)
    
    # Define iteration of posterior coefficients
    gamma_x_c <- gamma_c[1:(k_patients + 1)]
    gamma_q_c <- gamma_c[(k_patients + 2):(k_patients + 1 + k_organs)]
    # gamma_z_c <- gamma_c[(k_patients + 2 + k_organs):(k_patients + 1 + k_organs + k_instruments)]
    gamma_c[(k_patients + 2 + k_organs + k_instruments)] <- 1 # Coefficient of f_1 on decision outcome is normalized to 1
    gamma_c[(k_patients + 3 + k_organs + k_instruments)] <- 1 # Coefficient of f_3 on decision outcome is normalized to 1
    gamma_c[(k_patients + 4 + k_organs + k_instruments)] <- 1 # Coefficient of f_4 on decision outcome is normalized to 1
    
    # Compute matrix multiplications
    X_gamma_c <- X_obs %*% gamma_x_c
    Q_gamma_c <- Q_obs %*% gamma_q_c
    # Z_gamma_c <- Z_obs %*% gamma_z_c
    
    # Print that step 2 is done for gamma coefficients
    print(sprintf("On loop %i, step 2 for decision coefficients is done.", n_c))
    
    # --------------------------------------------------------------------------
    # Step 3: Draw posterior variances for unassigned outcomes
    # --------------------------------------------------------------------------
    
    # Draw posterior variance for unassigned outcomes
    error_unassigned_c <- (y_unassigned_c - RHS_unassigned_c %*% beta_c)
    error_unassigned_squared <- crossprod_armadillo_single(error_unassigned_c)
    posterior_scale_unassigned <- prior_scale + (error_unassigned_squared / 2)
    sigma_squared_tilde_0_c <- rinvgamma(1,
                                         shape = posterior_shape_unassigned,
                                         scale = posterior_scale_unassigned)
    
    # --------------------------------------------------------------------------
    # Step 3: Draw posterior variances for assigned outcomes
    # --------------------------------------------------------------------------
    
    # Draw posterior variance for assigned outcomes
    error_assigned_c <- (y_assigned_c - RHS_assigned_c %*% alpha_c)
    error_assigned_squared <- crossprod_armadillo_single(error_assigned_c)
    posterior_scale_assigned <- prior_scale + (error_assigned_squared / 2)
    sigma_squared_tilde_1_c <- rinvgamma(1,
                                         shape = posterior_shape_assigned,
                                         scale = posterior_scale_assigned)
    
    # Print that step 3 is done
    print(sprintf("On loop %i, step 3 is done.", n_c))
    
    # --------------------------------------------------------------------------
    # Step 4: Draw posterior factors f_1
    # --------------------------------------------------------------------------
    
    # Compute precision weights
    f_1_weight <- (1 / sigma_squared_1_c
                   + (weight_unassigned_patients * beta_f_1_c^2 / sigma_squared_tilde_0_c)
                   + (weight_assigned_patients * alpha_f_1_c^2 / sigma_squared_tilde_1_c)
                   + weight_decision_patients)
    
    # Residual error with respect to unassigned outcome
    f_1_unassigned <- y_unassigned_full_c - (X_beta_c + f_2_ext_c * beta_f_2_c) * indicator_unassigned
    f_1_unassigned <- f_1_unassigned * beta_f_1_c / sigma_squared_tilde_0_c
    f_1_unassigned <- as.numeric(crossprod(obs_decisions_patients, f_1_unassigned))
    
    # Residual error with respect to assigned outcome
    f_1_assigned <- (y_assigned_full_c -
                       (X_alpha_c + Q_alpha_c + f_2_ext_c + f_3_c * alpha_f_3_c + f_4_ext_c * alpha_f_4_c) * indicator_assigned)
    f_1_assigned <- f_1_assigned * alpha_f_1_c / sigma_squared_tilde_1_c
    f_1_assigned <- as.numeric(crossprod(obs_decisions_patients, f_1_assigned))
    
    # Residual error with respect to decision outcome
    f_1_decision <- (y_decision_c -
                       (X_gamma_c + Q_gamma_c + f_3_c + f_4_ext_c))
    f_1_decision <- as.numeric(crossprod(obs_decisions_patients, f_1_decision))
    
    # Compute posterior mean and standard deviations
    f_1_posterior_mean <- (f_1_unassigned + f_1_assigned + f_1_decision) / f_1_weight
    f_1_posterior_sd <- 1 / sqrt(f_1_weight)
    
    # Draw posterior factors
    f_1_c <- rnorm(n_patients, mean = f_1_posterior_mean, sd = f_1_posterior_sd)
    f_1_ext_c <- f_1_c[indices_decisions_patients]
    
    # Print that step 4 is done for factor f_1
    print(sprintf("On loop %i, step 4 for factor f_1 is done.", n_c))
    
    # --------------------------------------------------------------------------
    # Step 4: Draw posterior factor f_2
    # --------------------------------------------------------------------------
    
    # Compute precision weights
    f_2_weight <- (1 / sigma_squared_2_c
                   + (weight_unassigned_patients * beta_f_2_c^2 / sigma_squared_tilde_0_c)
                   + (weight_assigned_patients * alpha_f_2_c^2 / sigma_squared_tilde_1_c))
    
    # Residual error with respect to unassigned outcome
    f_2_unassigned <- y_unassigned_full_c - (X_beta_c + f_1_ext_c * beta_f_1_c) * indicator_unassigned
    f_2_unassigned <- f_2_unassigned * beta_f_2_c / sigma_squared_tilde_0_c
    f_2_unassigned <- as.numeric(crossprod(obs_decisions_patients, f_2_unassigned))
    
    # Residual error with respect to assigned outcome
    f_2_assigned <- (y_assigned_full_c -
                       (X_alpha_c + Q_alpha_c + f_1_ext_c * alpha_f_1_c + f_3_c * alpha_f_3_c + f_4_ext_c * alpha_f_4_c) * indicator_assigned)
    f_2_assigned <- f_2_assigned * alpha_f_2_c / sigma_squared_tilde_1_c
    f_2_assigned <- as.numeric(crossprod(obs_decisions_patients, f_2_assigned))
    
    # Compute posterior mean and standard deviations
    f_2_posterior_mean <- (f_2_unassigned + f_2_assigned) / f_2_weight
    f_2_posterior_sd <- 1 / sqrt(f_2_weight)
    
    # Draw posterior factors
    f_2_c <- rnorm(n_patients, mean = f_2_posterior_mean, sd = f_2_posterior_sd)
    f_2_ext_c <- f_2_c[indices_decisions_patients]
    RHS_assigned_c[, RHS_assigned_col_f2] <- f_2_ext_c[indicator_assigned == 1]
    
    # Print that step 4 is done for factor f_2
    print(sprintf("On loop %i, step 4 for factor f_2 is done.", n_c))
    
    # --------------------------------------------------------------------------
    # Step 4: Draw posterior factors f_3
    # --------------------------------------------------------------------------
    
    # Compute precision weights
    f_3_weight <- (1
                   + (indicator_assigned * alpha_f_3_c^2 / sigma_squared_tilde_0_c)
                   + 1)
    
    # Residual error with respect to assigned outcomes
    f_3_assigned <- (y_assigned_full_c -
                       (X_alpha_c + Q_alpha_c + f_1_ext_c * alpha_f_1_c + f_2_ext_c + f_4_ext_c * alpha_f_4_c) * indicator_assigned)
    f_3_assigned <- f_3_assigned * alpha_f_3_c / sigma_squared_tilde_1_c
    
    # Residual error with respect to decision outcomes
    f_3_decision <- (y_decision_c -
                       (X_gamma_c + Q_gamma_c + f_1_ext_c + f_4_ext_c))
    
    # Compute posterior mean and standard deviations
    f_3_posterior_mean <- (f_3_assigned + f_3_decision) / f_3_weight
    f_3_posterior_sd <- 1 / sqrt(f_3_weight)
    
    # Draw posterior factors
    f_3_c <- rnorm(n_obs, mean = f_3_posterior_mean, sd = f_3_posterior_sd)
    
    # Print that step 4 is done for factor f_3
    print(sprintf("On loop %i, step 4 for factor f_3 is done.", n_c))
    
    # --------------------------------------------------------------------------
    # Step 4: Draw posterior factors f_4
    # --------------------------------------------------------------------------
    
    # Compute precision weights directly as a vector
    f_4_weight <- (1 / sigma_squared_4_c
                   + weight_decision_organs
                   + (weight_assigned_organs * alpha_f_4_c^2 / sigma_squared_tilde_1_c))
    
    # Compute the assigned residual.
    f_4_assigned <- (y_assigned_full_c -
                       (X_alpha_c + Q_alpha_c + f_1_ext_c * alpha_f_1_c + f_2_ext_c + f_3_c * alpha_f_3_c) * indicator_assigned)
    f_4_assigned <- (alpha_f_4_c / sigma_squared_tilde_1_c) * f_4_assigned
    f_4_assigned <- as.numeric(crossprod(obs_decisions_organs, f_4_assigned))
    
    # Compute the decision residual
    f_4_decision <- (y_decision_c -
                       (X_gamma_c + Q_gamma_c + f_1_ext_c + f_3_c))
    f_4_decision <- as.numeric(crossprod(obs_decisions_organs, f_4_decision))
    
    # Compute the posterior mean and standard deviations
    f_4_posterior_mean <- (f_4_assigned + f_4_decision) / f_4_weight
    f_4_posterior_sd <- 1 / sqrt(f_4_weight)
    
    # Draw samples from the posterior distribution
    f_4_c <- rnorm(n_organs, mean = f_4_posterior_mean, sd = f_4_posterior_sd)
    f_4_ext_c <- f_4_c[indices_decisions_organs]
    
    # Print that step 4 is done for factor f_4
    print(sprintf("On loop %i, step 4 for factor f_4 is done.", n_c))
    
    # --------------------------------------------------------------------------
    # Step 5: Draw posterior variances of factors
    # --------------------------------------------------------------------------
    
    # Draw posterior variance of factor 1
    error_f_1_squared <- crossprod(f_1_c)
    posterior_scale_f_1 <- prior_scale + (error_f_1_squared / 2)
    sigma_squared_1_c <- rinvgamma(1,
                                   shape = posterior_shape_f_1,
                                   scale = posterior_scale_f_1)
    
    # Draw posterior variance of factor 2
    error_f_2_squared <- crossprod(f_2_c)
    posterior_scale_f_2 <- prior_scale + (error_f_2_squared / 2)
    sigma_squared_2_c <- rinvgamma(1,
                                   shape = posterior_shape_f_2,
                                   scale = posterior_scale_f_2)
    
    # Draw posterior variance of factor 4
    error_f_4_squared <- crossprod(f_4_c)
    posterior_scale_f_4 <- prior_scale + (error_f_4_squared / 2)
    sigma_squared_4_c <- rinvgamma(1,
                                   shape = posterior_shape_f_4,
                                   scale = posterior_scale_f_4)
    
    # Print that step 5 is done
    print(sprintf("On loop %i, step 5 is done.", n_c))
    
    # --------------------------------------------------------------------------
    # Iterate
    # --------------------------------------------------------------------------
    vector_c <- c(beta_x_c,
                  beta_f_1_c,
                  beta_f_2_c,
                  alpha_x_c,
                  alpha_q_c,
                  alpha_f_1_c,
                  alpha_f_2_c,
                  alpha_f_3_c,
                  alpha_f_4_c,
                  gamma_x_c,
                  gamma_q_c,
                  sigma_squared_1_c,
                  sigma_squared_2_c,
                  sigma_squared_4_c,
                  sigma_squared_tilde_0_c,
                  sigma_squared_tilde_1_c)
    
    # Save iterations
    if (n_c %% chain_partition == 0) {
      chain_c <- as.integer(n_c / chain_partition)
      chain[chain_c, ] <- vector_c
    }
  }
  
  # ----------------------------------------------------------------------------
  # Output
  # ----------------------------------------------------------------------------
  if (n_c == chain_total) {
    print("The Gibbs Sampler is done!")
  }
  return(chain)
}

################################################################################
# True Underlying Parameters
################################################################################

# Collect names of parameters
parameters <- c(rep("beta_x", k_patients + 1),
                "beta_f_1",
                "beta_f_2",
                rep("alpha_x", k_patients + 1),
                rep("alpha_q", k_organs),
                "alpha_f_1",
                "alpha_f_2",
                "alpha_f_3",
                "alpha_f_4",
                rep("gamma_x", k_patients + 1),
                rep("gamma_q", k_organs),
                "sigma_squared_1",
                "sigma_squared_2",
                "sigma_squared_4",
                "sigma_squared_tilde_0",
                "sigma_squared_tilde_1")

################################################################################
# Script
################################################################################

if (interactive()) {
  
  # Set seed for replication
  set.seed(42069)
  
  # Run Gibbs Sampler
  start_time <- Sys.time()
  chain <- gibbs_sampler()
  end_time <- Sys.time()
  execution_time <- end_time - start_time
  
  # Summarize parameters
  estimates_mean <- colMeans(chain[burn:n_chain, ])
  estimates_sd <- rowSd(t(chain[burn:n_chain, ]))
  chain_estimates <- cbind(parameters,
                           as.numeric(estimates_mean),
                           as.numeric(estimates_sd))
  print(chain_estimates)
  print(execution_time)
  
  # Save Output
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  file_path <- paste0("output/chain_nfe_niv_", timestamp, ".RData")
  save(chain, chain_estimates, execution_time, file = file_path)
  
}

# Plot chain
if (!interactive()) {
  
  # Plot function
  plot_chain <- function(est_index) {
    x_axis <- 1:(dim(chain)[1])
    y_axis <- chain[x_axis, est_index]
    title <- paste0("Chain for Posterior ", parameters[est_index])
    plot(x_axis, y_axis, main = title)
  }
  
  # Plot output
  plot_chain(40)
  
}

################################################################################
