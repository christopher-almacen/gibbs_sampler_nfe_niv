################################################################################
# Author(s)
################################################################################

# Chris Almacen

################################################################################
# Goal
################################################################################

# The goal of this script is to transform the data into an input file for the
# Gibbs Sampler. This input does not include fixed effects. Specifically, does
# not include year registration or listing center. This also does not include
# instrumental variables

################################################################################
# Personal Settings
################################################################################

# Clear variables from environment
rm(list = ls())

# Set (personal) working directory
setwd("~/Arizona_State_University/Research Assistant/gibbs_sampler_2025_08_25")

################################################################################
# Libraries
################################################################################

library(lubridate)
library(haven)

################################################################################
# Script
################################################################################

# Function to convert a date string ("DDmonYYYY") to a Date object.
convert_to_date <- function(date_str) {
  # Extract parts of the string:
  # Day: characters 1-2, Month: characters 3-5, Year: characters 6-9
  day <- substr(date_str, 1, 2)
  mon <- substr(date_str, 3, 5)
  yr  <- substr(date_str, 6, 9)
  
  # Convert month string to a format R can recognize.
  # R's %b specifier usually expects the first letter uppercase (e.g., "Jan").
  mon <- paste0(toupper(substr(mon, 1, 1)), tolower(substr(mon, 2, 3)))
  
  # Create a new date string with a separator (e.g., "01-Jan-2000")
  new_date_str <- paste(day, mon, yr, sep = "-")
  
  # Convert to Date using the appropriate format string
  as.Date(new_date_str, format = "%d-%b-%Y")
}

################################################################################
# Script
################################################################################

# Read csv file
data_full <- read.csv("data/PTR_DRAFT_testing_2019_2020.csv")

# Dimension and indices of observed patient-donor pairs
n_obs <- dim(data_full)[1] # Number of patient-donor observations
indices_decisions_patients_full <- data_full[["WL_ID_CODE"]] # Vector of observed patient IDs
indices_decisions_organs_full <- data_full[["DONOR_ID"]] # Vector of observed donor IDs

# Patients that "received a transplant"
indicator_assigned <- rep(0, n_obs)
indicator_assigned_index <- data_full[["TX_DATE"]] != ""
indicator_assigned[indicator_assigned_index] <- as.integer(!duplicated(indices_decisions_patients_full[indicator_assigned_index]))

# Patients that "did not receive a transplant"
indicator_unassigned <- as.numeric(!(indices_decisions_patients_full %in% indices_decisions_patients_full[indicator_assigned == 1]) & !duplicated(indices_decisions_patients_full))

# Observed binary decisions indicating acceptance
binary_decisions <- as.integer(data_full[["offer_accept"]] == "Y")

# Survival outcomes censored by ~ 30 Dec 2020
indicator_censor <- rep(1, n_obs)
indicator_censor[data_full[["COMPOSITE_DEATH_DATE"]] != ""] <- 0 # 0 corresponds to not censored

# Censor date
censor_date <- "30dec2020"

# ------------------------------------------------------------------------------
# Survival Outcomes for patients who "received transplants"
# ------------------------------------------------------------------------------

# Matrix of survival outcomes
survival_outcomes_mat <- cbind(data_full[["WL_ID_CODE"]],
                               data_full[["DONOR_ID"]],
                               data_full[["TX_DATE"]],
                               data_full[["COMPOSITE_DEATH_DATE"]])
colnames(survival_outcomes_mat) <- c("WL_ID_CODE",
                                     "DONOR_ID",
                                     "TX_DATE",
                                     "COMPOSITE_DEATH_DATE")

# Restrict to patients with transfers (conditioning by who accepted)
survival_outcomes_tx_mat <- survival_outcomes_mat[data_full[["TX_DATE"]] != "", ]

# Restrict to unique IDs
survival_outcomes_tx_mat <- unique(survival_outcomes_tx_mat)

# Input censor date
survival_outcomes_tx_mat[survival_outcomes_tx_mat[, 4] == "", 4] <- censor_date

# Convert date strings into date objects
tx_date_convert <- as.Date(unlist(lapply(survival_outcomes_tx_mat[, 3], convert_to_date)))
composite_death_date_tx_convert <- as.Date(unlist(lapply(survival_outcomes_tx_mat[, 4], convert_to_date)))

# Obtain survival outcomes for "assigned patients"
survival_outcomes_tx <- time_length(interval(tx_date_convert, composite_death_date_tx_convert), "year")

# Remark: The reason why I condition on which patients accepted is because I see
# empty transfer dates and yet patients said yes. The implicit assumption is that
# if a patient said yes, then I assumed they have received a transplant, but, 
# for whatever reason, we (as the econometrician) do not observe the transfer
# date.

# ------------------------------------------------------------------------------
# Survival Outcomes for patients who "did not receive transplants"
# ------------------------------------------------------------------------------

# Obtain unique IDs for patients that "received a transplant"
WL_ID_CODE_tx <- survival_outcomes_tx_mat[, 1]

# Construct matrix of survival outcomes for patients that "did not receive a transplant"
survival_outcomes_no_tx_mat <- cbind(data_full[["WL_ID_CODE"]],
                                     data_full[["INIT_DATE"]],
                                     data_full[["END_DATE"]],
                                     data_full[["COMPOSITE_DEATH_DATE"]])
colnames(survival_outcomes_no_tx_mat) <- c("WL_ID_CODE",
                                           "INIT_DATE",
                                           "END_DATE",
                                           "COMPOSITE_DEATH_DATE")

# Restrict to patient IDs that did not "receive a transplant"
survival_outcomes_no_tx_mat <- survival_outcomes_no_tx_mat[!(survival_outcomes_no_tx_mat[, 1] %in% WL_ID_CODE_tx), ]

# Unique patient IDs that did not "receive a transplant"
survival_outcomes_no_tx_mat <- unique(survival_outcomes_no_tx_mat)

# Convert date strings into date objects for INIT_DATE
init_date_convert <- as.Date(unlist(lapply(survival_outcomes_no_tx_mat[, 2], convert_to_date)))

# Convert date strings into date objects for upper censor (by END_DATE or COMPOSITIE_DEATH_DATE or upper_censor)
end_date_convert <- as.Date(unlist(lapply(survival_outcomes_no_tx_mat[, 3], convert_to_date)))
composite_death_date_no_tx_convert <- as.Date(unlist(lapply(survival_outcomes_no_tx_mat[, 4], convert_to_date)))
upper_date_no_tx <- end_date_convert
upper_date_no_tx[!(is.na(composite_death_date_no_tx_convert))] <- composite_death_date_no_tx_convert[!(is.na(composite_death_date_no_tx_convert))]

# Obtain survival outcomes for "unassigned patients"
survival_outcomes_no_tx <- time_length(interval(init_date_convert, upper_date_no_tx), "year")
survival_outcomes_no_tx[survival_outcomes_no_tx < 0] <- NA

# ------------------------------------------------------------------------------
# Survival Outcomes Full
# ------------------------------------------------------------------------------

# Construct and combine look up tables for outcomes for tx and no tx
lookup <- rbind(cbind(survival_outcomes_tx_mat[, 1], survival_outcomes_tx),
                cbind(survival_outcomes_no_tx_mat[, 1], survival_outcomes_no_tx))

# Turn column-2 into a named vector keyed by the IDs
outcome_lookup <- setNames(lookup[, 2], lookup[, 1])

# Pull the matching outcomes for every element of patient IDs
survival_outcomes <- as.numeric(outcome_lookup[as.character(indices_decisions_patients_full)])

# ------------------------------------------------------------------------------
# Patient Characteristics
# ------------------------------------------------------------------------------

# Intercept
intercept <- rep(1, n_obs)

# Diabetic
diab_y <- data_full[["diab_y"]]

# CPRA
INIT_CPRA <- as.numeric(data_full[["INIT_CPRA"]]) # CPRA
INIT_CPRA_gt80 <- as.numeric(INIT_CPRA >= 80) # Indicator if CPRA >= 80
INIT_CPRA_0 <- as.numeric(INIT_CPRA == 0) # Indicator if CPRA == 0
INIT_CPRA_g_80 <- pmax(INIT_CPRA - 80, 0) # CPRA - 0.8 if  CPRA >= 80

# Prior Transplant
prior_transplant <- rep(0, n_obs)
prior_transplant[data_full[["NUM_PREV_TX"]] != 0] <- 1 # Indicator if patient had prior transplant

# Dialysis
dial_y <- as.numeric(data_full[["dial_y"]]) # Indicator if patient is on dialysis at initial date

# Blood Types
ABO_1 <- data_full[["ABO_1"]] # Indicator for blood type A
ABO_7 <- data_full[["ABO_7"]] # Indicator for blood type B
ABO_8 <- data_full[["ABO_8"]] # Indicator for blood type O

# Age at registration
INIT_AGE <- as.numeric(data_full[["INIT_AGE"]])
INIT_AGE_gt_18 <- pmax(0, INIT_AGE - 18) # Age - 18 if Age >= 18
INIT_AGE_gt_35 <- pmax(0, INIT_AGE - 35) # Age - 35 if Age >= 35
INIT_AGE_gt_50 <- pmax(0, INIT_AGE - 50) # Age - 50 if Age >= 50
INIT_AGE_gt_65 <- pmax(0, INIT_AGE - 65) # Age - 65 if Age >= 65

# BMI at departure
END_BMI_CALC <- as.numeric(data_full[["END_BMI_CALC"]])
END_BMI_CALC_gt_18_5 <- pmax(0, END_BMI_CALC - 18.5) # BMI - 18.5 if BMI >= 18.5
END_BMI_CALC_gt_25 <- pmax(0, END_BMI_CALC - 25) # BMI - 25 if BMI >= 25
END_BMI_CALC_gt_30 <- pmax(0, END_BMI_CALC - 30) # BMI - 30 if BMI >= 30

# Serum Albumin
TOT_SERUM_ALBUM <- as.numeric(data_full[["TOT_SERUM_ALBUM"]])
TOT_SERUM_ALBUM_gt_3_7 <- pmax(0, TOT_SERUM_ALBUM - 3.7) # Serum Albumin - 3.7 if Serum Albumin >= 3.7
TOT_SERUM_ALBUM_gt_4_4 <- pmax(0, TOT_SERUM_ALBUM - 4.4) # Serum Albumin - 4.4 if Serum Albumin >= 4.4

# Dialysis time at registration
dial_time_at_reg <- as.numeric(data_full[["dial_time_at_reg"]])

# Log dialysis time at registration
log_dial_time_at_reg <- dial_time_at_reg
log_dial_time_at_reg[dial_time_at_reg <= 0] <- 0
log_dial_time_at_reg <- log_dial_time_at_reg + 1
log_dial_time_at_reg <- log(log_dial_time_at_reg)

# Log dialysis time at registration greater than 5 years
log_dial_time_at_reg_gt_1_6094 <- as.numeric(log_dial_time_at_reg >= 1.6094) 

# Initial CPRA missing
if (sum(is.na(data_full[["INIT_CPRA"]])) > 0) {
  INIT_CPRA_dum <- as.integer(is.na(INIT_CPRA))
  INIT_CPRA[is.na(data_full[["INIT_CPRA"]])] <- 0
}

# BMI missing
if (sum(is.na(data_full[["END_BMI_CALC"]])) > 0) {
  END_BMI_CALC_dum <- as.integer(is.na(END_BMI_CALC))
  END_BMI_CALC[is.na(data_full[["END_BMI_CALC"]])] <- 0
}

# Serum albumin missing
if (sum(is.na(data_full[["TOT_SERUM_ALBUM"]])) > 0) {
  TOT_SERUM_ALBUM_dum <- as.integer(is.na(TOT_SERUM_ALBUM))
  TOT_SERUM_ALBUM[is.na(data_full[["TOT_SERUM_ALBUM"]])] <- 0
}

# Collect patient characteristics
if (sum(is.na(data_full[["INIT_CPRA"]])) == 0) {
  X_full <- cbind(intercept,
                  diab_y,
                  INIT_CPRA,
                  INIT_CPRA_gt80,
                  INIT_CPRA_0,
                  INIT_CPRA_g_80,
                  prior_transplant,
                  dial_y,
                  ABO_1,
                  ABO_7,
                  ABO_8,
                  INIT_AGE,
                  INIT_AGE_gt_18,
                  INIT_AGE_gt_35,
                  INIT_AGE_gt_50,
                  INIT_AGE_gt_65,
                  END_BMI_CALC,
                  END_BMI_CALC_gt_18_5,
                  END_BMI_CALC_gt_25,
                  END_BMI_CALC_gt_30,
                  TOT_SERUM_ALBUM,
                  TOT_SERUM_ALBUM_gt_3_7,
                  TOT_SERUM_ALBUM_gt_4_4,
                  log_dial_time_at_reg,
                  log_dial_time_at_reg_gt_1_6094,
                  END_BMI_CALC_dum,
                  TOT_SERUM_ALBUM_dum)
}

if (sum(is.na(data_full[["INIT_CPRA"]])) > 0) {
  X_full <- cbind(intercept,
                  diab_y,
                  INIT_CPRA,
                  INIT_CPRA_gt80,
                  INIT_CPRA_0,
                  INIT_CPRA_g_80,
                  prior_transplant,
                  dial_y,
                  ABO_1,
                  ABO_7,
                  ABO_8,
                  INIT_AGE,
                  INIT_AGE_gt_18,
                  INIT_AGE_gt_35,
                  INIT_AGE_gt_50,
                  INIT_AGE_gt_65,
                  END_BMI_CALC,
                  END_BMI_CALC_gt_18_5,
                  END_BMI_CALC_gt_25,
                  END_BMI_CALC_gt_30,
                  TOT_SERUM_ALBUM,
                  TOT_SERUM_ALBUM_gt_3_7,
                  TOT_SERUM_ALBUM_gt_4_4,
                  log_dial_time_at_reg,
                  log_dial_time_at_reg_gt_1_6094,
                  INIT_CPRA_dum,
                  END_BMI_CALC_dum,
                  TOT_SERUM_ALBUM_dum)
}

# ------------------------------------------------------------------------------
# Donor characteristics
# ------------------------------------------------------------------------------

# Donor age
AGE_DON <- data_full[["AGE_DON"]]
AGE_DON_lt_18 <- as.numeric(AGE_DON < 18) # Age < 18
AGE_DON_18_35 <- as.numeric(AGE_DON >= 18 & AGE_DON <= 35) # Age 18-35
AGE_DON_ge_50 <- as.numeric(AGE_DON >= 50) # Age >= 50

# Cause of death: anoxia
cod_anox <- data_full[["cod_anox"]]

# Cause of death: stroke
cod_stroke <- data_full[["cod_stroke"]]

# Cause of death: CNS
cod_cns <- data_full[["cod_cns"]]

# Cause of death: head trauma
cod_ht <- data_full[["cod_ht"]]

# Creatinine
creat_low <- data_full[["creat_low"]]
creat_mid <- data_full[["creat_mid"]]
creat_high <- data_full[["creat_high"]] # let creat_high be the dummy

# Expanded criteria donor (ECD)
ECD_DONOR <- data_full[["ECD_DONOR"]]

# Donation After Cardiac Death (DCD)
NON_HRT_DON_dummy <- data_full[["NON_HRT_DON_dummy"]]

# Male donor
Male_Don <- data_full[["Male_Don"]]

# History of hypertension
hyp_y <- data_full[["hyp_y"]]

# Perfect Tissue Type Match
hla_perfect <- as.integer(data_full[["HLAMIS_calc"]] > 0)

# 2 A Mismatches
A2_HLA_mm <- data_full[["a2_miss"]]

# 2 B Mismatches
B2_HLA_mm <- data_full[["b2_miss"]]

# 2 DR Mismatches
DR2_HLA_mm <- data_full[["dr2_miss"]]

# Regional and local
is_regional <- data_full[["regional"]]
is_local <- data_full[["local_match"]]

# Log wait time years
log_wait_time_years <- data_full[["log_wait_time_years"]]
wait_time_years <- data_full[["wait_time_years"]]
log_wait_time_gt_0 <- ifelse(wait_time_years > 1, log_wait_time_years, 0)
log_wait_time_gt_0_6935 <- ifelse(wait_time_years > 2, log_wait_time_years, 0)

# Collect donor characteristics
Q_full <- cbind(AGE_DON_lt_18,
                AGE_DON_18_35,
                AGE_DON_ge_50,
                cod_anox,
                cod_stroke,
                cod_cns,
                cod_ht,
                creat_low,
                creat_mid,
                ECD_DONOR,
                NON_HRT_DON_dummy,
                Male_Don,
                hyp_y,
                hla_perfect,
                A2_HLA_mm,
                B2_HLA_mm,
                DR2_HLA_mm,
                is_regional,
                is_local,
                log_wait_time_years,
                log_wait_time_gt_0,
                log_wait_time_gt_0_6935)

# ------------------------------------------------------------------------------
# Interactions
# ------------------------------------------------------------------------------

hla_perfect_prior_transplant <- hla_perfect * prior_transplant
hla_perfect_diab_y <- hla_perfect * diab_y
hla_perfect_INIT_AGE <- hla_perfect * INIT_AGE
hla_perfect_INIT_CPRA <- hla_perfect * INIT_CPRA
hla_perfect_INIT_CPRA_gt80 <- hla_perfect * INIT_CPRA_gt80
hla_perfect_ECD_DONOR <- hla_perfect * ECD_DONOR
hla_perfect_NON_HRT_DON_dummy <- hla_perfect * NON_HRT_DON_dummy
hla_perfect_is_local <- hla_perfect * is_local
is_local_A2_HLA_mm <- is_local * A2_HLA_mm
is_local_B2_HLA_mm <- is_local * B2_HLA_mm
is_local_DR2_HLA_mm <- is_local * DR2_HLA_mm
is_local_AGE_DON_lt_18 <- is_local * AGE_DON_lt_18
is_local_AGE_DON_18_35 <- is_local * AGE_DON_18_35
is_local_AGE_DON_ge_50 <- is_local * AGE_DON_ge_50
INIT_AGE_AGE_DON_lt_18 <- INIT_AGE * AGE_DON_lt_18
INIT_AGE_AGE_DON_18_35 <- INIT_AGE * AGE_DON_18_35
INIT_AGE_AGE_DON_ge_50 <- INIT_AGE * AGE_DON_ge_50
INIT_AGE_gt_35_AGE_DON_18_35 <- INIT_AGE_gt_35 * AGE_DON_18_35
INIT_AGE_gt_35_AGE_DON_ge_50 <- INIT_AGE_gt_35 * AGE_DON_ge_50

# Collect interaction terms
interactions_full <- cbind(hla_perfect_prior_transplant,
                           hla_perfect_diab_y,
                           hla_perfect_INIT_AGE,
                           hla_perfect_INIT_CPRA,
                           hla_perfect_INIT_CPRA_gt80,
                           hla_perfect_ECD_DONOR,
                           hla_perfect_NON_HRT_DON_dummy,
                           hla_perfect_is_local,
                           is_local_A2_HLA_mm,
                           is_local_B2_HLA_mm,
                           is_local_DR2_HLA_mm,
                           is_local_AGE_DON_lt_18,
                           is_local_AGE_DON_18_35,
                           is_local_AGE_DON_ge_50,
                           INIT_AGE_AGE_DON_lt_18,
                           INIT_AGE_AGE_DON_18_35,
                           INIT_AGE_AGE_DON_ge_50,
                           INIT_AGE_gt_35_AGE_DON_18_35,
                           INIT_AGE_gt_35_AGE_DON_ge_50)

# ------------------------------------------------------------------------------
# Construct data
# ------------------------------------------------------------------------------

# Rid of incomplete rows
data_all <- cbind(indices_decisions_patients_full,
                  indices_decisions_organs_full,
                  survival_outcomes,
                  binary_decisions,
                  indicator_assigned,
                  indicator_unassigned,
                  indicator_censor,
                  X_full,
                  Q_full,
                  interactions_full)

# Set column names
colnames(data_all) <- c("patient_indices",
                        "donor_indices",
                        "survival_outcomes",
                        "binary_decisions",
                        "indicator_assigned",
                        "indicator_unassigned",
                        "indicator_censor",
                        colnames(X_full),
                        colnames(Q_full),
                        colnames(interactions_full))

# Subset of only complete cases
data_pre <- data_all[complete.cases(data_all), ]

# Relabel indices of observed patient-donor pairs
indices_decisions_patients <- as.integer(factor(data_pre[, "patient_indices"], levels = unique(data_pre[, "patient_indices"])))
indices_decisions_organs <- as.integer(factor(data_pre[, "donor_indices"], levels = unique(data_pre[, "donor_indices"])))

# Final labeling
survival_outcomes <- data_pre[, "survival_outcomes"]
binary_decisions <- data_pre[, "binary_decisions"]
indicator_assigned <- data_pre[, "indicator_assigned"]
indicator_unassigned <- data_pre[, "indicator_unassigned"]
indicator_censor <- data_pre[, "indicator_censor"]

# Final patient characteristics and ensure full rank
X_obs <- data_pre[, colnames(X_full)]
qr_X <- qr(X_obs) # by default uses pivoting
r_X <- qr_X$rank # numerical rank
pivot_X <- qr_X$pivot # permutation of 1:ncol(X)
indep_X <- pivot_X[seq_len(r_X)] # indices of independent columns
dep_X <- pivot_X[-seq_len(r_X)] # indices of dependent columns
if (length(dep_X) > 0) {
  X_obs <- X_obs[, indep_X]
}

# Final donor characteristics and ensure full rank
Q_obs <- data_pre[, colnames(Q_full)]
qr_Q <- qr(Q_obs) # by default uses pivoting
r_Q <- qr_Q$rank # numerical rank
pivot_Q <- qr_Q$pivot # permutation of 1:ncol(Q)
indep_Q <- pivot_Q[seq_len(r_Q)] # indices of independent columns
dep_Q <- pivot_Q[-seq_len(r_Q)] # indices of dependent columns
if (length(dep_Q) > 0) {
  Q_obs <- Q_obs[, indep_Q]
}

# Final interactions and ensure full rank
interactions_obs <- data_pre[, colnames(interactions_full)]
qr_inter <- qr(interactions_obs) # by default uses pivoting
r_inter <- qr_inter$rank # numerical rank
pivot_inter <- qr_inter$pivot # permutation of 1:ncol(Q)
indep_inter <- pivot_inter[seq_len(r_inter)] # indices of independent columns
dep_inter <- pivot_inter[-seq_len(r_inter)] # indices of dependent columns
if (length(dep_inter) > 0) {
  interactions_obs <- interactions_obs[, indep_inter]
}

################################################################################
# Save Data
################################################################################

# Save as RData
if (interactive()) {
  
  save(indices_decisions_patients, file = "input/indices_decisions_patients.RData")
  save(indices_decisions_organs, file = "input/indices_decisions_organs.RData")
  save(survival_outcomes, file = "input/survival_outcomes.RData")
  save(binary_decisions, file = "input/binary_decisions.RData")
  save(indicator_assigned, file = "input/indicator_assigned.RData")
  save(indicator_unassigned, file = "input/indicator_unassigned.RData")
  save(indicator_censor, file = "input/indicator_censor.RData")
  save(X_obs, file = "input/X_obs.RData")
  save(Q_obs, file = "input/Q_obs.RData")
  save(interactions_obs, file = "input/interactions_obs.RData")
  
}

################################################################################
# Sanity Checks
################################################################################

# Save as csv files
if (!interactive()) {
  
  # Output patient and donor indices
  write.table(data_pre[, "patient_indices"], file = "input/patient_indices.csv", row.names = FALSE, col.names = c("patient_indices"))
  write.table(data_pre[, "donor_indices"], file = "input/donor_indices.csv", row.names = FALSE, col.names = c("donor_indices"))
  
  # Output outcome variables
  write.table(data_pre[, "survival_outcomes"], file = "input/survival_outcomes_years.csv", row.names = FALSE, col.names = c("survival_outcomes"))
  write.table(data_pre[, "binary_decisions"], file = "input/binary_decisions.csv", row.names = FALSE, col.names = c("binary_decisions"))
  
  # Output indicators
  write.table(data_pre[, "indicator_assigned"], file = "input/indicator_assigned.csv", row.names = FALSE, col.names = c("indicator_assigned"))
  write.table(data_pre[, "indicator_unassigned"], file = "input/indicator_unassigned.csv", row.names = FALSE, col.names = c("indicator_unassigned"))
  write.table(data_pre[, "indicator_censor"], file = "input/indicator_censor.csv", row.names = FALSE, col.names = c("indicator_censor"))
  
  # Output RHS variables
  write.table(X_pre, file = "input/patient_data.csv", row.names = FALSE)
  write.table(Q_pre, file = "input/donor_data.csv", row.names = FALSE)
  write.table(interactions_pre, file = "input/interactions_data.csv", row.names = FALSE)
  
}

if (!interactive()) {
  
  # Current PRA (to test relative to INIT_CPRA)
  if (!interactive()) {
    current_cpra <- as.numeric(data_full[["CURRENT_PRA"]]) # CPRA
    current_cpra_gt80 <- as.numeric(current_cpra >= 80) # Indicator if CPRA >= 80
    current_cpra_0 <- as.numeric(current_cpra == 0) # Indicator if CPRA == 0
    current_cpra_g_80 <- pmax(current_cpra - 80, 0)
  }
  
  # Check for dependent columns in complete cases of patient characteristics
  X_pre <- data_pre[, colnames(X_full)]
  qr_X <- qr(X_pre) # by default uses pivoting
  r_X <- qr_X$rank # numerical rank
  pivot_X <- qr_X$pivot # permutation of 1:ncol(X)
  indep_X <- pivot_X[seq_len(r_X)] # indices of independent columns
  dep_X <- pivot_X[-seq_len(r_X)] # indices of dependent columns
  dep_X
  
  # Check for dependent columns in complete cases of donor characteristics
  Q <- data_pre[, colnames(Q_full)]
  qr_Q <- qr(Q) # by default uses pivoting
  r_Q <- qr_Q$rank # numerical rank
  pivot_Q <- qr_Q$pivot # permutation of 1:ncol(X)
  indep_Q <- pivot_Q[seq_len(r_Q)] # indices of independent columns
  dep_Q <- pivot_Q[-seq_len(r_Q)] # indices of dependent columns
  dep_Q
  
}
