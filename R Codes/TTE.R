
#STEP 1 ==================

library(TrialEmulation)

#STEP 2 ==================

# Load the dataset
data("data_censored")

# Save as CSV
write.csv(data_censored, "data_censored.csv", row.names = FALSE)

# View first few rows to confirm
head(data_censored)

# Define estimation methods
trial_pp  <- trial_sequence(estimand = "PP")  # Per-protocol
trial_itt <- trial_sequence(estimand = "ITT") # Intention-to-treat

# Create directories for saving models

# Assign observational dataset
trial_pp <- set_data(
  trial_pp,
  data      = data_censored,
  id        = "id",
  period    = "period",
  treatment = "treatment",
  outcome   = "outcome",
  eligible  = "eligible"
)

trial_itt <- set_data(
  trial_itt,
  data      = data_censored,
  id        = "id",
  period    = "period",
  treatment = "treatment",
  outcome   = "outcome",
  eligible  = "eligible"
)

# Print trial_itt summary
trial_itt

#STEP 3.1 ==================

# Treatment switching censoring (PP only)
trial_pp <- trial_pp |>
  set_switch_weight_model(
    numerator    = ~ age,
    denominator  = ~ age + x1 + x3,
    model_fitter = stats_glm_logit(save_path = file.path(trial_pp_dir, "switch_models"))
  )
trial_pp@switch_weights

#STEP 3.2 ==================

# Informative censoring (for both PP and ITT)
trial_pp <- trial_pp |>
  set_censor_weight_model(
    censor_event = "censored",
    numerator    = ~ x2,
    denominator  = ~ x2 + x1,
    pool_models  = "none",
    model_fitter = stats_glm_logit(save_path = file.path(trial_pp_dir, "switch_models"))
  )
trial_pp@censor_weights

trial_itt <- set_censor_weight_model(
  trial_itt,
  censor_event = "censored",
  numerator    = ~x2,
  denominator  = ~ x2 + x1,
  pool_models  = "numerator",
  model_fitter = stats_glm_logit(save_path = file.path(trial_itt_dir, "switch_models"))
)
trial_itt@censor_weights

#STEP 4 ==================

trial_pp  <- trial_pp |> calculate_weights()
trial_itt <- calculate_weights(trial_itt)

# Print weight models
show_weight_models(trial_itt)
show_weight_models(trial_pp)

#STEP 5 ==================

trial_pp  <- set_outcome_model(trial_pp)
trial_itt <- set_outcome_model(trial_itt, adjustment_terms = ~x2)

#STEP 6 ==================

trial_pp <- set_expansion_options(
  trial_pp,
  output     = save_to_datatable(),
  chunk_size = 500 # the number of patients to include in each expansion iteration
)
trial_itt <- set_expansion_options(
  trial_itt,
  output     = save_to_datatable(),
  chunk_size = 500
)

#STEP 6.1 ==================

trial_pp  <- expand_trials(trial_pp)
trial_itt <- expand_trials(trial_itt)

trial_pp@expansion

#STEP 7 ==================

trial_itt <- load_expanded_data(trial_itt, seed = 1234, p_control = 0.5)

#STEP 8 ==================

trial_itt <- fit_msm(
  trial_itt,
  weight_cols    = c("weight", "sample_weight"),
  modify_weights = function(w) { # winsorization of extreme weights
    q99 <- quantile(w, probs = 0.99)
    pmin(w, q99)
  }
)

trial_itt@outcome_model


trial_itt@outcome_model@fitted@model$model

trial_itt@outcome_model@fitted@model$vcov

trial_itt



#STEP 9 ==================

preds <- predict(
  trial_itt,
  newdata       = outcome_data(trial_itt)[trial_period == 1, ],
  predict_times = 0:10,
  type          = "survival"
)

plot(preds$difference$followup_time, preds$difference$survival_diff,
     type = "l", xlab = "Follow up", ylab = "Survival difference")
lines(preds$difference$followup_time, preds$difference$`2.5%`, type = "l", col = "red", lty = 2)
lines(preds$difference$followup_time, preds$difference$`97.5%`, type = "l", col = "red", lty = 2)
