# packages for this _targets.R file (not for the pipeline)
library(targets)
library(tidyverse)
library(tarchetypes)
library(furrr)
library(here)
library(ggplot2)
library(future)
library(cmdstanr)
library(stantargets)
library(clustermq)
library(lme4)

# Source the user-defined functions from functions.R and figs.R
source("R/functions.R")
source("R/figs.R")
# parallel computing on local or on the same node
plan(multicore)
options(clustermq.scheduler = "multicore")

# Set options for the targets pipeline, specifying required packages
tar_option_set(
  packages = c(
    "tidyverse",
    "lme4",
    "janitor",
    "yaml",
    "e1071",
    "car",
    "sjPlot",
    "gridExtra",
    "loo",
    "knitr",
    "kableExtra",
    "MASS",
    "patchwork",
    "cmdstanr",
    "cowplot",
    "MASS",
    "stringr",
    "grid",
    "gsubfn",
    "scales",
    "maps",
    "writexl",
    "openxlsx"
  )
)

# tar_option_set(
#   garbage_collection = TRUE,
#   memory = "transient"
# )


# mcmc (model) names -------------------------------------
lr_values <- expand_grid(
  lr = "lr",
  model = c("pl", "pl_biome"),
  div = c("ang", "gym"),
  pred = c("cr", "h", "dbh", "dbh1", "dbh2", "dbh3")
)

nlr_values <- expand_grid(
  lr = "nlr",
  model = c("gmm", "weibull", "gmm_biome", "weibull_biome"),
  div = c("ang", "gym"),
  pred = c("cr", "h")
)

values <- bind_rows(lr_values, nlr_values)

data_names <- values |>
  mutate(data_names = str_c(model, div, pred, sep = "_")) |>
  pull(data_names)

mcmc_names <- values |>
  mutate(mcmc_names = str_c("fit_", lr, "_nou_mcmc_",  data_names)) |>
  filter(mcmc_names != "fit_nlr_nou_mcmc_gmm_biome_ang_h") |>
  filter(mcmc_names != "fit_nlr_nou_mcmc_gmm_biome_ang_cr") |>
  pull(mcmc_names)

mcmc_names <- c(mcmc_names,
  "fit_nlr_nou_re_mcmc_weibull_biome_re",
  "fit_wd_gym_dbh_mcmc_pl_wd",
  "fit_wd_gym_cr_mcmc_gmm_wd",
  "fit_wd_gym_h_mcmc_weibull_wd",
  "fit_wd_ang_dbh_mcmc_pl_wd",
  "fit_wd_ang_cr_mcmc_pl_wd",
  "fit_wd_ang_h_mcmc_weibull_wd")

biomes <- c(
      "Boreal.montane.forest",
      "Tropical.rain.forest",
      "Temperate.broadleaf.forest",
      "Tropical.savanna",
      "Temperate.grassland",
      "Temperate.conifer.forest",
      "Mediterranean.woodland",
      "Tropical.dry.forest",
      "Dryland",
      "Mangrove")



# Define the tar_map for LOOIC calculation
loo_map <- tar_map(
  values = list(mcmc = rlang::syms(mcmc_names)),
  tar_target(
    loo,
    my_loo(mcmc)
  )
)



# Define the targets and their dependencies
# Data ------------------------------------------
data_ <- list(
  tar_target(
    tallo_csv,
    "data-raw/Tallo.csv",
    format = "file"
  ),
  tar_target(
    tallo_env_csv,
    "data-raw/Tallo_environment.csv",
    format = "file"
  ),
  tar_target(
    try_wd_txt,
    "data-raw/30953.txt",
    format = "file"
  ),
  # Define the target to prepare overlap.csv data directly
  tar_target(
    tallo_wd_df0,
    clean_tallo_try(tallo_csv, tallo_env_csv, try_wd_txt),
  ),
  NULL
)

# Stan and mapping -----------------------------------------
tar_map_reduced_data <- tar_map(
  values = list(
    variable = c("dbh", "h", "cr")),
  tar_target(
    tallo_wd_df, {
      if (variable == "cr") {
        tmp <- tallo_wd_df0 |>
          filter(!is.na(cr))
      } else if (variable == "h") {
        tmp <- tallo_wd_df0 |>
          filter(!is.na(h))
      } else if (variable == "dbh") {
        tmp <- tallo_wd_df0 |>
          filter(!is.na(dbh))
      }
      reduce_trees(tmp, 100, 50)
    }
  ),
  NULL
)

tar_map_loo_test <-
  tar_map(
    values = list(model = c("xy", "x", "none")),
    tar_target(
      stan_sim_data,
      generate_sim_stan_data(sim_data, model = model)
    ),
    tar_stan_mcmc(
      fit_test,
      c("stan/test.stan", "stan/test2.stan"),
      data = stan_sim_data,
      refresh = 0,
      chains = 3,
      parallel_chains = 3,
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.9,
      max_treedepth = 15,
      seed = 123,
      return_draws = FALSE,
      return_diagnostics = TRUE,
      return_summary = TRUE
   )
  )

# Stan
tar_map_lr_nou <-
  tar_map(
    values = expand_grid(
      div = c("ang", "gym"),
      variable = c("cr", "h", "dbh", "dbh1", "dbh2", "dbh3")),
    tar_target(
      tallo_reduced_lr_df, {
        if (variable == "cr") {
          tmp <- tallo_wd_df0 |>
            filter(!is.na(cr))
        } else if (variable == "h") {
          tmp <- tallo_wd_df0 |>
            filter(!is.na(h))
        } else if (str_detect(variable, "dbh")) {
          tmp <- tallo_wd_df0 |>
            filter(!is.na(dbh))
        }
        reduce_trees(tmp, 100, 50)
      }
    ),
    tar_target(
      stan_data_lr,
      generate_stan_data(tallo_reduced_lr_df, model = "lr", div = div, variable = variable)
    ),
    tar_stan_mcmc(
      fit_lr_nou,
      c("stan/pl.stan", "stan/pl_biome.stan"),
      data = stan_data_lr,
      refresh = 0,
      chains = 3,
      parallel_chains = 3,
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.9,
      seed = 123,
      return_draws = TRUE,
      return_diagnostics = TRUE,
      return_summary = TRUE,
      summaries = list(
        mean = ~mean(.x),
        sd = ~sd(.x),
        mad = ~mad(.x),
        ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
        posterior::default_convergence_measures()
      )
    )
  )

tar_map_nlr_nou <-
  tar_map(
    values = expand_grid(
      div = c("ang", "gym"),
      variable = c("cr", "h")),
    tar_target(
      tallo_reduced_nlr_df, {
        if (variable == "cr") {
          tmp <- tallo_wd_df0 |>
            filter(!is.na(cr))
        } else if (variable == "h") {
          tmp <- tallo_wd_df0 |>
            filter(!is.na(h))
        }
        reduce_trees(tmp, 100, 50)
      }
    ),
    tar_target(
      stan_data_nlr,
      generate_stan_data(tallo_reduced_nlr_df,
        model = "nlr", div = div, variable = variable)
    ),
    tar_stan_mcmc(
      fit_nlr_nou,
      c("stan/weibull_biome.stan", "stan/gmm_biome.stan", "stan/weibull.stan", "stan/gmm.stan"),
      data = stan_data_nlr,
      refresh = 0,
      chains = 3,
      parallel_chains = 3,
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.9,
      seed = 123,
      return_draws = TRUE,
      return_diagnostics = TRUE,
      return_summary = TRUE,
      summaries = list(
        mean = ~mean(.x),
        sd = ~sd(.x),
        mad = ~mad(.x),
        ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
        posterior::default_convergence_measures()
      )
     )
   )

tar_nlr_nou2 <-
    tar_stan_mcmc(
      fit_nlr_nou_re,
      "stan/weibull_biome_re.stan",
      data = stan_data_nlr_ang_h,
      refresh = 0,
      chains = 3,
      parallel_chains = 3,
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.9,
      seed = 123,
      return_draws = FALSE,
      return_diagnostics = TRUE,
      return_summary = TRUE,
      summaries = list(
        mean = ~mean(.x),
        sd = ~sd(.x),
        mad = ~mad(.x),
        ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
        posterior::default_convergence_measures()
      )
   )



main_ <- list(
  tar_target(
    sim_data,
    generate_sim_data()
  ),
  tar_map_loo_test,
  # tar_map_reduced_data,
  tar_map_lr_nou,
  tar_map_nlr_nou,
  tar_nlr_nou2,
  # model with wood density -----------------
  tar_stan_mcmc(
    fit_wd_ang_h,
    "stan/weibull_wd.stan",
    data = stan_data_nlr_ang_h,
    refresh = 0,
    chains = 3,
    parallel_chains = 3,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123,
    return_draws = TRUE,
    return_diagnostics = TRUE,
    return_summary = TRUE,
    summaries = list(
      mean = ~mean(.x),
      sd = ~sd(.x),
      mad = ~mad(.x),
      ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
      posterior::default_convergence_measures()
    )
  ),
  tar_stan_mcmc(
    fit_wd_gym_h,
    "stan/weibull_wd.stan",
    data = stan_data_nlr_gym_h,
    refresh = 0,
    chains = 3,
    parallel_chains = 3,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123,
    return_draws = TRUE,
    return_diagnostics = TRUE,
    return_summary = TRUE,
    summaries = list(
      mean = ~mean(.x),
      sd = ~sd(.x),
      mad = ~mad(.x),
      ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
      posterior::default_convergence_measures()
    )
  ),
  tar_stan_mcmc(
    fit_wd_ang_cr,
    "stan/pl_wd.stan",
    data = stan_data_lr_ang_cr,
    refresh = 0,
    chains = 3,
    parallel_chains = 3,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123,
    return_draws = TRUE,
    return_diagnostics = TRUE,
    return_summary = TRUE,
    summaries = list(
      mean = ~mean(.x),
      sd = ~sd(.x),
      mad = ~mad(.x),
      ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
      posterior::default_convergence_measures()
    )
  ),
  tar_stan_mcmc(
    fit_wd_gym_cr,
    "stan/gmm_wd.stan",
    data = stan_data_nlr_gym_cr,
    refresh = 0,
    chains = 3,
    parallel_chains = 3,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123,
    return_draws = TRUE,
    return_diagnostics = TRUE,
    return_summary = TRUE,
    summaries = list(
      mean = ~mean(.x),
      sd = ~sd(.x),
      mad = ~mad(.x),
      ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
      posterior::default_convergence_measures()
    )
  ),
  tar_stan_mcmc(
    fit_wd_ang_dbh,
    "stan/pl_wd.stan",
    data = stan_data_lr_ang_dbh,
    refresh = 0,
    chains = 3,
    parallel_chains = 3,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123,
    return_draws = TRUE,
    return_diagnostics = TRUE,
    return_summary = TRUE,
    summaries = list(
      mean = ~mean(.x),
      sd = ~sd(.x),
      mad = ~mad(.x),
      ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
      posterior::default_convergence_measures()
    )
  ),
  tar_stan_mcmc(
    fit_wd_gym_dbh,
    "stan/pl_wd.stan",
    data = stan_data_lr_gym_dbh,
    refresh = 0,
    chains = 3,
    parallel_chains = 3,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123,
    return_draws = TRUE,
    return_diagnostics = TRUE,
    return_summary = TRUE,
    summaries = list(
      mean = ~mean(.x),
      sd = ~sd(.x),
      mad = ~mad(.x),
      ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
      posterior::default_convergence_measures()
    )
  ),
 loo_map,
  tar_combine(
    loo_list,
    loo_map,
    command = list(!!!.x)
  ),
  tar_target(
    loo_tbl,
    generate_loo_tbl(loo_list)
  ),
  tar_target(
    test_loo,
    list(
      none = my_loo(fit_test_mcmc_test_none),
      xy = my_loo(fit_test_mcmc_test_xy),
      x = my_loo(fit_test_mcmc_test_x),
      none2 = my_loo(fit_test_mcmc_test2_none),
      xy2 = my_loo(fit_test_mcmc_test2_xy),
      x2 = my_loo(fit_test_mcmc_test2_x)
    )
  ),
  # BEST PREDICTIVE MODELS
  tar_target(
    best_pred_model,
    find_best_models(loo_tbl)
  ),
  # COMPARISON TABLE
  tar_target(
    com_tbl,
    generate_comparison_table(loo_tbl)
  ),
  # COMPUTING R-SQUARED
  ## H-DBH ANG
  tar_target(
    r2_pl_nou_h_ang,
    calculate_bayes_R2(
      stan_data_lr_ang_h,
      fit_lr_nou_draws_pl_ang_h,
      c("gamma[1]", "gamma[2]"),
      "lr"
    )
  ),
  tar_target(
    r2_gmm_nou_h_ang,
    calculate_bayes_R2(
      stan_data_nlr_ang_h,
      fit_nlr_nou_draws_gmm_ang_h,
      c("gamma[1]", "gamma[2]", "gamma[3]"),
      "gmm"
    )
  ),
  tar_target(
    r2_weibull_nou_h_ang,
    calculate_bayes_R2(
      stan_data_nlr_ang_h,
      fit_nlr_nou_draws_weibull_ang_h,
      c("gamma[1]", "gamma[2]", "gamma[3]"),
      "weibull"
    )
  ),
  tar_target(
    r2_weibull_h_ang,
    calculate_bayes_R2(
      stan_data_nlr_ang_h,
      fit_wd_ang_h_draws_weibull_wd,
      c("gamma[1,1]", "gamma[1,2]", "gamma[1,3]"),
      "weibull"
    )
  ),
  ## H-DBH GYM
  tar_target(
    r2_pl_nou_h_gym,
    calculate_bayes_R2(
      stan_data_lr_gym_h,
      fit_lr_nou_draws_pl_gym_h,
      c("gamma[1]", "gamma[2]"),
      "lr"
    )
  ),
  tar_target(
    r2_gmm_nou_h_gym,
    calculate_bayes_R2(
      stan_data_nlr_gym_h,
      fit_nlr_nou_draws_gmm_gym_h,
      c("gamma[1]", "gamma[2]", "gamma[3]"),
      "gmm"
    )
  ),
  tar_target(
    r2_weibull_nou_h_gym,
    calculate_bayes_R2(
      stan_data_nlr_gym_h,
      fit_nlr_nou_draws_weibull_gym_h,
      c("gamma[1]", "gamma[2]", "gamma[3]"),
      "weibull"
    )
  ),
  tar_target(
    r2_weibull_h_gym,
    calculate_bayes_R2(
      stan_data_nlr_gym_h,
      fit_wd_gym_h_draws_weibull_wd,
      c("gamma[1,1]", "gamma[1,2]", "gamma[1,3]"),
      "weibull"
    )
  ),
  ## CR-DBH ANG
  tar_target(
    r2_pl_nou_cr_ang,
    calculate_bayes_R2(
      stan_data_lr_ang_cr,
      fit_lr_nou_draws_pl_ang_cr,
      c("gamma[1]", "gamma[2]"),
      "lr"
    )
  ),
  tar_target(
    r2_pl_cr_ang,
    calculate_bayes_R2(
      stan_data_lr_ang_cr,
      fit_wd_ang_cr_draws_pl_wd,
      c("gamma[1,1]", "gamma[2,1]"),
      "lr"
    )
  ),
  tar_target(
    r2_gmm_nou_cr_ang,
    calculate_bayes_R2(
      stan_data_nlr_ang_cr,
      fit_nlr_nou_draws_gmm_ang_cr,
      c("gamma[1]", "gamma[2]", "gamma[3]"),
      "gmm"
    )
  ),
  tar_target(
    r2_weibull_nou_cr_ang,
    calculate_bayes_R2(
      stan_data_nlr_ang_cr,
      fit_nlr_nou_draws_weibull_ang_cr,
      c("gamma[1]", "gamma[2]", "gamma[3]"),
      "weibull"
    )
  ),

  ## CR-DBH GYM
  tar_target(
    r2_pl_nou_cr_gym,
    calculate_bayes_R2(
      stan_data_lr_gym_cr,
      fit_lr_nou_draws_pl_gym_cr,
      c("gamma[1]", "gamma[2]"),
      "lr"
    )
  ),
  tar_target(
    r2_gmm_nou_cr_gym,
    calculate_bayes_R2(
      stan_data_nlr_gym_cr,
      fit_nlr_nou_draws_gmm_gym_cr,
      c("gamma[1]", "gamma[2]", "gamma[3]"),
      "gmm"
    )
  ),
  tar_target(
    r2_gmm_cr_gym,
    calculate_bayes_R2(
      stan_data_nlr_gym_cr,
      fit_wd_gym_cr_draws_gmm_wd,
      c("gamma[1,1]", "gamma[1,2]", "gamma[1,3]"),
      "gmm"
    )
  ),
  tar_target(
    r2_weibull_nou_cr_gym,
    calculate_bayes_R2(
      stan_data_nlr_gym_cr,
      fit_nlr_nou_draws_weibull_gym_cr,
      c("gamma[1]", "gamma[2]", "gamma[3]"),
      "weibull"
    )
  ),
  ## DBH allometry
  ### ANG
  tar_target(
    r2_pl_nou_dbh_ang,
    calculate_bayes_R2(
      stan_data_lr_ang_dbh,
      fit_lr_nou_draws_pl_ang_dbh,
      c("gamma[1]", "gamma[2]", "gamma[2]"),
      "lr"
    )
  ),
  tar_target(
    r2_pl_dbh_ang,
    calculate_bayes_R2(
      stan_data_lr_ang_dbh,
      fit_wd_ang_dbh_draws_pl_wd,
      c("gamma[1,1]", "gamma[2,1]", "gamma[3,1]"),
      "lr"
    )
  ),
  tar_target(
    r2_pl_dbh1_ang,
    calculate_bayes_R2(
    stan_data_lr_ang_dbh1,
    fit_lr_nou_draws_pl_ang_dbh1,
    c("gamma[1]", "gamma[2]"),
    "lr"
    )
  ),
  tar_target(
    r2_pl_dbh2_ang,
    calculate_bayes_R2(
    stan_data_lr_ang_dbh2,
    fit_lr_nou_draws_pl_ang_dbh2,
    c("gamma[1]", "gamma[2]"),
    "lr"
    )
  ),
  tar_target(
    r2_pl_dbh3_ang,
    calculate_bayes_R2(
    stan_data_lr_ang_dbh3,
    fit_lr_nou_draws_pl_ang_dbh3,
    c("gamma[1]", "gamma[2]"),
    "lr"
    )
  ),
  ### GYM
  tar_target(
    r2_pl_nou_dbh_gym,
    calculate_bayes_R2(
      stan_data_lr_gym_dbh,
      fit_lr_nou_draws_pl_gym_dbh,
      c("gamma[1]", "gamma[2]", "gamma[2]"),
      "lr"
    )
  ),
  tar_target(
    r2_pl_dbh_gym,
    calculate_bayes_R2(
      stan_data_lr_gym_dbh,
      fit_wd_gym_dbh_draws_pl_wd,
      c("gamma[1,1]", "gamma[2,1]", "gamma[3,1]"),
      "lr"
    )
  ),
  tar_target(
    r2_pl_dbh1_gym,
    calculate_bayes_R2(
      stan_data_lr_gym_dbh1,
      fit_lr_nou_draws_pl_gym_dbh1,
      c("gamma[1]", "gamma[2]"),
      "lr"
    )
  ),
  tar_target(
    r2_pl_dbh2_gym,
    calculate_bayes_R2(
      stan_data_lr_gym_dbh2,
      fit_lr_nou_draws_pl_gym_dbh2,
      c("gamma[1]", "gamma[2]"),
      "lr"
    )
  ),
  tar_target(
    r2_pl_dbh3_gym,
    calculate_bayes_R2(
      stan_data_lr_gym_dbh3,
      fit_lr_nou_draws_pl_gym_dbh3,
      c("gamma[1]", "gamma[2]"),
      "lr"
    )
  ),
  ## R2 list
  tar_target(
    R2_list,
    list(
      # H ANG
      bayes_R2_pl_nou_h_ang = r2_pl_nou_h_ang,
      bayes_R2_gmm_nou_h_ang = r2_gmm_nou_h_ang,
      bayes_R2_weibull_nou_h_ang = r2_weibull_nou_h_ang,
      bayes_R2_weibull_h_ang = r2_weibull_h_ang,
      # H GYM
      bayes_R2_pl_nou_h_gym = r2_pl_nou_h_gym,
      bayes_R2_gmm_nou_h_gym = r2_gmm_nou_h_gym,
      bayes_R2_weibull_nou_h_gym = r2_weibull_nou_h_gym,
      bayes_R2_weibull_h_gym = r2_weibull_h_gym,
      # CR ANG
      bayes_R2_pl_nou_cr_ang = r2_pl_nou_cr_ang,
      bayes_R2_pl_cr_ang = r2_pl_cr_ang,
      bayes_R2_gmm_nou_cr_ang = r2_gmm_nou_cr_ang,
      bayes_R2_weibull_nou_cr_ang = r2_weibull_nou_cr_ang,
      # CR GYM
      bayes_R2_pl_nou_cr_gym = r2_pl_nou_cr_gym,
      bayes_R2_gmm_nou_cr_gym = r2_gmm_cr_gym,
      bayes_R2_gmm_cr_gym = r2_gmm_cr_gym,
      bayes_R2_weibull_nou_cr_gym = r2_weibull_nou_cr_gym,
      # DBH ANG
      bayes_R2_pl_nou_dbh_ang = r2_pl_nou_dbh_ang,
      bayes_R2_pl_dbh_ang = r2_pl_dbh_ang,
      bayes_R2_pl_nou_dbh1_ang = r2_pl_dbh1_ang,
      bayes_R2_pl_nou_dbh2_ang = r2_pl_dbh2_ang,
      bayes_R2_pl_nou_dbh3_ang = r2_pl_dbh3_ang,
      # DBH GYM
      bayes_R2_pl_nou_dbh_gym = r2_pl_nou_dbh_gym,
      bayes_R2_pl_dbh_gym = r2_pl_dbh_gym,
      bayes_R2_pl_nou_dbh1_gym = r2_pl_dbh1_gym,
      bayes_R2_pl_nou_dbh2_gym = r2_pl_dbh2_gym,
      bayes_R2_pl_nou_dbh3_gym = r2_pl_dbh3_gym
    )
  ),
  tar_target(
    R2_df,
    {
      R2_summary <- lapply(R2_list, function(x) {
        data.frame(
          median_R2 = median(x),
          lower_95CI = quantile(x, 0.025),
          upper_95CI = quantile(x, 0.975)
        )
      })
      R2_df <- do.call(rbind, R2_summary)
      R2_df <- R2_df[, c("median_R2","lower_95CI", "upper_95CI")]
      R2_df
    }
  ),
  tar_target(
    R2_tbl,
    generate_r2_table(
      R2_df
    )
  ),
  # LOO + R2 table
  tar_target(
    loo_r2_tbl,
    generate_loo_r2_table(
      R2_tbl,
      com_tbl
    )
  ),
  tar_target(
    loo_r2_tbl_csv,
    my_write_csv(loo_r2_tbl, "data/loo_r2_tbl.csv"),
    format = "file"
  ),
  # CURVE PLOT
  tar_target(
    curve_plot, {
      p <- generate_custom_curve_plot()
      my_ggsave(
        filename = "figs/curve",
        plot = p,
        dpi = 600,
        width = 82,
        height = 60,
        units = "mm"
      )
    },
    format = "file"
  ),
  # DATA MAP PLOT
  tar_target(
    geographical_coverage_map, {
      p <- data_map(tallo_wd_df0)
      my_ggsave(
        filename = "figs/map",
        plot = p,
        dpi = 600,
        width = 173,
        height = 140,
        units = "mm"
      )
    },
    format = "file"
  ),
  # Models' comparisons plots
  # ## H-DBH
  # tar_target(
  #   h_dbh, {
  #     p <- generate_h_dbh_plot(
  #       tallo_reduced_lr_df_ang_h,
  #       fit_lr_nou_summary_pl_ang_h,
  #       fit_lr_nou_summary_pl_gym_h,
  #       fit_nlr_nou_summary_gmm_ang_h,
  #       fit_nlr_nou_summary_gmm_gym_h,
  #       fit_nlr_nou_summary_weibull_ang_h,
  #       fit_nlr_nou_summary_weibull_gym_h
  #     )

  #     my_ggsave(
  #       plot = p,
  #       filename = "figs/h_dbh",
  #       dpi = 600,
  #       width = 82,
  #       height = 120,
  #       units = "mm"
  #     )
  #   },
  #   format = "file"
  # ),
  ## ALLOMETRY PLOT
  tar_target(
    h_cr_dbh_log, {
      p <- generate_allo_plot(
        tallo_reduced_lr_df_ang_h,
        tallo_reduced_lr_df_ang_cr,
        tallo_reduced_lr_df_gym_h,
        tallo_reduced_lr_df_gym_cr,
        fit_lr_nou_summary_pl_ang_h,
        fit_lr_nou_summary_pl_gym_h,
        fit_nlr_nou_summary_gmm_ang_h,
        fit_nlr_nou_summary_gmm_gym_h,
        fit_nlr_nou_summary_weibull_ang_h,
        fit_nlr_nou_summary_weibull_gym_h,
        fit_lr_nou_summary_pl_ang_cr,
        fit_lr_nou_summary_pl_gym_cr,
        fit_nlr_nou_summary_gmm_ang_cr,
        fit_nlr_nou_summary_gmm_gym_cr,
        fit_nlr_nou_summary_weibull_ang_cr,
        fit_nlr_nou_summary_weibull_gym_cr,
        log_scale = TRUE
      )
      my_ggsave(
        plot = p,
        filename = "figs/h_cr_dbh_log",
        dpi = 300,
        width = 173,
        height = 140,
        units = "mm"
      )
    },
    format = "file"
  ),

  tar_target(
    h_cr_dbh_non_log, {
      p <- generate_allo_plot(
        tallo_reduced_lr_df_ang_h,
        tallo_reduced_lr_df_ang_cr,
        tallo_reduced_lr_df_gym_h,
        tallo_reduced_lr_df_gym_cr,
        fit_lr_nou_summary_pl_ang_h,
        fit_lr_nou_summary_pl_gym_h,
        fit_nlr_nou_summary_gmm_ang_h,
        fit_nlr_nou_summary_gmm_gym_h,
        fit_nlr_nou_summary_weibull_ang_h,
        fit_nlr_nou_summary_weibull_gym_h,
        fit_lr_nou_summary_pl_ang_cr,
        fit_lr_nou_summary_pl_gym_cr,
        fit_nlr_nou_summary_gmm_ang_cr,
        fit_nlr_nou_summary_gmm_gym_cr,
        fit_nlr_nou_summary_weibull_ang_cr,
        fit_nlr_nou_summary_weibull_gym_cr,
        log_scale = FALSE
      )
      my_ggsave(
        plot = p,
        filename = "figs/h_cr_dbh_non_log",
        dpi = 300,
        width = 173,
        height = 140,
        units = "mm"
      )
    },
    format = "file"
  ),

  # Wood density vs. H-DBH plot
  tar_target(
    wd_h, {
      p <- generate_wd_para_h(
        tallo_reduced_nlr_df_ang_h,
        tallo_reduced_nlr_df_gym_h,
        fit_wd_ang_h_draws_weibull_wd,
        fit_wd_gym_h_draws_weibull_wd,
        stan_data_nlr_ang_h,
        stan_data_nlr_gym_h
      )

      my_ggsave(
        plot = p,
        filename = "figs/wd_h",
        dpi = 600,
        width = 173,
        height = 190,
        units = "mm"
      )
    },
    format = "file"
  ),
  # Wood density vs. CR-DBH plot
  tar_target(
    wd_cr, {
      p <- generate_wd_para_cr(
        tallo_reduced_lr_df_ang_cr,
        tallo_reduced_nlr_df_gym_cr,
        fit_wd_ang_cr_draws_pl_wd,
        fit_wd_gym_cr_draws_gmm_wd,
        stan_data_lr_ang_cr,
        stan_data_nlr_gym_cr
      )

      my_ggsave(
        plot = p,
        filename = "figs/wd_cr",
        dpi = 600,
        width = 173,
        height = 190,
        units = "mm"
      )
    },
    format = "file"
  ),
  # Combined WD-H and WD-CR plot
  tar_target(
    combined_wd_plot, {
      p <- generate_wd_para_com(
        tallo_reduced_nlr_df_ang_h,
        tallo_reduced_nlr_df_gym_h,
        fit_wd_ang_h_draws_weibull_wd,
        fit_wd_gym_h_draws_weibull_wd,
        stan_data_nlr_ang_h,
        stan_data_nlr_gym_h,
        tallo_reduced_lr_df_ang_cr,
        tallo_reduced_nlr_df_gym_cr,
        fit_wd_ang_cr_draws_pl_wd,
        fit_wd_gym_cr_draws_gmm_wd,
        stan_data_lr_ang_cr,
        stan_data_nlr_gym_cr
      )

      my_ggsave(
        plot = p,
        filename = "figs/wd_com",
        dpi = 600,
        width = 173,
        height = 190,
        units = "mm"
      )
    },
    format = "file"
  ),

  # WOOD DENSIRY's EFFECT
## SP POSTERIOR DF WITH WD
  # tar_target(
  #   sp_posterior_df_wd,
  #   generate_sp_posterior_df_wd(
  #     fit_wd_ang_h_summary_weibull_wd,
  #     fit_wd_ang_cr_summary_pl_wd,
  #     fit_wd_gym_cr_summary_gmm_wd,
  #     tallo_reduced_nlr_df_ang_h,
  #     tallo_reduced_lr_df_ang_cr,
  #     tallo_reduced_nlr_df_gym_cr,
  #     stan_data_nlr_ang_h,
  #     stan_data_lr_ang_cr,
  #     stan_data_nlr_gym_cr
  #   )
  # ),
  # tar_target(
  #   sp_posterior_h_ang_df_wd,
  #   sp_posterior_df_wd$sp_posterior_h_ang_df_wd
  # ),
  # tar_target(
  #   sp_posterior_cr_ang_df_wd,
  #   sp_posterior_df_wd$sp_posterior_cr_df_ang_wd
  # ),
  # tar_target(
  #   sp_posterior_cr_gym_df_wd,
  #   sp_posterior_df_wd$sp_posterior_cr_df_gym_wd
  # ),
  tar_target(
    wd_ef_plot, {
    p <- generate_wd_ef(
      tallo_reduced_nlr_df_ang_h,
      fit_wd_ang_h_draws_weibull_wd,
      tallo_reduced_lr_df_ang_cr,
      fit_wd_ang_cr_draws_pl_wd,
      tallo_reduced_nlr_df_gym_cr,
      fit_wd_gym_cr_draws_gmm_wd
    )
    my_ggsave(
        plot = p,
        filename = "figs/wd_ef",
        dpi = 600,
        width = 173,
        height = 47,
        units = "mm"
      )
    },
    format = "file"
  ),
  # POSTERIOR DATAFRAME
  tar_target(
    posterior_df,
    generate_posterior_df(
      pl_nou_ang_h <- fit_lr_nou_summary_pl_ang_h,
      gmm_nou_ang_h <- fit_nlr_nou_summary_gmm_ang_h,
      weibull_nou_ang_h <- fit_nlr_nou_summary_weibull_ang_h,
      weibull_ang_h <- fit_wd_ang_h_summary_weibull_wd,

      pl_nou_gym_h <- fit_lr_nou_summary_pl_gym_h,
      gmm_nou_gym_h <- fit_nlr_nou_summary_gmm_gym_h,
      weibull_nou_gym_h <- fit_nlr_nou_summary_weibull_gym_h,
      weibull_gym_h <- fit_wd_gym_h_summary_weibull_wd,

      pl_nou_ang_cr <- fit_lr_nou_summary_pl_ang_cr,
      pl_ang_cr <- fit_wd_ang_cr_summary_pl_wd,
      gmm_nou_ang_cr <- fit_nlr_nou_summary_gmm_ang_cr,
      weibull_nou_ang_cr <- fit_nlr_nou_summary_weibull_ang_cr,

      pl_nou_gym_cr <- fit_lr_nou_summary_pl_gym_cr,
      gmm_nou_gym_cr <- fit_nlr_nou_summary_gmm_gym_cr,
      gmm_gym_cr <- fit_wd_gym_cr_summary_gmm_wd,
      weibull_nou_gym_cr <- fit_nlr_nou_summary_weibull_gym_cr,

      pl_nou_ang_dbh <- fit_lr_nou_summary_pl_ang_dbh,
      pl_ang_dbh <- fit_wd_ang_dbh_summary_pl_wd,
      pl_nou_ang_dbh1 <- fit_lr_nou_summary_pl_ang_dbh1,
      pl_nou_ang_dbh2 <- fit_lr_nou_summary_pl_ang_dbh2,
      pl_nou_ang_dbh3 <- fit_lr_nou_summary_pl_ang_dbh3,

      pl_nou_gym_dbh <- fit_lr_nou_summary_pl_gym_dbh,
      pl_gym_dbh <- fit_wd_gym_dbh_summary_pl_wd,
      pl_nou_gym_dbh1 <- fit_lr_nou_summary_pl_gym_dbh1,
      pl_nou_gym_dbh2 <- fit_lr_nou_summary_pl_gym_dbh2,
      pl_nou_gym_dbh3 <- fit_lr_nou_summary_pl_gym_dbh3
    )
  ),
  tar_target(
    posterior_yaml,
    write_posterior_yaml(posterior_df),
    format = "file"
  ),
  # tar_target(
  #   posterior_rounded_df,
  #   round_posterior_df(posterior_df)
  # ),
  # tar_target(
  #   best_posterior_df,
  #   extract_best_models(posterior_rounded_df)
  # ),
  tar_target(
    best_posterior_df,
    extract_best_models(posterior_df)
  ),
  tar_target(
    best_posterior_df_csv,
    my_write_csv(best_posterior_df, "data/best_posterior_df.csv"),
    format = "file"
  ),
  # AGB ESTIMATION
  ## Species's parameters
  tar_target(
    sp_posterior_agb_df,
    generate_sp_posterior_agb_df(
      fit_nlr_nou_summary_weibull_ang_h,
      fit_lr_nou_summary_pl_ang_h,
      fit_lr_nou_summary_pl_ang_dbh,
      fit_lr_nou_summary_pl_ang_dbh1,
      fit_lr_nou_summary_pl_ang_dbh2,
      fit_lr_nou_summary_pl_ang_dbh3,
      tallo_reduced_nlr_df_ang_h,
      tallo_reduced_lr_df_ang_h,
      tallo_reduced_lr_df_ang_dbh,
      tallo_reduced_lr_df_ang_dbh1,
      tallo_reduced_lr_df_ang_dbh2,
      tallo_reduced_lr_df_ang_dbh3,
      stan_data_nlr_ang_h,
      stan_data_lr_ang_h,
      stan_data_lr_ang_dbh,
      stan_data_lr_ang_dbh1,
      stan_data_lr_ang_dbh2,
      stan_data_lr_ang_dbh3,

      fit_nlr_nou_summary_weibull_gym_h,
      fit_lr_nou_summary_pl_gym_h,
      fit_lr_nou_summary_pl_gym_dbh,
      fit_lr_nou_summary_pl_gym_dbh1,
      fit_lr_nou_summary_pl_gym_dbh2,
      fit_lr_nou_summary_pl_gym_dbh3,
      tallo_reduced_nlr_df_gym_h,
      tallo_reduced_lr_df_gym_h,
      tallo_reduced_lr_df_gym_dbh,
      tallo_reduced_lr_df_gym_dbh1,
      tallo_reduced_lr_df_gym_dbh2,
      tallo_reduced_lr_df_gym_dbh3,
      stan_data_nlr_gym_h,
      stan_data_lr_gym_h,
      stan_data_lr_gym_dbh,
      stan_data_lr_gym_dbh1,
      stan_data_lr_gym_dbh2,
      stan_data_lr_gym_dbh3
    )
  ),
# ALTERNATIVE AGB Eq1 ANGIOSPERM
  tar_target(
    agb,
    generate_agb_estimation(
      tallo_wd_df0,
      sp_posterior_agb_df,
      export_yaml = TRUE,
      yaml_file = "agb_metrics_ang.yaml"
    )
  ),
  tar_target(
    agb_ang_plot,
    {
      my_ggsave(
        filename = "figs/agb_ang",
        plot = agb$p,
        dpi = 600,
        width = 173,
        height = 115,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    agb_ang_metrics_yaml,
    {
      agb$metrics_agb_ang
      "agb_metrics_ang.yaml"
      },
    format = "file"
  ),
# DATA DESCRIPTION
  tar_target(
    data_yaml,
    write_data_yaml(
      tallo_wd_df0,
      tallo_reduced_lr_df_ang_h,
      tallo_reduced_lr_df_gym_h,
      tallo_reduced_lr_df_ang_cr,
      tallo_reduced_lr_df_gym_cr,
      tallo_reduced_lr_df_ang_dbh,
      tallo_reduced_lr_df_gym_dbh,
      tallo_reduced_lr_df_ang_dbh1,
      tallo_reduced_lr_df_gym_dbh1,
      tallo_reduced_lr_df_ang_dbh2,
      tallo_reduced_lr_df_gym_dbh2,
      tallo_reduced_lr_df_ang_dbh3,
      tallo_reduced_lr_df_gym_dbh3
    ),
    format = "file"
  ),
#SUPPLEMENT
## DATA
  tar_target(
    tbl_subdata,
    subset_charac(
      tallo_reduced_lr_df_ang_h,
      tallo_reduced_lr_df_ang_cr,
      tallo_reduced_lr_df_ang_dbh,
      tallo_reduced_lr_df_gym_dbh,
      tallo_reduced_lr_df_ang_dbh1,
      tallo_reduced_lr_df_gym_dbh1,
      tallo_reduced_lr_df_ang_dbh2,
      tallo_reduced_lr_df_gym_dbh2,
      tallo_reduced_lr_df_ang_dbh3,
      tallo_reduced_lr_df_gym_dbh3
    )
  ),
## SP POSTERIOR DF
  tar_target(
    sp_posterior_df,
    generate_sp_posterior_df(
      fit_nlr_nou_summary_weibull_ang_h,
      fit_nlr_nou_summary_weibull_gym_h,
      fit_lr_nou_summary_pl_ang_cr,
      fit_nlr_nou_summary_gmm_gym_cr,
      fit_lr_nou_summary_pl_ang_dbh,
      fit_lr_nou_summary_pl_gym_dbh,
      tallo_reduced_nlr_df_ang_h,
      tallo_reduced_nlr_df_gym_h,
      tallo_reduced_lr_df_ang_cr,
      tallo_reduced_nlr_df_gym_cr,
      tallo_reduced_lr_df_ang_dbh,
      tallo_reduced_lr_df_gym_dbh,
      stan_data_nlr_ang_h,
      stan_data_nlr_gym_h,
      stan_data_lr_ang_cr,
      stan_data_nlr_gym_cr,
      stan_data_lr_ang_dbh,
      stan_data_lr_gym_dbh

    )
  ),
  # 3 SIGNIFICANT FIGURES
  tar_target(
    sp_posterior_h_df,
    format_sp_posterior_df(sp_posterior_df$sp_posterior_h_df)
  ),
  tar_target(
    sp_posterior_cr_df,
    {
      df <- format_sp_posterior_df(sp_posterior_df$sp_posterior_cr_df) |>
        dplyr::mutate(
          Functional_form = ifelse(is.na(k), "Power-law", "gMM")
        ) |>
        dplyr::relocate(Functional_form, .after = Dependent_variable)
      df
    }
  ),
  tar_target(
    sp_posterior_dbh_df,
    format_sp_posterior_df(sp_posterior_df$sp_posterior_dbh_df)
  ),
  tar_target(
    save_sp_posterior_df_xlsx,
    {
      # Remove 'Dependent_variable' column
      sp_posterior_h_df_clean <- sp_posterior_h_df |> dplyr::select(-Dependent_variable)
      sp_posterior_cr_df_clean <- sp_posterior_cr_df |> dplyr::select(-Dependent_variable)
      sp_posterior_dbh_df_clean <- sp_posterior_dbh_df |> dplyr::select(-Dependent_variable)

      # Captions
      caption_height <- "Table S6. Posterior species-level parameter estimates (median and 95% Bayesian credible intervals) from the Weibull function (Eq. 4), the best-performing model for predicting tree height (m) from DBH (cm) across 1,290 species in two clades."
      caption_cr <- "Table S7. Posterior species-level parameter estimates (median and 95% Bayesian credible intervals) from the best-performing models for predicting crown radius (m) from DBH (cm) across 821 species in two clades. For angiosperms, the power-law function (Eq. 2) yielded the best predictions; for gymnosperms, the generalized Michaelis-Menten function (Eq. 3) was optimal."
      caption_dbh <- "Table S8. Posterior species-level parameter estimates (median and 95% Bayesian credible intervals) from the CR,H model, the best-performing model for predicting DBH (cm) from tree height (m) and crown radius (m) across 800 species in two clades."

      # Define bold & wrap style
      caption_style <- createStyle(textDecoration = "bold", wrapText = TRUE)

      # Save Table S6
      wb <- createWorkbook()
      addWorksheet(wb, "H_sp_parameter_estimates")
      writeData(wb, "H_sp_parameter_estimates", caption_height, startRow = 1, startCol = 1)
      mergeCells(wb, "H_sp_parameter_estimates", cols = 1:5, rows = 1)
      addStyle(wb, "H_sp_parameter_estimates", style = caption_style, rows = 1, cols = 1, gridExpand = TRUE)
      setRowHeights(wb, "H_sp_parameter_estimates", rows = 1, heights = 40)  # Adjust height for visibility
      writeData(wb, "H_sp_parameter_estimates", sp_posterior_h_df_clean, startRow = 2, startCol = 1)
      saveWorkbook(wb, "data/TableS6_height_sp_estimates.xlsx", overwrite = TRUE)

      # Save Table S7
      wb <- createWorkbook()
      addWorksheet(wb, "CR_sp_parameter_estimates")
      writeData(wb, "CR_sp_parameter_estimates", caption_cr, startRow = 1, startCol = 1)
      mergeCells(wb, "CR_sp_parameter_estimates", cols = 1:6, rows = 1)
      addStyle(wb, "CR_sp_parameter_estimates", style = caption_style, rows = 1, cols = 1, gridExpand = TRUE)
      setRowHeights(wb, "CR_sp_parameter_estimates", rows = 1, heights = 40)
      writeData(wb, "CR_sp_parameter_estimates", sp_posterior_cr_df_clean, startRow = 2, startCol = 1)
      saveWorkbook(wb, "data/TableS7_crown_radius_sp_estimates.xlsx", overwrite = TRUE)

      # Save Table S8
      wb <- createWorkbook()
      addWorksheet(wb, "DBH_sp_parameter_estimates")
      writeData(wb, "DBH_sp_parameter_estimates", caption_dbh, startRow = 1, startCol = 1)
      mergeCells(wb, "DBH_sp_parameter_estimates", cols = 1:5, rows = 1)
      addStyle(wb, "DBH_sp_parameter_estimates", style = caption_style, rows = 1, cols = 1, gridExpand = TRUE)
      setRowHeights(wb, "DBH_sp_parameter_estimates", rows = 1, heights = 40)
      writeData(wb, "DBH_sp_parameter_estimates", sp_posterior_dbh_df_clean, startRow = 2, startCol = 1)
      saveWorkbook(wb, "data/TableS8_dbh_sp_estimates.xlsx", overwrite = TRUE)

      c(
        "data/TableS6_height_sp_estimates.xlsx",
        "data/TableS7_crown_radius_sp_estimates.xlsx",
        "data/TableS8_dbh_sp_estimates.xlsx"
      )
    },
    format = "file"
  ),
## SPECIES-LEVEL PLOTTING
  tar_target(
    h_sp_plot, {
      p <- generate_h_sp_plot(
        tallo_reduced_lr_df_ang_h,
        sp_posterior_h_df
      )
      my_ggsave(
        plot = p,
        filename = "figs/h_sp",
        dpi = 600,
        width = 200,
        height = 300,
        units = "mm"
      )
    },
    format = "file"
  ),
## PLOT SPECIES WITH LARGE TREES
  tar_target(
    h_sp_large_plot, {
      p <- generate_h_sp_large_plot(
        tallo_reduced_lr_df_ang_h,
        sp_posterior_h_df
      )
      my_ggsave(
        plot = p,
        filename = "figs/h_sp_large",
        dpi = 600,
        width = 200,
        height = 240,
        units = "mm"
      )
    },
    format = "file"
  ),
NULL
)

list(data_, main_)
