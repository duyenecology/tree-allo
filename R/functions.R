#' @inheritParams readr::write_csv
my_write_csv <- function(x, path, append = FALSE, col_names = !append) {
    write_csv(x, path, append = FALSE, col_names = !append)
    paste(path)
}

#=============================
# CLEANING DATA
#=============================f

clean_tallo_try <- function(tallo_csv,  tallo_env_csv, try_wd_txt) {
  tallo <- read_csv(tallo_csv)
  tallo2 <- tallo |>
    dplyr::select(tree_id,
          division, family, genus, latitude, longitude, reference_id,
          sp = species, dbh = stem_diameter_cm, h = height_m, cr = crown_radius_m,
          height_outlier, crown_radius_outlier) |>
          filter(height_outlier != 'Y' & crown_radius_outlier != 'Y' & !is.na(sp) &
          sp != "unknown" & !is.na(dbh))

  env_df <- read_csv(tallo_env_csv)

  ai_df <- env_df |>
    group_by(biome) |>
    summarize(mean_aridity_index= mean(aridity_index, na.rm = TRUE))

  tallo_env <- left_join(tallo2, env_df, by = c("tree_id" = "tree_id"))

  try_data <- read_delim(try_wd_txt, delim = "\t") |>
    dplyr::select(sp = AccSpeciesName, wd = StdValue, unit = UnitName) |>
    filter(!is.na(sp) & sp != "unknown" & !is.na(wd) & !is.na(unit) & unit == "g/cm3")

  # Calculate mean wood density for each species
  mean_density <- try_data |>
    group_by(sp) |>
    summarize(wd = mean(as.numeric(wd), na.rm = TRUE))

  # Join selected_data and wd_filtered
  left_join(tallo_env, mean_density, by = "sp") |>
    left_join(ai_df) |>
    mutate(biome_sp = paste(biome, sp, sep = "_")) |>
    filter(!is.na(wd)) |>
    mutate(biome_div = paste(biome, division, sep = "_"))

}


#=============================
# STAN DATA
#=============================

generate_stan_data <- function(tallo_wd_df_200, model = c("lr", "nlr"), div = c("ang", "gym"), variable = c("cr", "h", "dbh", "dbh1", "dbh2", "dbh3")) {

   # Filter out rows with missing values in dbh, cr, h, and wd_s columns
  if (variable == "cr") {
    data <- tallo_wd_df_200 |>
      filter(!is.na(cr))
  } else if (variable == "h") {
    data <- tallo_wd_df_200 |>
      filter(!is.na(h))
  } else {
    data <- tallo_wd_df_200 |>
      filter(!is.na(cr)) |>
      filter(!is.na(h))
  }

  species_counts <- data |>
    group_by(sp) |>
    summarise(
      count = n()) |>
    filter(count >= 20)

  data2 <- data |>
    filter(sp %in% species_counts$sp) |>
    filter(division == if_else(div == "ang", "Angiosperm", "Gymnosperm"))

  dbh_max <- max(data2$dbh)

  data2 <- data2 |>
    mutate(log_dbh = log(dbh),
           log_h = log(h),
           log_cr = log(cr),
           log_cr_h = log(cr * h)) |>
    mutate(log_dbh_s = scale(log_dbh) |> as.numeric(),
           log_h_s = scale(log_h) |> as.numeric(),
           log_cr_s = scale(log_cr) |> as.numeric(),
           log_cr_h_s = scale(log_cr_h) |> as.numeric()) |>
    mutate(dbh_s = scale(dbh) |> as.numeric())

  # Calculate necessary statistics for the Stan model
  N <- nrow(data2)
  J <- nrow(distinct(data2, sp))  # Count the number of distinct species
  S <- nrow(distinct(data2, biome))

  # sp-level data
  wd_df <- data2 |>
    group_by(sp) |>
    summarize(wd = mean(wd)) |>
    mutate(wd_s = scale(wd) |> as.numeric())

  # Create species indicator array
  jj <- as.integer(as.factor(data2$sp))
  ss <- as.integer(as.factor(data2$biome))

  # Create species-level predictors matrix

  if (model == "lr") {
    x <- 1
    xs <- 1
    u <- rbind(1, wd_df$wd_s)
    L <- nrow(u)
    if (variable == "cr") {
      log_y <- data2$log_cr
      log_x <- cbind(1, data2$log_dbh)
    } else if (variable == "h") {
      log_y <- data2$log_h
      log_x <- cbind(1, data2$log_dbh)
    } else if (variable == "dbh") {
      # use scaled y for the fast convergence
      log_y <- data2$log_dbh_s
      log_x <- cbind(1, data2$log_cr_s, data2$log_h_s)
    } else if(variable == "dbh1") {
      log_y <- data2$log_dbh_s
      log_x <- cbind(1, data2$log_cr_h_s)
    } else if (variable == "dbh2") {
      log_y <- data2$log_dbh_s
      log_x <- cbind(1, data2$log_cr_s)
    } else if (variable == "dbh3") {
      log_y <- data2$log_dbh_s
      log_x <- cbind(1, data2$log_h_s)
    }
  } else if (model == "nlr") {
    x <- data2$dbh
    xs <- data2$dbh / dbh_max
    u <- cbind(1, wd_df$wd_s)
    L <- ncol(u)
    log_x <- 1
    if (variable == "cr") {
      log_y <- data2$log_cr
    } else if (variable == "h") {
      log_y <- data2$log_h
    }
  }

  list(
    N = N,
    J = J,
    S = S,
    K = ifelse(model == "lr", ncol(log_x), 3),
    L = L,
    log_y = log_y,
    x = x,
    xs = xs,
    dbh_max = dbh_max,
    log_x = log_x,
    u = u,
    jj = jj,
    ss = ss
  )
}


sp_name_match <- function(bb_sp_list_csv, abund_xtbg_csv, bb_allometry_csv, path) {
  sp_list <- read_csv(bb_sp_list_csv) |>
    janitor::clean_names()
  abund <- read_csv(abund_xtbg_csv) |>
    janitor::clean_names()
  abund <- abund |> dplyr::select(cname, sp, life)
  allo <- read_csv(bb_allometry_csv) |>
    janitor::clean_names()
  allo_re <- left_join(allo, abund, by = c("sp_cn" = "cname")) |>
    mutate(sp = ifelse(sp_cn == "橄榄", "CANATO", sp))
  my_write_csv(allo_re, path)

}

#=============================
# LOO COMPUTING
#=============================

my_loo <- function(x) x$loo(cores = parallel::detectCores())

#=============================
# LOO TABLE
#=============================

#' Generate a LOO Table from LOO Objects
#'
#' This function generates a tibble containing model names and their corresponding LOO (Leave-One-Out)
#' cross-validation statistics, including Expected Log Predictive Density (ELPD),
#' effective number of parameters (p_loo), and LOO Information Criterion (LOOIC).
#'
#' @param loo_list A named list of LOO objects, where each element corresponds to a model's LOO analysis.
#' @return A tibble with columns for model names, ELPD, p_loo, and LOOIC values.
#' @export
generate_loo_tbl <- function(loo_list, each = FALSE) {
  # Extract model names from the loo_list
  loo_names <- names(loo_list)
  # loo_names <- "fit_h_wd_mcmc_weibull_nc_wd_ang_h"
  # str_split(loo_names, "_")
  # Create a tibble with model names and corresponding LOO statistics
  loo_tbl <- tibble(model = loo_names) |>
    mutate(
      elpd = map_dbl(loo_list, ~ .x$estimates["elpd_loo", "Estimate"]),
      p_loo = map_dbl(loo_list, ~ .x$estimates["p_loo", "Estimate"]),
      looic = map_dbl(loo_list, ~ .x$estimates["looic", "Estimate"])
    ) |>
    rename(tmp = model) |>
    mutate(
      lr = case_when(
        str_detect(tmp, "pl") ~ "lr",
        TRUE  ~ "nlr"
      ),
      model = case_when(
        str_detect(tmp, "pl") ~ "pl",
        str_detect(tmp, "weibull") ~ "weibull",
        str_detect(tmp, "gmm") ~ "gmm",
        str_detect(tmp, "eg") ~ "eg"
      ),
      sp_biome = case_when(
        str_detect(tmp, "biome") ~ "biome",
        TRUE  ~ "sp-only"
      ),
      div = case_when(
        str_detect(tmp, "ang") ~ "ang",
        str_detect(tmp, "gym") ~ "gym"
      ),
      # pred = case_when(
      #   each ~ map_chr(str_split(tmp, "_"), ~ .x[length(.x) - 1]),
      #   TRUE ~ map_chr(str_split(tmp, "_"), ~ .x[length(.x)])
      # ),
      pred = case_when(
        str_detect(tmp, "_h") ~ "h",
        str_detect(tmp, "_cr") ~ "cr",
        str_detect(tmp, "_dbh1") ~ "dbh1",
        str_detect(tmp, "_dbh2") ~ "dbh2",
        str_detect(tmp, "_dbh3") ~ "dbh3",
        str_detect(tmp, "_dbh") ~ "dbh"
      ),
      biome = case_when(
        each ~ map_chr(str_split(tmp, "_"), ~ .x[length(.x)]),
        TRUE ~ NA_character_
      ),
      wd = case_when(
        str_detect(tmp, "wd") ~ "yes",
        TRUE ~ "no",
      )
    ) |>
    arrange(looic) |>
    arrange(pred)

  # Return the resulting table
  return(loo_tbl)
}


#=============================
# THE BEST PREDICTIVE MODELS
#=============================

#' Find the Best Model Based on LOOIC
#' Find the Best Models for 'sp-only', 'ang', and 'gym'
#'
#' This function identifies the best model for each predictor (cr, dbh, h) in the 'sp-only', 'ang', and 'gym' categories
#' based on the lowest LOOIC value.
#'
#' @param loo_tbl A tibble containing model names, LOOIC values, and other metadata.
#' @return A tibble with the best model for each predictor in the 'sp-only', 'ang', and 'gym' categories.
#' @export
find_best_models <- function(loo_tbl) {
  best_models <- loo_tbl |>
    filter(sp_biome == "sp-only" | div %in% c("ang", "gym")) |>  # Filter to include 'sp-only', 'ang', and 'gym'
    group_by(pred, div) |>                                       # Group by predictor and division (ang, gym)
    slice_min(looic, n = 1) |>                                   # Select the model with the lowest LOOIC in each group
    ungroup()                                                    # Ungroup after processing

  return(best_models)
}



#' Generate MCMC Summary and Extract Draws
#'
#' This function generates a summary of MCMC results and extracts the posterior draws as a data frame.
#'
#' @param mcmc An MCMC object containing the results of a Bayesian model.
#' @return A list containing the MCMC summary statistics and the posterior draws in a data frame format.
#' @export
generate_mcmc_summary <- function(mcmc) {
  list(
    summary = mcmc$summary(),
    draws = posterior::as_draws_df(mcmc$draws())
  )
}

#' @title Check divergence from draws
div_check <- function(diags) {
  diagnostics <- diags
  num_divergent <- sum(diagnostics$divergent__ == 1)
  total_samples <- nrow(diagnostics)
  percentage_divergent <- (num_divergent / total_samples) * 100
  cat("Number of divergent transitions:", num_divergent, "\n")
  cat("Percentage of divergent transitions:", round(percentage_divergent, 2), "%\n")

  list(
    diagnostics = diagnostics,
    num_divergent = num_divergent,
    percentage_divergent = percentage_divergent
  )
}

divergent_sum <- function(divergent_list) {
  tibble::tibble(
    model = names(divergent_list),
    num_divergent = sapply(divergent_list, function(x) x$num_divergent),
    percentage_divergent = sapply(divergent_list, function(x) x$percentage_divergent)
  )
}

#' Extract Posterior Summary for Gamma Parameters
#'
#' This function filters and extracts the posterior summary statistics for the `gamma` parameters from a model summary.
#'
#' @param sum A data frame containing the summary statistics of a model.
#' @return A filtered data frame with the `gamma` parameters, including their mean, standard deviation, and credible intervals.
#' @export
extract_posterior <- function(sum) {
  # Filter the summary data for the `gamma` parameters
  gamma <- sum |>
    filter(str_detect(variable, "gamma")) |>
    dplyr::select(variable, mean, sd, q2.5, q50, q97.5)
  return(gamma)
}

#=============================
# MODELS COMPARISON' TABLE
#=============================

#' Generate Comparison Table for Different Models
#'
#' This function processes LOO results from various models (Tree Height, Crown Radius, DBH)
#' and generates a formatted table with important metrics like ELPD, pLOO, LOOIC, and ΔLOOIC.
#'
#' @param loo_tbl A tibble containing the LOO results for multiple models.
#'
#' @return A list of formatted data frames for Tree Height, Crown Radius, and DBH models.
#'
#' @export
generate_comparison_table <- function(loo_tbl) {
  
  calculate_delta_looic <- function(df) {
    df <- df |>
      group_by(div) |>
      arrange(looic, .by_group = TRUE) |>
      mutate(ΔLOOIC = looic - min(looic)) |>
      ungroup()
    return(df)
  }
  
  format_table <- function(df) {
    formatted_df <- df |>
      mutate(
        `Dependent variable` = case_when(
          grepl("dbh", tmp) ~ "DBH",
          grepl("_h", tmp) ~ "Tree Height",
          grepl("_cr", tmp) ~ "Crown Radius",
          TRUE ~ NA_character_
        ),
        `Predictor variable` = case_when(
          grepl("dbh_mcmc_pl_wd", tmp) ~ "CR, H",
          grepl("dbh1", tmp) ~ "CR × H",
          grepl("dbh2", tmp) ~ "CR",
          grepl("dbh3", tmp) ~ "H",
          grepl("dbh$", tmp) ~ "CR, H",
          grepl("_h", tmp) ~ "DBH",
          grepl("_cr", tmp) ~ "DBH",
          TRUE ~ NA_character_
        ),
        `Functional form` = case_when(
          grepl("pl", model, ignore.case = TRUE) ~ "Power-law",
          grepl("gmm", model, ignore.case = TRUE) ~ "gMM",
          grepl("weibull", model, ignore.case = TRUE) ~ "Weibull",
          TRUE ~ model
        ),
        `Division` = case_when(
          div == "ang" ~ "Angiosperm",
          div == "gym" ~ "Gymnosperm",
          TRUE ~ div
        )
      ) |>
      rename(
        "Wood density" = Wood_Density,
        "ELPD" = elpd,
        "pLOO" = p_loo,
        "LOOIC" = looic,
        "ΔLOOIC" = ΔLOOIC
      ) |>
      mutate(across(c(ELPD, pLOO, LOOIC, `ΔLOOIC`), \(x) round(x, 2))) |>
      dplyr::select(`Dependent variable`, `Predictor variable`, Division, `Functional form`, `Wood density`, ELPD, pLOO, LOOIC, `ΔLOOIC`) |>
      as.data.frame()
    
    return(formatted_df)
  }

  # Filter and format the models
  h_models_tbl <- loo_tbl |>
    filter(str_detect(tmp, "_h"), !str_detect(sp_biome, "biome")) |>
    mutate(Wood_Density = ifelse(wd == "no", "Without", "With")) |>
    dplyr::select(tmp, model, div, Wood_Density, elpd, p_loo, looic) |>
    calculate_delta_looic() |>
    format_table()

  cr_models_tbl <- loo_tbl |>
    filter(str_detect(tmp, "_cr"), !str_detect(sp_biome, "biome")) |>
    mutate(Wood_Density = ifelse(wd == "no", "Without", "With")) |>
    dplyr::select(tmp, model, div, Wood_Density, elpd, p_loo, looic) |>
    calculate_delta_looic() |>
    format_table()

  dbh_models_tbl <- loo_tbl |>
    filter(str_detect(tmp, "dbh|dbh1|dbh2|dbh3"), !str_detect(sp_biome, "biome")) |>
    mutate(Wood_Density = ifelse(wd == "no", "Without", "With")) |>
    dplyr::select(tmp, model, div, Wood_Density, elpd, p_loo, looic) |>
    calculate_delta_looic() |>
    format_table()

  # Combine the three tables with Tree Height first, then Crown Radius, then DBH
  com_tbl <- dplyr::bind_rows(h_models_tbl, cr_models_tbl, dbh_models_tbl)

  return(com_tbl)
}

#' Calculate Bayesian R-squared for 26 tested Models (3 functions)
#'
#' @description
#' This function calculates Bayesian R-squared values for various models, including Power-law,
#' generalized Michaelis-Menten (gMM), and Weibull models.
#' @param stan_data A list containing the input data used in the Stan models. Each element in the list represents a specific dataset for a model.
#' @param fit_draws A list of Stan fit objects containing posterior samples for each model.
#' @param predictors A character vector specifying the names of the predictors (e.g., "gamma[1]", "gamma[2]") used in the model.
#' @param model_type A character string indicating the type of model ("lr" for linear regression, "gmm" for generalized Michaelis-Menten, "weibull" for Weibull).
#'
#' @return A list of Bayesian R-squared values for each model. Each element in the list contains the R-squared values calculated from the posterior samples.
#' @export
calculate_bayes_R2 <- function(stan_data, fit_draws, predictors, model_type) {

  # Extract the gamma and sigma draws
  gamma_draws <- lapply(predictors, function(p) as.matrix(fit_draws[, p]))
  sigma_draws <- as.matrix(fit_draws[, "sigma"])

  # Determine the number of posterior samples and observations
  num_sam <- nrow(gamma_draws[[1]])

  if (model_type == "lr") {
    # Linear Regression Case
    num_obs <- nrow(stan_data$log_x)
    log_y_pred_draws <- matrix(NA, nrow = num_sam, ncol = num_obs)

    # Calculate predicted log_y values for each posterior sample
    for (i in 1:num_sam) {
      log_y_pred_draws[i, ] <- gamma_draws[[1]][i] + gamma_draws[[2]][i] * stan_data$log_x[, 2]
      if (length(predictors) == 3) {
        log_y_pred_draws[i, ] <- log_y_pred_draws[i, ] + gamma_draws[[3]][i] * stan_data$log_x[, 3]
      }
    }

  } else if (model_type == "gmm" || model_type == "weibull") {
    # Nonlinear Regression Case (GMM or Weibull)
    num_obs <- length(stan_data$x)
    log_y_pred_draws <- matrix(NA, nrow = num_sam, ncol = num_obs)

    # Calculate predicted log_y values for each posterior sample based on model type
    for (i in 1:num_sam) {
      if (model_type == "gmm") {
        if (length(predictors) == 3) {
          log_y_pred_draws[i, ] <- gamma_draws[[1]][i] +
                                   (gamma_draws[[2]][i]) * log(stan_data$x) -
                                   log((gamma_draws[[3]][i]) + (stan_data$x)^(gamma_draws[[2]][i]))
        }
      } else if (model_type == "weibull") {
        if (length(predictors) == 3) {
          log_y_pred_draws[i, ] <- gamma_draws[[1]][i] +
                                   log(1 - exp(-(gamma_draws[[2]][i]) *
                                   (stan_data$x)^(gamma_draws[[3]][i])))
        }
      }
    }
  } else {
    stop("Unknown model type")
  }

  # Calculate variance of predicted log_y values and residual variance
  var_log_y_pred <- apply(log_y_pred_draws, 1, var, na.rm = TRUE)
  var_res <- sigma_draws^2

  # Calculate Bayesian R^2 for each posterior sample
  bayes_R2_draws <- var_log_y_pred / (var_log_y_pred + var_res)

  return(bayes_R2_draws)
}


#' Generate R2 Table
#'
#' This function processes a data frame containing R2 values for different models.
#' It renames the 'median_R2' column to 'R2', rounds it to three digits, and adds
#' relevant columns describing the dependent variable, predictor variable,
#' functional form, and whether wood density is included. The resulting data frame
#' is ordered by dependent variable.
#'
#' @param r2_df A data frame containing the R2 values and related information.
#'
#' @return A data frame with renamed and additional columns, ordered by dependent variable.
#' @export
generate_r2_table <- function(r2_df) {
  # Convert row names to a column for easier manipulation
  r2_df <- tibble::rownames_to_column(r2_df, var = "Model")

  # Add the relevant columns without rounding the R2 values
  r2_df <- r2_df |>
    dplyr::mutate(
      Dependent_variable = dplyr::case_when(
        grepl("h", Model) & !grepl("dbh", Model) ~ "Tree height",
        grepl("cr", Model) ~ "Crown radius",
        grepl("dbh", Model) ~ "DBH",
        TRUE ~ "Unknown"
      ),
      Division = dplyr::case_when(
        grepl("_ang", Model) ~ "Angiosperm",
        grepl("_gym", Model) ~ "Gymnosperm",
        TRUE ~ "Unknown"
      ),
      Predictor_variable = dplyr::case_when(
        Dependent_variable == "Tree height" ~ "DBH",
        Dependent_variable == "Crown radius" ~ "DBH",
        Dependent_variable == "DBH" & grepl("dbh1", Model) ~ "CR × H",
        Dependent_variable == "DBH" & grepl("dbh2", Model) ~ "CR",
        Dependent_variable == "DBH" & grepl("dbh3", Model) ~ "H",
        Dependent_variable == "DBH" & !grepl("dbh1|dbh2|dbh3", Model) ~ "CR, H",
        TRUE ~ "Unknown"
      ),
      Functional_form = dplyr::case_when(
        grepl("pl", Model) ~ "Power-law",
        grepl("gmm", Model) ~ "gMM",
        grepl("weibull", Model) ~ "Weibull",
        TRUE ~ "Unknown"
      ),
      Wood_density = dplyr::case_when(
        grepl("nou", Model) ~ "Without",
        TRUE ~ "With"
      )
    )

  # Arrange by Dependent_variable and Division
  r2_df <- r2_df |>
    dplyr::arrange(
      factor(Dependent_variable, levels = c("Tree height", "Crown radius", "DBH")),
      factor(Division, levels = c("Angiosperm", "Gymnosperm"))
    )

  # Select and order the columns with Division as the third column
  r2_df <- r2_df |>
    dplyr::select(
      Dependent_variable, Predictor_variable, Division, Functional_form, Wood_density, median_R2, lower_95CI, upper_95CI, Model
    )

  return(r2_df)
}
#=============================
# LOO + R2
#=============================
#' Combine LOO and R2 Tables
#'
#' This function combines the LOO information criterion table (`com_tbl`) with the R2 table (`r2_tbl`),
#' ensuring that the R2 values are aligned with the corresponding model information.
#'
#' @param r2_tbl A data frame containing R2 values
#' @param com_tbl A data frame containing LOOIC results
#' @return A combined data frame with the relevant columns from both `r2_tbl` and `com_tbl`.
#' @export
#'
generate_loo_r2_table <- function(r2_tbl, com_tbl) {
  # Rename columns in R2_tbl to match com_tbl
  r2_tbl <- r2_tbl |>
    dplyr::rename(
      `Dependent variable` = Dependent_variable,
      `Predictor variable` = Predictor_variable,
      `Functional form` = Functional_form,
      `Wood density` = Wood_density
    )
  
  # Standardize the format of key columns for matching
  r2_tbl <- r2_tbl |>
    dplyr::mutate(
      `Dependent variable` = stringr::str_replace_all(
        stringr::str_to_title(stringr::str_trim(`Dependent variable`)), 
        c("\\bDbh\\b" = "DBH", "\\bCr\\b" = "CR", "\\bH\\b" = "H")
      ),
      `Predictor variable` = stringr::str_replace_all(
        stringr::str_to_title(stringr::str_trim(`Predictor variable`)), 
        c("\\bDbh\\b" = "DBH", "\\bCr\\b" = "CR", "\\bH\\b" = "H")
      ),
      Division = stringr::str_to_title(stringr::str_trim(Division)),
      `Functional form` = stringr::str_replace_all(
        stringr::str_to_title(stringr::str_trim(`Functional form`)), 
        "Gmm", "gMM"
      ),
      `Wood density` = stringr::str_to_title(stringr::str_trim(`Wood density`))
    )
  
  com_tbl <- com_tbl |>
    dplyr::mutate(
      `Dependent variable` = stringr::str_replace_all(
        stringr::str_to_title(stringr::str_trim(`Dependent variable`)), 
        c("\\bDbh\\b" = "DBH", "\\bCr\\b" = "CR", "\\bH\\b" = "H")
      ),
      `Predictor variable` = stringr::str_replace_all(
        stringr::str_to_title(stringr::str_trim(`Predictor variable`)), 
        c("\\bDbh\\b" = "DBH", "\\bCr\\b" = "CR", "\\bH\\b" = "H")
      ),
      Division = stringr::str_to_title(stringr::str_trim(Division)),
      `Functional form` = stringr::str_replace_all(
        stringr::str_to_title(stringr::str_trim(`Functional form`)), 
        "Gmm", "gMM"
      ),
      `Wood density` = stringr::str_to_title(stringr::str_trim(`Wood density`))
    )
  
  # Perform the join to include R2 values from R2_tbl into com_tbl
  loo_r2_tbl <- com_tbl |>
    dplyr::left_join(
      r2_tbl |> dplyr::select(
        `Dependent variable`, `Predictor variable`, Division, `Functional form`, 
        `Wood density`, median_R2
      ),
      by = c("Dependent variable", "Predictor variable", "Division", "Functional form", "Wood density")
    )
  
  # Adjust R2 column and drop original median_R2 column
  loo_r2_tbl <- loo_r2_tbl |>
    dplyr::mutate(
      R2 = round(median_R2, 3)
    ) |>
    dplyr::select(-median_R2)
  
  return(loo_r2_tbl)
}
#=============================
# POSTERIOR DATAFRAME
#=============================
#' Combine Posterior Estimates from 28 Models
#'
#' This function merges posterior estimates from various models related to tree height, crown radius, and DBH for both angiosperms and gymnosperms into a single data frame.
#'
#' @param pl_nou_ang_h, gmm_nou_ang_h, weibull_nou_ang_h, weibull_ang_h Posterior estimates for angiosperm tree height models.
#' @param pl_nou_gym_h, gmm_nou_gym_h, weibull_nou_gym_h, weibull_gym_h Posterior estimates for gymnosperm tree height models.
#' @param pl_nou_ang_cr, pl_ang_cr, gmm_nou_ang_cr, weibull_nou_ang_cr Posterior estimates for angiosperm crown radius models.
#' @param pl_nou_gym_cr, gmm_nou_gym_cr, gmm_gym_cr, weibull_nou_gym_cr Posterior estimates for gymnosperm crown radius models.
#' @param pl_nou_ang_dbh, pl_ang_dbh, pl_nou_ang_dbh1, pl_nou_ang_dbh2, pl_nou_ang_dbh3 Posterior estimates for angiosperm DBH models.
#' @param pl_nou_gym_dbh, pl_gym_dbh, pl_nou_gym_dbh1, pl_nou_gym_dbh2, pl_nou_gym_dbh3 Posterior estimates for gymnosperm DBH models.
#' @return A combined data frame with all posterior estimates.
#' @export
#' 
generate_posterior_df <- function(
  pl_nou_ang_h, gmm_nou_ang_h, weibull_nou_ang_h, weibull_ang_h,
  pl_nou_gym_h, gmm_nou_gym_h, weibull_nou_gym_h, weibull_gym_h,
  pl_nou_ang_cr, pl_ang_cr, gmm_nou_ang_cr, weibull_nou_ang_cr,
  pl_nou_gym_cr, gmm_nou_gym_cr, gmm_gym_cr, weibull_nou_gym_cr,
  pl_nou_ang_dbh, pl_ang_dbh, pl_nou_ang_dbh1, pl_nou_ang_dbh2, pl_nou_ang_dbh3,
  pl_nou_gym_dbh, pl_gym_dbh, pl_nou_gym_dbh1, pl_nou_gym_dbh2, pl_nou_gym_dbh3) {

  create_model_df <- function(gamma_data, dependent_variable, functional_form, wood_density, division) {
    tau_summary <- gamma_data |> filter(str_detect(variable, "tau"))

    if (dependent_variable == "DBH" && functional_form == "Power-Law") {
      if (wood_density == "With") {
        a_median <- exp(gamma_data |> filter(variable == "gamma[1,1]") |> pull(q50))
        a_upper <- exp(gamma_data |> filter(variable == "gamma[1,1]") |> pull(q97.5))
        a_lower <- exp(gamma_data |> filter(variable == "gamma[1,1]") |> pull(q2.5))

        b_median <- gamma_data |> filter(variable == "gamma[2,1]") |> pull(q50)
        b_upper <- gamma_data |> filter(variable == "gamma[2,1]") |> pull(q97.5)
        b_lower <- gamma_data |> filter(variable == "gamma[2,1]") |> pull(q2.5)

        c_median <- gamma_data |> filter(variable == "gamma[3,1]") |> pull(q50)
        c_upper <- gamma_data |> filter(variable == "gamma[3,1]") |> pull(q97.5)
        c_lower <- gamma_data |> filter(variable == "gamma[3,1]") |> pull(q2.5)

        tau_a <- tau_summary |> filter(variable == "tau[1]") 
        tau_b <- tau_summary |> filter(variable == "tau[2]") 
        tau_c <- tau_summary |> filter(variable == "tau[3]") 

        a_slope <- gamma_data |> filter(variable == "gamma[1,2]") |> pull(q50)
        a_slope_upper <- gamma_data |> filter(variable == "gamma[1,2]") |> pull(q97.5)
        a_slope_lower <- gamma_data |> filter(variable == "gamma[1,2]") |> pull(q2.5)

        b_slope <- gamma_data |> filter(variable == "gamma[2,2]") |> pull(q50)
        b_slope_upper <- gamma_data |> filter(variable == "gamma[2,2]") |> pull(q97.5)
        b_slope_lower <- gamma_data |> filter(variable == "gamma[2,2]") |> pull(q2.5)

        c_slope <- gamma_data |> filter(variable == "gamma[3,2]") |> pull(q50)
        c_slope_upper <- gamma_data |> filter(variable == "gamma[3,2]") |> pull(q97.5)
        c_slope_lower <- gamma_data |> filter(variable == "gamma[3,2]") |> pull(q2.5)

        data.frame(
          Dependent_variable = dependent_variable,
          Division = division,
          Functional_form = functional_form,
          Wood_density = wood_density,
          Parameter = c("a", "b", "c"),
          Intercept_CI = c(
            paste0(a_median, " (", a_lower, ", ", a_upper, ")"),
            paste0(b_median, " (", b_lower, ", ", b_upper, ")"),
            paste0(c_median, " (", c_lower, ", ", c_upper, ")")
          ),
          Slope_CI = c(
            paste0(a_slope, " (", a_slope_lower, ", ", a_slope_upper, ")"),
            paste0(b_slope, " (", b_slope_lower, ", ", b_slope_upper, ")"),
            paste0(c_slope, " (", c_slope_lower, ", ", c_slope_upper, ")")
          ),
          Tau_CI = c(
            paste0(tau_a$mean, " (", tau_a$q2.5, ", ", tau_a$q97.5, ")"),
            paste0(tau_b$mean, " (", tau_b$q2.5, ", ", tau_b$q97.5, ")"),
            paste0(tau_c$mean, " (", tau_c$q2.5, ", ", tau_c$q97.5, ")")
          )
        )
      } else {
        # Special case for gamma_pl_nou_dbh without wood density
        a_median <- exp(gamma_data |> filter(variable == "gamma[1]") |> pull(q50))
        a_upper <- exp(gamma_data |> filter(variable == "gamma[1]") |> pull(q97.5))
        a_lower <- exp(gamma_data |> filter(variable == "gamma[1]") |> pull(q2.5))

        b_median <- gamma_data |> filter(variable == "gamma[2]") |> pull(q50)
        b_upper <- gamma_data |> filter(variable == "gamma[2]") |> pull(q97.5)
        b_lower <- gamma_data |> filter(variable == "gamma[2]") |> pull(q2.5)

        c_median <- gamma_data |> filter(variable == "gamma[3]") |> pull(q50)
        c_upper <- gamma_data |> filter(variable == "gamma[3]") |> pull(q97.5)
        c_lower <- gamma_data |> filter(variable == "gamma[3]") |> pull(q2.5)

        tau_a <- tau_summary |> filter(variable == "tau[1]") 
        tau_b <- tau_summary |> filter(variable == "tau[2]") 
        tau_c <- tau_summary |> filter(variable == "tau[3]") 

        data.frame(
          Dependent_variable = dependent_variable,
          Division = division,
          Functional_form = functional_form,
          Wood_density = wood_density,
          Parameter = c("a", "b", "c"),
          Intercept_CI = c(
            paste0(a_median, " (", a_lower, ", ", a_upper, ")"),
            paste0(b_median, " (", b_lower, ", ", b_upper, ")"),
            paste0(c_median, " (", c_lower, ", ", c_upper, ")")
          ),
          Slope_CI = c("-", "-", "-"),  # No slopes without wood density
          Tau_CI = c(
            paste0(tau_a$mean, " (", tau_a$q2.5, ", ", tau_a$q97.5, ")"),
            paste0(tau_b$mean, " (", tau_b$q2.5, ", ", tau_b$q97.5, ")"),
            paste0(tau_c$mean, " (", tau_c$q2.5, ", ", tau_c$q97.5, ")")
          )
        )
      }
    } else {
      # Standard cases for other models without 'c'
      if (wood_density == "With") {
        if (functional_form %in% c("Weibull", "gMM")) {
          a_median <- exp(gamma_data |> filter(variable == "gamma[1,1]") |> pull(q50))
          a_upper <- exp(gamma_data |> filter(variable == "gamma[1,1]") |> pull(q97.5))
          a_lower <- exp(gamma_data |> filter(variable == "gamma[1,1]") |> pull(q2.5))

          b_median <- gamma_data |> filter(variable == "gamma[1,2]") |> pull(q50)
          b_upper <- gamma_data |> filter(variable == "gamma[1,2]") |> pull(q97.5)
          b_lower <- gamma_data |> filter(variable == "gamma[1,2]") |> pull(q2.5)

          k_median <- gamma_data |> filter(variable == "gamma[1,3]") |> pull(q50)
          k_upper <- gamma_data |> filter(variable == "gamma[1,3]") |> pull(q97.5)
          k_lower <- gamma_data |> filter(variable == "gamma[1,3]") |> pull(q2.5)

          tau_a <- tau_summary |> filter(variable == "tau[1]") 
          tau_b <- tau_summary |> filter(variable == "tau[2]") 
          tau_k <- tau_summary |> filter(variable == "tau[3]") 

          a_slope <- exp(gamma_data |> filter(variable == "gamma[2,1]") |> pull(q50))
          a_slope_upper <- exp(gamma_data |> filter(variable == "gamma[2,1]") |> pull(q97.5))
          a_slope_lower <- exp(gamma_data |> filter(variable == "gamma[2,1]") |> pull(q2.5))

          b_slope <- gamma_data |> filter(variable == "gamma[2,2]") |> pull(q50)
          b_slope_upper <- gamma_data |> filter(variable == "gamma[2,2]") |> pull(q97.5)
          b_slope_lower <- gamma_data |> filter(variable == "gamma[2,2]") |> pull(q2.5)

          k_slope <- gamma_data |> filter(variable == "gamma[2,3]") |> pull(q50)
          k_slope_upper <- gamma_data |> filter(variable == "gamma[2,3]") |> pull(q97.5)
          k_slope_lower <- gamma_data |> filter(variable == "gamma[2,3]") |> pull(q2.5)

          data.frame(
            Dependent_variable = dependent_variable,
            Division = division,
            Functional_form = functional_form,
            Wood_density = wood_density,
            Parameter = c("a", "b", "k"),
            Intercept_CI = c(
              paste0(a_median, " (", a_lower, ", ", a_upper, ")"),
              paste0(b_median, " (", b_lower, ", ", b_upper, ")"),
              paste0(k_median, " (", k_lower, ", ", k_upper, ")")
            ),
            Slope_CI = c(
              paste0(a_slope, " (", a_slope_lower, ", ", a_slope_upper, ")"),
              paste0(b_slope, " (", b_slope_lower, ", ", b_slope_upper, ")"),
              paste0(k_slope, " (", k_slope_lower, ", ", k_slope_upper, ")")
            ),
            Tau_CI = c(
              paste0(tau_a$mean, " (", tau_a$q2.5, ", ", tau_a$q97.5, ")"),
              paste0(tau_b$mean, " (", tau_b$q2.5, ", ", tau_b$q97.5, ")"),
              paste0(tau_k$mean, " (", tau_k$q2.5, ", ", tau_k$q97.5, ")")
            )
          )

        } else if (functional_form == "Power-Law") {
          a_median <- exp(gamma_data |> filter(variable == "gamma[1,1]") |> pull(q50))
          a_upper <- exp(gamma_data |> filter(variable == "gamma[1,1]") |> pull(q97.5))
          a_lower <- exp(gamma_data |> filter(variable == "gamma[1,1]") |> pull(q2.5))

          b_median <- gamma_data |> filter(variable == "gamma[2,1]") |> pull(q50)
          b_upper <- gamma_data |> filter(variable == "gamma[2,1]") |> pull(q97.5)
          b_lower <- gamma_data |> filter(variable == "gamma[2,1]") |> pull(q2.5)

          tau_a <- tau_summary |> filter(variable == "tau[1]") 
          tau_b <- tau_summary |> filter(variable == "tau[2]") 

          a_slope <- gamma_data |> filter(variable == "gamma[1,2]") |> pull(q50)
          a_slope_upper <- gamma_data |> filter(variable == "gamma[1,2]") |> pull(q97.5)
          a_slope_lower <- gamma_data |> filter(variable == "gamma[1,2]") |> pull(q2.5)

          b_slope <- gamma_data |> filter(variable == "gamma[2,2]") |> pull(q50)
          b_slope_upper <- gamma_data |> filter(variable == "gamma[2,2]") |> pull(q97.5)
          b_slope_lower <- gamma_data |> filter(variable == "gamma[2,2]") |> pull(q2.5)

          data.frame(
            Dependent_variable = dependent_variable,
            Division = division,
            Functional_form = functional_form,
            Wood_density = wood_density,
            Parameter = c("a", "b"),
            Intercept_CI = c(
              paste0(a_median, " (", a_lower, ", ", a_upper, ")"),
              paste0(b_median, " (", b_lower, ", ", b_upper, ")")
            ),
            Slope_CI = c(
              paste0(a_slope, " (", a_slope_lower, ", ", a_slope_upper, ")"),
              paste0(b_slope, " (", b_slope_lower, ", ", b_slope_upper, ")")
            ),
            Tau_CI = c(
              paste0(tau_a$mean, " (", tau_a$q2.5, ", ", tau_a$q97.5, ")"),
              paste0(tau_b$mean, " (", tau_b$q2.5, ", ", tau_b$q97.5, ")")
            )
          )
        }

      } else {  # Without wood density
        if (functional_form %in% c("Weibull", "gMM")) {
          a_median <- exp(gamma_data |> filter(variable == "gamma[1]") |> pull(q50))
          a_upper <- exp(gamma_data |> filter(variable == "gamma[1]") |> pull(q97.5))
          a_lower <- exp(gamma_data |> filter(variable == "gamma[1]") |> pull(q2.5))

          b_median <- gamma_data |> filter(variable == "gamma[2]") |> pull(q50)
          b_upper <- gamma_data |> filter(variable == "gamma[2]") |> pull(q97.5)
          b_lower <- gamma_data |> filter(variable == "gamma[2]") |> pull(q2.5)

          k_median <- gamma_data |> filter(variable == "gamma[3]") |> pull(q50)
          k_upper <- gamma_data |> filter(variable == "gamma[3]") |> pull(q97.5)
          k_lower <- gamma_data |> filter(variable == "gamma[3]") |> pull(q2.5)

          tau_a <- tau_summary |> filter(variable == "tau[1]") 
          tau_b <- tau_summary |> filter(variable == "tau[2]") 
          tau_k <- tau_summary |> filter(variable == "tau[3]") 

          data.frame(
            Dependent_variable = dependent_variable,
            Division = division,
            Functional_form = functional_form,
            Wood_density = wood_density,
            Parameter = c("a", "b", "k"),
            Intercept_CI = c(
              paste0(a_median, " (", a_lower, ", ", a_upper, ")"),
              paste0(b_median, " (", b_lower, ", ", b_upper, ")"),
              paste0(k_median, " (", k_lower, ", ", k_upper, ")")
            ),
            Slope_CI = c("-", "-", "-"),  # No slopes without wood density
            Tau_CI = c(
              paste0(tau_a$mean, " (", tau_a$q2.5, ", ", tau_a$q97.5, ")"),
              paste0(tau_b$mean, " (", tau_b$q2.5, ", ", tau_b$q97.5, ")"),
              paste0(tau_k$mean, " (", tau_k$q2.5, ", ", tau_k$q97.5, ")")
            )
          )

        } else if (functional_form == "Power-Law") {
          a_median <- exp(gamma_data |> filter(variable == "gamma[1]") |> pull(q50))
          a_upper <- exp(gamma_data |> filter(variable == "gamma[1]") |> pull(q97.5))
          a_lower <- exp(gamma_data |> filter(variable == "gamma[1]") |> pull(q2.5))

          b_median <- gamma_data |> filter(variable == "gamma[2]") |> pull(q50)
          b_upper <- gamma_data |> filter(variable == "gamma[2]") |> pull(q97.5)
          b_lower <- gamma_data |> filter(variable == "gamma[2]") |> pull(q2.5)

          tau_a <- tau_summary |> filter(variable == "tau[1]") 
          tau_b <- tau_summary |> filter(variable == "tau[2]") 

          data.frame(
            Dependent_variable = dependent_variable,
            Division = division,
            Functional_form = functional_form,
            Wood_density = wood_density,
            Parameter = c("a", "b"),
            Intercept_CI = c(
              paste0(a_median, " (", a_lower, ", ", a_upper, ")"),
              paste0(b_median, " (", b_lower, ", ", b_upper, ")")
            ),
            Slope_CI = c("-", "-"),  # No slopes without wood density
            Tau_CI = c(
              paste0(tau_a$mean, " (", tau_a$q2.5, ", ", tau_a$q97.5, ")"),
              paste0(tau_b$mean, " (", tau_b$q2.5, ", ", tau_b$q97.5, ")")
            )
          )
        }
      }
    }
  }
# Create individual data frames for all model types
  df_list <- list(
   # Angiosperm models
    create_model_df(pl_nou_ang_h, "Tree Height", "Power-Law", "Without", "Angiosperm"),
    create_model_df(gmm_nou_ang_h, "Tree Height", "gMM", "Without", "Angiosperm"),
    create_model_df(weibull_nou_ang_h, "Tree Height", "Weibull", "Without", "Angiosperm"),
    create_model_df(weibull_ang_h, "Tree Height", "Weibull", "With", "Angiosperm"),
    create_model_df(pl_nou_ang_cr, "Crown Radius", "Power-Law", "Without", "Angiosperm"),
    create_model_df(pl_ang_cr, "Crown Radius", "Power-Law", "With", "Angiosperm"),
    create_model_df(gmm_nou_ang_cr, "Crown Radius", "gMM", "Without", "Angiosperm"),
    create_model_df(weibull_nou_ang_cr, "Crown Radius", "Weibull", "Without", "Angiosperm"),
    create_model_df(pl_nou_ang_dbh, "DBH", "Power-Law", "Without", "Angiosperm"),
    create_model_df(pl_ang_dbh, "DBH", "Power-Law", "With", "Angiosperm"),
    create_model_df(pl_nou_ang_dbh1, "DBH1", "Power-Law", "Without", "Angiosperm"),
    create_model_df(pl_nou_ang_dbh2, "DBH2", "Power-Law", "Without", "Angiosperm"),
    create_model_df(pl_nou_ang_dbh3, "DBH3", "Power-Law", "Without", "Angiosperm"),

    # Gymnosperm models
    create_model_df(pl_nou_gym_h, "Tree Height", "Power-Law", "Without", "Gymnosperm"),
    create_model_df(gmm_nou_gym_h, "Tree Height", "gMM", "Without", "Gymnosperm"),
    create_model_df(weibull_nou_gym_h, "Tree Height", "Weibull", "Without", "Gymnosperm"),
    create_model_df(weibull_gym_h, "Tree Height", "Weibull", "With", "Gymnosperm"),
    create_model_df(pl_nou_gym_cr, "Crown Radius", "Power-Law", "Without", "Gymnosperm"),
    create_model_df(gmm_nou_gym_cr, "Crown Radius", "gMM", "Without", "Gymnosperm"),
    create_model_df(gmm_gym_cr, "Crown Radius", "gMM", "With", "Gymnosperm"),
    create_model_df(weibull_nou_gym_cr, "Crown Radius", "Weibull", "Without", "Gymnosperm"),
    create_model_df(pl_nou_gym_dbh, "DBH", "Power-Law", "Without", "Gymnosperm"),
    create_model_df(pl_gym_dbh, "DBH", "Power-Law", "With", "Gymnosperm"),
    create_model_df(pl_nou_gym_dbh1, "DBH1", "Power-Law", "Without", "Gymnosperm"),
    create_model_df(pl_nou_gym_dbh2, "DBH2", "Power-Law", "Without", "Gymnosperm"),
    create_model_df(pl_nou_gym_dbh3, "DBH3", "Power-Law", "Without", "Gymnosperm")
  )
  # Combine all data frames into one
  combined_df <- bind_rows(df_list)

  # Order by Dependent_variable, Functional_form, and Wood_density
  combined_df <- combined_df |>
    arrange(Dependent_variable, Division, Functional_form, Wood_density)

  return(combined_df)
}

#=============================
# ROUNDING POSTERIOR 
#=============================

#' round_posterior_df: Rounding Numerical Values in Posterior Summary Data
#'
#' This function takes a posterior summary data frame and rounds all numerical values in the `Median_CI`
#' and `Slope_CI` columns, as well as the `SD` column, to three decimal places. It ensures consistency
#' and readability in the presentation of posterior estimates.
#'
#' @param posterior_df A data frame containing posterior summary statistics, including confidence intervals
#' (in character format) and standard deviations (in numeric format).
#'
#' @return A data frame with rounded numerical values in the `Median_CI`, `Slope_CI`, and `SD` columns.
# round_posterior_df <- function(posterior_df) {
#   posterior_df |>
#     dplyr::mutate(
#       Tau_CI = purrr::map_chr(Tau_CI, function(x) {
#         gsubfn::gsubfn("(\\d+\\.\\d+)", ~ format(round(as.numeric(.x), 3), nsmall = 3), x)
#       }),
#       Intercept_CI = purrr::map_chr(Intercept_CI, function(x) {
#         gsubfn::gsubfn("(\\d+\\.\\d+)", ~ format(round(as.numeric(.x), 3), nsmall = 3), x)
#       }),
#       Slope_CI = purrr::map_chr(Slope_CI, function(x) {
#         if (x == "-") return(x) # Skip rounding for non-numeric entries
#         gsubfn::gsubfn("(\\d+\\.\\d+)", ~ format(round(as.numeric(.x), 3), nsmall = 3), x)
#       })
#     )
# }

round_posterior_df <- function(posterior_df) {
  round_number <- function(x) {
    num <- as.numeric(x)
    if (is.na(num)) return(x)
    if (abs(num) < 0.01) {
      format(round(num, 3), nsmall = 3)
    } else {
      format(round(num, 2), nsmall = 2)
    }
  }
  
  posterior_df |>
    dplyr::mutate(
      Tau_CI = purrr::map_chr(Tau_CI, function(x) {
        gsubfn::gsubfn("(\\d+\\.\\d+)", ~ round_number(.x), x)
      }),
      Intercept_CI = purrr::map_chr(Intercept_CI, function(x) {
        gsubfn::gsubfn("(\\d+\\.\\d+)", ~ round_number(.x), x)
      }),
      Slope_CI = purrr::map_chr(Slope_CI, function(x) {
        if (x == "-") return(x)
        gsubfn::gsubfn("(\\d+\\.\\d+)", ~ round_number(.x), x)
      })
    )
}


#==============================
# BEST MODELS' POSTERIOR TABLE 
#==============================

#' Extract Best Models Based on Specific Criteria
#'
#' This function filters and extracts the best models for tree height, crown radius, 
#' and DBH based on the specified criteria for angiosperms and gymnosperms.
#'
#' @param posterior_df A dataframe containing the posterior summary of models.
#' @return A dataframe of the best models filtered by the criteria.
#' @export
# extract_best_models <- function(posterior_rounded_df) {
#   # Filter for Tree Height with Weibull Model for Angiosperm and Wood Density "Without"
#     wb_h_ang_nou_wd_df <- posterior_rounded_df |>
#     dplyr::filter(Dependent_variable == "Tree Height" &
#                   Division == "Angiosperm" &
#                   Functional_form == "Weibull" &
#                   Wood_density == "Without")

#   # Filter for Tree Height with Weibull Model for Angiosperm and Wood Density "With"
#   wb_h_ang_wd_df <- posterior_rounded_df |>
#     dplyr::filter(Dependent_variable == "Tree Height" &
#                   Division == "Angiosperm" &
#                   Functional_form == "Weibull" &
#                   Wood_density == "With")
  
#   # Filter for Tree Height with Weibull Model for Gymnosperm and Wood Density "Without"
#   wb_h_gym_nou_wd_df <- posterior_rounded_df |>
#     dplyr::filter(Dependent_variable == "Tree Height" &
#                   Division == "Gymnosperm" &
#                   Functional_form == "Weibull" &
#                   Wood_density == "Without")

#   # Filter for Tree Height with Weibull Model for Gymnosperm and Wood Density "With"
#   wb_h_gym_wd_df <- posterior_rounded_df |>
#     dplyr::filter(Dependent_variable == "Tree Height" &
#                   Division == "Gymnosperm" &
#                   Functional_form == "Weibull" &
#                   Wood_density == "With")
  
#   # Filter for Crown Radius with Power-Law Model for Angiosperm and Wood Density "Without"
#   pl_cr_ang_nou_wd_df <- posterior_rounded_df |>
#     dplyr::filter(Dependent_variable == "Crown Radius" &
#                   Division == "Angiosperm" &
#                   Functional_form == "Power-Law" &
#                   Wood_density == "Without")

#   # Filter for Crown Radius with Power-Law Model for Angiosperm and Wood Density "With"
#   pl_cr_ang_wd_df <- posterior_rounded_df |>
#     dplyr::filter(Dependent_variable == "Crown Radius" &
#                   Division == "Angiosperm" &
#                   Functional_form == "Power-Law" &
#                   Wood_density == "With")
  
#   # Filter for Crown Radius with gMM Model for Gymnosperm and Wood Density "Without"
#   gmm_cr_gym_nou_wd_df <- posterior_rounded_df |>
#     dplyr::filter(Dependent_variable == "Crown Radius" &
#                   Division == "Gymnosperm" &
#                   Functional_form == "gMM" &
#                   Wood_density == "Without")

#   # Filter for Crown Radius with gMM Model for Gymnosperm and Wood Density "With"
#   gmm_cr_gym_wd_df <- posterior_rounded_df |>
#     dplyr::filter(Dependent_variable == "Crown Radius" &
#                   Division == "Gymnosperm" &
#                   Functional_form == "gMM" &
#                   Wood_density == "With")
  
#   # Filter for DBH with Power-Law Model for Angiosperm and Wood Density "Without"
#   pl_dbh_ang_nou_wd_df <- posterior_rounded_df |>
#     dplyr::filter(Dependent_variable == "DBH" &
#                   Division == "Angiosperm" &
#                   Functional_form == "Power-Law" &
#                   Wood_density == "Without")
  
#   # Filter for DBH with Power-Law Model for Angiosperm and Wood Density "With"
#   pl_dbh_ang_wd_df <- posterior_rounded_df |>
#     dplyr::filter(Dependent_variable == "DBH" &
#                   Division == "Angiosperm" &
#                   Functional_form == "Power-Law" &
#                   Wood_density == "With")
  
#   # Filter for DBH with Power-Law Model for Gymnosperm and Wood Density "Without"
#   pl_dbh_gym_nou_wd_df <- posterior_rounded_df |>
#     dplyr::filter(Dependent_variable == "DBH" &
#                   Division == "Gymnosperm" &
#                   Functional_form == "Power-Law" &
#                   Wood_density == "Without")

#   # Filter for DBH with Power-Law Model for Gymnosperm and Wood Density "With"
#   pl_dbh_gym_wd_df <- posterior_rounded_df |>
#     dplyr::filter(Dependent_variable == "DBH" &
#                   Division == "Gymnosperm" &
#                   Functional_form == "Power-Law" &
#                   Wood_density == "With")
  

#   # Combine all filtered data frames
#   best_model_tbl <- dplyr::bind_rows(
#       wb_h_ang_nou_wd_df,
#       wb_h_ang_wd_df,
#       wb_h_gym_nou_wd_df,
#       wb_h_gym_wd_df,
#       pl_cr_ang_nou_wd_df,
#       pl_cr_ang_wd_df,
#       gmm_cr_gym_nou_wd_df,
#       gmm_cr_gym_wd_df,
#       pl_dbh_ang_nou_wd_df,
#       pl_dbh_ang_wd_df,
#       pl_dbh_gym_nou_wd_df,
#       pl_dbh_gym_wd_df
#     )
  
#   return(best_model_tbl)
# }

extract_best_models <- function(posterior_df) {
  # Filter for Tree Height with Weibull Model for Angiosperm and Wood Density "Without"
    wb_h_ang_nou_wd_df <- posterior_df |>
    dplyr::filter(Dependent_variable == "Tree Height" &
                  Division == "Angiosperm" &
                  Functional_form == "Weibull" &
                  Wood_density == "Without")

  # Filter for Tree Height with Weibull Model for Angiosperm and Wood Density "With"
  wb_h_ang_wd_df <- posterior_df |>
    dplyr::filter(Dependent_variable == "Tree Height" &
                  Division == "Angiosperm" &
                  Functional_form == "Weibull" &
                  Wood_density == "With")
  
  # Filter for Tree Height with Weibull Model for Gymnosperm and Wood Density "Without"
  wb_h_gym_nou_wd_df <- posterior_df |>
    dplyr::filter(Dependent_variable == "Tree Height" &
                  Division == "Gymnosperm" &
                  Functional_form == "Weibull" &
                  Wood_density == "Without")

  # Filter for Tree Height with Weibull Model for Gymnosperm and Wood Density "With"
  wb_h_gym_wd_df <- posterior_df |>
    dplyr::filter(Dependent_variable == "Tree Height" &
                  Division == "Gymnosperm" &
                  Functional_form == "Weibull" &
                  Wood_density == "With")
  
  # Filter for Crown Radius with Power-Law Model for Angiosperm and Wood Density "Without"
  pl_cr_ang_nou_wd_df <- posterior_df |>
    dplyr::filter(Dependent_variable == "Crown Radius" &
                  Division == "Angiosperm" &
                  Functional_form == "Power-Law" &
                  Wood_density == "Without")

  # Filter for Crown Radius with Power-Law Model for Angiosperm and Wood Density "With"
  pl_cr_ang_wd_df <- posterior_df |>
    dplyr::filter(Dependent_variable == "Crown Radius" &
                  Division == "Angiosperm" &
                  Functional_form == "Power-Law" &
                  Wood_density == "With")
  
  # Filter for Crown Radius with gMM Model for Gymnosperm and Wood Density "Without"
  gmm_cr_gym_nou_wd_df <- posterior_df |>
    dplyr::filter(Dependent_variable == "Crown Radius" &
                  Division == "Gymnosperm" &
                  Functional_form == "gMM" &
                  Wood_density == "Without")

  # Filter for Crown Radius with gMM Model for Gymnosperm and Wood Density "With"
  gmm_cr_gym_wd_df <- posterior_df |>
    dplyr::filter(Dependent_variable == "Crown Radius" &
                  Division == "Gymnosperm" &
                  Functional_form == "gMM" &
                  Wood_density == "With")
  
  # Filter for DBH with Power-Law Model for Angiosperm and Wood Density "Without"
  pl_dbh_ang_nou_wd_df <- posterior_df |>
    dplyr::filter(Dependent_variable == "DBH" &
                  Division == "Angiosperm" &
                  Functional_form == "Power-Law" &
                  Wood_density == "Without")
  
  # Filter for DBH with Power-Law Model for Angiosperm and Wood Density "With"
  pl_dbh_ang_wd_df <- posterior_df |>
    dplyr::filter(Dependent_variable == "DBH" &
                  Division == "Angiosperm" &
                  Functional_form == "Power-Law" &
                  Wood_density == "With")
  
  # Filter for DBH with Power-Law Model for Gymnosperm and Wood Density "Without"
  pl_dbh_gym_nou_wd_df <- posterior_df |>
    dplyr::filter(Dependent_variable == "DBH" &
                  Division == "Gymnosperm" &
                  Functional_form == "Power-Law" &
                  Wood_density == "Without")

  # Filter for DBH with Power-Law Model for Gymnosperm and Wood Density "With"
  pl_dbh_gym_wd_df <- posterior_df |>
    dplyr::filter(Dependent_variable == "DBH" &
                  Division == "Gymnosperm" &
                  Functional_form == "Power-Law" &
                  Wood_density == "With")
  

  # Combine all filtered data frames
  best_model_tbl <- dplyr::bind_rows(
      wb_h_ang_nou_wd_df,
      wb_h_ang_wd_df,
      wb_h_gym_nou_wd_df,
      wb_h_gym_wd_df,
      pl_cr_ang_nou_wd_df,
      pl_cr_ang_wd_df,
      gmm_cr_gym_nou_wd_df,
      gmm_cr_gym_wd_df,
      pl_dbh_ang_nou_wd_df,
      pl_dbh_ang_wd_df,
      pl_dbh_gym_nou_wd_df,
      pl_dbh_gym_wd_df
    )
  
  return(best_model_tbl)
}

format_ci <- function(ci_string, sig_fig = 3) {
  if (is.na(ci_string) || ci_string == "-") return(ci_string)
  
  # Extract numbers
  numbers <- str_extract_all(ci_string, "\\-?\\d+\\.?\\d*")[[1]]
  numbers <- as.numeric(numbers)
  
  # Format each number to avoid scientific notation and use significant figures
  formatted <- sapply(numbers, function(x) {
    if (x == 0) return("0")
    digits <- sig_fig - floor(log10(abs(x))) - 1
    sprintf(paste0("%.", max(0, digits), "f"), x)
  })
  
  # Rebuild string
  if (length(formatted) == 1) {
    return(formatted[1])
  } else if (length(formatted) == 3) {
    return(paste0(formatted[1], " (", formatted[2], ", ", formatted[3], ")"))
  } else {
    return(ci_string)  # fallback
  }
}


format_ci_sp <- function(ci_string, sig_fig = 3) {
  if (is.na(ci_string) || ci_string == "-") return(ci_string)

  # Extract numeric parts
  numbers <- stringr::str_extract_all(ci_string, "\\-?\\d+\\.?\\d*")[[1]]
  numbers <- as.numeric(numbers)

  # Format to significant figures
  formatted <- sapply(numbers, function(x) {
    if (x == 0) return("0")
    digits <- sig_fig - floor(log10(abs(x))) - 1
    sprintf(paste0("%.", max(0, digits), "f"), x)
  })

  # Return formatted string
  if (length(formatted) == 1) {
    return(formatted[1])
  } else if (length(formatted) == 3) {
    return(paste0(formatted[1], " (", formatted[2], ", ", formatted[3], ")"))
  } else {
    return(ci_string)
  }
}

format_sp_posterior_df <- function(df) {
  # Format 'a', 'b', 'k' or 'c' columns based on their availability in the data frame
  if ("k" %in% colnames(df)) {
    df <- df %>%
      dplyr::mutate(
        a = sapply(a, function(x) ifelse(is.na(x), NA, format_ci_sp(x))),
        b = sapply(b, function(x) ifelse(is.na(x), NA, format_ci_sp(x))),
        k = sapply(k, function(x) ifelse(is.na(x), NA, format_ci_sp(x)))
      )
  }
  if ("c" %in% colnames(df)) {
    df <- df %>%
      dplyr::mutate(
        a = sapply(a, function(x) ifelse(is.na(x), NA, format_ci_sp(x))),
        b = sapply(b, function(x) ifelse(is.na(x), NA, format_ci_sp(x))),
        c = sapply(c, function(x) ifelse(is.na(x), NA, format_ci_sp(x)))
      )
  }
  return(df)
}

#=============================
# TABLE 1, POSTERIOR DF
#=============================



#===============================================
# SPECIES POSTERIOR DATAFRAME FOR AGB ESTIMATION
#===============================================
generate_sp_posterior_agb_df <- function(
  # Angiosperm datasets
  fit_nlr_nou_summary_weibull_ang_h, fit_lr_nou_summary_pl_ang_h,
  fit_lr_nou_summary_pl_ang_dbh, fit_lr_nou_summary_pl_ang_dbh1, 
  fit_lr_nou_summary_pl_ang_dbh2, fit_lr_nou_summary_pl_ang_dbh3,
  tallo_reduced_nlr_df_ang_h, tallo_reduced_lr_df_ang_h,
  tallo_reduced_lr_df_ang_dbh, tallo_reduced_lr_df_ang_dbh1,
  tallo_reduced_lr_df_ang_dbh2, tallo_reduced_lr_df_ang_dbh3,
  stan_data_nlr_ang_h, stan_data_lr_ang_h,
  stan_data_lr_ang_dbh, stan_data_lr_ang_dbh1,
  stan_data_lr_ang_dbh2, stan_data_lr_ang_dbh3,
  
  # Gymnosperm datasets
  fit_nlr_nou_summary_weibull_gym_h, fit_lr_nou_summary_pl_gym_h,
  fit_lr_nou_summary_pl_gym_dbh, fit_lr_nou_summary_pl_gym_dbh1, 
  fit_lr_nou_summary_pl_gym_dbh2, fit_lr_nou_summary_pl_gym_dbh3,
  tallo_reduced_nlr_df_gym_h, tallo_reduced_lr_df_gym_h,
  tallo_reduced_lr_df_gym_dbh, tallo_reduced_lr_df_gym_dbh1,
  tallo_reduced_lr_df_gym_dbh2, tallo_reduced_lr_df_gym_dbh3,
  stan_data_nlr_gym_h, stan_data_lr_gym_h,
  stan_data_lr_gym_dbh, stan_data_lr_gym_dbh1,
  stan_data_lr_gym_dbh2, stan_data_lr_gym_dbh3
) {
 # Helper function to filter species for DBH datset
  filter_species <- function(dataset, division) {
      dataset |> 
        filter(division == division, !is.na(dbh), !is.na(cr), !is.na(h)) |> 
        group_by(sp) |> 
        summarise(count = n()) |> 
        filter(count >= 20) |> 
        pull(sp)
    }
  # Helper function to extract posterior summaries and prepare species-specific output
    create_species_model_df <- function(beta_data, dependent_variable, functional_form, sp_list) {
      if (functional_form == "Weibull" || functional_form == "gMM") {
        a_median <- exp(beta_data |> filter(str_detect(variable, "beta\\[\\d+,1\\]")) |> pull(q50))
        a_upper <- exp(beta_data |> filter(str_detect(variable, "beta\\[\\d+,1\\]")) |> pull(q97.5))
        a_lower <- exp(beta_data |> filter(str_detect(variable, "beta\\[\\d+,1\\]")) |> pull(q2.5))

        b_median <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,2\\]")) |> pull(q50)
        b_upper <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,2\\]")) |> pull(q97.5)
        b_lower <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,2\\]")) |> pull(q2.5)

        k_median <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,3\\]")) |> pull(q50)
        k_upper <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,3\\]")) |> pull(q97.5)
        k_lower <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,3\\]")) |> pull(q2.5)

        ci_95 <- data.frame(
          sp = sp_list,
          a = paste0(format(a_median, scientific = FALSE), " (", format(a_lower, scientific = FALSE), ", ", format(a_upper, scientific = FALSE), ")"),
          b = paste0(format(b_median, scientific = FALSE), " (", format(b_lower, scientific = FALSE), ", ", format(b_upper, scientific = FALSE), ")"),
          k = paste0(format(k_median, scientific = FALSE), " (", format(k_lower, scientific = FALSE), ", ", format(k_upper, scientific = FALSE), ")"),
          Dependent_variable = dependent_variable,
          stringsAsFactors = FALSE
        )
      } else if (functional_form == "Power-Law") {
        a_median <- exp(beta_data |> filter(str_detect(variable, "beta\\[1,\\d+\\]")) |> pull(q50))
        a_upper <- exp(beta_data |> filter(str_detect(variable, "beta\\[1,\\d+\\]")) |> pull(q97.5))
        a_lower <- exp(beta_data |> filter(str_detect(variable, "beta\\[1,\\d+\\]")) |> pull(q2.5))

        b_median <- beta_data |> filter(str_detect(variable, "beta\\[2,\\d+\\]")) |> pull(q50)
        b_upper <- beta_data |> filter(str_detect(variable, "beta\\[2,\\d+\\]")) |> pull(q97.5)
        b_lower <- beta_data |> filter(str_detect(variable, "beta\\[2,\\d+\\]")) |> pull(q2.5)

        c_median <- beta_data |> filter(str_detect(variable, "beta\\[3,\\d+\\]")) |> pull(q50)
        c_upper <- beta_data |> filter(str_detect(variable, "beta\\[3,\\d+\\]")) |> pull(q97.5)
        c_lower <- beta_data |> filter(str_detect(variable, "beta\\[3,\\d+\\]")) |> pull(q2.5)

      ci_95 <- data.frame(
            sp = sp_list,
            a = paste0(format(a_median, scientific = FALSE), " (", format(a_lower, scientific = FALSE), ", ", format(a_upper, scientific = FALSE), ")"),
            b = paste0(format(b_median, scientific = FALSE), " (", format(b_lower, scientific = FALSE), ", ", format(b_upper, scientific = FALSE), ")"),
            c = paste0(format(c_median, scientific = FALSE), " (", format(c_lower, scientific = FALSE), ", ", format(c_upper, scientific = FALSE), ")"),
            Dependent_variable = dependent_variable,
            stringsAsFactors = FALSE
        )
      }

      return(ci_95)
    }
  # Function to process datasets and align species IDs with Stan data
  process_dataset <- function(dataset, stan_data, division_filter) {
  # Filter the dataset by division
  dataset <- dataset |> filter(division == division_filter)
      
   # Prepare species data
    sp_data <- dataset |>
      group_by(sp) |>
      summarise(wd = mean(wd, na.rm = TRUE)) |>
      arrange(sp) |>
      mutate(sp_id = row_number()) |>
      mutate(wd_s = scale(wd) |> as.numeric()) |>
      dplyr::select(sp_id, sp, wd, wd_s)
    
    if (!all(c("jj", "u") %in% names(stan_data))) {
      stop("Stan data must contain 'jj' and 'u'.")
    }
    
    # Handle `nlr` structure: `u` has rows indexed by `jj`
    if (nrow(stan_data$u) >= max(stan_data$jj, na.rm = TRUE)) {
      stan_processed <- data.frame(
        sp_id = as.factor(stan_data$jj),
        wd_s = stan_data$u[stan_data$jj, 2]  # Extract `wd_s` for species indices
      )
    } 
    # Handle `lr` structure: `u` has two rows (variables) and columns for species
    else if (nrow(stan_data$u) == 2) {
      stan_processed <- data.frame(
        sp_id = 1:ncol(stan_data$u),
        wd_s = stan_data$u[2, ]  # Use the second row for `wd_s`
      )
    } 
    # Error handling for unexpected structures
    else {
      stop("Unexpected structure for 'u' in Stan data.")
    }
    
    # Add species names and calculate tree counts
    stan_processed <- stan_processed |>
      group_by(sp_id, wd_s) |>
      summarise(tree_count = n(), .groups = "drop") |>
      mutate(sp_id = as.integer(as.character(sp_id))) |>
      left_join(sp_data, by = "sp_id")
    
    return(list(sp_data = sp_data, stan_processed = stan_processed))
  }

  
  # Filter species for all datasets
  ang_dbh_species <- filter_species(tallo_reduced_lr_df_ang_dbh, "Angiosperm")
  gym_dbh_species <- filter_species(tallo_reduced_lr_df_gym_dbh, "Gymnosperm")
  
  ang_dbh1_species <- filter_species(tallo_reduced_lr_df_ang_dbh1, "Angiosperm")
  gym_dbh1_species <- filter_species(tallo_reduced_lr_df_gym_dbh1, "Gymnosperm")
  
  ang_dbh2_species <- filter_species(tallo_reduced_lr_df_ang_dbh2, "Angiosperm")
  gym_dbh2_species <- filter_species(tallo_reduced_lr_df_gym_dbh2, "Gymnosperm")
  
  ang_dbh3_species <- filter_species(tallo_reduced_lr_df_ang_dbh3, "Angiosperm")
  gym_dbh3_species <- filter_species(tallo_reduced_lr_df_gym_dbh3, "Gymnosperm")
  
  # Process datasets for Angiosperm
  h_wb_ang <- process_dataset(tallo_reduced_nlr_df_ang_h, stan_data_nlr_ang_h, "Angiosperm")
  h_pl_ang <- process_dataset(tallo_reduced_lr_df_ang_h, stan_data_lr_ang_h, "Angiosperm")
  dbh_ang <- process_dataset(tallo_reduced_lr_df_ang_dbh |> filter(sp %in% ang_dbh_species), stan_data_lr_ang_dbh, "Angiosperm")
  dbh1_ang <- process_dataset(tallo_reduced_lr_df_ang_dbh1 |> filter(sp %in% ang_dbh1_species), stan_data_lr_ang_dbh1, "Angiosperm")
  dbh2_ang <- process_dataset(tallo_reduced_lr_df_ang_dbh2 |> filter(sp %in% ang_dbh2_species), stan_data_lr_ang_dbh2, "Angiosperm")
  dbh3_ang <- process_dataset(tallo_reduced_lr_df_ang_dbh3 |> filter(sp %in% ang_dbh3_species), stan_data_lr_ang_dbh3, "Angiosperm")
  
  # Process datasets for Gymnosperm
  h_wb_gym <- process_dataset(tallo_reduced_nlr_df_gym_h, stan_data_nlr_gym_h, "Gymnosperm")
  h_pl_gym <- process_dataset(tallo_reduced_lr_df_gym_h, stan_data_lr_gym_h, "Gymnosperm")
  dbh_gym <- process_dataset(tallo_reduced_lr_df_gym_dbh |> filter(sp %in% gym_dbh_species), stan_data_lr_gym_dbh, "Gymnosperm")
  dbh1_gym <- process_dataset(tallo_reduced_lr_df_gym_dbh1 |> filter(sp %in% gym_dbh1_species), stan_data_lr_gym_dbh1, "Gymnosperm")
  dbh2_gym <- process_dataset(tallo_reduced_lr_df_gym_dbh2 |> filter(sp %in% gym_dbh2_species), stan_data_lr_gym_dbh2, "Gymnosperm")
  dbh3_gym <- process_dataset(tallo_reduced_lr_df_gym_dbh3 |> filter(sp %in% gym_dbh3_species), stan_data_lr_gym_dbh3, "Gymnosperm")
  
  # Generate posterior dataframes
  h_wb_df_ang <- create_species_model_df(fit_nlr_nou_summary_weibull_ang_h, "Tree Height", "Weibull", h_wb_ang$sp_data$sp)
  h_wb_df_gym <- create_species_model_df(fit_nlr_nou_summary_weibull_gym_h, "Tree Height", "Weibull", h_wb_gym$sp_data$sp)
  h_pl_df_ang <- create_species_model_df(fit_lr_nou_summary_pl_ang_h, "Tree Height", "Power-Law", h_pl_ang$sp_data$sp)
  h_pl_df_gym <- create_species_model_df(fit_lr_nou_summary_pl_gym_h, "Tree Height", "Power-Law", h_pl_gym$sp_data$sp)
  dbh_df_ang <- create_species_model_df(fit_lr_nou_summary_pl_ang_dbh, "DBH", "Power-Law", dbh_ang$sp_data$sp)
  dbh_df_gym <- create_species_model_df(fit_lr_nou_summary_pl_gym_dbh, "DBH", "Power-Law", dbh_gym$sp_data$sp)
  dbh1_df_ang <- create_species_model_df(fit_lr_nou_summary_pl_ang_dbh1, "DBH1", "Power-Law", dbh1_ang$sp_data$sp)
  dbh1_df_gym <- create_species_model_df(fit_lr_nou_summary_pl_gym_dbh1, "DBH1", "Power-Law", dbh1_gym$sp_data$sp)
  dbh2_df_ang <- create_species_model_df(fit_lr_nou_summary_pl_ang_dbh2, "DBH2", "Power-Law", dbh2_ang$sp_data$sp)
  dbh2_df_gym <- create_species_model_df(fit_lr_nou_summary_pl_gym_dbh2, "DBH2", "Power-Law", dbh2_gym$sp_data$sp)
  dbh3_df_ang <- create_species_model_df(fit_lr_nou_summary_pl_ang_dbh3, "DBH3", "Power-Law", dbh3_ang$sp_data$sp)
  dbh3_df_gym <- create_species_model_df(fit_lr_nou_summary_pl_gym_dbh3, "DBH3", "Power-Law", dbh3_gym$sp_data$sp)
  
  # Combine results
  sp_posterior_agb_df <- bind_rows(
    h_wb_df_ang |> mutate(Division = "Angiosperm"),
    h_wb_df_gym |> mutate(Division = "Gymnosperm"),
    h_pl_df_ang |> mutate(Division = "Angiosperm"),
    h_pl_df_gym |> mutate(Division = "Gymnosperm"),
    dbh_df_ang |> mutate(Division = "Angiosperm"),
    dbh_df_gym |> mutate(Division = "Gymnosperm"),
    dbh1_df_ang |> mutate(Division = "Angiosperm"),
    dbh1_df_gym |> mutate(Division = "Gymnosperm"),
    dbh2_df_ang |> mutate(Division = "Angiosperm"),
    dbh2_df_gym |> mutate(Division = "Gymnosperm"),
    dbh3_df_ang |> mutate(Division = "Angiosperm"),
    dbh3_df_gym |> mutate(Division = "Gymnosperm")
  )
  
  # Rename columns to remove `sp.` prefix
  colnames(sp_posterior_agb_df) <- gsub("^sp\\.", "", colnames(sp_posterior_agb_df))
  
  return(sp_posterior_agb_df)
}


#==========================================================
# SPECIES POSTERIOR DATAFRAME OF THE BEST PREDICTIVE MODELS
#==========================================================
generate_sp_posterior_df<- function(
  fit_nlr_nou_summary_weibull_ang_h, fit_nlr_nou_summary_weibull_gym_h,
  fit_lr_nou_summary_pl_ang_cr, fit_nlr_nou_summary_gmm_gym_cr,
  fit_lr_nou_summary_pl_ang_dbh, fit_lr_nou_summary_pl_gym_dbh,
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
) {

 # Helper function to filter species for DBH datset
  filter_species <- function(dataset, division) {
      dataset |> 
        filter(division == division, !is.na(dbh), !is.na(cr), !is.na(h)) |> 
        group_by(sp) |> 
        summarise(count = n()) |> 
        filter(count >= 20) |> 
        pull(sp)
    }
  # Helper function to extract posterior summaries and prepare species-specific output
    create_species_model_df <- function(beta_data, dependent_variable, functional_form, sp_list) {
      if (functional_form == "Weibull" || functional_form == "gMM") {
        a_median <- exp(beta_data |> filter(str_detect(variable, "beta\\[\\d+,1\\]")) |> pull(q50))
        a_upper <- exp(beta_data |> filter(str_detect(variable, "beta\\[\\d+,1\\]")) |> pull(q97.5))
        a_lower <- exp(beta_data |> filter(str_detect(variable, "beta\\[\\d+,1\\]")) |> pull(q2.5))

        b_median <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,2\\]")) |> pull(q50)
        b_upper <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,2\\]")) |> pull(q97.5)
        b_lower <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,2\\]")) |> pull(q2.5)

        k_median <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,3\\]")) |> pull(q50)
        k_upper <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,3\\]")) |> pull(q97.5)
        k_lower <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,3\\]")) |> pull(q2.5)

        ci_95 <- data.frame(
          sp = sp_list,
          a = paste0(a_median, " (", a_lower, ", ", a_upper, ")"),
          b = paste0(b_median, " (", b_lower, ", ", b_upper, ")"),
          k = paste0(k_median, " (", k_lower, ", ", k_upper, ")"),
          Dependent_variable = dependent_variable,
          stringsAsFactors = FALSE
        )
      } else if (functional_form == "Power-Law") {
        a_median <- exp(beta_data |> filter(str_detect(variable, "beta\\[1,\\d+\\]")) |> pull(q50))
        a_upper <- exp(beta_data |> filter(str_detect(variable, "beta\\[1,\\d+\\]")) |> pull(q97.5))
        a_lower <- exp(beta_data |> filter(str_detect(variable, "beta\\[1,\\d+\\]")) |> pull(q2.5))

        b_median <- beta_data |> filter(str_detect(variable, "beta\\[2,\\d+\\]")) |> pull(q50)
        b_upper <- beta_data |> filter(str_detect(variable, "beta\\[2,\\d+\\]")) |> pull(q97.5)
        b_lower <- beta_data |> filter(str_detect(variable, "beta\\[2,\\d+\\]")) |> pull(q2.5)

        c_median <- beta_data |> filter(str_detect(variable, "beta\\[3,\\d+\\]")) |> pull(q50)
        c_upper <- beta_data |> filter(str_detect(variable, "beta\\[3,\\d+\\]")) |> pull(q97.5)
        c_lower <- beta_data |> filter(str_detect(variable, "beta\\[3,\\d+\\]")) |> pull(q2.5)

        ci_95 <- data.frame(
          sp = sp_list,
          a = paste0(a_median, " (", a_lower, ", ", a_upper, ")"),
          b = paste0(b_median, " (", b_lower, ", ", b_upper, ")"),
          c = paste0(c_median, " (", c_lower, ", ", c_upper, ")"),
          Dependent_variable = dependent_variable,
          stringsAsFactors = FALSE
        )
      }

      return(ci_95)
    }
  # Function to process datasets and align species IDs with Stan data
  process_dataset <- function(dataset, stan_data, division_filter) {
  # Filter the dataset by division
  dataset <- dataset |> filter(division == division_filter)
      
   # Prepare species data
    sp_data <- dataset |>
      group_by(sp) |>
      summarise(wd = mean(wd, na.rm = TRUE)) |>
      arrange(sp) |>
      mutate(sp_id = row_number()) |>
      mutate(wd_s = scale(wd) |> as.numeric()) |>
      dplyr::select(sp_id, sp, wd, wd_s)
    
    if (!all(c("jj", "u") %in% names(stan_data))) {
      stop("Stan data must contain 'jj' and 'u'.")
    }
    
    # Handle `nlr` structure: `u` has rows indexed by `jj`
    if (nrow(stan_data$u) >= max(stan_data$jj, na.rm = TRUE)) {
      stan_processed <- data.frame(
        sp_id = as.factor(stan_data$jj),
        wd_s = stan_data$u[stan_data$jj, 2]  # Extract `wd_s` for species indices
      )
    } 
    # Handle `lr` structure: `u` has two rows (variables) and columns for species
    else if (nrow(stan_data$u) == 2) {
      stan_processed <- data.frame(
        sp_id = 1:ncol(stan_data$u),
        wd_s = stan_data$u[2, ]  # Use the second row for `wd_s`
      )
    } 
    # Error handling for unexpected structures
    else {
      stop("Unexpected structure for 'u' in Stan data.")
    }
    
    # Add species names and calculate tree counts
    stan_processed <- stan_processed |>
      group_by(sp_id, wd_s) |>
      summarise(tree_count = n(), .groups = "drop") |>
      mutate(sp_id = as.integer(as.character(sp_id))) |>
      left_join(sp_data, by = "sp_id")
    
    return(list(sp_data = sp_data, stan_processed = stan_processed))
  }

  # DBH allometry data preparing
  ang_dbh_species <- filter_species(tallo_reduced_lr_df_ang_dbh, "Angiosperm")
  gym_dbh_species <- filter_species(tallo_reduced_lr_df_ang_dbh, "Gymnosperm")
  
  # Process datasets for Angiosperm and Gymnosperm
  h_ang <- process_dataset(tallo_reduced_nlr_df_ang_h, stan_data_nlr_ang_h, "Angiosperm")
  h_gym <- process_dataset(tallo_reduced_nlr_df_gym_h, stan_data_nlr_gym_h, "Gymnosperm")
  cr_ang <- process_dataset(tallo_reduced_lr_df_ang_cr, stan_data_lr_ang_cr, "Angiosperm")
  cr_gym <- process_dataset(tallo_reduced_nlr_df_gym_cr, stan_data_nlr_gym_cr, "Gymnosperm")
  dbh_ang <- process_dataset(tallo_reduced_lr_df_ang_dbh |> filter(sp %in% ang_dbh_species), stan_data_lr_ang_dbh, "Angiosperm")
  dbh_gym <- process_dataset(tallo_reduced_lr_df_gym_dbh |> filter(sp %in% gym_dbh_species), stan_data_lr_gym_dbh, "Gymnosperm")

  # Generate posterior dataframes
  h_df_ang <- create_species_model_df(fit_nlr_nou_summary_weibull_ang_h, "Tree Height", "Weibull", h_ang$sp)
  h_df_gym <- create_species_model_df(fit_nlr_nou_summary_weibull_gym_h, "Tree Height", "Weibull", h_gym$sp)
  cr_df_ang <- create_species_model_df(fit_lr_nou_summary_pl_ang_cr, "Crown Radius", "Power-Law", cr_ang$sp)
  cr_df_gym <- create_species_model_df(fit_nlr_nou_summary_gmm_gym_cr, "Crown Radius", "gMM", cr_gym$sp)
  dbh_df_ang <- create_species_model_df(fit_lr_nou_summary_pl_ang_dbh, "DBH", "Power-Law", dbh_ang$sp)
  dbh_df_gym <- create_species_model_df(fit_lr_nou_summary_pl_gym_dbh, "DBH", "Power-Law", dbh_gym$sp)

  # Combine results
  sp_posterior_h_df <- bind_rows(
    h_df_ang |> mutate(Division = "Angiosperm"),
    h_df_gym |> mutate(Division = "Gymnosperm")
  )
  sp_posterior_cr_df <- bind_rows(
    cr_df_ang |> mutate(Division = "Angiosperm"),
    cr_df_gym |> mutate(Division = "Gymnosperm")
  )
  sp_posterior_dbh_df <- bind_rows(
    dbh_df_ang |> mutate(Division = "Angiosperm"),
    dbh_df_gym |> mutate(Division = "Gymnosperm")
  )
  # Rename columns to remove `sp.` prefix
  colnames(sp_posterior_h_df) <- gsub("^sp\\.", "", colnames(sp_posterior_h_df))
  sp_posterior_h_df <- sp_posterior_h_df |> 
    dplyr::rename(Species = sp) |> 
    dplyr::select(
      Dependent_variable, Species, Division, a, b, k)

  colnames(sp_posterior_cr_df) <- gsub("^sp\\.", "", colnames(sp_posterior_cr_df))
  sp_posterior_cr_df <- sp_posterior_cr_df |>
    dplyr::rename(Species = sp) |> 
    dplyr::select(
      Dependent_variable, Species, Division, a, b, k)
      
  colnames(sp_posterior_dbh_df) <- gsub("^sp\\.", "", colnames(sp_posterior_dbh_df))
  sp_posterior_dbh_df <- sp_posterior_dbh_df |>
    dplyr::rename(Species = sp) |> 
    dplyr::select(
      Dependent_variable, Species, Division, a, b, c)

  return(list(
    sp_posterior_h_df = sp_posterior_h_df,
    sp_posterior_cr_df = sp_posterior_cr_df,
    sp_posterior_dbh_df = sp_posterior_dbh_df
  ))
}

# generate_sp_posterior_df <- function(
#   fit_nlr_nou_summary_weibull_ang_h, fit_nlr_nou_summary_weibull_gym_h,
#   fit_lr_nou_summary_pl_ang_cr, fit_nlr_nou_summary_gmm_gym_cr,
#   fit_lr_nou_summary_pl_ang_dbh, fit_lr_nou_summary_pl_gym_dbh,
#   tallo_reduced_nlr_df_ang_h,
#   tallo_reduced_nlr_df_gym_h,
#   tallo_reduced_lr_df_ang_cr,
#   tallo_reduced_nlr_df_gym_cr,
#   tallo_reduced_lr_df_ang_dbh,
#   tallo_reduced_lr_df_gym_dbh,
#   stan_data_nlr_ang_h,
#   stan_data_nlr_gym_h,
#   stan_data_lr_ang_cr,
#   stan_data_nlr_gym_cr,
#   stan_data_lr_ang_dbh,
#   stan_data_lr_gym_dbh
# ) {

#   # Helper function to filter valid species
#   filter_valid_species <- function(data, variable) {
#     if (variable == "dbh") {
#       data <- data |> filter(!is.na(dbh), !is.na(cr), !is.na(h))
#     } else if (variable == "cr") {
#       data <- data |> filter(!is.na(cr))
#     } else if (variable == "h") {
#       data <- data |> filter(!is.na(h))
#     }
#     species_counts <- data |> 
#       group_by(sp) |> 
#       summarise(count = n(), .groups = "drop") |> 
#       filter(count >= 20)
#     return(species_counts$sp)
#   }

#   # Helper function to extract posterior summaries and prepare species-specific output
#     create_species_model_df <- function(beta_data, dependent_variable, functional_form, sp_list) {
#       if (functional_form == "Weibull" || functional_form == "gMM") {
#         a_median <- exp(beta_data |> filter(str_detect(variable, "beta\\[\\d+,1\\]")) |> pull(q50))
#         a_upper <- exp(beta_data |> filter(str_detect(variable, "beta\\[\\d+,1\\]")) |> pull(q97.5))
#         a_lower <- exp(beta_data |> filter(str_detect(variable, "beta\\[\\d+,1\\]")) |> pull(q2.5))

#         b_median <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,2\\]")) |> pull(q50)
#         b_upper <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,2\\]")) |> pull(q97.5)
#         b_lower <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,2\\]")) |> pull(q2.5)

#         k_median <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,3\\]")) |> pull(q50)
#         k_upper <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,3\\]")) |> pull(q97.5)
#         k_lower <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,3\\]")) |> pull(q2.5)

#         ci_95 <- data.frame(
#           sp = sp_list,
#           a = paste0(a_median, " (", a_lower, ", ", a_upper, ")"),
#           b = paste0(b_median, " (", b_lower, ", ", b_upper, ")"),
#           k = paste0(k_median, " (", k_lower, ", ", k_upper, ")"),
#           Dependent_variable = dependent_variable,
#           stringsAsFactors = FALSE
#         )
#       } else if (functional_form == "Power-Law") {
#         a_median <- exp(beta_data |> filter(str_detect(variable, "beta\\[1,\\d+\\]")) |> pull(q50))
#         a_upper <- exp(beta_data |> filter(str_detect(variable, "beta\\[1,\\d+\\]")) |> pull(q97.5))
#         a_lower <- exp(beta_data |> filter(str_detect(variable, "beta\\[1,\\d+\\]")) |> pull(q2.5))

#         b_median <- beta_data |> filter(str_detect(variable, "beta\\[2,\\d+\\]")) |> pull(q50)
#         b_upper <- beta_data |> filter(str_detect(variable, "beta\\[2,\\d+\\]")) |> pull(q97.5)
#         b_lower <- beta_data |> filter(str_detect(variable, "beta\\[2,\\d+\\]")) |> pull(q2.5)

#         c_median <- beta_data |> filter(str_detect(variable, "beta\\[3,\\d+\\]")) |> pull(q50)
#         c_upper <- beta_data |> filter(str_detect(variable, "beta\\[3,\\d+\\]")) |> pull(q97.5)
#         c_lower <- beta_data |> filter(str_detect(variable, "beta\\[3,\\d+\\]")) |> pull(q2.5)

#         ci_95 <- data.frame(
#           sp = sp_list,
#           a = paste0(a_median, " (", a_lower, ", ", a_upper, ")"),
#           b = paste0(b_median, " (", b_lower, ", ", b_upper, ")"),
#           c = paste0(c_median, " (", c_lower, ", ", c_upper, ")"),
#           Dependent_variable = dependent_variable,
#           stringsAsFactors = FALSE
#         )
#       }

#       return(ci_95)
#     }

#   # Process datasets and align species with Stan data
#   process_dataset <- function(dataset, stan_data, variable, division_filter) {
#     valid_species <- filter_valid_species(dataset, variable)
#     dataset <- dataset |> 
#       filter(sp %in% valid_species, division == division_filter) |> 
#       group_by(sp) |> 
#       summarise(wd = mean(wd, na.rm = TRUE), .groups = "drop") |> 
#       arrange(sp) |> 
#       mutate(sp_id = as.integer(as.factor(sp))) |> 
#       mutate(wd_s = scale(wd) |> as.numeric())
    
#     if (!all(c("jj", "u") %in% names(stan_data))) {
#       stop("Stan data must contain 'jj' and 'u'.")
#     }
    
#     if (nrow(stan_data$u) >= max(stan_data$jj, na.rm = TRUE)) {
#       stan_processed <- data.frame(
#         sp_id = as.factor(stan_data$jj),
#         wd_s = stan_data$u[stan_data$jj, 2]
#       )
#     } else if (nrow(stan_data$u) == 2) {
#       stan_processed <- data.frame(
#         sp_id = 1:ncol(stan_data$u),
#         wd_s = stan_data$u[2, ]
#       )
#     } else {
#       stop("Unexpected structure in Stan data.")
#     }
    
#     return(list(sp_data = dataset, stan_processed = stan_processed))
#   }

#   # Process datasets for each division and variable
#   h_ang <- process_dataset(tallo_reduced_nlr_df_ang_h, stan_data_nlr_ang_h, "h", "Angiosperm")
#   h_gym <- process_dataset(tallo_reduced_nlr_df_gym_h, stan_data_nlr_gym_h, "h", "Gymnosperm")
#   cr_ang <- process_dataset(tallo_reduced_lr_df_ang_cr, stan_data_lr_ang_cr, "cr", "Angiosperm")
#   cr_gym <- process_dataset(tallo_reduced_nlr_df_gym_cr, stan_data_nlr_gym_cr, "cr", "Gymnosperm")
#   dbh_ang <- process_dataset(tallo_reduced_lr_df_ang_dbh, stan_data_lr_ang_dbh, "dbh", "Angiosperm")
#   dbh_gym <- process_dataset(tallo_reduced_lr_df_gym_dbh, stan_data_lr_gym_dbh, "dbh", "Gymnosperm")

#   # Create posterior summaries
#   h_df_ang <- create_species_model_df(fit_nlr_nou_summary_weibull_ang_h, "Tree Height", "Weibull", h_ang$sp_data$sp)
#   h_df_gym <- create_species_model_df(fit_nlr_nou_summary_weibull_gym_h, "Tree Height", "Weibull", h_gym$sp_data$sp)
#   cr_df_ang <- create_species_model_df(fit_lr_nou_summary_pl_ang_cr, "Crown Radius", "Power-Law", cr_ang$sp_data$sp)
#   cr_df_gym <- create_species_model_df(fit_nlr_nou_summary_gmm_gym_cr, "Crown Radius", "gMM", cr_gym$sp_data$sp)
#   dbh_df_ang <- create_species_model_df(fit_lr_nou_summary_pl_ang_dbh, "DBH", "Power-Law", dbh_ang$sp_data$sp)
#   dbh_df_gym <- create_species_model_df(fit_lr_nou_summary_pl_gym_dbh, "DBH", "Power-Law", dbh_gym$sp_data$sp)

#   # Combine results
#   sp_posterior_h_df <- bind_rows(
#     h_df_ang |> mutate(Division = "Angiosperm"),
#     h_df_gym |> mutate(Division = "Gymnosperm")
#   )
#   sp_posterior_cr_df <- bind_rows(
#     cr_df_ang |> mutate(Division = "Angiosperm"),
#     cr_df_gym |> mutate(Division = "Gymnosperm")
#   )
#   sp_posterior_dbh_df <- bind_rows(
#     dbh_df_ang |> mutate(Division = "Angiosperm"),
#     dbh_df_gym |> mutate(Division = "Gymnosperm")
#   )
#   # Rename columns to remove `sp.` prefix
#   colnames(sp_posterior_h_df) <- gsub("^sp\\.", "", colnames(sp_posterior_h_df))
#   sp_posterior_h_df <- sp_posterior_h_df |> 
#     dplyr::rename(Species = sp) |> 
#     dplyr::select(
#       Dependent_variable, Species, Division, a, b, k)

#   colnames(sp_posterior_cr_df) <- gsub("^sp\\.", "", colnames(sp_posterior_cr_df))
#   sp_posterior_cr_df <- sp_posterior_cr_df |>
#     dplyr::rename(Species = sp) |> 
#     dplyr::select(
#       Dependent_variable, Species, Division, a, b, k)
      
#   colnames(sp_posterior_dbh_df) <- gsub("^sp\\.", "", colnames(sp_posterior_dbh_df))
#   sp_posterior_dbh_df <- sp_posterior_dbh_df |>
#     dplyr::rename(Species = sp) |> 
#     dplyr::select(
#       Dependent_variable, Species, Division, a, b, c)

#   return(list(
#     sp_posterior_h_df = sp_posterior_h_df,
#     sp_posterior_cr_df = sp_posterior_cr_df,
#     sp_posterior_dbh_df = sp_posterior_dbh_df
#   ))
# }


# round_sp_posterior_df <- function(sp_posterior_df) {
#   sp_posterior_df |>
#     dplyr::mutate(
#       a = purrr::map_chr(a, function(x) {
#         gsubfn::gsubfn("(\\d+\\.\\d+)", ~ format(round(as.numeric(.x), 3), nsmall = 3), x)
#       }),
#       b = purrr::map_chr(b, function(x) {
#         gsubfn::gsubfn("(\\d+\\.\\d+)", ~ format(round(as.numeric(.x), 3), nsmall = 3), x)
#       }),
#       k = purrr::map_chr(k, function(x) {
#         if (is.na(x)) return(NA) # Handle NA values
#         gsubfn::gsubfn("(\\d+\\.\\d+)", ~ format(round(as.numeric(.x), 3), nsmall = 3), x)
#       }),
#       c = purrr::map_chr(c, function(x) {
#         if (is.na(x)) return(NA) # Handle NA values
#         gsubfn::gsubfn("(\\d+\\.\\d+)", ~ format(round(as.numeric(.x), 3), nsmall = 3), x)
#       })
#     )
# }

#=============================

#=============================

generate_sim_data <- function() {
  set.seed(123)
  N <- 100
  beta <- c(-2, 1)
  xs <- rnorm(N)
  ys <- rnorm(N, beta[1] + beta[2] * xs, 0.3)
  mu1 <- 2
  mu2 <- -3
  sig1 <- 0.8
  sig2 <- 0.8

  x <- mu1 + xs * sig1
  y <- mu2 + ys * sig2

  list(
    x = x,
    xs = xs,
    y = y,
    ys = ys,
    N = N
  )
}

generate_sim_stan_data <- function(sim_data, model = c("xy", "x", "none")) {
  if (model == "xy") {
    list(x = sim_data$xs, y = sim_data$ys, N = sim_data$N)
  } else if (model == "x") {
    list(x = sim_data$xs, y = sim_data$y, N = sim_data$N)
  } else {
    list(x = sim_data$x, y = sim_data$y, N = sim_data$N)
  }
}


reduce_trees_single <- function(tallo_wd_df, n_single_biome = 100) {
  set.seed(123)

  species_counts <- tallo_wd_df |>
    group_by(sp) |>
    summarise(
      count = n()) |>
    filter(count >= 20)

  tallo_wd_df |>
    filter(sp %in% species_counts$sp) |>
    group_by(sp) |>
    group_modify(~ slice_sample(.x, n = min(nrow(.x), n_single_biome))) |>
    ungroup()

}


sample_100sp <- function(data) {
  # data <- tar_read(tallo_wd_df_100)
  sp_name <- data |>
    pull(sp) |>
    unique() |>
    sample(100)

  data |>
    filter(sp %in% sp_name)

}

reduce_trees <- function(data, n_single_biome = 100, n_multi_biome = 50) {
  filtered_data <- data |>
    group_by(sp, biome) |>
    summarize(n = n()) |>
    filter(n >= 20) |>
    ungroup()

  unique_biome_sp <- filtered_data |>
    dplyr::select(sp, biome) |>
    distinct() |>
    mutate(biome_sp = paste(biome, sp, sep = "_"))

  sp_name <- unique_biome_sp |>
    pull(sp)
  biome_sp_name <- unique_biome_sp |>
    pull(biome_sp)

  multi_sp <- sp_name[duplicated(sp_name)] |> unique()
  single_sp <- sp_name[!sp_name %in% sp_name[duplicated(sp_name)]]

  data2 <- data |>
    filter(biome_sp %in% biome_sp_name)

  multi_df <- data2 |>
    filter(sp %in% multi_sp) |>
    group_by(biome_sp) |>
    group_modify(~ slice_sample(.x, n = min(nrow(.x), n_multi_biome))) |>
    ungroup()

  single_df <- data2 |>
    filter(sp %in% single_sp) |>
    group_by(sp) |>
    group_modify(~ slice_sample(.x, n = min(nrow(.x), n_single_biome))) |>
    ungroup()

  bind_rows(multi_df, single_df)

}

reduce_trees_simple <- function(data, n_max = 100) {
  filtered_data <- data |>
    group_by(sp) |>
    summarize(n = n()) |>
    filter(n >= 20) |>
    ungroup()

  sp_name <- filtered_data |>
    pull(sp) |>
    unique()

  data |>
    filter(sp %in% sp_name) |>
    group_by(sp) |>
    group_modify(~ slice_sample(.x, n = min(nrow(.x), n_max))) |>
    ungroup()

}

#==========================
# DATA PARAMS
#==========================
write_data_yaml <- function(
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
  tallo_reduced_lr_df_gym_dbh3,

  file_path = "data.yaml"
) {
  # Dataset
  d0 <- tallo_wd_df0
  ori_sample_size <- comma(nrow(d0))
  ori_species_number <- comma(d0 |> distinct(sp) |> nrow())
  
  # Process H allometry data
  ang_data_h <- tallo_reduced_lr_df_ang_h |> filter(division == "Angiosperm")
  gym_data_h <- tallo_reduced_lr_df_gym_h |> filter(division == "Gymnosperm")
  
  h_params <- list(
    angiosperm = list(
      sample_size = comma(nrow(ang_data_h)),
      species_number = comma(ang_data_h |> distinct(sp) |> nrow()),
      trait_ranges = list(
        dbh = list(min = min(ang_data_h$dbh), max = max(ang_data_h$dbh)),
        h = list(min = min(ang_data_h$h), max = max(ang_data_h$h)),
        wd = list(min = round(min(ang_data_h$wd), 3), max = round(max(ang_data_h$wd), 3))
      )
    ),
    gymnosperm = list(
      sample_size = comma(nrow(gym_data_h)),
      species_number = comma(gym_data_h |> distinct(sp) |> nrow()),
      trait_ranges = list(
        dbh = list(min = min(gym_data_h$dbh), max = max(gym_data_h$dbh)),
        h = list(min = min(gym_data_h$h), max = max(gym_data_h$h)),
        wd = list(min = round(min(gym_data_h$wd), 3), max = round(max(gym_data_h$wd), 3))
      )
    )
  )
  
  # Process CR allometry data
  ang_data_cr <- tallo_reduced_lr_df_ang_cr |> filter(division == "Angiosperm")
  gym_data_cr <- tallo_reduced_lr_df_gym_cr |> filter(division == "Gymnosperm")
  
  cr_params <- list(
    angiosperm = list(
      sample_size = comma(nrow(ang_data_cr)),
      species_number = comma(ang_data_cr |> distinct(sp) |> nrow()),
      trait_ranges = list(
        dbh = list(min = min(ang_data_cr$dbh), max = max(ang_data_cr$dbh)),
        cr = list(min = min(ang_data_cr$cr), max = max(ang_data_cr$cr)),
        wd = list(min = round(min(ang_data_cr$wd), 3), max = round(max(ang_data_cr$wd), 3))
      )
    ),
    gymnosperm = list(
      sample_size = comma(nrow(gym_data_cr)),
      species_number = comma(gym_data_cr |> distinct(sp) |> nrow()),
      trait_ranges = list(
        dbh = list(min = min(gym_data_cr$dbh), max = max(gym_data_cr$dbh)),
        cr = list(min = min(gym_data_cr$cr), max = max(gym_data_cr$cr)),
        wd = list(min = round(min(gym_data_cr$wd), 3), max = round(max(gym_data_cr$wd), 3))
      )
    )
  )
  
  # Process DBH allometry data
  ang_data_dbh <- tallo_reduced_lr_df_ang_dbh |> 
    filter(division == "Angiosperm", 
           !is.na(dbh), 
           !is.na(cr), 
           !is.na(h)) |> 
    group_by(sp) |> 
    filter(n() >= 20) |> 
    ungroup()
  
  gym_data_dbh <- tallo_reduced_lr_df_gym_dbh |> 
    filter(division == "Gymnosperm", 
           !is.na(dbh), 
           !is.na(cr), 
           !is.na(h)) |> 
    group_by(sp) |> 
    filter(n() >= 20) |> 
    ungroup()
  
  dbh_params <- list(
    angiosperm = list(
      sample_size = comma(nrow(ang_data_dbh)),
      species_number = comma(ang_data_dbh |> distinct(sp) |> nrow()),
      trait_ranges = list(
        dbh = list(min = min(ang_data_dbh$dbh), max = max(ang_data_dbh$dbh)),
        cr = list(min = min(ang_data_dbh$cr), max = max(ang_data_dbh$cr)),
        h = list(min = min(ang_data_dbh$h), max = max(ang_data_dbh$h)),
        wd = list(min = round(min(ang_data_dbh$wd), 3), max = round(max(ang_data_dbh$wd), 3))
      )
    ),
    gymnosperm = list(
      sample_size = comma(nrow(gym_data_dbh)),
      species_number = comma(gym_data_dbh |> distinct(sp) |> nrow()),
      trait_ranges = list(
        dbh = list(min = min(gym_data_dbh$dbh), max = max(gym_data_dbh$dbh)),
        cr = list(min = min(gym_data_dbh$cr), max = max(gym_data_dbh$cr)),
        h = list(min = min(gym_data_dbh$h), max = max(gym_data_dbh$h)),
        wd = list(min = round(min(gym_data_dbh$wd), 3), max = round(max(gym_data_dbh$wd), 3))
      )
    )
  )
  # Process DBH1 allometry data
  ang_data_dbh1 <- tallo_reduced_lr_df_ang_dbh1 |> 
    filter(division == "Angiosperm", 
           !is.na(dbh), 
           !is.na(cr), 
           !is.na(h)) |> 
    group_by(sp) |> 
    filter(n() >= 20) |> 
    ungroup()
  
  gym_data_dbh1 <- tallo_reduced_lr_df_gym_dbh1 |> 
    filter(division == "Gymnosperm", 
           !is.na(dbh), 
           !is.na(cr), 
           !is.na(h)) |> 
    group_by(sp) |> 
    filter(n() >= 20) |> 
    ungroup()
    
  dbh1_params <- list(
      angiosperm = list(
        sample_size = comma(nrow(ang_data_dbh1)),
        species_number = comma(ang_data_dbh1 |> distinct(sp) |> nrow()),
        trait_ranges = list(
          dbh = list(min = min(ang_data_dbh1$dbh), max = max(ang_data_dbh1$dbh)),
          cr = list(min = min(ang_data_dbh1$cr), max = max(ang_data_dbh1$cr)),
          h = list(min = min(ang_data_dbh1$h), max = max(ang_data_dbh1$h)),
          wd = list(min = round(min(ang_data_dbh1$wd), 3), max = round(max(ang_data_dbh1$wd), 3))
        )
      ),
      gymnosperm = list(
        sample_size = comma(nrow(gym_data_dbh1)),
        species_number = comma(gym_data_dbh1 |> distinct(sp) |> nrow()),
        trait_ranges = list(
          dbh = list(min = min(gym_data_dbh1$dbh), max = max(gym_data_dbh1$dbh)),
          cr = list(min = min(gym_data_dbh1$cr), max = max(gym_data_dbh1$cr)),
          h = list(min = min(gym_data_dbh1$h), max = max(gym_data_dbh1$h)),
          wd = list(min = round(min(gym_data_dbh1$wd), 3), max = round(max(gym_data_dbh1$wd), 3))
        )
      )
    )

 # Process DBH2
  ang_data_dbh2 <- tallo_reduced_lr_df_ang_dbh2 |> 
    filter(division == "Angiosperm", 
           !is.na(dbh), 
           !is.na(cr), 
           !is.na(h)) |> 
    group_by(sp) |> 
    filter(n() >= 20) |> 
    ungroup()
  
  gym_data_dbh2 <- tallo_reduced_lr_df_gym_dbh2 |> 
    filter(division == "Gymnosperm", 
           !is.na(dbh), 
           !is.na(cr), 
           !is.na(h)) |> 
    group_by(sp) |> 
    filter(n() >= 20) |> 
    ungroup()

  dbh2_params <- list(
    angiosperm = list(
      sample_size = comma(nrow(ang_data_dbh2)),
      species_number = comma(ang_data_dbh2 |> distinct(sp) |> nrow()),
      trait_ranges = list(
        dbh = list(min = min(ang_data_dbh2$dbh), max = max(ang_data_dbh2$dbh)),
        cr = list(min = min(ang_data_dbh2$cr), max = max(ang_data_dbh2$cr)),
        wd = list(min = round(min(ang_data_dbh2$wd), 3), max = round(max(ang_data_dbh2$wd), 3))
      )
    ),
    gymnosperm = list(
      sample_size = comma(nrow(gym_data_dbh2)),
      species_number = comma(gym_data_dbh2 |> distinct(sp) |> nrow()),
      trait_ranges = list(
        dbh = list(min = min(gym_data_dbh2$dbh), max = max(gym_data_dbh2$dbh)),
        cr = list(min = min(gym_data_dbh2$cr), max = max(gym_data_dbh2$cr)),
        wd = list(min = round(min(gym_data_dbh2$wd), 3), max = round(max(gym_data_dbh2$wd), 3))
      )
    )
  )
  
  ang_data_dbh3 <- tallo_reduced_lr_df_ang_dbh3 |> 
    filter(division == "Angiosperm", 
           !is.na(dbh), 
           !is.na(cr), 
           !is.na(h)) |> 
    group_by(sp) |> 
    filter(n() >= 20) |> 
    ungroup()
  
  gym_data_dbh3 <- tallo_reduced_lr_df_gym_dbh3 |> 
    filter(division == "Gymnosperm", 
           !is.na(dbh), 
           !is.na(cr), 
           !is.na(h)) |> 
    group_by(sp) |> 
    filter(n() >= 20) |> 
    ungroup()

  dbh3_params <- list(
    angiosperm = list(
      sample_size = comma(nrow(ang_data_dbh3)),
      species_number = comma(ang_data_dbh3 |> distinct(sp) |> nrow()),
      trait_ranges = list(
        dbh = list(min = min(ang_data_dbh3$dbh), max = max(ang_data_dbh3$dbh)),
        h = list(min = min(ang_data_dbh3$h), max = max(ang_data_dbh3$h)),
        wd = list(min = round(min(ang_data_dbh3$wd), 3), max = round(max(ang_data_dbh3$wd), 3))
      )
    ),
    gymnosperm = list(
      sample_size = comma(nrow(gym_data_dbh3)),
      species_number = comma(gym_data_dbh3 |> distinct(sp) |> nrow()),
      trait_ranges = list(
        dbh = list(min = min(gym_data_dbh3$dbh), max = max(gym_data_dbh3$dbh)),
        h = list(min = min(gym_data_dbh3$h), max = max(gym_data_dbh3$h)),
        wd = list(min = round(min(gym_data_dbh3$wd), 3), max = round(max(gym_data_dbh3$wd), 3))
      )
    )
  )
  
  # Combine all parameters into a single list
  params <- list(
    ori_sample_size = ori_sample_size,
    ori_species_number = ori_species_number,
    h_allometry = h_params,
    cr_allometry = cr_params,
    dbh_allometry = dbh_params,
    dbh1_allometry = dbh1_params,
    dbh2_allometry = dbh2_params,
    dbh3_allometry = dbh3_params
  )
  
  # Write to YAML
  write_yaml(params, file = file_path)
}

#====================
# POSTERIOR YAML
#====================
# write_posterior_yaml <- function(posterior_df) {
#   params_list <- list()
  
#   round_for_small_values <- function(x) {
#     if (!is.na(x) && abs(as.numeric(x)) < 0.001) {  # Check if the value is less than 0.001
#       round(as.numeric(x), 4)  # Round to 4 decimal places for small values
#     } else {
#       round(as.numeric(x), 3)  # Round to 3 decimal places for other values
#     }
#   }
  
#   extract_ci_values <- function(ci_string) {
#     if (is.na(ci_string) || ci_string == "-") return(list(Median = '-', Lower = '-', Upper = '-'))
#     median_value <- str_extract(ci_string, "^[0-9.eE-]+")
#     lower_bound <- str_extract(ci_string, "(?<=\\()[0-9.eE-]+")
#     upper_bound <- str_extract(ci_string, "[0-9.eE-]+(?=\\))")
#     list(
#       Median = ifelse(is.na(median_value), '-', round_for_small_values(median_value)),
#       Lower = ifelse(is.na(lower_bound), '-', round_for_small_values(lower_bound)),
#       Upper = ifelse(is.na(upper_bound), '-', round_for_small_values(upper_bound))
#     )
#   }
  
#   for (i in seq_len(nrow(posterior_df))) {
#     row <- posterior_df[i, ]
#     dep_var <- row$Dependent_variable
#     division <- row$Division  # Adding division to the hierarchy
#     func_form <- row$Functional_form
#     wood_density <- row$Wood_density
#     param <- row$Parameter
    
#     if (!is.list(params_list[[dep_var]])) params_list[[dep_var]] <- list()
#     if (!is.list(params_list[[dep_var]][[division]])) params_list[[dep_var]][[division]] <- list()  # Handle division
#     if (!is.list(params_list[[dep_var]][[division]][[func_form]])) params_list[[dep_var]][[division]][[func_form]] <- list()
#     if (!is.list(params_list[[dep_var]][[division]][[func_form]][[wood_density]])) {
#       params_list[[dep_var]][[division]][[func_form]][[wood_density]] <- list()
#     }
    
#     intercept_ci <- extract_ci_values(row$Intercept_CI)
#     slope_ci <- extract_ci_values(row$Slope_CI)
#     tau_ci <- extract_ci_values(row$Tau_CI)
    
#     params_list[[dep_var]][[division]][[func_form]][[wood_density]][[param]] <- list(
#       Intercept_CI = intercept_ci,
#       Slope_CI = slope_ci,
#       Tau_CI = tau_ci
#     )
#   }

#   write_yaml(params_list, file = "posterior.yaml")
# }

write_posterior_yaml <- function(posterior_df) {
  params_list <- list()
  
  extract_ci_values <- function(ci_string) {
    if (is.na(ci_string) || ci_string == "-") return(list(Median = '-', Lower = '-', Upper = '-'))
    median_value <- str_extract(ci_string, "^[0-9.eE-]+")
    lower_bound <- str_extract(ci_string, "(?<=\\()[0-9.eE-]+")
    upper_bound <- str_extract(ci_string, "[0-9.eE-]+(?=\\))")
    list(
      Median = ifelse(is.na(median_value), '-', as.numeric(median_value)),
      Lower = ifelse(is.na(lower_bound), '-', as.numeric(lower_bound)),
      Upper = ifelse(is.na(upper_bound), '-', as.numeric(upper_bound))
    )
  }
  
  for (i in seq_len(nrow(posterior_df))) {
    row <- posterior_df[i, ]
    dep_var <- row$Dependent_variable
    division <- row$Division  # Adding division to the hierarchy
    func_form <- row$Functional_form
    wood_density <- row$Wood_density
    param <- row$Parameter
    
    if (!is.list(params_list[[dep_var]])) params_list[[dep_var]] <- list()
    if (!is.list(params_list[[dep_var]][[division]])) params_list[[dep_var]][[division]] <- list()  # Handle division
    if (!is.list(params_list[[dep_var]][[division]][[func_form]])) params_list[[dep_var]][[division]][[func_form]] <- list()
    if (!is.list(params_list[[dep_var]][[division]][[func_form]][[wood_density]])) {
      params_list[[dep_var]][[division]][[func_form]][[wood_density]] <- list()
    }
    
    intercept_ci <- extract_ci_values(row$Intercept_CI)
    slope_ci <- extract_ci_values(row$Slope_CI)
    tau_ci <- extract_ci_values(row$Tau_CI)
    
    params_list[[dep_var]][[division]][[func_form]][[wood_density]][[param]] <- list(
      Intercept_CI = intercept_ci,
      Slope_CI = slope_ci,
      Tau_CI = tau_ci
    )
  }

  write_yaml(params_list, file = "posterior.yaml")
}

#================
#SUPPLEMENT
#================
subset_charac <- function(tallo_reduced_lr_df_ang_h,
                          tallo_reduced_lr_df_ang_cr,
                          tallo_reduced_lr_df_ang_dbh,
                          tallo_reduced_lr_df_gym_dbh,
                          tallo_reduced_lr_df_ang_dbh1,
                          tallo_reduced_lr_df_gym_dbh1,
                          tallo_reduced_lr_df_ang_dbh2,
                          tallo_reduced_lr_df_gym_dbh2,
                          tallo_reduced_lr_df_ang_dbh3,
                          tallo_reduced_lr_df_gym_dbh3) {
  
  # Helper function to format numbers with commas
  format_number <- function(x) {
    format(x, big.mark = ",", scientific = FALSE)
  }

  ang_data_h <- tallo_reduced_lr_df_ang_h |> filter(division == "Angiosperm")
  gym_data_h <- tallo_reduced_lr_df_ang_h |> filter(division == "Gymnosperm")

  # Process CR allometry dataset
  ang_data_cr <- tallo_reduced_lr_df_ang_cr |> filter(division == "Angiosperm")
  gym_data_cr <- tallo_reduced_lr_df_ang_cr |> filter(division == "Gymnosperm")

  # Process DBH allometry datasets with specific filters
  ang_data_dbh <- tallo_reduced_lr_df_ang_dbh |>
    filter(division == "Angiosperm", !is.na(dbh), !is.na(cr), !is.na(h)) |>
    group_by(sp) |> filter(n() >= 20) |>
    ungroup()

  gym_data_dbh <- tallo_reduced_lr_df_gym_dbh |>
    filter(division == "Gymnosperm", !is.na(dbh), !is.na(cr), !is.na(h)) |>
    group_by(sp) |> filter(n() >= 20) |>
    ungroup()

  ang_data_dbh1 <- tallo_reduced_lr_df_ang_dbh1 |>
    filter(division == "Angiosperm", !is.na(dbh), !is.na(cr), !is.na(h)) |>
    group_by(sp) |> filter(n() >= 20) |>
    ungroup()

  gym_data_dbh1 <- tallo_reduced_lr_df_gym_dbh1 |>
    filter(division == "Gymnosperm", !is.na(dbh), !is.na(cr), !is.na(h)) |>
    group_by(sp) |> filter(n() >= 20) |>
    ungroup()

  ang_data_dbh2 <- tallo_reduced_lr_df_ang_dbh2 |>
    filter(division == "Angiosperm", !is.na(dbh), !is.na(cr), !is.na(h)) |>
    group_by(sp) |> filter(n() >= 20) |>
    ungroup()

  gym_data_dbh2 <- tallo_reduced_lr_df_gym_dbh2 |>
    filter(division == "Gymnosperm", !is.na(dbh), !is.na(cr), !is.na(h)) |>
    group_by(sp) |> filter(n() >= 20) |>
    ungroup()

  ang_data_dbh3 <- tallo_reduced_lr_df_ang_dbh3 |>
    filter(division == "Angiosperm", !is.na(dbh), !is.na(cr), !is.na(h)) |>
    group_by(sp) |> filter(n() >= 20) |>
    ungroup()

  gym_data_dbh3 <- tallo_reduced_lr_df_gym_dbh3 |>
    filter(division == "Gymnosperm", !is.na(dbh), !is.na(cr), !is.na(h)) |>
    group_by(sp) |> filter(n() >= 20) |>
    ungroup()

  # Recreate `combined_data` with corrected Predictor variable range logic
  combined_data <- tibble(
    `Dependent variables` = rep(c("H", "CR", "DBH", "DBH", "DBH", "DBH"), each = 2),
    `Predictor variable` = c(
      "DBH", "DBH", 
      "DBH", "DBH",
      "CR, H", "CR, H", 
      "CR × H", "CR × H", 
      "CR", "CR", 
      "H", "H"
    ),
    Division = rep(c("Angiosperm", "Gymnosperm"), times = 6),
    # `Number of trees` = as.character(c(
    #   nrow(ang_data_h), nrow(gym_data_h),
    #   nrow(ang_data_cr), nrow(gym_data_cr),
    #   nrow(ang_data_dbh), nrow(gym_data_dbh),
    #   nrow(ang_data_dbh1), nrow(gym_data_dbh1),
    #   nrow(ang_data_dbh2), nrow(gym_data_dbh2),
    #   nrow(ang_data_dbh3), nrow(gym_data_dbh3)
    # )),
    # `Number of species` = as.character(c(
    #   ang_data_h |> distinct(sp) |> nrow(),
    #   gym_data_h |> distinct(sp) |> nrow(),
    #   ang_data_cr |> distinct(sp) |> nrow(),
    #   gym_data_cr |> distinct(sp) |> nrow(),
    #   ang_data_dbh |> distinct(sp) |> nrow(),
    #   gym_data_dbh |> distinct(sp) |> nrow(),
    #   ang_data_dbh1 |> distinct(sp) |> nrow(),
    #   gym_data_dbh1 |> distinct(sp) |> nrow(),
    #   ang_data_dbh2 |> distinct(sp) |> nrow(),
    #   gym_data_dbh2 |> distinct(sp) |> nrow(),
    #   ang_data_dbh3 |> distinct(sp) |> nrow(),
    #   gym_data_dbh3 |> distinct(sp) |> nrow()
    # )),
   `Number of trees` = c(
      format_number(nrow(ang_data_h)), format_number(nrow(gym_data_h)),
      format_number(nrow(ang_data_cr)), format_number(nrow(gym_data_cr)),
      format_number(nrow(ang_data_dbh)), format_number(nrow(gym_data_dbh)),
      format_number(nrow(ang_data_dbh1)), format_number(nrow(gym_data_dbh1)),
      format_number(nrow(ang_data_dbh2)), format_number(nrow(gym_data_dbh2)),
      format_number(nrow(ang_data_dbh3)), format_number(nrow(gym_data_dbh3))
    ),
    `Number of species` = c(
      format_number(ang_data_h |> distinct(sp) |> nrow()),
      format_number(gym_data_h |> distinct(sp) |> nrow()),
      format_number(ang_data_cr |> distinct(sp) |> nrow()),
      format_number(gym_data_cr |> distinct(sp) |> nrow()),
      format_number(ang_data_dbh |> distinct(sp) |> nrow()),
      format_number(gym_data_dbh |> distinct(sp) |> nrow()),
      format_number(ang_data_dbh1 |> distinct(sp) |> nrow()),
      format_number(gym_data_dbh1 |> distinct(sp) |> nrow()),
      format_number(ang_data_dbh2 |> distinct(sp) |> nrow()),
      format_number(gym_data_dbh2 |> distinct(sp) |> nrow()),
      format_number(ang_data_dbh3 |> distinct(sp) |> nrow()),
      format_number(gym_data_dbh3 |> distinct(sp) |> nrow())
    ),
    `Dependent variable range` = c(
      range(ang_data_h$h) |> paste(collapse = " - "), range(gym_data_h$h) |> paste(collapse = " - "),
      range(ang_data_cr$cr) |> paste(collapse = " - "), range(gym_data_cr$cr) |> paste(collapse = " - "),
      range(ang_data_dbh$dbh) |> paste(collapse = " - "), range(gym_data_dbh$dbh) |> paste(collapse = " - "),
      range(ang_data_dbh1$dbh) |> paste(collapse = " - "), range(gym_data_dbh1$dbh) |> paste(collapse = " - "),
      range(ang_data_dbh2$dbh) |> paste(collapse = " - "), range(gym_data_dbh2$dbh) |> paste(collapse = " - "),
      range(ang_data_dbh3$dbh) |> paste(collapse = " - "), range(gym_data_dbh3$dbh) |> paste(collapse = " - ")
    ),
    `Predictor variable range` = c(
      range(ang_data_h$dbh) |> paste(collapse = " - "), range(gym_data_h$dbh) |> paste(collapse = " - "),
      range(ang_data_cr$dbh) |> paste(collapse = " - "), range(gym_data_cr$dbh) |> paste(collapse = " - "),
      paste("CR:", range(ang_data_dbh$cr) |> paste(collapse = " - "), 
            "H:", range(ang_data_dbh$h) |> paste(collapse = " - "), sep = " "),
      paste("CR:", range(gym_data_dbh$cr) |> paste(collapse = " - "), 
            "H:", range(gym_data_dbh$h) |> paste(collapse = " - "), sep = " "),
      paste("CR:", range(ang_data_dbh1$cr) |> paste(collapse = " - "),
            "H:", range(ang_data_dbh1$h) |> paste(collapse = " - "), sep = " "),
      paste("CR:", range(gym_data_dbh1$cr) |> paste(collapse = " - "),
            "H:", range(gym_data_dbh1$h) |> paste(collapse = " - "), sep = " "),
      range(ang_data_dbh2$cr) |> paste(collapse = " - "),
      range(gym_data_dbh2$cr) |> paste(collapse = " - "),
      range(ang_data_dbh3$h) |> paste(collapse = " - "),
      range(gym_data_dbh3$h) |> paste(collapse = " - ")
    ),
    `WD range` = c(
      ang_data_h$wd |> range(na.rm = TRUE) |> round(3) |> paste(collapse = " - "),
      gym_data_h$wd |> range(na.rm = TRUE) |> round(3) |> paste(collapse = " - "),
      ang_data_cr$wd |> range(na.rm = TRUE) |> round(3) |> paste(collapse = " - "),
      gym_data_cr$wd |> range(na.rm = TRUE) |> round(3) |> paste(collapse = " - "),
      ang_data_dbh$wd |> range(na.rm = TRUE) |> round(3) |> paste(collapse = " - "),
      gym_data_dbh$wd |> range(na.rm = TRUE) |> round(3) |> paste(collapse = " - "),
      ang_data_dbh1$wd |> range(na.rm = TRUE) |> round(3) |> paste(collapse = " - "),
      gym_data_dbh1$wd |> range(na.rm = TRUE) |> round(3) |> paste(collapse = " - "),
      ang_data_dbh2$wd |> range(na.rm = TRUE) |> round(3) |> paste(collapse = " - "),
      gym_data_dbh2$wd |> range(na.rm = TRUE) |> round(3) |> paste(collapse = " - "),
      ang_data_dbh3$wd |> range(na.rm = TRUE) |> round(3) |> paste(collapse = " - "),
      gym_data_dbh3$wd |> range(na.rm = TRUE) |> round(3) |> paste(collapse = " - ")
    )
  )

  tbl_sub <- combined_data |>
    pivot_longer(
      cols = c(`Number of trees`, `Number of species`, `Dependent variable range`, `Predictor variable range`, `WD range`),
      names_to = "Characteristics",
      values_to = "Values"
    ) |>
    pivot_wider(
      names_from = Division,
      values_from = Values
    ) |>
    dplyr::select(`Dependent variables`, `Predictor variable`, Characteristics, Angiosperm, Gymnosperm)
    
  return(tbl_sub)
}


#====================================
# SPECIES POSTERIOR DATAFRAME WITH WD
#====================================
generate_sp_posterior_df_wd <- function(
  fit_wd_ang_h_summary_weibull_wd,
  fit_wd_ang_cr_summary_pl_wd,
  fit_wd_gym_cr_summary_gmm_wd,
  tallo_reduced_nlr_df_ang_h,
  tallo_reduced_lr_df_ang_cr,
  tallo_reduced_nlr_df_gym_cr,
  stan_data_nlr_ang_h,
  stan_data_lr_ang_cr,
  stan_data_nlr_gym_cr
) {
  
  # Helper function to extract posterior summaries
  create_species_model_df <- function(beta_data, dependent_variable, functional_form, sp_list) {
    if (functional_form == "Weibull" || functional_form == "gMM") {
      a_median <- exp(beta_data |> filter(str_detect(variable, "beta\\[\\d+,1\\]")) |> pull(q50))
      a_upper <- exp(beta_data |> filter(str_detect(variable, "beta\\[\\d+,1\\]")) |> pull(q97.5))
      a_lower <- exp(beta_data |> filter(str_detect(variable, "beta\\[\\d+,1\\]")) |> pull(q2.5))

      b_median <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,2\\]")) |> pull(q50)
      b_upper <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,2\\]")) |> pull(q97.5)
      b_lower <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,2\\]")) |> pull(q2.5)

      k_median <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,3\\]")) |> pull(q50)
      k_upper <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,3\\]")) |> pull(q97.5)
      k_lower <- beta_data |> filter(str_detect(variable, "beta\\[\\d+,3\\]")) |> pull(q2.5)

      ci_95 <- data.frame(
        sp = sp_list,
        a = paste0(a_median, " (", a_lower, ", ", a_upper, ")"),
        b = paste0(b_median, " (", b_lower, ", ", b_upper, ")"),
        k = paste0(k_median, " (", k_lower, ", ", k_upper, ")"),
        Dependent_variable = dependent_variable,
        stringsAsFactors = FALSE
      )
    } else if (functional_form == "Power-Law") {
      a_median <- exp(beta_data |> filter(str_detect(variable, "beta\\[1,\\d+\\]")) |> pull(q50))
      a_upper <- exp(beta_data |> filter(str_detect(variable, "beta\\[1,\\d+\\]")) |> pull(q97.5))
      a_lower <- exp(beta_data |> filter(str_detect(variable, "beta\\[1,\\d+\\]")) |> pull(q2.5))

      b_median <- beta_data |> filter(str_detect(variable, "beta\\[2,\\d+\\]")) |> pull(q50)
      b_upper <- beta_data |> filter(str_detect(variable, "beta\\[2,\\d+\\]")) |> pull(q97.5)
      b_lower <- beta_data |> filter(str_detect(variable, "beta\\[2,\\d+\\]")) |> pull(q2.5)

      c_median <- beta_data |> filter(str_detect(variable, "beta\\[3,\\d+\\]")) |> pull(q50)
      c_upper <- beta_data |> filter(str_detect(variable, "beta\\[3,\\d+\\]")) |> pull(q97.5)
      c_lower <- beta_data |> filter(str_detect(variable, "beta\\[3,\\d+\\]")) |> pull(q2.5)

      ci_95 <- data.frame(
        sp = sp_list,
        a = paste0(a_median, " (", a_lower, ", ", a_upper, ")"),
        b = paste0(b_median, " (", b_lower, ", ", b_upper, ")"),
        c = paste0(c_median, " (", c_lower, ", ", c_upper, ")"),
        Dependent_variable = dependent_variable,
        stringsAsFactors = FALSE
      )
    }

    return(ci_95)
  }
  # Function to process datasets and align species with Stan data
  process_dataset <- function(dataset, stan_data, division_filter) {
    # Filter the dataset by division
    dataset <- dataset |> filter(division == division_filter)
    
    # Prepare species data
    sp_data <- dataset |>
      group_by(sp) |>
      summarise(wd = mean(wd, na.rm = TRUE)) |>
      arrange(sp) |>
      mutate(sp_id = row_number()) |>
      mutate(wd_s = scale(wd) |> as.numeric()) |>
      dplyr::select(sp_id, sp, wd, wd_s)
    
    # Prepare Stan processed data
    stan_processed <- if (nrow(stan_data$u) >= max(stan_data$jj, na.rm = TRUE)) {
      data.frame(sp_id = as.integer(stan_data$jj), wd_s = stan_data$u[stan_data$jj, 2])
    } else if (nrow(stan_data$u) == 2) {
      data.frame(sp_id = 1:ncol(stan_data$u), wd_s = stan_data$u[2, ])
    } else {
      stop("Unexpected structure for 'u' in Stan data.")
    }
    
    # Convert sp_id to the same type
    sp_data <- sp_data |> mutate(sp_id = as.integer(sp_id))
    
    # Add species names and calculate tree counts
    stan_processed |>
      group_by(sp_id, wd_s) |>
      summarise(tree_count = n(), .groups = "drop") |>
      left_join(sp_data, by = "sp_id")
  }

  
  # Process datasets
  h_ang <- process_dataset(tallo_reduced_nlr_df_ang_h, stan_data_nlr_ang_h, "Angiosperm")
  cr_ang <- process_dataset(tallo_reduced_lr_df_ang_cr, stan_data_lr_ang_cr, "Angiosperm")
  cr_gym <- process_dataset(tallo_reduced_nlr_df_gym_cr, stan_data_nlr_gym_cr, "Gymnosperm")
  
  # Generate posterior dataframes
  h_df_ang <- create_species_model_df(fit_wd_ang_h_summary_weibull_wd, "Tree Height", "Weibull", h_ang$sp)
  cr_df_ang <- create_species_model_df(fit_wd_ang_cr_summary_pl_wd, "Crown Radius", "Power-Law", cr_ang$sp)
  cr_df_gym <- create_species_model_df(fit_wd_gym_cr_summary_gmm_wd, "Crown Radius", "gMM", cr_gym$sp)
  
  # Combine results
  sp_posterior_h_ang_df_wd <- h_df_ang |> mutate(Division = "Angiosperm")
  sp_posterior_cr_df_ang_wd <- cr_df_ang |> mutate(Division = "Angiosperm")
  sp_posterior_cr_df_gym_wd <- cr_df_gym |> mutate(Division = "Gymnosperm")

  return(list(
    sp_posterior_h_ang_df_wd = sp_posterior_h_ang_df_wd,
    sp_posterior_cr_df_ang_wd = sp_posterior_cr_df_ang_wd,
    sp_posterior_cr_df_gym_wd = sp_posterior_cr_df_gym_wd)
  )
}