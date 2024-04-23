# Fit Bayesian logistic regression model

library(tidyverse)
library(brms)
library(ggnewscale)
library(ggbreak)
library(patchwork)

RUN_BRMS <- TRUE # run entire analysis (T) or just plotting functions (F)
run_grid <- read.csv("control.csv", row.names = 1)

## Run all analyses
for (i in 1:4) {
  control <- run_grid[i, ]
  set.seed(1234 + control$inc_sitrep)

  data <- read.csv("data.csv") %>%
    dplyr::mutate(
      vaccinated_est2 = factor(vaccinated_est2,
                               levels = c("No", "Doubtful", "Probable", "Yes"),
                               ordered = TRUE),
      is_vaccinated = vaccinated_est == "yes",
      w_post_70 = year_post_70 / (2024 - 1970)) %>%
    dplyr::filter(!is.na(age_est))


  if (control$inc_sitrep) {
    data <- dplyr::filter(data, str_detect(study, "WHODRC2024", negate = TRUE) |
                            study == control$sitrep)
    if (control$drc_age_sens) {
      data <- data %>%
        dplyr::mutate(age_est = if_else(study == control$sitrep & age_est > 15,
                                        true = 21, false = age_est))
      results_path <- "outputs_inc_drc_sens/"
    } else {
  
      results_path <- "outputs_inc_drc/"
    }
  } else if (control$only_sitrep) {
    data <- dplyr::filter(data, study == control$sitrep)
    results_path <- "outputs_only_drc/"
  } else if (!control$inc_sitrep) {
  
    data <- dplyr::filter(data, str_detect(study, "WHODRC2024", negate = TRUE))
    results_path <- "outputs_exc_drc/"
  
    if (control$exclude_probable) {
      data <- dplyr::filter(data, vaccinated_est2 != "Probable")
      results_path <- "outputs_exc_drc_no_probable/"
    }
  }


  dir.create(results_path, FALSE, TRUE)

  data %>%
    dplyr::select(study, country, age_est, vaccinated_est2, year_post_70,
                  died, recovered, total) %>%
    write_csv(paste0(results_path, "data.csv"))

  if (RUN_BRMS) {
  
    mod <- list()
    fam <- brms::brmsfamily("binomial", "logit")
  
    prior <-  prior(uniform(0, 3), nlpar = "b", lb = 0, ub = 3) +
      prior(uniform(0, 70), nlpar = "a50", lb = 0, ub = 70)
  
    prior_theta <- prior(uniform(0, 1), nlpar = "theta", lb = 0, ub = 1)
  
    prior_ve <-  prior(uniform(-1, 1), nlpar = "ve", lb = -1, ub = 1)
  
    prior_year <-  prior(uniform(0, 1), nlpar = "theta1970", lb = 0, ub = 1) +
      prior(uniform(0, 1), nlpar = "theta2024", lb = 0, ub = 1)
  
    prior_y0 <- prior + prior_theta
    prior_y1 <- prior + prior_year
  
    #' Models
    #' 1: Logistic, age only
    #' 2. Logistic, vacc only
    #' 3. Logistic, year only
    #' 4. Logistic, age + vacc
    #' 5. Logistic, age + year
    #' 6. Logistic, age + vacc + year
    #' 7. Hill, age only
    #' 8. Hill, age + vacc
    #' 9. Hill, age + year
    #' 10. Hill, age + vacc + year
  
    # M1: age only
    mod$m01_age <-
      brms::brm(bf(died | trials(died + recovered) ~
                     log(theta / (1 - theta + exp(age_est - a50) ^ b)),
                   a50 + b + theta ~ 1, nl = TRUE),
                prior = prior_y0, control = list(adapt_delta = 0.98),
                data = data, family = fam, chains = 4, iter = 5e3,
                save_pars = brms::save_pars(all = TRUE)) %>%
      brms::add_criterion("loo", moment_match = TRUE, reloo = TRUE)
    
    
    if (!control$only_sitrep) {  
      
      # M2: vacc only
      mod$m02_vacc <-
        brms::brm(bf(died | trials(died + recovered) ~
                       log(theta * (1 - ve * is_vaccinated) /
                             (1 - theta * (1 - ve * is_vaccinated))),
                     theta + ve ~ 1, nl = TRUE),
                  prior = prior_theta + prior_ve,
                  control = list(adapt_delta = 0.98),
                  data = data, family = fam, chains = 4, iter = 5e3,
                  save_pars = brms::save_pars(all = TRUE)) %>%
        brms::add_criterion("loo", moment_match = TRUE, reloo = TRUE)
      
      
      # M3 age only
      
      mod$m03_year <-
        brms::brm(bf(died | trials(died + recovered) ~
                       log((theta1970 + (theta2024 - theta1970) * w_post_70) /
                             (1 - (theta1970 + (theta2024 - theta1970) * w_post_70))),
                     theta1970 + theta2024 ~ 1, nl = TRUE),
                  prior = prior_year,
                  control = list(adapt_delta = 0.98),
                  data = data, family = fam, chains = 4, iter = 5e3,
                  save_pars = brms::save_pars(all = TRUE)) %>%
        brms::add_criterion("loo", moment_match = TRUE, reloo = TRUE)
      
      ## M4: age + vacc
      mod$m04_age_vacc <-
        brms::brm(bf(died | trials(died + recovered) ~
                       log(theta * (1 - ve * is_vaccinated) /
                             (1 - theta * (1 - ve * is_vaccinated) +
                                exp(age_est - a50) ^ b)),
                     a50 + b + theta + ve ~ 1, nl = TRUE),
                  prior = prior_y0 + prior_ve,
                  control = list(adapt_delta = 0.98),
                  data = data, family = fam, chains = 4, iter = 5e3,
                  save_pars = brms::save_pars(all = TRUE)) %>%
        brms::add_criterion("loo", moment_match = TRUE, reloo = TRUE)
      
      # M5: age + year post 70
      
      mod$m05_age_year <-
        brms::brm(bf(died | trials(died + recovered) ~
                       log((theta1970 + (theta2024 - theta1970) * w_post_70) /
                             (1 - (theta1970 + (theta2024 - theta1970) * w_post_70) +
                                exp(age_est - a50) ^ b)),
                     a50 + b + theta1970 + theta2024 ~ 1, nl = TRUE),
                  prior = prior_y1, control = list(adapt_delta = 0.98),
                  data = data, family = fam, chains = 4, iter = 5e3,
                  save_pars = brms::save_pars(all = TRUE)) %>%
        brms::add_criterion("loo", moment_match = TRUE, reloo = TRUE)
      
      
      
      ## M6: age + vacc + year
      mod$m06_age_vacc_year <-
        brms::brm(bf(died | trials(died + recovered) ~
                       log((theta1970 + (theta2024 - theta1970) * w_post_70) *
                             (1 - ve * is_vaccinated) /
                             (1 - (theta1970 + (theta2024 - theta1970) * w_post_70) *
                                (1 - ve * is_vaccinated) + exp(age_est - a50) ^ b)),
                     a50 + b + ve + theta1970 + theta2024 ~ 1,
                     nl = TRUE),
                  prior = prior_y1 + prior_ve,
                  control = list(adapt_delta = 0.98),
                  data = data, family = fam, chains = 4,
                  iter = 5e3,
                  save_pars = brms::save_pars(all = TRUE)) %>%
        brms::add_criterion("loo", moment_match = TRUE, reloo = TRUE)
    }


    ## M7: log(age) + theta


    mod$m07_log_age_theta <-
      brms::brm(bf(died | trials(died + recovered) ~
                     log(theta / (1 - theta + (age_est / a50) ^ b)),
                   a50 + b + theta ~ 1, nl = TRUE),
                prior = prior_y0,
                control = list(adapt_delta = 0.98),
                data = data, family = fam, chains = 4, iter = 5e3,
                save_pars = brms::save_pars(all = TRUE)) %>%
      brms::add_criterion("loo", moment_match = TRUE, reloo = TRUE)


    ## M8: log(age) + theta*vacc
    if (!control$only_sitrep) {
      
      mod$m08_log_age_theta_vacc <-
        brms::brm(bf(died | trials(died + recovered) ~
                       log(theta * (1 - ve * is_vaccinated) /
                             (1 - theta * (1 - ve * is_vaccinated) +
                                (age_est / a50) ^ b)),
                     a50 + b + theta + ve ~ 1, nl = TRUE),
                  prior = prior_y0 + prior_ve,
                  control = list(adapt_delta = 0.98),
                  data = data, family = fam, chains = 4, iter = 5e3,
                  save_pars = brms::save_pars(all = TRUE)) %>%
        brms::add_criterion("loo", moment_match = TRUE, reloo = TRUE)
      
      ## M9:  log(age) + theta period

      
      mod$m09_log_age_theta_year_post_70 <-
        brms::brm(bf(died | trials(died + recovered) ~
                       log((theta1970 + (theta2024 - theta1970) * w_post_70) /
                             (1 - (theta1970 + (theta2024 - theta1970) * w_post_70) +
                                (age_est / a50) ^ b)),
                     a50 + b + theta1970 + theta2024 ~ 1, nl = TRUE),
                  prior = prior_y1, control = list(adapt_delta = 0.98),
                  data = data, family = fam, chains = 4, iter = 5e3,
                  save_pars = brms::save_pars(all = TRUE)) %>%
        brms::add_criterion("loo", moment_match = TRUE, reloo = TRUE)



      # M10: log(age) + year + vacc
      
      mod$m10_log_age_theta_vacc_year_post_70 <-
        brms::brm(bf(died | trials(died + recovered) ~
                       log((theta1970 + (theta2024 - theta1970) * w_post_70) *
                             (1 - ve * is_vaccinated) /
                             (1 - (theta1970 + (theta2024 - theta1970) * w_post_70) *
                                (1 - ve * is_vaccinated) + (age_est / a50) ^ b)),
                     a50 + b + ve + theta1970 + theta2024 ~ 1, nl = TRUE),
                  prior = prior_y1 + prior_ve,
                  control = list(adapt_delta = 0.98),
                  data = data, family = fam, chains = 4, iter = 5e3,
                  save_pars = brms::save_pars(all = TRUE)) %>%
        brms::add_criterion("loo", moment_match = TRUE, reloo = TRUE)


    }
    ## Model comparison with LOO-CV
    # LOO provides an estimate of expected log posterior predictive distribution (ELPD)
    
    dir.create(results_path, FALSE, TRUE)
    saveRDS(mod, paste0(results_path, "models.rds"))
  }

  mod <- readRDS(paste0(results_path, "models.rds"))
  
  
  if (control$only_sitrep) {
    elpd <- brms::loo_compare(mod$m01_age, mod$m07_log_age_theta) %>%
      as.data.frame()
  } else {
    elpd <- brms::loo_compare(mod$m01_age, 
                              mod$m02_vacc,
                              mod$m03_year,
                              mod$m04_age_vacc, 
                              mod$m05_age_year,
                              mod$m06_age_vacc_year,
                              mod$m07_log_age_theta, 
                              mod$m08_log_age_theta_vacc, 
                              mod$m09_log_age_theta_year_post_70,
                              mod$m10_log_age_theta_vacc_year_post_70) %>%
      as.data.frame()
  }
  
  elpd$model <- gsub("mod\\$", "", rownames(elpd))
  
  
  log_odds <- lapply(mod, function(x) {
    posterior <- brms::posterior_summary(x)
    posterior <- data.frame(parameter = rownames(posterior), posterior)
    rownames(posterior) <- NULL
    posterior }) %>%
    dplyr::bind_rows(.id = "model")
  
  log_odds %>%
    dplyr::mutate(beta = sprintf("%.3f (%.3f,%.3f)", Estimate, Q2.5, Q97.5)) %>%
    dplyr::select(model, parameter, beta) %>%
    dplyr::filter(!startsWith(parameter, "lp")) %>%
    tidyr::pivot_longer(beta) %>%
    dplyr::mutate(parameter = gsub("1970", "", parameter)) %>%
    tidyr::pivot_wider(names_from = parameter) %>%
    dplyr::arrange(name, model) %>%
    dplyr::select(name, model, starts_with("b_theta"), b_ve_Intercept,
                  b_a50_Intercept, b_b_Intercept) %>%
    dplyr::left_join(dplyr::group_by(elpd, model) %>%
                       dplyr::transmute(elpd = sprintf("%.1f (%.1f)",
                                                       elpd_diff, se_diff))) %>%
    writexl::write_xlsx(sprintf("%slog_odds_%s_drc.xlsx", results_path,
                                ifelse(control$inc_sitrep, "inc", "exc")))
  
  
  pred_age <- seq(0, max(data$age_est, na.rm = TRUE), by = 0.2)
  posteriors <- lapply(mod, function(x) {
    brms::as_draws_df(x) %>%
      dplyr::mutate(parameter_set = seq_len(nrow(.))) %>%
      as.data.frame()
  })
  
  saveRDS(posteriors, paste0(results_path, "posteriors.rds"))
  
  ## Plotting results
  vaccine_cols <- RColorBrewer::brewer.pal(5, "PuOr")[-3]
  results <- list()
  
  ## Posteriors
  
  posteriors_tidy <- posteriors %>%
    dplyr::bind_rows(.id = "model") %>%
    tidyr::pivot_longer(starts_with("b_"), names_to = "parameter") %>%
    dplyr::mutate(model = as.numeric(str_extract(model, "[0-9]{2}")),
                  theta_year = str_extract(parameter, "[0-9]{4}"),
                  theta_year = replace_na(theta_year, "All"),
                  parameter = str_replace(parameter, "[0-9]{4}", ""),
                  parameter = str_replace(parameter, "b_", ""),
                  parameter = str_replace(parameter, "_Intercept", ""))
  
  estimate_mode <- function(x) {
    if (any(is.na(x))) return(NA)
    dens <- density(x)
    dens$x[which.max(dens$y)]
  }
  
  posterior_mode <- posteriors_tidy %>%
    dplyr::group_by(model, parameter, theta_year) %>%
    dplyr::summarise(mode = estimate_mode(value))
  
  param_labels <- c(a50 = "a[50]", b = "b", theta = "theta", ve = "phi")
  g <- posteriors_tidy %>%
    dplyr::filter(model > 6) %>%
    dplyr::mutate(
      param_label = factor(parameter, levels = names(param_labels),
                           labels = param_labels, ordered = TRUE)) %>%
    ggplot(aes(x = value, group = theta_year, fill = theta_year)) +
    geom_density(alpha = 0.7) +
    ggh4x::facet_grid2(vars(model), vars(param_label), scales = "free",
                       independent = "y",
                       labeller = labeller(model = function(x) paste0("M", x),
                                           param_label = label_parsed)) +
    theme_bw()  +
    labs(fill = "Year", x = "Posterior parameter value", y = "Density") +
    theme(legend.position = "bottom", legend.direction = "horizontal")
  
  ggsave(paste0(results_path, "posterior_density.png"), g,
         width = 17, height = 10, unit = "cm", scale = 1.3)
  
  posterior_mode %>%
    tidyr::pivot_wider(names_from = c(parameter, theta_year),
                       values_from = mode) %>%
    writexl::write_xlsx(sprintf("%sposterior_mode_%s_drc.xlsx", results_path,
                                ifelse(control$inc_sitrep, "inc", "exc")))
  
  
  
  cfr <- function(x_a, theta, a50, b) theta / (1 + (x_a / a50) ^ b)
  

    # Plot M7 (log age) theta
    results$log_age_theta <- tidyr::expand_grid(x_a = pred_age, posteriors$m07_log_age_theta) %>%
      dplyr::mutate(
        p = cfr(x_a, b_theta_Intercept, b_a50_Intercept, b_b_Intercept)) %>%
      dplyr::group_by(x_a) %>%
      dplyr::summarise(mean = mean(p),
                       median = median(p),
                       q2.5 = quantile(p, 0.025),
                       q25 = quantile(p, 0.25),
                       q75 = quantile(p, 0.75),
                       q97.5 = quantile(p, 0.975))
    
    g_log_age_theta <- results$log_age_theta %>%
      ggplot(aes(x = x_a, y = median)) +
      geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.2) +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3) +
      geom_line() +
      geom_point(data = data,
                 aes(x = age_est, y = died / (died + recovered),
                     size = died + recovered, colour = study),
                 position = position_jitter(width = 0.5, height = 0, seed = 1),
                 inherit.aes = FALSE, alpha = 0.5) +
      scale_size_continuous(
        breaks = c(1, 50, ifelse(control$inc_sitrep, 500, 150))) +
      theme_bw() +
      theme(legend.text = element_text(size = 8),
            legend.title = element_text(size = 10)) +
      labs(y = "Case Fatality Ratio (CFR)", x = "Age (years)",
           colour = "Author / Year", size = "Study size (N)")
    
    
    ggsave(paste0(results_path, "cfr_by_age_M7_log_age_theta.png"),
           g_log_age_theta, width = 17, height = 15, unit = "cm")
    
    
    ## Plot M8: log(age) + theta vacc
    results$log_age_vacc <- tidyr::expand_grid(
      x_a = pred_age,
      x_v = c(TRUE, FALSE),
      posteriors$m08_log_age_theta_vacc) %>%
      dplyr::mutate(
        theta = b_theta_Intercept * (1 - b_ve_Intercept * x_v),
        p = cfr(x_a, theta, b_a50_Intercept, b_b_Intercept)) %>%
      dplyr::group_by(x_a, x_v) %>%
      dplyr::summarise(mean = mean(p),
                       median = median(p),
                       q2.5 = quantile(p, 0.025),
                       q25 = quantile(p, 0.25),
                       q75 = quantile(p, 0.75),
                       q97.5 = quantile(p, 0.975))
    
    g_log_age_vacc <- results$log_age_vacc %>%
      ggplot(aes(x = x_a, y = median, group = x_v, fill = x_v)) +
      geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.2) +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3) +
      geom_line(aes(colour = x_v)) +
      scale_fill_manual(values = vaccine_cols[c(1, 4)],
                        labels = c("FALSE" = "Doubtful / No",
                                   "TRUE" = "Probable / Yes"),
                        aesthetics = c("fill", "colour")) +
      labs(y = "Case Fatality Ratio (CFR)", x = "Age (years)",
           colour = "Vaccinated (Model)", fill = "Vaccinated (Model)") +
      ggnewscale::new_scale_colour() +
      ggnewscale::new_scale_fill() +
      geom_point(data = data,
                 aes(x = age_est, y = died / total, size = total,
                     fill = vaccinated_est2), 
                 position = position_jitter(width = 0.5, height = 0, seed = 1),
                 pch = 21, inherit.aes = FALSE, alpha = 0.5) +
      scale_colour_manual(values = vaccine_cols,
                          aesthetics = c("colour", "fill")) +
      scale_size_continuous(
        breaks = c(1, 50, ifelse(control$inc_sitrep, 500, 150))) +
      guides(colour = guide_legend(order = 2),
             fill = guide_legend(order = 2),
             size = guide_legend(order = 1)) +
      labs(y = "Case Fatality Ratio (CFR)", x = "Age (years)",
           colour = "Vaccinated (Data)", fill = "Vaccinated (Data)",
           size = "Study size (N)") +
      theme_bw()
    
    ggsave(paste0(results_path, "cfr_by_age_M8_log_age_vacc_data.png"),
           g_log_age_vacc, width = 17, height = 15, unit = "cm")
    
    
   
    ## M10 Plot models with year + vacc
    
    results$log_age_vacc_year <- tidyr::expand_grid(
      x_a = pred_age,
      x_v = c(TRUE, FALSE),
      year = seq(1970, 2024, length.out = 2),
      posteriors$m10_log_age_theta_vacc_year_post_70) %>%
      dplyr::mutate(
        x_y = (year - 1970) / (2024 - 1970),
        theta = b_theta1970_Intercept +
          x_y * (b_theta2024_Intercept - b_theta1970_Intercept),
        theta = theta * (1 - b_ve_Intercept * x_v),
        p = cfr(x_a, theta, b_a50_Intercept, b_b_Intercept)) %>%
      dplyr::group_by(x_a, x_v, year) %>%
      dplyr::summarise(mean = mean(p),
                       median = median(p),
                       q2.5 = quantile(p, 0.025),
                       q25 = quantile(p, 0.25),
                       q75 = quantile(p, 0.75),
                       q97.5 = quantile(p, 0.975))
    
    
    g_M10 <- results$log_age_vacc_year %>%
      ggplot(aes(x = x_a, y = median, group = year, fill = year)) +
      geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.3) +
      geom_line(aes(colour = year)) +
      labs(y = "Case Fatality Ratio (CFR)", x = "Age (years)", colour = "Year",
           fill = "Year", size = "Study size (N)") +
      geom_point(
        data = dplyr::rename(data, x_v = is_vaccinated),
        aes(x = age_est, y = died / (died + recovered),
            size = died + recovered, fill = year, colour = year),
        pch = 16, inherit.aes = FALSE, alpha = 0.7,
        position = position_jitter(width = 0.5, height = 0, seed = 1)) +
      scale_fill_distiller(palette = "Spectral",
                           aesthetics = c("fill", "colour")) +
      scale_fill_viridis_c(aesthetics = c("fill", "colour"),
                           end = 0.85, begin = 0.2, direction = -1,
                           option = "mako") +
      scale_size_continuous(breaks = c(1, 100, 1000)) +
      theme_bw() +
      ggbreak::scale_y_break(breaks = c(0.35, 0.9),
                             ticklabels = c(seq(0, 0.4, 0.1), c(0.9, 1)),
                             expand = c(0, 0.01)) +
      facet_grid(
        cols = vars(x_v),
        labeller = labeller(.cols = c("TRUE" = "Vaccinated / Probable",
                                      "FALSE" = "Unvaccinated / Doubtful")))
    
    
    ggsave(paste0(results_path, "cfr_by_age_M10.png"), g_M10,
           width = 15, height = 10, unit = "cm")
    
    
    
    ## Save CFR results
    results %>%
      dplyr::bind_rows(.id = "model") %>%
      dplyr::filter(x_a %in% c(1, seq(5, 40, 5))) %>%
      dplyr::group_by(model, x_a, x_v, year) %>%
      dplyr::transmute(
        value = sprintf("%.1f (%.1f - %.1f)",
                        median * 100, q2.5 * 100, q97.5 * 100)) %>%
      tidyr::pivot_wider(names_from = c(x_v, year), names_prefix = "x_v=") %>%
      writexl::write_xlsx(paste0(results_path, "cfr.xlsx"))
    
    
    m07_obs_vs_exp <- data %>%
      dplyr::group_by(age_est) %>%
      dplyr::summarise(died = sum(died),
                       total = sum(total)) %>%
      tidyr::expand_grid(dplyr::select(posteriors$m07_log_age_theta,
                                       c(.draw, starts_with("b_")))) %>%
      dplyr::mutate(
        p = cfr(age_est, b_theta_Intercept, b_a50_Intercept, b_b_Intercept),
        age_group = cut(age_est, breaks = c(0, 5, 15, 100),
                        include.lowest = TRUE, right = FALSE,
                        labels = c("<5", "5-14", "15+"))) %>%
      dplyr::group_by(.draw, age_group) %>%
      dplyr::summarise(died_observed = sum(died),
                       cfr_expected = sum(p * total) / sum(total),
                       total = sum(total)) %>%
      dplyr::group_by(age_group, died_observed, total) %>%
      dplyr::summarise(cfr_expected_mean = mean(cfr_expected),
                       cfr_expected_q2.5 = quantile(cfr_expected, 0.025),
                       cfr_expected_q97.5 = quantile(cfr_expected, 0.975)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        cfr_observed = died_observed / total,
        cfr_observed_q2.5 = binom.test(died_observed, total)$conf.int[1],
        cfr_observed_q97.5 = binom.test(died_observed, total)$conf.int[2])
    
    
    g_m07_obs_vs_exp <- m07_obs_vs_exp %>%
      ggplot(aes(x = age_group)) +
      geom_boxplot(aes(lower = cfr_expected_q2.5, upper = cfr_expected_q97.5,
                       middle = cfr_expected_mean, ymin = cfr_expected_q2.5,
                       ymax = cfr_expected_q97.5),
                   stat = "identity", fill = "steelblue", alpha = 0.5) +
      geom_point(aes(y = cfr_observed)) +
      geom_segment(aes(y = cfr_observed_q2.5,
                       yend = cfr_observed_q97.5,
                       xend = age_group)) +
      theme_bw() +
      labs(x = "Age group", y = "Case fatality ratio (CFR)")
    
    ggsave(paste0(results_path, "cfr_obs_vs_exp_M07.png"), g_m07_obs_vs_exp,
           width = 15, height = 10, unit = "cm")
    
    
    m08_obs_vs_exp <- data %>%
      dplyr::group_by(age_est, is_vaccinated) %>%
      dplyr::summarise(died = sum(died),
                       total = sum(total)) %>%
      tidyr::expand_grid(dplyr::select(posteriors$m08_log_age_theta_vacc,
                                       c(.draw, starts_with("b_")))) %>%
      dplyr::mutate(
        theta = b_theta_Intercept * (1 - b_ve_Intercept * is_vaccinated),
        p = cfr(age_est, theta, b_a50_Intercept, b_b_Intercept),
        age_group = cut(age_est, breaks = c(0, 5, 15, 100),
                        include.lowest = TRUE, right = FALSE,
                        labels = c("<5", "5-14", "15+"))) %>%
      dplyr::group_by(.draw, age_group, is_vaccinated) %>%
      dplyr::summarise(died_observed = sum(died),
                       cfr_expected = sum(p * total) / sum(total),
                       total = sum(total)) %>%
      dplyr::group_by(age_group, died_observed, total, is_vaccinated) %>%
      dplyr::summarise(cfr_expected_mean = mean(cfr_expected),
                       cfr_expected_q2.5 = quantile(cfr_expected, 0.025),
                       cfr_expected_q97.5 = quantile(cfr_expected, 0.975)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        cfr_observed = died_observed / total,
        cfr_observed_q2.5 = binom.test(died_observed, total)$conf.int[1],
        cfr_observed_q97.5 = binom.test(died_observed, total)$conf.int[2])
    
    
    g_m08_obs_vs_exp <- m08_obs_vs_exp %>%
      ggplot(aes(x = age_group, fill = is_vaccinated)) +
      geom_boxplot(aes(lower = cfr_expected_q2.5, upper = cfr_expected_q97.5,
                       middle = cfr_expected_mean, ymin = cfr_expected_q2.5,
                       ymax = cfr_expected_q97.5),
                   stat = "identity", alpha = 0.5) +
      geom_point(aes(y = cfr_observed)) +
      geom_segment(aes(y = cfr_observed_q2.5,
                       yend = cfr_observed_q97.5,
                       xend = age_group)) +
      theme_bw() +
      labs(x = "Age group", y = "Case fatality ratio (CFR)") +
      facet_wrap(
        vars(is_vaccinated),
        labeller = labeller(.cols = c("TRUE" = "Vaccinated/Probable",
                                      "FALSE" = "Unvaccinated/Doubtful"))) +
      coord_cartesian(ylim = c(0, 0.25)) +
      guides(fill = "none")
    
    ggsave(paste0(results_path, "cfr_obs_vs_exp_M08.png"), g_m08_obs_vs_exp,
           width = 15, height = 6, unit = "cm")

    
    
    results$log_age_vacc_year <- tidyr::expand_grid(
      x_a = pred_age,
      x_v = c(TRUE, FALSE),
      year = c(1970, 2024),
      posteriors$m10_log_age_theta_vacc_year_post_70) %>%
      dplyr::mutate(
        x_y = (year - 1970) / (2024 - 1970),
        theta = b_theta1970_Intercept + x_y *
          (b_theta2024_Intercept - b_theta1970_Intercept),
        theta = theta * (1 - b_ve_Intercept * x_v),
        p = cfr(x_a, theta, b_a50_Intercept, b_b_Intercept)) %>%
      dplyr::group_by(x_a, x_v, year) %>%
      dplyr::summarise(mean = mean(p),
                       median = median(p),
                       q2.5 = quantile(p, 0.025),
                       q25 = quantile(p, 0.25),
                       q75 = quantile(p, 0.75),
                       q97.5 = quantile(p, 0.975))
    
    m10_obs_vs_exp <- data %>%
      dplyr::group_by(age_est, is_vaccinated, w_post_70) %>%
      dplyr::summarise(died = sum(died),
                       total = sum(total)) %>%
      tidyr::expand_grid(
        dplyr::select(posteriors$m10_log_age_theta_vacc_year_post_70,
                      c(.draw, starts_with("b_")))) %>%
      dplyr::mutate(
        theta = b_theta1970_Intercept +
          w_post_70 * (b_theta2024_Intercept - b_theta1970_Intercept),
        theta = theta * (1 - b_ve_Intercept * is_vaccinated),
        p = cfr(age_est, theta, b_a50_Intercept, b_b_Intercept),
        age_group = cut(age_est, breaks = c(0, 5, 15, 100),
                        include.lowest = TRUE, right = FALSE,
                        labels = c("<5", "5-14", "15+"))) %>%
      dplyr::group_by(.draw, age_group, is_vaccinated) %>%
      dplyr::summarise(died_observed = sum(died),
                       cfr_expected = sum(p * total) / sum(total),
                       total = sum(total)) %>%
      dplyr::group_by(age_group, died_observed, total, is_vaccinated) %>%
      dplyr::summarise(cfr_expected_mean = mean(cfr_expected),
                       cfr_expected_q2.5 = quantile(cfr_expected, 0.025),
                       cfr_expected_q97.5 = quantile(cfr_expected, 0.975)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        cfr_observed = died_observed / total,
        cfr_observed_q2.5 = binom.test(died_observed, total)$conf.int[1],
        cfr_observed_q97.5 = binom.test(died_observed, total)$conf.int[2])
    
    
    g_m10_obs_vs_exp <- m10_obs_vs_exp %>%
      ggplot(aes(x = age_group, fill = is_vaccinated)) +
      geom_boxplot(aes(lower = cfr_expected_q2.5, upper = cfr_expected_q97.5,
                       middle = cfr_expected_mean, ymin = cfr_expected_q2.5,
                       ymax = cfr_expected_q97.5),
                   stat = "identity", alpha = 0.5) +
      geom_point(aes(y = cfr_observed)) +
      geom_segment(aes(y = cfr_observed_q2.5,
                       yend = cfr_observed_q97.5,
                       xend = age_group)) +
      theme_bw() +
      labs(x = "Age group", y = "Case fatality ratio (CFR)") +
      facet_grid(
        cols = vars(is_vaccinated),
        labeller = labeller(.cols = c("TRUE" = "Vaccinated/Probable",
                                      "FALSE" = "Unvaccinated/Doubtful"))) +
      coord_cartesian(ylim = c(0, 0.25)) +
      guides(fill = "none")
    
    ggsave(paste0(results_path, "cfr_obs_vs_exp_M10.png"), g_m10_obs_vs_exp,
           width = 15, height = 7, unit = "cm")
  }


