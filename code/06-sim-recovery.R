library(dplyr)
library(ggplot2)
library(drmr)
library(cowplot)

##--- simulation design & true parameters ----

n_patches <- 10
n_years   <- 25
n_ages    <- 5
n_tot     <- n_patches * n_years

true_pars <- list(
  beta_r = c(1.0, .6),   ## log-recruitment: intercept and slope on x_rec
  beta_s = c(.8, -.7),   ## logit-survival: intercept and slope on x_surv
  beta_t = c(-2.5),      ## logit-absence: intercept
  phi    = 10            ## gamma precision parameter
)

sim_toggles <- list(est_surv = 1,
                    minit = 1,
                    rho_mu = 0)

truth_df <-
  tibble(
    variable = c("beta_r[1]", "beta_r[2]",
                 "beta_s[1,1]", "beta_s[1,2]",
                 "beta_t[1]", "phi[1]"),
    param_label = c("beta[r*','*0]~'(Rec. intercept)'",
                    "beta[r*','*1]~'(Rec. slope)'",
                    "beta[s*','*0]~'(Surv. intercept)'",
                    "beta[s*','*1]~'(Surv. slope)'",
                    "beta[t*','*0]~'(Zero-infl. intercept)'",
                    "phi~'(Gamma precision)'"),
    true_val = c(true_pars$beta_r,
                 true_pars$beta_s,
                 true_pars$beta_t,
                 true_pars$phi)
  )

simulate_one <- function(seed, x_r, x_s) {
  set.seed(seed)
  grid <-
    expand.grid(year = seq_len(n_years),
                patch = seq_len(n_patches)) |>
    mutate(x_rec  = x_r,
           x_surv = x_s)
  sim_drm(.data = grid,
          time_col = "year",
          site_col = "patch",
          y_col = "dens",
          family = "gamma",
          formula_zero = ~ 1,
          formula_rec = ~ 1 + x_rec,
          formula_surv = ~ 1 + x_surv,
          n_ages = n_ages,
          .toggles = sim_toggles,
          pars = true_pars)
}

fit_one <- function(sim_obj, seed,
                    chains = 4,
                    parallel_chains = 4,
                    iter_warmup = 1000,
                    iter_sampling = 1000) {
  fit_drm(.data = sim_obj$data,
          y_col = "dens",
          time_col = "year",
          site_col = "patch",
          family = "gamma",
          formula_zero = ~ 1,
          formula_rec = ~ 1 + x_rec,
          formula_surv = ~ 1 + x_surv,
          n_ages = n_ages,
          .toggles = sim_toggles,
          seed = seed,
          algo_args = list(chains = chains,
                           parallel_chains = parallel_chains,
                           iter_warmup = iter_warmup,
                           iter_sampling = iter_sampling,
                           show_messages = FALSE,
                           show_exceptions = FALSE))
}

##--- single synthetic case ----

single_seed <- 202610
set.seed(single_seed)
x_rec  <- rnorm(n_tot)
x_surv <- rnorm(n_tot)
sim_single  <- simulate_one(seed = single_seed,
                            x_r = x_rec, x_s = x_surv)
fit_single  <- fit_one(sim_single, seed = single_seed)

##--- * convergence check ----

mcmc_diag(fit_single) |>
  print() |>
  summary()

##--- * parameter estimates ----

est_single <-
  summary(fit_single,
          probs = c(.025, .05, .5, .95, .975))$estimates |>
  left_join(truth_df, by = "variable")

print(est_single)

##--- * Table for supplement ----

est_single |>
  mutate(CI95 = sprintf("[%.2f, %.2f]", q2.5, q97.5),
         CI90 = sprintf("[%.2f, %.2f]", q5, q95)) |>
  select(variable, description, true_val,
         q50, CI90, CI95) |>
  print() |>
  xtable::xtable(caption = "Parameter estimates for the single synthetic dataset",
                 digits = 2) |>
  print(include.rownames = FALSE)

##--- * relationships with the environment ----

eff_rec <-
  effects_drm(fit_single,
              process = "rec",
              variable = "x_rec",
              prob = .9) |>
  mutate(true_curve = exp(true_pars$beta_r[1] +
                            true_pars$beta_r[2] * x_rec))

eff_surv <-
  effects_drm(fit_single,
              process = "surv",
              variable = "x_surv",
              prob = .9) |>
  mutate(true_curve = plogis(true_pars$beta_s[1] +
                               true_pars$beta_s[2] * x_surv))

fig_rec_curve <-
  ggplot(data = eff_rec,
         aes(x = x_rec)) +
  geom_ribbon(aes(ymin = q5, ymax = q95),
              fill = "steelblue",
              alpha = .3) +
  geom_line(aes(y = q50,
                color = "Estimated (median & 90% CI)"),
            linewidth = .9) +
  geom_line(aes(y = true_curve,
                color = "True relationship"),
            lty = 2,
            linewidth = .9) +
  scale_color_manual(name = NULL,
                     values = c("Estimated (median & 90% CI)" = "steelblue4",
                                "True relationship" = "firebrick")) +
  labs(x = expression(x^{(r)}~"(Recruitment covariate)"),
       y = "Recruitment density") +
  theme_bw() +
  theme(legend.position = "top")

fig_surv_curve <-
  ggplot(data = eff_surv,
         aes(x = x_surv)) +
  geom_ribbon(aes(ymin = q5, ymax = q95),
              fill = "steelblue",
              alpha = .3) +
  geom_line(aes(y = q50,
                color = "Estimated (median & 90% CI)"),
            linewidth = .9) +
  geom_line(aes(y = true_curve,
                color = "True relationship"),
            lty = 2,
            linewidth = .9) +
  scale_color_manual(name = NULL,
                     values = c("Estimated (median & 90% CI)" = "steelblue4",
                                "True relationship" = "firebrick")) +
  labs(x = expression(x^{(s)}~"(Survival covariate)"),
       y = "Survival probability") +
  theme_bw() +
  theme(legend.position = "top")

leg <- cowplot::get_legend(fig_rec_curve)

##--- * panel for single synthetic case ----

gratio <- 0.5 * (1 + sqrt(5))

panels <-
  plot_grid(fig_rec_curve + theme(legend.position = "none"),
            fig_surv_curve + theme(legend.position = "none"),
            labels = "AUTO",
            label_size = 10,
            ncol = 2)

fig_single_all <-
  plot_grid(leg, panels,
            ncol = 1,
            rel_heights = c(0.08, 1))

ggsave(filename = "overleaf/img/sim_single_recovery.pdf",
       plot = fig_single_all,
       width = 7,
       height = .75 * 7 / gratio)

##--- monte carlo study (100 replicates) ----

n_reps <- 100
mc_results_file <- "data/sim_recovery_100reps.rds"

if (file.exists(mc_results_file)) {
  mc_draws_df <- readRDS(mc_results_file)
} else {
  mc_list <- vector("list", n_reps)

  for (r in seq_len(n_reps)) {
    rep_seed <- 10000 + r
    sim_r <- simulate_one(seed = rep_seed,
                          x_r = x_rec, x_s = x_surv)
    fit_r <- fit_one(sim_r,
                     seed = rep_seed,
                     chains = 4,
                     parallel_chains = 4,
                     iter_warmup = 500,
                     iter_sampling = 500)

    diag_r <- mcmc_diag(fit_r)
    est_r  <- summary(fit_r,
                      probs = c(.025, .05, .5, .95, .975))$estimates

    mc_list[[r]] <-
      est_r |>
      left_join(select(diag_r, variable, rhat, ess_bulk),
                by = "variable") |>
      left_join(truth_df, by = "variable") |>
      mutate(rep = r,
             seed = rep_seed)

    if (r %% 10 == 0) {
      message(sprintf("Completed replicate %d / %d", r, n_reps))
    }
  }

  mc_draws_df <- bind_rows(mc_list)
  saveRDS(mc_draws_df, file = mc_results_file)
}

##--- * table summary ----

mc_summary_table <-
  mc_draws_df |>
  mutate(err = q50 - true_val,
         rel_bias = (q50 - true_val) / abs(true_val),
         cov_90 = (true_val >= q5) & (true_val <= q95),
         cov_95 = (true_val >= q2.5) & (true_val <= q97.5),
         width_90 = q95 - q5) |>
  group_by(variable, description, true_val) |>
  summarise(mean_est = mean(q50),
            bias = mean(err),
            rmse = sqrt(mean(err * err)),
            avg_sd = mean(sd),
            emp_se = sd(q50),
            coverage_90 = mean(cov_90) * 100,
            coverage_95 = mean(cov_95) * 100,
            max_rhat = max(rhat, na.rm = TRUE),
            rhat_g_1 = 100 * mean(rhat > 1.01, na.rm = TRUE),
            .groups = "drop") |>
  print()

mc_summary_table |>
  xtable::xtable(caption = "Monte Carlo parameter-recovery summary",
                 digits = 3) |>
  print(include.rownames = FALSE)

##--- * coverage plot ----

fig_mc_coverage <-
  mc_draws_df |>
  filter(grepl("^beta_[rs]", variable)) |>
  mutate(covered_90 = (true_val >= q5) & (true_val <= q95)) |>
  group_by(variable) |>
  arrange(q50, .by_group = TRUE) |>
  mutate(rank_id = row_number()) |>
  ungroup() |>
  ggplot(data = _,
         aes(x = rank_id,
             y = q50,
             color = covered_90)) +
  geom_hline(aes(yintercept = true_val),
             lty = 2,
             color = "black",
             linewidth = .6) +
  geom_linerange(aes(ymin = q5, ymax = q95),
                 alpha = .7,
                 linewidth = .5) +
  geom_point(size = .8) +
  facet_wrap(~ param_label,
             scales = "free_y",
             labeller = label_parsed) +
  scale_color_manual(name = "90% CI covers truth",
                     values = c("TRUE" = "steelblue",
                                "FALSE" = "firebrick")) +
  labs(x = "Simulation replicate (ordered by posterior median)",
       y = "Posterior median & 90% CI") +
  theme_bw() +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "white"))

ggsave(filename = "overleaf/img/sim_mc_coverage.pdf",
       plot = fig_mc_coverage,
       width = 7,
       height = 5)
