library(dplyr)
library(ggplot2)
library(bayesplot)
library(sf)
library(drmr)
library(cowplot)
library(cmdstanr)

bayesplot::color_scheme_set(scheme = "mix-pink-teal")

## loading data
data(sum_fl)

## loading map
map_name <- system.file("maps/sum_fl.shp", package = "drmr")

polygons <- st_read(map_name)

polygons |>
  st_area() |>
  units::set_units("km^2") |>
  summary()

##--- splitting data for validation ----

## reserving 5 years for forecast assessment
first_year_forecast <- max(sum_fl$year) - 4

## "year to id"
first_id_forecast <-
  first_year_forecast - min(sum_fl$year) + 1

years_all <- order(unique(sum_fl$year))
years_train <- years_all[years_all < first_id_forecast]
years_test <- years_all[years_all >= first_id_forecast]

dat_test <- sum_fl |>
  filter(year >= first_year_forecast)

dat_train <- sum_fl |>
  filter(year < first_year_forecast)

##--- centering covariates (for improved mcmc efficiency) ---

avgs <- c("stemp" = mean(dat_train$stemp),
          "btemp" = mean(dat_train$btemp),
          "depth" = mean(dat_train$depth),
          "n_hauls" = mean(dat_train$n_hauls),
          "lat" = mean(dat_train$lat),
          "lon" = mean(dat_train$lon))

min_year <- dat_train$year |>
  min()

## centering covariates
dat_train <- dat_train |>
  mutate(c_stemp = stemp - avgs["stemp"],
         c_btemp = btemp - avgs["btemp"],
         c_hauls = n_hauls - avgs["n_hauls"],
         ## depth = depth - avgs["depth"],
         c_lat   = lat - avgs["lat"],
         c_lon   = lon - avgs["lon"],
         time  = year - min_year)

dat_test <- dat_test |>
  mutate(c_stemp = stemp - avgs["stemp"],
         c_btemp = btemp - avgs["btemp"],
         c_hauls = n_hauls - avgs["n_hauls"],
         ## depth = depth - avgs["depth"],
         c_lat   = lat - avgs["lat"],
         c_lon   = lon - avgs["lon"],
         time  = year - min_year)

##--- turning response into density: 1k individuals per km2 ----

dat_train <- dat_train |>
  mutate(dens = 100 * y / area_km2,
         .before = y)

dat_test <- dat_test |>
  mutate(dens = 100 * y / area_km2,
         .before = y)

chains <- 4
cores <- 4

##--- fitting DRMs ----

adj_mat <- gen_adj(st_buffer(st_geometry(polygons),
                             dist = 2500))

## instantaneous fishing mortality rates
fmat <-
  system.file("fmat.rds", package = "drmr") |>
  readRDS()

f_train <- fmat[, years_train]
f_test  <- fmat[, years_test]

## vizualizing different beta priors
mode_zeta <- .4
conc_zeta1 <- 10
conc_zeta2 <- 3.5

alpha1 <- mode_zeta * (conc_zeta1 - 2) + 1
alpha2 <- mode_zeta * (conc_zeta2 - 2) + 1

beta1 <- (1 - mode_zeta) * (conc_zeta1 - 2) + 1
beta2 <- (1 - mode_zeta) * (conc_zeta2 - 2) + 1

ggplot() +
  stat_function(fun =
                  \(x) dbeta(x, shape1 = alpha1, shape2 = beta1),
                color = 2) +
  stat_function(fun =
                  \(x) dbeta(x, shape1 = alpha2, shape2 = beta2),
                color = 4) +
  theme_bw()

mode_zeta <- .75
conc_zeta1 <- 10
conc_zeta2 <- 3.5

alpha1 <- mode_zeta * (conc_zeta1 - 2) + 1
alpha2 <- mode_zeta * (conc_zeta2 - 2) + 1

beta1 <- (1 - mode_zeta) * (conc_zeta1 - 2) + 1
beta2 <- (1 - mode_zeta) * (conc_zeta2 - 2) + 1

ggplot() +
  stat_function(fun =
                  \(x) dbeta(x, shape1 = alpha1, shape2 = beta1),
                color = 2) +
  stat_function(fun =
                  \(x) dbeta(x, shape1 = alpha2, shape2 = beta2),
                color = 4) +
  theme_bw()

##--- DRM recruitment ----

drm_rec <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "patch",
          family = "gamma",
          seed = 202505,
          formula_zero = ~ 1 + c_hauls,
          formula_rec = ~ 1 + c_stemp + I(c_stemp * c_stemp),
          formula_surv = ~ 1,
          f_mort = f_train,
          n_ages = NROW(f_train),
          adj_mat = adj_mat, ## A matrix for movement routine
          ages_movement = c(0, 0,
                            rep(1, 12),
                            0, 0), ## ages allowed to move
          .toggles = list(ar_re = "rec",
                          movement = 1,
                          est_surv = 1,
                          est_init = 0,
                          minit = 1),
          .priors = list(pr_alpha_a = 4.2, pr_alpha_b = 5.8,
                         pr_zeta_a = 7, pr_zeta_b = 3),
          algo_args = list(parallel_chains = 4))

##--- Convergence check ----

mcmc_diag(drm_rec) |>
  print() |>
  summary()

par(mfrow = c(4, 6))
plot(drm_rec)

par(mfrow = c(3, 4))
plot(drm_rec, type = "density")

par(mfrow = c(3, 4))
plot(drm_rec, type = "trace")

##--- parameter estimates ----

summary(drm_rec)
## specific quantiles
summary(drm_rec, probs = c(.1, .9))

##--- comparing some priors and posteriors ----

par(mfrow = c(1, 1))

##--- * Base R ----
plot(drm_rec, variables = "phi", type = "density")
curve(dgamma(x, 
             shape = drm_rec$data$pr_phi_a, 
             rate = drm_rec$data$pr_phi_b), 
      add = TRUE, 
      lty = 2,
      lwd = 2)

plot(drm_rec, variables = "alpha", type = "density")
curve(dbeta(x,
            shape1 = drm_rec$data$pr_alpha_a,
            shape2 = drm_rec$data$pr_alpha_b), 
      add = TRUE, 
      lty = 2,
      lwd = 2)

##--- * ggplot ----

draws(drm_rec, variables = c("zeta")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dbeta(x,
                                 shape1 = drm_rec$data$pr_zeta_a,
                                 shape2 = drm_rec$data$pr_zeta_b),
                xlim = c(0, 1),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

draws(drm_rec, variables = c("sigma_t")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dlnorm(x,
                                  meanlog = drm_rec$data$pr_lsigma_t_mu,
                                  sdlog = drm_rec$data$pr_lsigma_t_sd),
                xlim = c(0, .5),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2)

draws(drm_rec, variables = c("xi")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) {
    y <- - x
    dnorm(log(y),
          mean = drm_rec$data$pr_lmxi_mu,
          sd = drm_rec$data$pr_lmxi_sd) / y
  },
  xlim = c(-5, -1e-16),
  n = 501,
  inherit.aes = FALSE,
  color = 2,
  lwd = 1.2)

##--- DRM Survival ----

drm_surv <-
  update(drm_rec,
         formula_rec = ~ 1,
         formula_surv = ~ 1 + c_btemp + I(c_btemp * c_btemp))

##--- Convergence check ----

mcmc_diag(drm_surv) |>
  print() |>
  summary()

par(mfrow = c(3, 4))
plot(drm_surv, type = "density")

par(mfrow = c(3, 4))
plot(drm_surv, type = "trace")

##--- parameter estimates ----

summary(drm_surv)
## specific quantiles
summary(drm_surv, probs = c(.1, .9))

##--- An even more complex model ----

drm_rs <-
  update(drm_rec,
         formula_surv = ~ 1 + c_btemp + I(c_btemp * c_btemp))

##--- Convergence check ----

mcmc_diag(drm_rs) |>
  print() |>
  summary()

par(mfrow = c(3, 5))
plot(drm_rs, type = "density")

par(mfrow = c(3, 5))
plot(drm_rs, type = "trace")

##--- parameter estimates ----

summary(drm_rs)
## specific quantiles
summary(drm_rs, probs = c(.1, .9))

##--- Projections ----

##--- * DRM ----

proj_rec <- predict(drm_rec,
                    new_data = dat_test,
                    past_data = filter(dat_train,
                                       year == max(year)),
                    seed = 125,
                    f_test = f_test,
                    cores = 4) |>
  summary(probs = c(.05, .5, .95))

proj_surv <- predict(drm_surv,
                     new_data = dat_test,
                     past_data = filter(dat_train,
                                        year == max(year)),
                     seed = 125,
                     f_test = f_test,
                     cores = 4) |>
  summary(probs = c(.05, .5, .95))

proj_rs <- predict(drm_rs,
                   new_data = dat_test,
                   past_data = filter(dat_train,
                                      year == max(year)),
                   seed = 125,
                   f_test = f_test,
                   cores = 4) |>
  summary(probs = c(.05, .5, .95))

##--- Viz predicted and observed ----

ci_mass <- .8
tails <- 1 - ci_mass
l_t <- tails * .5
u_t <- 1 - .5 * tails

fitted_rec <-
  fitted(drm_rec) |>
  summary(probs = c(.05, .5, .95))

fitted_surv <-
  fitted(drm_surv) |>
  summary(probs = c(.05, .5, .95))

fitted_rs <-
  fitted(drm_rs) |>
  summary(probs = c(.05, .5, .95))

##--- Figure 2 ----

aux_fig <-
  bind_rows(dat_train, dat_test) |>
  mutate(patch = factor(as.integer(patch),
                        levels = rev(unique(as.integer(patch)))))
  

bind_rows(
    bind_rows(fitted_rec, proj_rec) |>
    mutate(model = "DRM (rec)"),
    bind_rows(fitted_surv, proj_surv) |>
    mutate(model = "DRM (surv)"),
    bind_rows(fitted_rs, proj_rs) |>
    mutate(model = "DRM (rec-surv)")
) |>
  mutate(patch = factor(as.integer(patch),
                        levels = rev(unique(as.integer(patch))))) |>
  ## filter(model != "DRM (surv)") |>
  ggplot(data = _) +
  geom_vline(xintercept = first_year_forecast,
             lty = 2) +
  geom_ribbon(aes(x = year,
                  ymin = q5, ymax = q95,
                  fill = model,
                  color = model),
              alpha = .4) +
  geom_line(aes(x = year, y = q50, color = model)) +
  geom_point(data = aux_fig,
             aes(x = year, y = dens), size = .5) +
  facet_grid(rows = patch ~ model, scales = "free_y") +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3),
                     trans = "log1p") +
  theme_bw() +
  guides(color = "none",
         fill = "none") +
  labs(y = "Density (in hundreds of individuals per square-km)",
       x = "Year") +
  theme(strip.background = element_rect(fill = "white"))

ggsave(filename = "overleaf/img/forecast_sf.pdf",
       width = 6,
       height = 7)

##--- * Table 3 ----

bind_rows(
    bind_rows(fitted_rec, proj_rec) |>
    mutate(model = "DRM (rec)"),
    bind_rows(fitted_surv, proj_surv) |>
    mutate(model = "DRM (surv)"),
    bind_rows(fitted_rs, proj_rs) |>
    mutate(model = "DRM (rec-surv)")
) |>
  mutate(patch = as.integer(patch)) |>
  left_join(bind_rows(dat_train, dat_test),
            by = c("patch", "year")) |>
  mutate(type = ifelse(year < first_year_forecast, "in-sample",
                       "out-of-sample")) |>
  mutate(bias = dens - q50) |>
  mutate(rmse = bias * bias) |>
  mutate(is = int_score(dens, l = q5, u = q95, alpha = .1)) |>
  ungroup() |>
  group_by(type, model) |>
  summarise(across(rmse:is, mean)) |>
  ungroup() |>
  rename_all(toupper) |>
  rename("Model" = "MODEL",
         "IS (90%)" = "IS") |>
  arrange(RMSE) |>
  print() |>
  xtable::xtable(caption = "Forecasting skill according to different metrics",
                 digits = 2) |>
  print(include.rownames = FALSE)

##--- Relationships with the environment ----

rec_fig <-
  effects_drm(drm_rs,
              process = "rec",
              variable = "c_stemp") |>
  plot() +
  scale_x_continuous(labels = \(x) round(x + avgs["stemp"], 1),
                     breaks = c(10, 15, 20, 25) - avgs["stemp"]) +
  theme_bw() +
  labs(x = "SST (in Celsius)",
       y = "Est. Recruitment (per km2)")

rec_fig

gratio <- 0.5 * (1 + sqrt(5))

ggsave(filename = "overleaf/img/recruitment.pdf",
       plot = rec_fig,
       width = 6,
       height = 6 / gratio)

##--- surv and environment ----

surv_fig <-
  effects_drm(drm_rs,
              process = "surv",
              variable = "c_btemp") |>
  plot() +
  scale_x_continuous(labels = \(x) round(x + avgs["btemp"], 1),
                     breaks = c(10, 15, 20, 25) - avgs["btemp"]) +
  theme_bw() +
  labs(x = "SST (in Celsius)",
       y = "Est. Survival Rates")

surv_fig

gratio <- 0.5 * (1 + sqrt(5))

ggsave(filename = "overleaf/img/surv.pdf",
       plot = surv_fig,
       width = 6,
       height = 6 / gratio)

##--- Figure 3 ----
##--- Panel for recruitment and survival ----

plot_grid(rec_fig, surv_fig,
          labels = "AUTO",
          label_size = 10)

ggsave(filename = "overleaf/img/rec_surv.pdf",
       width = 7,
       height = .75 * 7 / gratio)

##--- ** SST that optimizes recruitment ----

betas_rec <-
  draws(drm_rs, variables = "beta_r",
        format = "matrix")

max_quad_x(betas_rec[, 2], betas_rec[, 3],
           offset = avgs["stemp"]) |>
  apply(2, quantile, probs = c(.05, .5, .95))

##--- ** SBT that maximizes surival ----

betas_surv <-
  draws(drm_rs, variables = "beta_s",
        format = "matrix")

max_quad_x(betas_surv[, 2], betas_surv[, 3],
           offset = avgs["btemp"]) |>
  apply(2, quantile, probs = c(.05, .5, .95))
