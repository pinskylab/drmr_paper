library(dplyr)
library(ggplot2)
library(bayesplot)
library(sf)
library(drmr)
library(cowplot)
library(cmdstanr)

bayesplot::color_scheme_set(scheme = "mix-pink-teal")

## loading data
my_dt <- readRDS("data/birds/processed.rds")

my_map <- st_read("data/birds/shape/grid.shp") |>
  mutate(lon = st_coordinates(st_centroid(geometry))[, 1],
         lat = st_coordinates(st_centroid(geometry))[, 2],
         .before = "geometry")

my_dt <- my_dt |>
  mutate(id = as.integer(factor(id)),) |>
  arrange(id, year) |>
  left_join(select(st_drop_geometry(my_map),
                 - area), by = "id")

polygons <- my_map |>
  st_geometry()

polygons |>
  st_area() |>
  units::set_units("km^2") |>
  summary()

##--- splitting data for validation ----

## reserving 5 years for forecast assessment
first_year_forecast <- max(my_dt$year) - 5

## "year to id"
first_id_forecast <-
  first_year_forecast - min(my_dt$year) + 1

years_all <- order(unique(my_dt$year))
years_train <- years_all[years_all < first_id_forecast]
years_test <- years_all[years_all >= first_id_forecast]

dat_test <- my_dt |>
  filter(year >= first_year_forecast)

dat_train <- my_dt |>
  filter(year < first_year_forecast)

##--- centering covariates (for improved mcmc efficiency) ---

avgs <- c("tavg" = mean(dat_train$tavg),
          "lon" = mean(dat_train$lon),
          "lat" = mean(dat_train$lat))

min_year <- dat_train$year |>
  min()

## centering covariates
dat_train <- dat_train |>
  mutate(c_tavg = tavg - avgs["tavg"],
         c_lat  = lat - avgs["lat"],
         c_lon  = lon - avgs["lon"],
         time   = year - min_year)

dat_test <- dat_test |>
  mutate(c_tavg = tavg - avgs["tavg"],
         c_lat  = lat - avgs["lat"],
         c_lon  = lon - avgs["lon"],
         time   = year - min_year)

##--- turning response into density: 1k individuals per km2 ----

dat_train <- dat_train |>
  mutate(dens = 100 * y,
         .before = y)

dat_test <- dat_test |>
  mutate(dens = 100 * y,
         .before = y)

chains <- 4
cores <- 4

##--- viz ap ----

my_map |>
  left_join(filter(dat_train, year == 2011),
            by = "id") |>
  ggplot(data = _,
         aes(fill = dens)) +
  geom_sf(alpha = .9) +
  scale_fill_viridis_c(option = "H") +
  theme_bw()

##--- fitting DRMs ----

adj_mat <- gen_adj(st_geometry(polygons))

## row-standardized matrix
adj_mat <-
  t(apply(adj_mat, 1, \(x) x / (sum(x))))

##--- * Recruitment ----

n_ages <- 12

## algo <- "laplace" ## change the inference method
## algo_args <- list() ## change inference algorithm defaults

algo <- "nuts" ## change the inference method
algo_args <- list(parallel_chains = 4,
                  chains = 4) ## change inference algorithm defaults

drm_rec <-
  fit_drm(.data = dat_train,
          y_col = "dens", ## response variable: density
          time_col = "year", ## vector of time points
          site_col = "id",
          family = "gamma",
          seed = 202505,
          formula_zero = ~ 1 + n_routes,
          formula_rec = ~ 1 + c_tavg + I(c_tavg * c_tavg),
          formula_surv = ~ 1,
          n_ages = 12,
          adj_mat = adj_mat, ## A matrix for movement routine
          ## init_data = init_ldens,
          ages_movement = c(0, 0,
                            rep(1, n_ages - 4),
                            0, 0), ## ages allowed to move
          .toggles = list(ar_re = "rec",
                          sp_re = "surv",
                          movement = 1,
                          est_surv = 1,
                          est_init = 0,
                          minit = 1),
          .priors = list(pr_phi_a = 1, pr_phi_b = .1,
                         pr_alpha_a = 4.2, pr_alpha_b = 5.8,
                         pr_zeta_a = 7, pr_zeta_b = 3),
          algo = algo,
          algo_args = algo_args)
##--- Convergence & estimates ----

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
                linewidth = 1.2)

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

draws(drm_rec, variables = c("sigma_s")) |>
  mcmc_dens_overlay() +
  stat_function(fun = \(x) dlnorm(x,
                                  meanlog = drm_rec$data$pr_lsigma_s_mu,
                                  sdlog = drm_rec$data$pr_lsigma_s_sd),
                xlim = c(0, .5),
                n = 501,
                inherit.aes = FALSE,
                color = 2,
                lwd = 1.2
)

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
         formula_surv = ~ 1 + c_tavg + I(c_tavg * c_tavg))

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

##--- Projections ----

##--- * DRM ----

proj_rec <- predict(drm_rec,
                    new_data = dat_test,
                    past_data = filter(dat_train,
                                       year == max(year)),
                    seed = 125,
                    cores = 4) |>
  summary(probs = c(.05, .5, .95))

proj_surv <- predict(drm_surv,
                     new_data = dat_test,
                     past_data = filter(dat_train,
                                        year == max(year)),
                     seed = 125,
                     cores = 4) |>
  summary(probs = c(.05, .5, .95))

##--- Viz predicted and observed ----

ci_mass <- .9
tails <- 1 - ci_mass
l_t <- tails * .5
u_t <- 1 - .5 * tails

fitted_rec <-
  fitted(drm_rec) |>
  summary(probs = c(.05, .5, .95))

fitted_surv <-
  fitted(drm_surv) |>
  summary(probs = c(.05, .5, .95))

##--- Figure not included in the manuscript ----

set.seed(2026)
sample_ids <- sample(unique(my_dt$id), 10)

bind_rows(
    bind_rows(fitted_rec, proj_rec) |>
    mutate(model = "DRM (rec)"),
    bind_rows(fitted_surv, proj_surv) |>
    mutate(model = "DRM (surv)")
) |>
  mutate(id = as.integer(id)) |>
  filter(id %in% sample_ids) |>
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
  geom_point(data = filter(bind_rows(dat_train, dat_test),
                           id %in% sample_ids),
             aes(x = year, y = dens), size = .5) +
  facet_grid(id ~ model, scales = "free_y") +
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity,
                                                   n = 3),
                     trans = "log1p") +
  theme_bw() +
  guides(color = "none",
         fill = "none") +
  labs(y = "Density (in hundreds of individuals per square-km)",
       x = "Year") +
  theme(strip.background = element_rect(fill = "white"))

ggsave(filename = "overleaf/img/forecast_rbw.pdf",
       width = 6,
       height = 7)

##--- * Table 3 ----

bind_rows(
    bind_rows(fitted_rec, proj_rec) |>
    mutate(model = "DRM (rec)"),
    bind_rows(fitted_surv, proj_surv) |>
    mutate(model = "DRM (surv)")
) |>
  mutate(id = as.integer(id)) |>
  left_join(bind_rows(dat_train, dat_test),
            by = c("id", "year")) |>
  mutate(type = ifelse(year < first_year_forecast, "in-sample",
                       "out-of-sample")) |>
  mutate(bias = dens - q50) |>
  mutate(rmse = bias * bias) |>
  mutate(is = int_score(dens, l = q5, u = q95, alpha = .2)) |>
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
  effects_drm(drm_rec,
              process = "rec",
              variable = "c_tavg") |>
  plot() +
  scale_x_continuous(labels = \(x) round(x + avgs["tavg"], 1),
                     breaks = c(10, 15, 20, 25) - avgs["tavg"]) +
  theme_bw() +
  labs(x = "Averate air-temperature (in C)",
       y = "Est. Recruitment (per km2)")

rec_fig

gratio <- 0.5 * (1 + sqrt(5))

ggsave(filename = "overleaf/img/recruitment_rbw.pdf",
       plot = rec_fig,
       width = 6,
       height = 6 / gratio)

##--- surv and environment ----

surv_fig <-
  effects_drm(drm_surv,
              process = "surv",
              variable = "c_tavg") |>
  plot() +
  scale_x_continuous(labels = \(x) round(x + avgs["tavg"], 1),
                     breaks = c(10, 15, 20, 25) - avgs["tavg"]) +
  theme_bw() +
  labs(x = "Averate air-temperature (in C)",
       y = "Est. survival rates")

surv_fig

gratio <- 0.5 * (1 + sqrt(5))

ggsave(filename = "overleaf/img/surv_rbw.pdf",
       plot = surv_fig,
       width = 6,
       height = 6 / gratio)

##--- Figure 3 ----
##--- Panel for recruitment and survival ----

plot_grid(rec_fig, surv_fig,
          labels = "AUTO",
          label_size = 10)

ggsave(filename = "overleaf/img/rec_surv_rbw.pdf",
       width = 7,
       height = .75 * 7 / gratio)

##--- ** SST that optimizes recruitment ----

betas_rec <-
  draws(drm_rec, variables = "beta_r",
        format = "matrix")

max_quad_x(betas_rec[, 2], betas_rec[, 3],
           offset = avgs["tavg"]) |>
  apply(2, quantile, probs = c(.05, .5, .95))

##--- ** SBT that maximizes surival ----

betas_surv <-
  draws(drm_surv, variables = "beta_s",
        format = "matrix")

max_quad_x(betas_surv[, 2], betas_surv[, 3],
           offset = avgs["tavg"]) |>
  apply(2, quantile, probs = c(.05, .5, .95))
