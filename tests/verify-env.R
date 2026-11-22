cat("=== Environment info ===\n")
cat("R:", R.version.string, "\n")
cat("OS:", Sys.info()[["sysname"]], Sys.info()[["machine"]], "\n")

cat("\n=== Checking renv lockfile status ===\n")
st <- renv::status()
if (!isTRUE(st$synchronized)) {
  stop("renv::status() reported an inconsistent environment.")
}

cat("\n=== Loading locked manuscript packages ===\n")
required_pkgs <- c(
  "drmr",
  "cmdstanr",
  "instantiate",
  "sf",
  "spdep",
  "units",
  "dplyr",
  "tidyr",
  "ggplot2",
  "bayesplot",
  "posterior",
  "cowplot",
  "scales",
  "xtable"
)
for (pkg in required_pkgs) {
  library(pkg, character.only = TRUE)
}

cat("\n=== Checking CmdStan and drmr pre-compiled Stan models ===\n")
if (!isTRUE(instantiate::stan_cmdstan_exists())) {
  stop("instantiate::stan_cmdstan_exists() returned FALSE.")
}
cmdstan_ver <- cmdstanr::cmdstan_version()
cat("CmdStan version:", cmdstan_ver, "\n")
cat("CmdStan path:", cmdstanr::cmdstan_path(), "\n")
if (package_version(cmdstan_ver) < "2.36") {
  stop("CmdStan version must be >= 2.36, found: ", cmdstan_ver)
}
expected_cmdstan <- Sys.getenv("CMDSTAN_VERSION", unset = "")
if (nzchar(expected_cmdstan) && !identical(as.character(cmdstan_ver), expected_cmdstan)) {
  stop("Expected CmdStan ", expected_cmdstan, ", found: ", cmdstan_ver)
}

stan_dir <- system.file("bin/stan", package = "drmr")
if (!nzchar(stan_dir) || !dir.exists(stan_dir)) {
  stop("drmr bin/stan directory not found.")
}
ext <- if (.Platform$OS.type == "windows") ".exe" else ""
model_names <- c(
  "adults_drm",
  "drm",
  "fitted_drm",
  "fitted_sdm",
  "lambda_gq",
  "loglik_drm",
  "loglik_sdm",
  "predict",
  "predict_sdm",
  "sdm"
)
expected_bins <- paste0(model_names, ext)
missing_bins <- expected_bins[!file.exists(file.path(stan_dir, expected_bins))]
if (length(missing_bins) > 0L) {
  stop("Missing pre-compiled drmr Stan executables: ", paste(missing_bins, collapse = ", "))
}
cat("Verified all", length(expected_bins), "pre-compiled drmr Stan executables.\n")

cat("\n=== Checking spatial & manuscript data I/O ===\n")
data(sum_fl, package = "drmr")
stopifnot(is.data.frame(sum_fl), nrow(sum_fl) > 0L)

map_path <- system.file("maps/sum_fl.shp", package = "drmr")
stopifnot(nzchar(map_path), file.exists(map_path))
polygons <- sf::st_read(map_path, quiet = TRUE)
areas <- polygons |>
  sf::st_area() |>
  units::set_units("km^2")
stopifnot(length(areas) == nrow(polygons), all(as.numeric(areas) > 0))
adj_mat <- drmr::gen_adj(sf::st_buffer(sf::st_geometry(polygons), dist = 2500))
stopifnot(is.matrix(adj_mat), nrow(adj_mat) == nrow(polygons))

birds_dt <- readRDS("data/birds/processed.rds")
stopifnot(is.data.frame(birds_dt), nrow(birds_dt) > 0L)
birds_grid <- sf::st_read("data/birds/shape/grid.shp", quiet = TRUE)
stopifnot(nrow(birds_grid) > 0L)
sim_mc <- readRDS("data/sim_recovery_100reps.rds")
stopifnot(is.data.frame(sim_mc), nrow(sim_mc) > 0L)

cat("\n=== Running synthetic drmr + Stan smoke test ===\n")
set.seed(2026)
n_patches <- 4L
n_years <- 10L
n_ages <- 3L
grid <- expand.grid(year = seq_len(n_years), patch = seq_len(n_patches)) |>
  dplyr::mutate(
    x_rec = rnorm(n_patches * n_years),
    x_surv = rnorm(n_patches * n_years)
  )

sim_toggles <- list(est_surv = 1, minit = 1, rho_mu = 0)
true_pars <- list(
  beta_r = c(1.0, 0.5),
  beta_s = c(0.8, -0.5),
  beta_t = c(-2.0),
  phi = 10
)

sim_obj <- drmr::sim_drm(
  .data = grid,
  time_col = "year",
  site_col = "patch",
  y_col = "dens",
  family = "gamma",
  formula_zero = ~ 1,
  formula_rec = ~ 1 + x_rec,
  formula_surv = ~ 1 + x_surv,
  n_ages = n_ages,
  .toggles = sim_toggles,
  pars = true_pars
)

train_df <- dplyr::filter(sim_obj$data, year <= 8L)
test_df <- dplyr::filter(sim_obj$data, year > 8L)
past_df <- dplyr::filter(train_df, year == 8L)

smoke_algo <- list(
  chains = 2L,
  parallel_chains = 2L,
  iter_warmup = 50L,
  iter_sampling = 50L,
  show_messages = FALSE,
  show_exceptions = FALSE
)

drm_fit <- drmr::fit_drm(
  .data = train_df,
  y_col = "dens",
  time_col = "year",
  site_col = "patch",
  family = "gamma",
  formula_zero = ~ 1,
  formula_rec = ~ 1 + x_rec,
  formula_surv = ~ 1 + x_surv,
  n_ages = n_ages,
  .toggles = sim_toggles,
  seed = 2026L,
  algo_args = smoke_algo
)
stopifnot(inherits(drm_fit, "adrm"))

diag_df <- drmr::mcmc_diag(drm_fit)
stopifnot(is.data.frame(diag_df), nrow(diag_df) > 0L)
summ_drm <- summary(drm_fit)
stopifnot(is.data.frame(summ_drm$estimates), nrow(summ_drm$estimates) > 0L)

fitted_drm <- fitted(drm_fit) |>
  summary()
stopifnot(is.data.frame(fitted_drm), nrow(fitted_drm) == nrow(train_df))

pred_drm <- predict(
  drm_fit,
  newdata = test_df,
  past_data = past_df,
  seed = 2026L
) |>
  summary()
stopifnot(is.data.frame(pred_drm), nrow(pred_drm) == nrow(test_df))

sdm_fit <- drmr::fit_sdm(
  .data = train_df,
  y_col = "dens",
  time_col = "year",
  site_col = "patch",
  family = "gamma",
  formula_zero = ~ 1,
  formula_dens = ~ 1 + x_rec,
  seed = 2026L,
  algo_args = smoke_algo
)
stopifnot(inherits(sdm_fit, "sdm"))

fitted_sdm <- fitted(sdm_fit) |>
  summary()
stopifnot(is.data.frame(fitted_sdm), nrow(fitted_sdm) == nrow(train_df))

pred_sdm <- predict(sdm_fit, newdata = test_df, seed = 2026L) |>
  summary()
stopifnot(is.data.frame(pred_sdm), nrow(pred_sdm) == nrow(test_df))

eff_rec <- drmr::effects_drm(drm_fit, process = "rec", variable = "x_rec")
tmp_pdf <- tempfile(fileext = ".pdf")
on.exit(unlink(tmp_pdf), add = TRUE)
p <- plot(eff_rec) + ggplot2::theme_bw()
ggplot2::ggsave(filename = tmp_pdf, plot = p, width = 5, height = 4)
stopifnot(file.exists(tmp_pdf), file.info(tmp_pdf)$size > 0)

cat("\nAll environment and Stan smoke tests passed!\n")
