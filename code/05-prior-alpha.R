library(ggplot2)

alpha <- seq(from = -.99, to = .99, len = 201)

normal_tanh <- function(x, mean = 0, sd = 1, log = FALSE) {
  stopifnot(all(abs(x) <= 1))
  theta <- atanh(x)
  
  lp_theta <- dnorm(theta, mean = mean, sd = sd, log = TRUE) - log1p(- x * x)
  
  if (!log) return(exp(lp_theta))
  return(lp_theta)
}

pc_base0 <- function(x, u, alpha_param, log = FALSE) {
  stopifnot(all(abs(x) <= 1))
  theta <- - log(alpha_param) / sqrt(- log1p(- u * u))
  log1mx <- log1p(- x * x)
  lp_tx <- rep(log(theta) - log(2), length(x))
  idx <- x != 0
  if (any(idx)) {
    x_sub <- x[idx]
    log1mx_sub <- log1mx[idx] 
    lp_tx[idx] <- 
      log(theta) - log(2) -
      theta * sqrt(-log1mx_sub) + 
      log(abs(x_sub)) -
      log1mx_sub - 0.5 * log(-log1mx_sub)
  } 
  if (!log) return(exp(lp_tx))
  return(lp_tx)
}

trans_beta <- function(x, a = 1, b = 1, log = FALSE) {
  theta <- 0.5 * (x + 1)
  lp_theta <- dbeta(theta, shape1 = a, shape2 = b, log = TRUE) - log(2)
  
  if (!log) return(exp(lp_theta))
  return(lp_theta)
}

dt_aux <-
  data.frame(alpha = alpha) |>
  transform(atanh = normal_tanh(alpha, sd = 0.25),
            pc    = pc_base0(alpha, 0.9, 0.05),
            tbeta = trans_beta(alpha, 0.5, 0.5),
            beta  = ifelse(alpha >= 0, dbeta(alpha, 0.5, 0.5), NA_real_)) # FIX: Typo corrected

tidyr::pivot_longer(dt_aux,
                    cols = atanh:beta,
                    names_to = "prior",
                    values_to = "density") |>
  ggplot(data = _,
         aes(x = alpha,
             y = density,
             color = prior)) +
  geom_line(lwd = 1)
