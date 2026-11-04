library(ggplot2)

alpha <- seq(from = -.99, to = .99, len = 201)

normal_tanh <- function(x, mean = 0, sd = 1, log = FALSE) {
  stopifnot(all(abs(x) <= 1))
  theta <- atanh(x)
  lp_theta <- dnorm(theta, mean = mean, sd = sd, log = TRUE) - log1p(- x * x)  
  if (!log) return(exp(lp_theta))
  return(lp_theta)
}

tidyr::expand_grid(alpha = alpha,
                   mean  = c(-1, 0, 1),
                   sd    = c(0.25, 0.75, 1)) |>
  transform(density =
              normal_tanh(alpha, mean = mean, sd = sd)) |>
  dplyr::rename("prior mean" = "mean",
                "prior sd"   = "sd") |>
  ggplot(aes(x = alpha, y = density)) +
  geom_line(lwd = 1) +
  facet_grid(`prior sd` ~ `prior mean`,
             labeller = label_both) +
  labs(x = expression(alpha),
       y = expression(pi(alpha))) +
  theme_bw()


