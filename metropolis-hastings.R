# =================================================================================================
# METROPOLIS-HASTINGS ALGORITHM
# =================================================================================================

# Metropolis-Hastings algorithm:
#   - if the posterior is given, sample directly from it
#   - if no posterior is given, compute it using the likelihood and the prior
metropolis_hastings <- function(theta_init, posterior, proposal, max_iter, data, log_lik, prior) {

  theta <- theta_init
  dims <- length(theta_init)

  results <- list(
    theta = matrix(NA_real_, nrow = max_iter, ncol = dims),
    proposed_theta = matrix(NA_real_, nrow = max_iter, ncol = dims),
    acceptance_prob = rep(NA_real_, max_iter),
    decision = rep("rejected", max_iter)
  )

  for (i in seq_len(max_iter)) {

    results$proposed_theta[i, ] <- theta_proposed <- proposal(theta)

    if (!missing(posterior)) {
      results$acceptance_prob[i] <- acceptance_prob <- min(1, posterior(theta_proposed) / posterior(theta))
    } else {
      log_posterior <- log_lik(theta, data) + log(prior(theta))
      log_posterior_proposed <- log_lik(theta_proposed, data) + log(prior(theta_proposed))
      results$acceptance_prob[i] <- acceptance_prob <- min(1, exp(log_posterior_proposed - log_posterior))
    }

    if (runif(1) <= acceptance_prob) {
      theta <- theta_proposed
      results$decision[i] <- "accepted"
    }

    results$theta[i, ] <- theta
  }

  results
}

# =================================================================================================
# HELPER FUNCTIONS FOR GAUSSIAN MIXTURE MODEL EXAMPLE
# =================================================================================================

gmm_posterior <- Vectorize(function(x, weights = c(0.4, 0.6), mean = c(-1, 2), sd = c(0.5, 2)) {
  weights[1] * dnorm(x, mean[1], sd[1]) + weights[2] * dnorm(x, mean[2], sd[2])
})

gaussian_proposal_1d <- function(x, sd = 4) {
  rnorm(1, x, sd)
}

# =================================================================================================
# HELPER FUNCTIONS FOR SUNSPOT EXAMPLE
# =================================================================================================

unif_prior <- function(x) {
  alpha <- x[1]
  beta <- x[2]
  ifelse(alpha <= 0 || beta <= 0, 0, 1)
}

gamma_log_lik <- function(x, data) {
  alpha <- x[1]
  beta <- x[2]
  n <- length(data)
  if (alpha <= 0 || beta <= 0) {
    -Inf
  } else {
    (alpha - 1) * sum(log(data)) - (1 / beta) * sum(data) - n * log(gamma(alpha)) - n * alpha * log(beta)
  }
}

gaussian_proposal_2d <- function(x, sd = c(0.05, 5)) {
  rnorm(2, mean = x, sd = sd)
}
