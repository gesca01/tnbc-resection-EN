# Load libraries
library(brms)
library(dplyr)
library(posterior)

# Manually input data
df <- data.frame(
  time = c(
    5.714, 12, 3.429, 2.571, 12, 12, 12, 12, 12, 5.571, 12,     # Combination
    2, 4.857, 4.857, 12, 5, 2.429, 12, 4.857, 12,               # Untreated
    8.143, 8.143, 8.143, 6.714, 4.714, 8.143, 8.143, 3.571, 6, 8.143, # Combination
    7.857, 8.143, 8.143, 7, 4, 2.143, 8.143, 8.143              # Untreated
  ),
  event = c(
    1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0,
    1, 1, 1, 0, 1, 1, 0, 1, 0,
    0, 0, 0, 1, 1, 0, 0, 1, 1, 0,
    1, 0, 0, 1, 1, 1, 0, 0
  ),
  group = c(
    rep("Combination", 11),
    rep("Untreated", 9),
    rep("Combination", 10),
    rep("Untreated", 8)
  )
)

# Replace survival times of censored mice with 30
df$time[df$event == 0] <- 30

# Fit Bayesian Cox model using brms
fit_bayes <- brm(
  formula = bf(time | cens(1 - event) ~ group),
  data = df,
  family = brmsfamily("cox"),
  chains = 4,
  iter = 4000,
  seed = 42,
  refresh = 0
)

# Extract posterior draws
post_draws <- as_draws_df(fit_bayes)

# Compute posterior probability that combination improves survival (HR < 1)
posterior_prob <- mean(post_draws$b_groupUntreated > 0)

# Print result
cat("\nPosterior probability that combination improves survival (HR < 1):\n")
print(posterior_prob)
