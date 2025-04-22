# Load necessary libraries
library(glmmTMB)
library(dplyr)

# Example data - simulated scRNA-seq counts
set.seed(123)
n_cells <- 100
n_genes <- 10

# Create example data frame
data <- expand.grid(cell_id = 1:n_cells, gene_id = 1:n_genes)
data$treatment <- rep(c("control", "treated"), each = n_cells*n_genes/2)
data$batch <- rep(rep(c("batch1", "batch2"), each = n_cells/2), n_genes)

# Simulate counts with excess zeros
# True mean depends on treatment
data$true_mean <- ifelse(data$treatment == "treated", 5, 2)
# Add some batch effect
data$true_mean <- data$true_mean * ifelse(data$batch == "batch2", 1.5, 1)
# Zero-inflation probability - higher for low expression genes
data$pi <- 0.5 - 0.03 * data$true_mean
data$pi <- pmax(0.1, pmin(0.9, data$pi))
# Generate counts
data$is_zero <- rbinom(nrow(data), 1, data$pi)
data$count <- ifelse(data$is_zero == 1, 0, 
                     rnbinom(nrow(data), mu = data$true_mean, size = 1))

# Fit ZINB model
zinb_model <- glmmTMB(
  count ~ treatment + batch + (1|gene_id) + (1|cell_id),
  ziformula = ~ treatment + batch,
  family = nbinom2,
  data = data
)

# Extract model parameters for each observation
# 1. Get linear predictors for the count model (on log scale)
eta_count <- predict(zinb_model, type = "link")
# 2. Transform to get mu (negative binomial mean)
mu_hat <- exp(eta_count)

# 3. Get linear predictors for zero-inflation model (on logit scale)
eta_zero <- predict(zinb_model, type = "zlink")
# 4. Transform to get pi (zero-inflation probability)
pi_hat <- plogis(eta_zero)

# 5. Calculate expected count
expected_y <- (1 - pi_hat) * mu_hat

# Add these to our data frame
data$mu_hat <- mu_hat
data$pi_hat <- pi_hat
data$expected_y <- expected_y

# Examine results for the first few observations
head(data[, c("treatment", "batch", "count", "mu_hat", "pi_hat", "expected_y")])

# Compare average expected values with observed averages by treatment
data %>%
  group_by(treatment) %>%
  summarize(
    observed_mean = mean(count),
    expected_mean = mean(expected_y)
  )