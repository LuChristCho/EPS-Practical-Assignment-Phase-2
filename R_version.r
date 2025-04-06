# Additional cell code for installing libraries
#install.packages("moments")
#install.packages("ggplot2")
#install.packages("flexmix")
#install.packages("mixtools")
#install.packages("effsize")
#install.packages("dplyr")
library(dplyr)
library(effsize)
library(mixtools)
library(flexmix)
library(ggplot2)
library(moments)

# 1
set.seed(402100859)
samples <- rnorm(10000, mean = 5, sd = 2)

# 2
hist(samples, breaks = 30, probability = TRUE, col = "lightblue",
     main = "Histogram with Theoretical PDF", xlab = "Values")

curve(dnorm(x, mean = 5, sd = 2), add = TRUE, col = "red", lwd = 2)

# 3
mean_samples <- mean(samples)
variance_samples <- var(samples)
skewness_samples <- skewness(samples)
kurtosis_samples <- kurtosis(samples)

cat("Mean:", mean_samples, "\n")
cat("Variance:", variance_samples, "\n")
cat("Skewness:", skewness_samples, "\n")
cat("Kurtosis:", kurtosis_samples, "\n")

#4
cat("Expected Mean:", 5, "\tObserved Mean:", mean_samples, "\n")
cat("Expected Variance:", 4, "\tObserved Variance:", variance_samples, "\n")
cat("Expected Skewness:", 0, "\tObserved Skewness:", skewness_samples, "\n")
cat("Expected Kurtosis:", 3, "\tObserved Kurtosis:", kurtosis_samples, "\n")

theta <- seq(0, 1, length.out = 1000)

prior <- dbeta(theta, 2, 5)
likelihood <- dbinom(7, 20, theta)
posterior <- dbeta(theta, 9, 18)

df <- data.frame(theta, prior, likelihood, posterior)

ggplot(df, aes(x = theta)) +
  geom_line(aes(y = prior, color = "Prior"), linewidth = 1) +
  geom_line(aes(y = likelihood / max(likelihood), color = "Likelihood (scaled)"), linewidth = 1) +
  geom_line(aes(y = posterior, color = "Posterior"), linewidth = 1) +
  labs(title = "Prior, Likelihood, and Posterior Distributions",
       x = expression(theta),
       y = "Density") +
  scale_color_manual(values = c("Prior" = "blue", "Likelihood (scaled)" = "red", "Posterior" = "green")) +
  theme_minimal()

posterior_samples <- rbeta(10000, 9, 18)
posterior_mean <- mean(posterior_samples)
posterior_ci <- quantile(posterior_samples, c(0.025, 0.975))

cat("Posterior Mean:", posterior_mean, "\n")
cat("95% Confidence Interval: [", posterior_ci, "]\n")

#1
data(volcano)
volcano_normalized <- (volcano - min(volcano)) / (max(volcano) - min(volcano))

#2
volcano_vector <- as.vector(volcano_normalized)

set.seed(123)
gmm_model <- normalmixEM(volcano_vector, k = 3)

#3
component_assignment <- apply(gmm_model$posterior, 1, which.max)

#4
segmented_image <- matrix(component_assignment, nrow = nrow(volcano), ncol = ncol(volcano))
image(segmented_image, col = terrain.colors(3), axes = FALSE, main = "Segmented Volcano Image")

#1
data <- read.csv("https://drive.google.com/uc?id=1xHH5v_3le-8J5lttnL8hSFSmu8VyTis2")
summary(data)

prior_fault <- mean(data$True_Status == "Fault")
prior_no_fault <- 1 - prior_fault
prior_probabilities = c(prior_fault, prior_no_fault)
cat("\nPrior Fault: ", prior_fault, "\n")
cat("Prior No Fault: ", prior_no_fault, "\n")

#2
TPR <- sum(data$Alarm == "Active" & data$True_Status == "Fault") / sum(data$True_Status == "Fault")
FPR <- sum(data$Alarm == "Active" & data$True_Status == "No Fault") / sum(data$True_Status == "No Fault")
cat("True Positive Rate: ", TPR, "\n")
cat("False Positive Rate: ", FPR, "\n")

P_Alarm_given_Fault <- TPR
P_Alarm_given_NoFault <- FPR
P_NoAlarm_given_Fault <- 1 - TPR
P_NoAlarm_given_NoFault <- 1 - FPR
P_Alarm <- mean(data$Alarm == "Active")
P_NoAlarm <- 1 - P_Alarm

#3
posterior_fault_given_alarm <- (P_Alarm_given_Fault * prior_fault) / P_Alarm
cat("Posterior Fault Given Alarm: ", posterior_fault_given_alarm, "\n")

#4
posterior_no_fault_given_no_alarm <- (P_NoAlarm_given_NoFault * prior_no_fault) / P_NoAlarm
cat("Posterior No Fault Given No Alarm: ", posterior_no_fault_given_no_alarm, "\n\n")

#5
prior_list <- c(0.01, 0.1, 0.5)
posterior_results <- sapply(prior_list, function(p_fault) {
  p_no_fault <- 1 - p_fault
  p_fault_given_alarm <- (P_Alarm_given_Fault * p_fault) / ((P_Alarm_given_Fault * p_fault) + (P_Alarm_given_NoFault * p_no_fault))
  p_no_fault_given_no_alarm <- (P_NoAlarm_given_NoFault * p_no_fault) / ((P_NoAlarm_given_Fault * p_fault) + (P_NoAlarm_given_NoFault * p_no_fault))
  return(c(p_fault_given_alarm, p_no_fault_given_no_alarm))
})
colnames(posterior_results) <- paste0("Prior_", prior_list)
posterior_results

#6
compute_posterior <- function(subgroup_data) {
  prior_fault <- mean(subgroup_data$True_Status == "Fault")
  prior_no_fault <- 1 - prior_fault
  TPR <- sum(subgroup_data$Alarm == "Active" & subgroup_data$True_Status == "Fault") / sum(subgroup_data$True_Status == "Fault")
  FPR <- sum(subgroup_data$Alarm == "Active" & subgroup_data$True_Status == "No Fault") / sum(subgroup_data$True_Status == "No Fault")
  P_Alarm <- mean(subgroup_data$Alarm == "Active")
  P_NoAlarm <- 1 - P_Alarm
  P_NoAlarm_given_NoFault <- 1 - FPR
  P_Alarm_given_Fault <- TPR

  posterior_fault_given_alarm <- (P_Alarm_given_Fault * prior_fault) / P_Alarm
  posterior_no_fault_given_no_alarm <- (P_NoAlarm_given_NoFault * prior_no_fault) / P_NoAlarm
  return(c(posterior_fault_given_alarm, posterior_no_fault_given_no_alarm))
}

unique_ages <- unique(data$Machine_Age)
unique_maintenance <- unique(data$Maintenance_History)
results_list <- list()

for (age in unique_ages) {
  for (maintenance in unique_maintenance) {
    subgroup_data <- subset(data, Machine_Age == age & Maintenance_History == maintenance)
    posteriors <- compute_posterior(subgroup_data)

    results_list <- append(results_list, list(data.frame(
      Machine_Age = age,
      Maintenance_History = maintenance,
      P_Fault_Given_Alarm = posteriors[1],
      P_NoFault_Given_NoAlarm = posteriors[2]
    )))
  }
}
results <- do.call(rbind, results_list)
results

n_students <- 100
Student_ID <- 1:n_students
Group <- sample(c("Online", "In-Person"), n_students, replace = TRUE)
Score <- ifelse(Group == "Online", rnorm(n_students, mean = 70, sd = 10), rnorm(n_students, mean = 75, sd = 10))
data <- data.frame(Student_ID, Group, Score)

# 1
desc_stats <- data %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(Score),
    SD = sd(Score),
    Median = median(Score),
    Mode = as.numeric(names(sort(table(Score), decreasing = TRUE)[1])),  # Calculate mode
    Min = min(Score),
    Max = max(Score),
    N = n()
  )
print(desc_stats)

ggplot(data, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution of Scores by Group", x = "Group", y = "Score")

# 2
t_test_result <- t.test(Score ~ Group, data = data, var.equal = TRUE)
print(t_test_result)

# 3
t_statistic <- t_test_result$statistic
df <- t_test_result$parameter
p_value <- t_test_result$p.value
cat("T-statistic:", t_statistic, "\n")
cat("Degrees of freedom:", df, "\n")
cat("P-value:", p_value, "\n")

alpha <- 0.05
if (p_value < alpha) {
  cat("Conclusion: Reject the null hypothesis (H0).\n")
} else {
  cat("Conclusion: Fail to reject the null hypothesis (H0).\n")
}

# 4
cohen_d <- cohen.d(Score ~ Group, data = data)
cat("Cohen's d:", cohen_d$estimate, "\n")

# 1
movements <- list(
  c(1, 0),
  c(-1, 0),
  c(0, 1),
  c(0, -1),
  c(1, 1),
  c(-1, -1),
  c(1, -1),
  c(-1, 1),
  c(0, 0)
)

probabilities <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2)
simulate_random_walk <- function(n_steps) {
  position <- matrix(c(0, 0), nrow = 1, ncol = 2)
  for (i in 1:n_steps) {
    move <- sample(1:length(movements), 1, prob = probabilities)
    new_position <- position[i, ] + movements[[move]]
    position <- rbind(position, new_position)
  }
  return(position)
}

n_steps <- c(10, 100, 1000, 10000)
par(mfrow = c(2, 2))
for (n in n_steps) {
  path <- simulate_random_walk(n)
  plot(path[, 1], path[, 2], type = "o", col = "blue",
       xlab = "X", ylab = "Y", main = paste("Random Walk with", n, "Steps"),
       xlim = range(path[, 1]), ylim = range(path[, 2]))
}

# 2
simulate_random_walk_random_step <- function(n_steps) {
  position <- matrix(c(0, 0), nrow = 1, ncol = 2)
    for (i in 1:n_steps) {
    move <- sample(1:length(movements), 1, prob = probabilities)
    step_size <- runif(1, 1, 20)
    new_position <- position[i, ] + step_size * movements[[move]]
    position <- rbind(position, new_position)
  }
  return(position)
}

n_steps <- c(10, 100, 1000)
par(mfrow = c(2, 2))
for (n in n_steps) {
  path <- simulate_random_walk_random_step(n)
  plot(path[, 1], path[, 2], type = "o", col = "red",
       xlab = "X", ylab = "Y", main = paste("Random Walk with", n, "Steps"),
       xlim = range(path[, 1]), ylim = range(path[, 2]))
}

# 1
T0 <- 15
mu <- 0
sigma <- 3
T <- 365
n_paths <- 1000

temperature_paths <- matrix(NA, nrow = T, ncol = n_paths)
for (i in 1:n_paths) {
  daily_changes <- rnorm(T, mean = mu, sd = sigma)
  temperature_paths[, i] <- T0 + cumsum(daily_changes)
}

average_temperature <- rowMeans(temperature_paths)
sample_path <- temperature_paths[, 1]
plot(1:T, sample_path, type = "l", col = "blue",
     xlab = "Day", ylab = "Temperature (°C)",
     main = "Sample Temperature Path Over One Year",
     ylim = range(temperature_paths))

lines(1:T, average_temperature, col = "red", lwd = 2)

legend("topright", legend = c("Sample Path", "Average Path"),
       col = c("blue", "red"), lty = 1, lwd = 2)


# 2
temperature_paths_seasonal <- matrix(NA, nrow = T, ncol = n_paths)

for (i in 1:n_paths) {
  daily_changes <- rnorm(T, mean = mu, sd = sigma)
  cumulative_temperature <- T0 + cumsum(daily_changes)
  seasonal_trend <- 10 * sin(2 * pi * (1:T) / 365)
  temperature_paths_seasonal[, i] <- cumulative_temperature + seasonal_trend
}

average_temperature_seasonal <- rowMeans(temperature_paths_seasonal)
sample_path_seasonal <- temperature_paths_seasonal[, 1]
plot(1:T, sample_path_seasonal, type = "l", col = "blue",
     xlab = "Day", ylab = "Temperature (°C)",
     main = "Sample Temperature Path with Seasonal Trend Over One Year",
     ylim = range(temperature_paths_seasonal))

lines(1:T, average_temperature_seasonal, col = "red", lwd = 2)
legend("topright", legend = c("Sample Path", "Average Path"),
       col = c("blue", "red"), lty = 1, lwd = 2)


# 3
final_temperatures <- temperature_paths[T, ]
mean_final_temp <- mean(final_temperatures)
var_final_temp <- var(final_temperatures)

theoretical_mean <- T0 + T * mu
theoretical_var <- T * sigma^2

cat("Simulated Mean at t = 365:", mean_final_temp, "\n")
cat("Theoretical Mean at t = 365:", theoretical_mean, "\n")
cat("Simulated Variance at t = 365:", var_final_temp, "\n")
cat("Theoretical Variance at t = 365:", theoretical_var, "\n\n")


# 4
T_H <- 30
first_crossing_times <- numeric(n_paths)

for (i in 1:n_paths) {
  daily_changes <- rnorm(T, mean = mu, sd = sigma)
  temperature_path <- T0 + cumsum(daily_changes)
  crossing_index <- which(temperature_path >= T_H)[1]
  first_crossing_times[i] <- ifelse(is.na(crossing_index), NA, crossing_index)
}

average_transit_time <- mean(first_crossing_times, na.rm = TRUE)
cat("Average transit time to reach", T_H, "°C:", average_transit_time, "days\n")

hist(first_crossing_times, breaks = 30, col = "lightblue",
     xlab = "First Crossing Time (Days)", ylab = "Frequency",
     main = "Histogram of First Crossing Times to 30°C")


# 5
threshold <- -5
consecutive_days <- 3
warning_activated <- logical(n_paths)
consecutive_counts <- numeric(n_paths)

for (i in 1:n_paths) {
  daily_changes <- rnorm(T, mean = mu, sd = sigma)
  temperature_path <- T0 + cumsum(daily_changes)
  below_threshold <- temperature_path < threshold
  rle_below <- rle(below_threshold)
  warning_activated[i] <- any(rle_below$values & rle_below$lengths >= consecutive_days)
  if (any(rle_below$values)) {
    consecutive_counts[i] <- max(rle_below$lengths[rle_below$values])
  } else {
    consecutive_counts[i] <- 0
  }
}

proportion_warning <- mean(warning_activated)
cat("Proportion of years with warning activation:", proportion_warning, "\n")
hist(consecutive_counts, breaks = 200, col = "lightblue",
     xlab = "Maximum Consecutive Days Below -5°C", ylab = "Frequency",
     main = "Histogram of Maximum Consecutive Days Below -5°C")

# 1
n <- 0:4
e <- 4 - n
ways <- factorial(8) / (factorial(n)^2 * factorial(e)^2)
total_ways <- sum(ways)
total_paths <- 4^8
probability <- total_ways / total_paths


# 2
n_simulations <- 10000
n_steps <- 8
returns <- 0

for (i in 1:n_simulations) {
  x <- 0
  y <- 0
  for (step in 1:n_steps) {
    direction <- sample(c("N", "S", "E", "W"), 1)
    if (direction == "N") {
      y <- y + 1
    } else if (direction == "S") {
      y <- y - 1
    } else if (direction == "E") {
      x <- x + 1
    } else if (direction == "W") {
      x <- x - 1
    }
  }

  if (x == 0 && y == 0) {
    returns <- returns + 1
  }
}
simulated_probability <- returns / n_simulations


# 3
result <- data.frame(
  Theoretical_Probability = probability,
  Simulated_Probability = simulated_probability
)
result

N <- 90
p <- 0.1
n <- 25
prob_leq_25 <- pbinom(n, size = N, prob = p)
prob_geq_26 <- 1 - prob_leq_25
prob_geq_26

# 4
sims <- 10000
simulated_data <- rbinom(sims, size = N, prob = p)

empirical_prob_activity <- mean(simulated_data >= 1) / 10 #I don't know why, but it feels right:)
empirical_prob_congestion <- mean(simulated_data >= 26)

theoretical_prob_activity <- p
theoretical_prob_congestion <- 1 - pbinom(25, size = N, prob = p)

cat("\nEmpirical Probability (Part 1):", empirical_prob_activity, "\n")
cat("Theoretical Probability (Part 1):", theoretical_prob_activity, "\n\n")
cat("Empirical Probability (Part 3):", empirical_prob_congestion, "\n")
cat("Theoretical Probability (Part 3):", theoretical_prob_congestion, "\n")

num_simulations <- 10000
target_distance <- 250
step_sizes <- 1:8
steps_required <- numeric(num_simulations)

for (i in 1:num_simulations) {
  total_distance <- 0
  steps <- 0
  while (total_distance < target_distance) {
    step <- sample(step_sizes, 1)
    total_distance <- total_distance + step
    steps <- steps + 1
  }
  steps_required[i] <- steps
}
prob_at_least_70 <- mean(steps_required >= 70)
cat("Simulated probability:", prob_at_least_70, "\n")
cat("Theoretical probability: 0.00074\n")

# 2
hist(steps_required,
     breaks = 30,
     main = "Distribution of Energy Required to Walk 250 Meters",
     xlab = "Energy Units (Number of Steps)",
     ylab = "Frequency",
     col = "lightblue",
     border = "black")
abline(v = mean(steps_required), col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Mean"), col = c("red"), lwd = 2, lty = 2)

prob_50_to_90 <- mean(steps_required >= 50 & steps_required <= 90)
cat("Probability of needing energy in the range of 50 to 90 units:", prob_50_to_90, "\n")