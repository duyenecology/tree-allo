# Load necessary libraries
library(ggplot2)
library(dplyr)
library(stringr)

# Load MCMC results and Stan data
gmm_h <- tar_read(fit_nlr3_mcmc_gmm_h)
stan_nlr_data_h <- tar_read(stan_data_nlr_h)

# Extract relevant data
log_y <- stan_nlr_data_h$log_y
x <- stan_nlr_data_h$x
sp <- as.factor(stan_nlr_data_h$jj)

# Extract beta parameters
beta <- gmm_h$summary(variable = "beta")

# Ensure dplyr and stringr are loaded and then extract the mean values of beta parameters
log_a2 <- beta %>% filter(str_detect(variable, "beta\\[\\d+,1\\]")) %>% pull(mean)
b2 <- beta %>% filter(str_detect(variable, "beta\\[\\d+,2\\]")) %>% pull(mean)
k1 <- beta %>% filter(str_detect(variable, "beta\\[\\d+,3\\]")) %>% pull(mean)

# Verify the lengths of the extracted parameters
print(length(log_a2))
print(length(b2))
print(length(k1))
print(length(sp))

# Compute fitted values for H using the gMM formula
fitted_log_y_h <- sapply(1:length(log_a2), function(i) {
  log_a2[i] + b2[i] * log(x) - log(k1[i] + x^b2[i])
})

# Calculate community-level fitted line by averaging the species-level coefficients
community_log_a2 <- mean(log_a2)
community_b2 <- mean(b2)
community_k1 <- mean(k1)
community_fitted_log_y_h <- community_log_a2 + community_b2 * log(x) - log(community_k1 + x^community_b2)
community_fitted_H <- exp(community_fitted_log_y_h)

# Create a data frame for plotting
data_h <- data.frame(H = exp(log_y), DBH = x, sp = sp)

# Add species-specific fitted values to the data frame
for (i in 1:length(log_a2)) {
  data_h[[paste0("fitted_H_sp_", i)]] <- exp(fitted_log_y_h[, i])
}

# Plot observed H vs DBH with fitted lines for each species and community-level line
gmm_h_plot <- ggplot(data_h, aes(x = DBH, y = H)) +
  geom_point(alpha = 0.5, color = "grey") +
  geom_line(aes(y = community_fitted_H), color = "white", linewidth = 1.2) +
  geom_line(aes(y = rowMeans(fitted_log_y_h)), color = "blue", alpha = 0.3) +
  labs(title = "Scatter Plot of H vs DBH with Fitted Lines for Each Species and Community-Level Line",
       x = "DBH (cm)",
       y = "H (m)") +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "black"))

# Save the plot
ggsave("figs/beta_plot_H_vs_DBH.png", plot = gmm_h_plot, width = 8, height = 6, units = "in", dpi = 300, bg = "white")

gmm_h_plot


# Plot observed H vs DBH with fitted lines for each species and community-level line
gmm_h_plot1 <- ggplot(data_h, aes(x = DBH, y = H, group = sp)) +
  geom_point(alpha = 0.5, color = "grey") +
  # geom_line(aes(y = community_fitted_H), color = "#4682B4", linewidth = 0.3) +
  geom_line(aes(y = rowMeans(fitted_log_y_h), group = sp), color = "#4682B4", alpha = 0.3) + # SteelBlue color
  labs(title = "Scatter Plot of H vs DBH with Fitted Lines for Each Species and Community-Level Line",
       x = "DBH (cm)",
       y = "H (m)") +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"))

# Save the plot
ggsave("figs/beta1_plot_H_vs_DBH.png", plot = gmm_h_plot1, width = 8, height = 6, units = "in", dpi = 300, bg = "white")

# Display the plot
print(gmm_h_plot1)



# Melt data_h to long format for easier ggplot handling
library(reshape2)
data_long <- melt(data_h, id.vars = c("H", "DBH", "sp"), measure.vars = grep("fitted_H_sp_", names(data_h), value = TRUE))

# Plot observed H vs DBH with fitted lines for each species and community-level line
gmm_h_plot2 <- ggplot(data_h, aes(x = DBH, y = H)) +
  geom_point(alpha = 0.5, color = "grey") +
  geom_line(data = data_long, aes(x = DBH, y = value, group = variable), color = "#4682B4", alpha = 0.3) + # SteelBlue color for species fits
  geom_line(aes(y = community_fitted_H), color = "white", linewidth = 1.2) +
  labs(title = "Scatter Plot of H vs DBH with Fitted Lines for Each Species and Community-Level Line",
       x = "DBH (cm)",
       y = "H (m)") +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "black"))

# Save the plot
ggsave("figs/beta2_plot_H_vs_DBH.png", plot = gmm_h_plot2, width = 8, height = 6, units = "in", dpi = 300, bg = "white")

# Display the plot
print(gmm_h_plot2)

# Plot observed H vs DBH with fitted lines for each species and community-level line
gmm_h_plot3 <- ggplot(data_h, aes(x = DBH, y = H)) +
  geom_point(alpha = 0.5, color = "grey") +
  geom_line(data = data_long, aes(x = DBH, y = value, group = variable), color = "#4682B4", alpha = 0.3) + # SteelBlue color for species fits
  geom_line(aes(y = community_fitted_H), color = "white", linewidth = 1.2) +
  labs(title = "Scatter Plot of H vs DBH with Fitted Lines for Each Species and Community-Level Line",
       x = "DBH (cm)",
       y = "H (m)") +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal()

# Save the plot
ggsave("figs/beta3_plot_H_vs_DBH.png", plot = gmm_h_plot3, width = 8, height = 6, units = "in", dpi = 300, bg = "white")

# Display the plot
print(gmm_h_plot3)
