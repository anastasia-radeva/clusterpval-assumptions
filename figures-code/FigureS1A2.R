library(ggplot2)
library(patchwork)
library(dplyr)
ev_cat <- NULL

# for evaluation plotting
sim <- load_simulation("A2-type1-est-n100-q2-pt1")
draws <- draws(sim)
ev_cat <- rbind(ev_cat, ev)

ev_cat$len <- rep(NA, nrow(ev_cat))
ev_cat$len[grep("len_2", ev_cat$Model)] <- 2
ev_cat$len[grep("len_4", ev_cat$Model)] <- 4
ev_cat$len[grep("len_6", ev_cat$Model)] <- 6
ev_cat$len <- as.factor(ev_cat$len)

# Draws where we can't reject H0
average_h0 <- ev_cat[ev_cat$Method == "average-iso-est-test-K-3" & ev_cat$effect != 0, ]
average_h0_filter <- as.data.frame(average_h0 %>% group_by(len))

# Draws where we can reject H0
average_h1 <- ev_cat[ev_cat$Method == "average-iso-est-test-K-3" & ev_cat$effect == 0, ]
average_h1_filter <- as.data.frame(average_h1 %>% group_by(len))

#### Data Plot ####

# WHY ONLY LEN 2 MODELS??
# First draw where effect is 0 (H0)
exemplary_h0_draw <- average_h0$Draw[1]
exemplary_h0_model <- average_h0$Model[1]
# Last draw where effect is not 0 (H1)
exemplary_h1_draw <- tail(average_h1$Draw, n=1)
exemplary_h1_model <- tail(average_h1$Model, n=1)
# Get this draw
index_0 <- which(sapply(draws, function(x) x@name == exemplary_h0_model))
draws_0 <- draws[[index_0]]@draws[[exemplary_h0_draw]]
index_1 <- which(sapply(draws, function(x) x@name == exemplary_h1_model))
draws_1 <- draws[[index_1]]@draws[[exemplary_h1_draw]]
# Extract data
h0_data <- draws_0$data
h1_data <- draws_1$data
# Extract cluster labels
h0_clusters <- draws_0$clusters
h1_clusters <- draws_1$clusters

# Get the corresponding len value
h0_len <- as.numeric(as.character(average_h0_filter$len)[1])
h1_len <- as.numeric(as.character(average_h1_filter$len)[1])

# Data into df
h0_df <- as.data.frame(h0_data)
h0_df$Cluster <- as.factor(h0_clusters)
colnames(h0_df)[1:2] <- c("Feature1", "Feature2")
h1_df <- as.data.frame(h1_data)
h1_df$Cluster <- as.factor(h1_clusters)
colnames(h1_df)[1:2] <- c("Feature1", "Feature2")

library(ggplot2)
library(patchwork)
library(dplyr)

# Determine global x and y limits for consistent scaling
x_limits <- range(c(h0_df$Feature1, h1_df$Feature1))
y_limits <- range(c(h0_df$Feature2, h1_df$Feature2))

# Scatter Plot for H0 (With Legend)
p1 <- ggplot(h0_df, aes(x = Feature1, y = Feature2, color = Cluster)) +
  geom_point(size = 3) +
  labs(
    title = bquote("Exemplary Draw (H0 holds, " ~ delta ~ " = " ~ .(h0_len) ~ ")"),
    x = "Feature 1",
    y = "Feature 2",
    color = "Cluster"
  ) +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 10)
  )

# Scatter Plot for H1 (With Legend)
p2 <- ggplot(h1_df, aes(x = Feature1, y = Feature2, color = Cluster)) +
  geom_point(size = 3) +
  labs(
    title = bquote("Exemplary Draw (H0 rejected, " ~ delta ~ " = " ~ .(h1_len) ~ ")"),
    x = "Feature 1",
    y = "Feature 2",
    color = "Cluster"
  ) +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 10)
  )

# QQ Plot Limits (Ensuring Consistency)
qq_x_limits <- c(0, 1)
qq_y_limits <- c(0, 1)

# QQ Plot for H0 (Distance Between Clusters Legend)
p3 <- ggplot(average_h0_filter) + 
  geom_qq(aes(sample = pval, group = len, colour = len), size = 0.8, distribution = qunif) + 
  geom_abline(slope = 1, intercept = 0, col = "black") + 
  xlab("Uniform(0, 1) Quantiles") + 
  ylab("Empirical Quantiles") +
  scale_x_continuous(limits = qq_x_limits) +
  scale_y_continuous(limits = qq_y_limits) +
  labs(title = "QQ-plot of p-values - Average Linkage u.H0", colour = expression("Distance between clusters ("~delta~")")) +
  scale_colour_manual(values = c("#9ecae1", "#4292c6", "#08519c")) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom")

# QQ Plot for H1 (Distance Between Clusters Legend)
p4 <- ggplot(average_h1_filter) + 
  geom_qq(aes(sample = pval, group = len, colour = len), size = 0.8, distribution = qunif) + 
  geom_abline(slope = 1, intercept = 0, col = "black") + 
  xlab("Uniform(0, 1) Quantiles") + 
  ylab("Empirical Quantiles") +
  scale_x_continuous(limits = qq_x_limits) +
  scale_y_continuous(limits = qq_y_limits) +
  labs(title = "QQ-plot of p-values - Average Linkage u.H1") + 
  scale_colour_manual(values = c("#9ecae1", "#4292c6", "#08519c")) +
  guides(colour = "none") +  # ✅ This removes the legend
  theme_bw(base_size = 15)

# Combine Plots
top_row <- (p1 + p2) + 
  plot_layout(guides = "collect", widths = c(1, 1)) & 
  theme(legend.position = "bottom")

bottom_row <- (p3 + p4) + 
  plot_layout(guides = "collect", widths = c(1, 1)) & 
  theme(legend.position = "bottom")

# Final Combined Plot
combined <- (top_row) / (bottom_row) + 
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Different Variance Across Clusters",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10))
    )
  )

# Save the Combined Plot
ggsave("~/clusterpval-assumptions/figures/FigureS1A2.pdf", 
       plot = combined, 
       height = 8, width = 12.5)