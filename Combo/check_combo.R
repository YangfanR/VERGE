library(R.matlab)
library(ggplot2)
library(dplyr)
library(gridExtra)

results_path <- "results_combo.mat"
results <- readMat(results_path)

### check selected predictors
sel_pred = results$PPI.pred > 0.5

### check selected covariates
sel_cov = results$PPI.cov > 0.5

### check graph selection
P <- 69
indmx <- matrix(1:(P^2), nrow = P, ncol = P)
upperind <- indmx[upper.tri(indmx)]
G_est <- results$adj.avg > 0.5
G_est_up <- G_est[upperind]
sum(G_est_up)

####### Plots for covariate effects
Z = read.csv("Z.csv")
Z = Z[,c(2, 28, 30)]
Beta <- results$Beta.s
pred_inds <- which(results$PPI.pred > 0.5)
inds <- which(results$PPI.cov > 0.5, arr.ind = TRUE)
inds <- inds[order(inds[,1]),]
Beta <- Beta[,pred_inds]

Beta_df <- as.data.frame(Beta)
Z_df <- as.data.frame(Z)

gg_color_hue <- function(n) {
     hues = seq(15, 375, length = n + 1)
     hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)

point_color <- "#F8766D"
point_size <- 0.5

subtitles <- c("UCG-004","Megasphaera", "Megasphaera","Megasphaera",
               "Unclassified", "Unclassified",
               "UCG-002","UCG-002", "Bacteroides", "Bacteroides", 
               "AD3011 group","AD3011 group",
               "Catenibacterium", "Catenibacterium", "Anaerostipes")
z_labels <- c("Sex", "Total fat", "AOAC fiber") 

plots <- list()
for(i in 1:nrow(inds)) {
  # Prepare data for the current plot
  if(inds[i,2] == 1){
    plot_data <- data.frame(Z = as.factor(Z_df[, inds[i, 2]]), 
                            Beta = Beta_df[, inds[i, 1]])
    p <- ggplot(plot_data, aes(x = Z, y = Beta)) +
      geom_point(color = point_color, size = point_size) +
      scale_x_discrete(expand = c(0.03, 0.03), labels = c("M","F"))  # Remove space for the first plot
    y_min <- min(plot_data$Beta) - 0.1
    y_max <- max(plot_data$Beta) + 0.1
  } else {
    plot_data <- data.frame(Z = Z_df[, inds[i, 2]], 
                            Beta = Beta_df[, inds[i, 1]])
    p <- ggplot(plot_data, aes(x = Z, y = Beta)) +
      geom_point(color = point_color, size = point_size)
    y_min <- min(plot_data$Beta) - 0.2
    y_max <- max(plot_data$Beta) + 0.2
  }
  # Set common y-axis label, custom x-axis labels, and subtitles
  p <- p + labs(x = z_labels[inds[i, 2]], y = expression(beta), subtitle = subtitles[i]) +
    theme_minimal() +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    ylim(y_min, y_max)  # Set y-axis limits
  
  plots[[i]] <- p
}

# Arrange the plots in a 3x2 grid using grid.arrange
grid_plot <- do.call(grid.arrange, c(plots, ncol = 3))

# Print the grid plot
print(grid_plot)

ggsave("Combo_cov.pdf", plot = grid_plot, width = 9, height = 12, device = "pdf")
