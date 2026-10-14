#################################
## Empirical Orthogonal Functions
#################################
## Requires: 00_init.R, 01_load_data.R (S_x_df, loc_x, time_step_df,
##           gridpolygon_sf, nT, nS, center_data(), plot_factors(),
##           plot_loadings())
## Reproduces manuscript Figures 2 (right panel), 4, 5, 6 and Tables 1-2,
## and Supplementary Material Figures S1.4, S1.5 (six-dimension real-data
## SVD, AMMI + spatial centering only) and S1.6 (AMMI residuals diagnostic).
## Also writes res/*.RData, consumed by 03_simulation_metrics_figures.R.

dir.create(paste0("images/", spp_to_plot, "/additionnal_dim"),
           showWarnings = FALSE, recursive = TRUE)
dir.create("res", showWarnings = FALSE)

## Data matrix Zs: n locations (rows, ordered by cell) x T time steps
## (columns, in chronological order, consistent with time_step_df)
Zs <- S_x_df %>%
  filter(!is.na(S_x), !is.na(x), !is.na(y), !is.na(cell)) %>%
  dplyr::select(cell, date, S_x) %>%
  arrange(cell, date) %>%
  pivot_wider(names_from = date, values_from = S_x) %>%
  dplyr::select(-cell) %>%
  as.matrix()

# Mean spatial pattern and mean temporal trend (Figure 2, right panel) ------

## Mean spatial pattern: average over time at each location
spat_average <- apply(Zs, 1, mean, na.rm = TRUE)

spatial_pattern_df <- data.frame(loc_x, S_x = spat_average - mean(Zs)) %>%
  inner_join(gridpolygon_sf, by = c("cell" = "layer")) %>%
  st_as_sf()

spatial_pattern_plot <- ggplot() +
  geom_sf(data = spatial_pattern_df, aes(col = S_x, fill = S_x)) +
  scale_color_viridis_c(
    option = "viridis",
    labels = scales::label_number(accuracy = 0.01),
    name = NULL
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    labels = scales::label_number(accuracy = 0.01),
    name = NULL
  ) +
  geom_sf(data = mapBase) +
  coord_sf(xlim = c(grid_xmin, grid_xmax),
           ylim = c(grid_ymin, grid_ymax),
           expand = FALSE) +
  theme_bw() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, face = "bold")) +
  xlab("") + ylab("") +
  ggtitle(label = "Mean spatial pattern",
          subtitle = expression("(" * mu + bold(alpha) * " )"))
plot(spatial_pattern_plot)

## Mean temporal trend: average over locations at each time step
temp_average_plot <- apply(Zs, 2, mean)
temporal_trend_df <- data.frame(time_step_df, S_x = temp_average_plot)
temporal_trend_plot <- ggplot() +
  geom_line(data = temporal_trend_df, aes(x = date, y = S_x)) +
  theme_classic() +
  theme(aspect.ratio = 1/4,
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, face = "bold")) +
  xlab("") + ylab("") + ylim(0, NA) +
  ggtitle(label = "Mean temporal trend",
          subtitle = expression("(" * mu + bold(beta) * " )"))
plot(temporal_trend_plot)

spatial_temporal_trend_plot <- cowplot::plot_grid(spatial_pattern_plot, temporal_trend_plot,
                                                  ncol = 1, rel_heights = c(1, 0.75))

ggsave(filename = paste0("images/", spp_to_plot, "/spatial_temporal_trend_plot.png"),
       plot = spatial_temporal_trend_plot,
       width = 6, height = 6, bg = "white")

# SVD of the matrix under each centering model ----------------------------
n_dim <- 3                     # number of dimensions to plot
n_rank <- min(dim(Zs))
switch_sign <- TRUE            # switch the sign of the EOFs (arbitrary in an SVD)
which_centering_to_switch <- c(1, 3)
uv_plot_list <- list()
type_center_df <- data.frame(center_type = c(1, 2, 3),
                             center_type_name = c("Spatial centering", "Temporal centering", "AMMI model"))
type_center_df$center_type_name <- factor(type_center_df$center_type_name,
                                          levels = c("AMMI model", "Spatial centering", "Temporal centering"))

# Single colour for all EOF time series
dim_colors <- c("Dim.1" = "steelblue", "Dim.2" = "steelblue", "Dim.3" = "steelblue",
                "Dim.4" = "steelblue", "Dim.5" = "steelblue", "Dim.6" = "steelblue")

# Per centering model (1 = spatial, 2 = temporal, 3 = AMMI):
n_dim_plot   <- c(3, 1, 3)  # number of EOFs shown in the main figures
n_dim_filter <- c(6, 1, 6)  # number of EOFs retained (the rest = residuals)

for (center_i in 1:3) {
  
  ## Center the data
  #-----------------
  res_center_data <- center_data(Zs, center_type = center_i)
  Zs_c <- res_center_data$data_matrix_c
  
  ## Perform SVD
  #-------------
  res_svd_i <- svd(Zs_c)
  
  ## Switch sign if necessary
  if (switch_sign & center_i %in% which_centering_to_switch) {
    res_svd_i$u <- -res_svd_i$u
    res_svd_i$v <- -res_svd_i$v
  }
  
  ## Plot factors (spatial EOFs) and loadings (temporal EOFs)
  #----------------------------------------------------------
  u_plot <- plot_factors(support = loc_x, res_svd = res_svd_i, n_dim = n_dim_plot[center_i]) +
    ggtitle(label = type_center_df$center_type_name[center_i], subtitle = "") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkgray"))
  v_plot <- plot_loadings(time_step_df, res_svd_i, n_dim = n_dim_plot[center_i]) +
    scale_color_manual(values = dim_colors)
  
  ## Temporal centering has a single EOF: add two empty facets so that the
  ## panels align with the two other centering models
  if (center_i == 2) {
    v_df <- data.frame(time_step_df, Dim.1 = res_svd_i$v[, 1]) %>%
      pivot_longer(cols = Dim.1, names_to = "dim", values_to = "v") %>%
      bind_rows(tibble(dim = c("Dim.2", "Dim.3")))
    v_plot <- ggplot(v_df) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2) +
      geom_line(linewidth = 1, aes(x = date, y = v, col = dim, group = dim), na.rm = TRUE) +
      scale_color_manual(values = dim_colors) +
      theme_bw() +
      xlab("") + ylab("") +
      theme(legend.title = element_blank(),
            aspect.ratio = 1/4,
            legend.position = "none",
            plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "pt")) +
      coord_cartesian(clip = "off") +
      facet_wrap(. ~ dim, ncol = 1) +
      geom_vline(xintercept = v_df$date[which(str_detect(v_df$date, "01-01"))],
                 linetype = "dashed", color = "skyblue")
  }
  
  v_plot <- v_plot + ylim(-0.4, 0.4)
  uv_plot <- cowplot::plot_grid(u_plot, v_plot, ncol = 1, rel_heights = c(0.35, 0.65))
  ggsave(paste0("images/", spp_to_plot, "/eofs_c", center_i, ".png"),
         plot = uv_plot, width = 9 * 3, height = 9, bg = "white")
  uv_plot_list[[center_i]] <- uv_plot
  
  ## Six-dimension figures: SM Figures S1.4 (AMMI) and S1.5 (spatial centering).
  ## No such figure for temporal centering (single EOF).
  if (center_i %in% c(1, 3)) {
    plot_fact <- plot_factors(support = loc_x, res_svd = res_svd_i, n_dim = 6)
    plot_load <- plot_loadings(time_step_df, res_svd_i, n_dim = 6) +
      scale_color_manual(values = dim_colors)
    plot_fact_load <- cowplot::plot_grid(plot_fact, plot_load, nrow = 1)
    ggsave(paste0("images/", spp_to_plot, "/additionnal_dim/eofs_c", center_i, "_additional_dim.png"),
           plot = plot_fact_load, width = 9 * 1.25, height = 6, bg = "white")
  }
  
  ## Sums of squares of each term (Tables 1-2, Figure 6)
  #-----------------------------------------------------
  total_ss    <- sum((Zs - mean(Zs))^2)
  spatial_ss  <- ncol(Zs) * sum((res_center_data$row_mean)^2)
  temporal_ss <- nrow(Zs) * sum((res_center_data$col_mean)^2)
  interac_ss  <- (res_svd_i$d)^2
  
  ss_df <- data.frame(ss = c(spatial_ss, temporal_ss, interac_ss),
                      dim = c("Spatial", "Temporal", paste0("SVD.", seq_along(interac_ss))),
                      center_type = center_i,
                      SVD_dim_in = c(NA, NA,
                                     rep("EOF", n_dim_filter[center_i]),
                                     rep("Residuals", min(dim(Zs_c)) - n_dim_filter[center_i])))
  if (center_i == 1) ss_df_full <- ss_df
  if (center_i != 1) ss_df_full <- rbind(ss_df, ss_df_full)
  
  ## Residuals diagnostic: SM Figure S1.6 (AMMI model only)
  #--------------------------------------------------------
  if (analyse_residuals & center_i == 3) {
    
    ## Residuals = matrix reconstructed from the non-retained SVD dimensions
    ncp <- n_dim_filter[center_i]
    residuals_mat <- res_svd_i$u[, (ncp + 1):n_rank] %*%
      diag(res_svd_i$d[(ncp + 1):n_rank]) %*%
      t(res_svd_i$v[, (ncp + 1):n_rank])
    x_resid <- as.vector(residuals_mat)
    sd_res  <- sd(x_resid)
    
    # Student-t: df estimated from the excess kurtosis (method of moments)
    kurt_res <- mean(((x_resid - mean(x_resid)) / sd_res)^4) - 3
    df_t     <- if (kurt_res > 0) max(4 + 6 / kurt_res, 4.01) else 30
    scale_t  <- sd_res * sqrt((df_t - 2) / df_t)
    
    # Candidate distributions, all with mean 0 and variance sd_res^2
    families <- list(
      normal = list(
        label = "Normal",
        qtheo = function(p) as.numeric(qnorm(p, mean = 0, sd = sd_res)),
        dtheo = function(x) dnorm(x, mean = 0, sd = sd_res)
      ),
      student = list(
        label = "Student-t",
        qtheo = function(p) as.numeric(scale_t * qt(p, df = df_t)),
        dtheo = function(x) dt(x / scale_t, df = df_t) / scale_t
      ),
      laplace = list(
        label = "Laplace",
        qtheo = function(p) as.numeric(qlaplace(p, scale = sd_res / sqrt(2))),
        dtheo = function(x) dlaplace(x, scale = sd_res / sqrt(2))
      )
    )
    
    ## One row per distribution: label | QQ-plot | histogram
    n_fam <- length(families)
    clip_lim <- c(-3, 3)   # displayed range for both plots
    
    png(paste0("images/", spp_to_plot, "/residuals_qq_hist_c", center_i, ".png"),
        width = 9 / 1.5, height = 4.5 / 1.5 * n_fam, units = "in", res = 300)
    
    layout(matrix(1:(n_fam * 3), nrow = n_fam, byrow = TRUE), widths = c(0.15, 0.425, 0.425))
    
    for (fam in families) {
      
      ## Row label panel (rotated 90 degrees)
      par(mar = c(0, 0, 0, 0), pty = "m")
      plot.new()
      text(0.5, 0.5, fam$label, srt = 90, cex = 1.3, font = 2)
      
      ## QQ-plot (square panel, clipped to [-3, 3] on both axes)
      par(mar = c(4, 4, 2, 1), pty = "s")
      qtheo_vals <- fam$qtheo(ppoints(length(x_resid)))
      qqplot(x = qtheo_vals,
             y = x_resid,
             main = "",
             xlab = "Theoretical quantiles",
             ylab = "Sample quantiles",
             xlim = clip_lim,
             ylim = clip_lim)
      qqline(x_resid, distribution = fam$qtheo)
      
      ## Histogram + empirical density (blue) + fitted density (red)
      x_resid_clip <- x_resid[x_resid >= clip_lim[1] & x_resid <= clip_lim[2]]
      grid_x <- seq(clip_lim[1], clip_lim[2], length.out = 1000)
      hist_breaks <- seq(clip_lim[1], clip_lim[2], length.out = 30)
      dens_max <- max(
        hist(x_resid_clip, plot = FALSE, breaks = hist_breaks)$density,
        density(x_resid_clip)$y,
        fam$dtheo(grid_x)
      )
      
      par(mar = c(4, 4, 2, 1), pty = "s")
      hist(x_resid_clip, freq = FALSE, ylim = c(0, dens_max * 1.05), xlim = clip_lim,
           breaks = hist_breaks,
           main = "",
           xlab = "Residuals")
      lines(density(x_resid_clip), col = "blue", lwd = 3)
      curve(fam$dtheo(x), add = TRUE, col = "red", lwd = 2, from = clip_lim[1], to = clip_lim[2])
    }
    
    dev.off()
  }
}

# Proportion of variance explained (Figure 6, Tables 1-2) ------------------

ss_df_full$dim <- factor(ss_df_full$dim, levels = unique(ss_df_full$dim))
ss_df_full$center_type_name <- type_center_df$center_type_name[ss_df_full$center_type]

## Group SVD dimensions into "Interaction Dim.k" (retained) and "Residuals"
prop_var_plot <- ss_df_full |>
  mutate(dim = as.character(dim),
         dim_num = as.numeric(str_extract(dim, "(?<=^SVD\\.)\\d+")),
         dim = case_when(
           !str_detect(dim, "^SVD\\.") ~ dim,
           SVD_dim_in == "EOF"         ~ paste0("Interaction Dim.", dim_num),
           TRUE                        ~ SVD_dim_in
         )) |>
  group_by(dim, center_type_name) |>
  dplyr::summarise(ss = sum(ss), .groups = "drop")

prop_var_plot$dim <- factor(prop_var_plot$dim,
                            levels = c("Spatial",
                                       "Temporal",
                                       paste0("Interaction Dim.", 1:max(n_dim_filter)),
                                       "Residuals"))

skyblue_palette <- rev(c(
  "#E3F6FF", "#CFEFFF", "#9BD7F0", "#87CEEB", "#5FA8D3", "#2E8BC0"
))

pattern_lookup <- c(
  "Spatial"           = "stripe",
  "Temporal"          = "crosshatch",
  "Interaction Dim.1" = "circle",
  "Interaction Dim.2" = "circle",
  "Interaction Dim.3" = "circle",
  "Interaction Dim.4" = "circle",
  "Interaction Dim.5" = "circle",
  "Interaction Dim.6" = "circle",
  "Residuals"         = "none"
)

lm_var_expl_plot <- ggplot(data = prop_var_plot,
                           aes(x = factor(center_type_name), y = ss / total_ss,
                               fill = dim, pattern = dim)) +
  geom_col_pattern(position = "stack",
                   color = "black",
                   pattern_fill = "black",
                   pattern_colour = "black",
                   pattern_density = 0.15,
                   pattern_spacing = 0.03,
                   pattern_angle = 45,
                   pattern_size = 0.3) +
  scale_pattern_manual(values = pattern_lookup) +
  scale_fill_manual(values = c("orange", "mediumseagreen", skyblue_palette, "lightgrey")) +
  labs(x = "Center Type", y = "Proportion of variance", fill = "Dimension", pattern = "Dimension") +
  theme_bw() +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(size = 10.5)) +
  ggtitle("Proportion of variance explained \n by each term of the linear models")

## Table of variance explained by each term of the AMMI model
## (converted to LaTeX for the manuscript)
table_perc_var <- prop_var_plot |>
  filter(center_type_name == "AMMI model") |>
  dplyr::select(-center_type_name) |>
  mutate(perc_var = ss / sum(ss))
print(table_perc_var)

## Variance explained by each SVD dimension (first 20)
vline_df <- data.frame(dim = n_dim_filter + 1,
                       center_type_name = type_center_df$center_type_name)
svd_var_expl_plot <- ss_df_full |>
  filter(str_detect(dim, "^SVD\\.")) |>
  dplyr::mutate(dim = as.numeric(str_remove(dim, "^SVD\\."))) |>
  filter(dim <= 20) |>
  ggplot() +
  geom_line(aes(x = dim, y = ss / total_ss)) +
  geom_point(aes(x = dim, y = ss / total_ss), size = 0.75) +
  facet_wrap(. ~ center_type_name, scales = "free_y") +
  theme_bw() +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_blank()) +
  ggtitle("Proportion of variance explained by the SVD dimensions \n conducted on the matrix R (interaction + residuals)") +
  ylab("Proportion of variance") +
  xlab("Dimensions") +
  geom_vline(data = vline_df, aes(xintercept = dim), color = "skyblue")

lm_var_expl_plot2 <- cowplot::plot_grid(lm_var_expl_plot, NULL, nrow = 1, rel_widths = c(0.7, 0.3))
prop_var_expl_plot <- cowplot::plot_grid(lm_var_expl_plot2, svd_var_expl_plot, ncol = 1)

ggsave(paste0("images/", spp_to_plot, "/prop_var_expl.png"),
       plot = prop_var_expl_plot, width = 11, height = 11, bg = "white")

# EOFs of the three centering models side by side (Figures 4-5) -----------
gray_line <- ggdraw() +   # vertical separator
  draw_line(
    x = c(0.5, 0.5),
    y = c(0, 1),
    color = "darkgray", linetype = "dashed",
    linewidth = 0.5
  )
uv_plot_full <- cowplot::plot_grid(uv_plot_list[[3]], gray_line,
                                   uv_plot_list[[1]], gray_line,
                                   uv_plot_list[[2]], nrow = 1,
                                   rel_widths = c(1, 0.05, 1, 0.05, 1))
ggsave(paste0("images/", spp_to_plot, "/eofs_c_full.png"),
       plot = uv_plot_full, width = 6 * 3, height = 9, bg = "white")

# ---------------------------------------------------------------------------
# Save intermediate objects used by 03_simulation_metrics_figures.R.
# The simulation draws its spatial and temporal supports, the mean spatial
# pattern, the mean temporal trend and the two leading spatial EOFs of the
# AMMI model (used as the "movement" pattern) from this real-data fit, so
# that it is calibrated on the actual sole dataset.
# ---------------------------------------------------------------------------
save(loc_x, file = "res/loc_x.RData")
save(time_step_df, file = "res/time_step_df.RData")
save(spatial_pattern_df, file = "res/spatial_pattern_df.RData")
save(temporal_trend_df, file = "res/temporal_trend_df.RData")

movement_df <- as.data.frame(res_svd_i$u[, 1:2])   # res_svd_i = last loop iteration (AMMI)
save(movement_df, file = "res/movement_df.RData")