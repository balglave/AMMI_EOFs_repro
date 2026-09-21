#####################################################################
## Simulation study: recovery of the true spatial/temporal/seasonal
## components under the three centering models
#####################################################################
## Requires: 02_eof_real_data.R to have been run first (writes
##           res/loc_x.RData, res/time_step_df.RData,
##           res/spatial_pattern_df.RData, res/movement_df.RData,
##           res/temporal_trend_df.RData, used below as the calibration
##           basis for the synthetic spatial/temporal/seasonal signals).
## Reproduces:
##   - manuscript Figure 3               (p_dim_full)
##   - SM Figure S1.1 (monthly maps of one simulated replicate)
##   - SM Figure S1.2 (main-effect recovery violin, p_mean)
##   - SM Figure S1.3 (eofs_c_full-style panel, one simulated replicate)
## A colorblind-accessible hatch/pattern layer (ggpattern) is applied on
## top of the color palette throughout, and every figure that encodes
## "Spatial centering" / "Temporal centering" / "AMMI model" shares one
## color+pattern lookup (center_type_colors / center_type_patterns) so the
## encoding is consistent across figures.
#####################################################################

library(dplyr); library(ggplot2); library(tidyr); library(Matrix);
library(cowplot); library(stringr);  library(ggh4x);  library(cowplot);
library(ggpattern);

set.seed(42)

# ---- fixed simulation settings (as in simulate_statistical.R) --------------
sigma2_spatial  <- 2
sigma2_seasonal <- 1
sigma2_temporal <- 2
sigma2_noise    <- 1
grf_rho <- 1; grf_m <- 20
movement_scale <- 0.15
n_reps <- 200

# ---- load data --------------------------------------------------------------
load("res/loc_x.RData")
load("res/time_step_df.RData")
load("res/spatial_pattern_df.RData")
load("res/movement_df.RData")
load("res/temporal_trend_df.RData")

n <- nrow(loc_x)
T <- nrow(time_step_df)
dates_t <- as.Date(time_step_df$date)

# ---- static ingredients (deterministic, computed once) ----------------------
w_x <- pmax(spatial_pattern_df$S_x, 0); w_x <- w_x / sum(w_x)
mu_raw <- w_x - mean(w_x)
mu <- mu_raw / sd(mu_raw) * sqrt(sigma2_spatial)                      # omega_i

anom_scaled <- movement_df$V1 * movement_scale *
  (sd(spatial_pattern_df$S_x) / sd(movement_df$V1))

## xi_ij = eta_i * phi_j : rank-1 spatio-seasonal interaction term ------------
eta_raw <- -anom_scaled
eta <- eta_raw / sd(eta_raw)

seasonal_alpha <- function(month) {
  dplyr::case_when(
    month %in% c(1,2,3,12) ~ 0,  month == 4 ~ 1/3,  month == 5 ~ 2/3,
    month %in% c(6,7,8,9)  ~ 1,  month == 10 ~ 2/3, month == 11 ~ 1/3)
}

temporal_trend_df <- temporal_trend_df[order(temporal_trend_df$t), ]
R_raw <- temporal_trend_df$S_x / mean(temporal_trend_df$S_x) - 1
R_vec <- R_raw / sd(R_raw)
trend <- R_vec * sqrt(sigma2_temporal)                                 # psi_j

month_t <- as.integer(format(dates_t, "%m"))
alpha_c <- seasonal_alpha(month_t)   # length T, deterministic

phi_raw <- alpha_c - mean(alpha_c)
phi <- phi_raw / sd(phi_raw)                                           # unit-variance seasonal scalar

xi_mat <- sqrt(sigma2_seasonal) * outer(eta, phi)                       # n x T, rank-1 interaction

# ---- GRF spatial covariance (noise term) -------------------------------------
x_coarse <- seq(min(loc_x$x), max(loc_x$x), length.out = grf_m)
y_coarse <- seq(min(loc_x$y), max(loc_x$y), length.out = grf_m)
nodes <- as.matrix(expand.grid(x = x_coarse, y = y_coarse))
m2 <- nrow(nodes)
D_coarse <- as.matrix(dist(nodes))
Sigma_c <- exp(-D_coarse / grf_rho)
L_chol <- t(chol(Sigma_c + diag(1e-8, m2)))
dx <- x_coarse[2] - x_coarse[1]; dy <- y_coarse[2] - y_coarse[1]

build_bilinear_A <- function(fx, fy, xc, yc, dx, dy) {
  nn <- length(fx); mx <- length(xc); my <- length(yc)
  ix <- pmax(1, pmin(mx-1, floor((fx-xc[1])/dx)+1))
  iy <- pmax(1, pmin(my-1, floor((fy-yc[1])/dy)+1))
  tx <- pmax(0, pmin(1, (fx-xc[ix])/dx))
  ty <- pmax(0, pmin(1, (fy-yc[iy])/dy))
  node <- function(i,j) (j-1L)*mx+i
  c00 <- node(ix,iy); c10 <- node(ix+1L,iy); c01 <- node(ix,iy+1L); c11 <- node(ix+1L,iy+1L)
  Matrix::sparseMatrix(i=rep(seq_len(nn),4), j=c(c00,c10,c01,c11),
    x=c((1-tx)*(1-ty), tx*(1-ty), (1-tx)*ty, tx*ty), dims=c(nn, mx*my))
}
A <- build_bilinear_A(loc_x$x, loc_x$y, x_coarse, y_coarse, dx, dy)
AL <- A %*% L_chol
set.seed(1); grf_sd_ref <- sd(as.numeric(AL %*% rnorm(m2)))
noise_scale <- sqrt(sigma2_noise) / grf_sd_ref

jitter_frac <- 0.05   # small perturbation, as a fraction of each term's sd

# ---- one simulated replicate (noise + small jitter on mu/seas/trend) --------
simulate_once <- function(seed) {
  set.seed(seed)
  mu_r    <- mu    + rnorm(n, sd = jitter_frac * sqrt(sigma2_spatial))
  trend_r <- trend + rnorm(T, sd = jitter_frac * sqrt(sigma2_temporal))

  total_B <- matrix(0, nrow = T, ncol = n)
  xi_last <- rep(0, n)
  for (step in seq_len(T)) {

    xi_ij <- xi_mat[, step] + rnorm(n, sd = jitter_frac * sqrt(sigma2_seasonal))
    xi_last <- xi_ij

    noise_raw <- as.numeric(AL %*% rnorm(m2)); noise_raw <- noise_raw - mean(noise_raw)
    noise <- noise_raw * noise_scale

    x <- mu_r + xi_ij + trend_r[step] + noise
    total_B[step, ] <- ifelse(x >= 0, x + log1p(exp(-x)), log1p(exp(x)))
  }
  list(Zs = t(total_B), seas = xi_last, mu = mu_r, trend = trend_r)   # Zs: n x T
}

# ---- centering (spatial / temporal / AMMI double centering) -----------------
center_data <- function(data_matrix, center_type = 3){
  rowcol_mean <- mean(data_matrix)
  col_mean <- rep(0, ncol(data_matrix))   # starts at zero...
  row_mean <- rep(0, nrow(data_matrix))   # ...and stays zero unless updated below
  data_matrix <- data_matrix - rowcol_mean
  if(center_type %in% c(1,3)){
    row_mean <- rowMeans(data_matrix)
    data_matrix <- data_matrix - outer(row_mean, rep(1, ncol(data_matrix)))
  }
  if(center_type %in% c(2,3)){
    col_mean <- colMeans(data_matrix)     # computed on the ALREADY row-centered residual for type=3
    data_matrix <- t(t(data_matrix) - outer(col_mean, rep(1, nrow(data_matrix))))
  }
  list(data_matrix_c = data_matrix, row_mean = row_mean, col_mean = col_mean, rowcol_mean = rowcol_mean)
}

# ---- evaluation metrics (correlation + angle) --------------------------------
angle_deg <- function(r) acos(pmin(1, abs(r))) * 180 / pi

eval_metrics <- function(Zs, seas, mu_true, trend_true) {
  purrr::map_dfr(1:3, function(center_i) {
    res_center <- center_data(Zs, center_type = center_i)
    res_svd <- svd(res_center$data_matrix_c)

    cor_meancol_meantemp <- if (center_i %in% c(2,3)) cor(res_center$col_mean, trend_true) else NA
    cor_meanrow_meanspat <- if (center_i %in% c(1,3)) cor(res_center$row_mean, mu_true)    else NA

    dim_metrics <- purrr::map_dfr(1:3, function(i) data.frame(
      dim = i,
      cor_u_meanspat = cor(res_svd$u[,i], mu_true),
      cor_v_meantemp = cor(res_svd$v[,i], trend_true),
      cor_u_season   = cor(res_svd$d[i] * res_svd$u[,i], seas),
      cor_v_season   = cor(res_svd$d[i] * res_svd$v[,i], alpha_c)
    ))

    dim_metrics$center_type <- center_i
    dim_metrics$cor_meancol_meantemp <- cor_meancol_meantemp
    dim_metrics$cor_meanrow_meanspat <- cor_meanrow_meanspat

    # angle counterpart for every correlation metric above
    dim_metrics |> dplyr::mutate(
      angle_u_meanspat       = angle_deg(cor_u_meanspat),
      angle_v_meantemp       = angle_deg(cor_v_meantemp),
      angle_u_season         = angle_deg(cor_u_season),
      angle_v_season         = angle_deg(cor_v_season),
      angle_meancol_meantemp = angle_deg(cor_meancol_meantemp),
      angle_meanrow_meanspat = angle_deg(cor_meanrow_meanspat)
    )
  })
}

# ---- replicate loop -----------------------------------------------------------
cat(sprintf("Running %d replicates ...\n", n_reps))
pb <- txtProgressBar(min = 0, max = n_reps, style = 3)

results_list <- vector("list", n_reps)
for (r in seq_len(n_reps)) {
  sim <- simulate_once(seed = r)
  res <- eval_metrics(sim$Zs, sim$seas, sim$mu, sim$trend)
  res$rep <- r
  results_list[[r]] <- res
  setTxtProgressBar(pb, r)
}
close(pb)

results_df <- dplyr::bind_rows(results_list)
results_df$center_type_name <- factor(
  results_df$center_type,
  levels = 1:3, labels = c("Spatial centering", "Temporal centering", "AMMI model")
)

saveRDS(results_df, "res/simulation_metrics_results.rds")

# ---- summary ------------------------------------------------------------------
summary_df <- results_df |>
  group_by(center_type_name, dim) |>
  summarise(across(starts_with("cor_") | starts_with("angle_"), ~ median(., na.rm = TRUE)),
            .groups = "drop")
print(summary_df)

# =================================================================================
# ---- SHARED COLOR PALETTE across ALL plots (bar charts + violin plots) --------
# =================================================================================
center_type_colors <- c(
  "Spatial centering"  = "orange",
  "Temporal centering" = "mediumseagreen",
  "AMMI model"          = "#2E8BC0"   # matches Interaction Dim.1 blue in the bar chart
)

# hatch/texture overlay for the violin plots too, same rationale as the bar
# chart's ggpattern layer: an independent visual channel so colorblind
# readers can tell "Spatial centering" / "Temporal centering" / "AMMI model"
# apart without relying on color. Mapped directly to center_type_name so the
# fill and pattern legends merge into one.
center_type_patterns <- c(
  "Spatial centering"  = "stripe",
  "Temporal centering" = "crosshatch",
  "AMMI model"          = "circle"
)

# ---- comparison plot: main-effect recovery (SM Figure S1.2) -------------------

# mean spatial pattern (alpha) vs spatial forcing (omega), and
#     mean temporal trend (beta) vs temporal forcing (psi) -- only
#     shown for the models where each term is actually estimated
# |cos(angle)| = |correlation| (angle = acos(|r|)): 1 = perfect match, 0 = orthogonal
metric_labels_mean <- c(
  meanrow_meanspat = "Spatial main effect\nvs. true spatial mean",
  meancol_meantemp = "Temporal main effect\nvs. true temporal mean"
)

mean_plot_df <- results_df |>
  distinct(center_type_name, rep, .keep_all = TRUE) |>
  select(center_type_name, rep, cor_meanrow_meanspat, cor_meancol_meantemp) |>
  pivot_longer(-c(center_type_name, rep), names_to = "metric", values_to = "value") |>
  mutate(metric = sub("^cor_", "", metric), value = abs(value)) |>
  filter(!is.na(value))

mean_plot_df$metric <- factor(mean_plot_df$metric, levels = names(metric_labels_mean))

make_mean_plot <- function(metric_to_keep, title_text, show_y_title = TRUE) {
  df_m <- mean_plot_df |> filter(metric == metric_to_keep)
  df_m$metric <- droplevels(df_m$metric)

  ggplot(df_m, aes(x = metric, y = value, fill = center_type_name, pattern = center_type_name)) +
    geom_violin_pattern(position = position_dodge(0.8), scale = "width",
                        color = "black", pattern_fill = "black", pattern_colour = "black",
                        pattern_density = 0.15, pattern_spacing = 0.03,
                        pattern_angle = 45, pattern_size = 0.3) +
    scale_x_discrete(labels = metric_labels_mean[metric_to_keep]) +
    scale_y_continuous(limits = c(0.75, 1), breaks = seq(0, 1, 0.25)) +
    scale_fill_manual(values = center_type_colors) +
    scale_pattern_manual(values = center_type_patterns) +
    labs(x = NULL, y = if (show_y_title) "|cos(angle)|" else NULL,
         fill = "Centering approach", pattern = "Centering approach", title = title_text) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          legend.position = "none")
}

p_mean_spatial  <- make_mean_plot("meanrow_meanspat", "Which model recovers\nthe true spatial\nmain effect?", show_y_title = TRUE)
p_mean_temporal <- make_mean_plot("meancol_meantemp", "Which model recovers\nthe true temporal\nmain effect?", show_y_title = FALSE)

legend_plot_mean <- ggplot(mean_plot_df, aes(x = metric, y = value, fill = center_type_name, pattern = center_type_name)) +
  geom_violin_pattern(color = "black", pattern_fill = "black", pattern_colour = "black",
                      pattern_density = 0.15, pattern_spacing = 0.03, pattern_angle = 45, pattern_size = 0.3) +
  scale_fill_manual(values = center_type_colors) +
  scale_pattern_manual(values = center_type_patterns) +
  labs(fill = "Centering approach", pattern = "Centering approach") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", legend.title = element_blank())
# see the note below (shared_legend) on why get_plot_component(..., return_all
# = TRUE) is used instead of cowplot::get_legend() here.
shared_legend_mean <- cowplot::get_plot_component(legend_plot_mean, "guide-box-right", return_all = TRUE)
if (is.list(shared_legend_mean) && length(shared_legend_mean) == 1) shared_legend_mean <- shared_legend_mean[[1]]

p_mean_row  <- cowplot::plot_grid(p_mean_spatial, p_mean_temporal, nrow = 1, align = "h", axis = "tb")
p_mean_body <- cowplot::plot_grid(p_mean_row, shared_legend_mean, nrow = 1, rel_widths = c(1, 0.22))

title_grob_mean <- cowplot::ggdraw() +
  cowplot::draw_label("Cosine of the angle between estimated and true main effects", fontface = "bold", size = 13)

p_mean <- cowplot::plot_grid(title_grob_mean, p_mean_body, ncol = 1, rel_heights = c(0.1, 1))

ggsave("images/sim/simulation_mean_effects_comparison.png", p_mean, width = 9, height = 5, bg = "white")
cat("Done.\n")

#---

# ---- shared data prep -----------------------------------------------------------
# NOTE: renamed to metric_labels_dim_q (not metric_labels_dim) so it can never
# clobber the theta(...) version used by p_dim/p_mean above -- reassigning the
# same name to plain text broke scales::parse_format() there if this block
# ran (or was re-run) after those plots were built.
metric_labels_dim_q <- c(
  u_meanspat = "Spatial eigenvector\nvs. true spatial mean",
  v_meantemp = "Temporal eigenvector\nvs. true temporal mean",
  u_season   = "Spatial eigenvector\nvs. true seasonal pattern",
  v_season   = "Temporal eigenvector\nvs. true seasonal cycle"
)

# |cos(angle)| = |correlation|  (angle = acos(|r|)), so we plot the absolute
# correlations directly: 1 = perfect match, 0 = orthogonal.
dim_plot_df <- results_df |>
  select(center_type_name, dim,
         cor_u_meanspat, cor_v_meantemp, cor_u_season, cor_v_season) |>
  pivot_longer(-c(center_type_name, dim), names_to = "metric", values_to = "value") |>
  mutate(metric = sub("^cor_", "", metric), value = abs(value)) |>
  filter(dim %in% c(1, 2))   # drop dimension 3

dim_plot_df$metric <- factor(dim_plot_df$metric, levels = names(metric_labels_dim_q))
dim_plot_df$dim_label <- factor(paste0("Dimension ", dim_plot_df$dim),
                                levels = paste0("Dimension ", sort(unique(dim_plot_df$dim))))

# ---- common building block, with fixed panel size across all plots ---------------
panel_w <- unit(3, "cm")
panel_h <- unit(3, "cm")

make_question_plot <- function(metrics_to_keep, title_text, show_y_title = TRUE) {
  df_q <- dim_plot_df |> filter(metric %in% metrics_to_keep)
  df_q$metric <- droplevels(df_q$metric)

  ggplot(df_q, aes(x = metric, y = value, fill = center_type_name, pattern = center_type_name)) +
    geom_violin_pattern(position = position_dodge(0.8), scale = "width",
                         color = "black",
                         pattern_fill = "black",
                         pattern_colour = "black",
                         pattern_density = 0.15,
                         pattern_spacing = 0.03,
                         pattern_angle = 45,
                         pattern_size = 0.3) +
    facet_grid(dim_label ~ .) +
    scale_x_discrete(labels = metric_labels_dim_q[levels(df_q$metric)]) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +  # 1 (perfect match) at top, 0 (orthogonal) at bottom
    scale_fill_manual(values = center_type_colors) +
    scale_pattern_manual(values = center_type_patterns) +
    labs(x = NULL, y = if (show_y_title) "|cos(angle)|" else NULL,
         fill = "Centering approach", pattern = "Centering approach", title = title_text) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          strip.text.y = element_text(angle = 0),
          legend.position = "none") +
    ggh4x::force_panelsizes(rows = panel_h, cols = panel_w, respect = TRUE)
}

# ---- 3 individual sub-question plots, each with identical panel cell size --------
# movement first, then spatial, then temporal
p_movement <- make_question_plot(c("u_season", "v_season"), "Which dimensions \n captures the true \nseasonal movement?", show_y_title = TRUE)
p_spatial  <- make_question_plot("u_meanspat", "Which dimensions \n captures the true\nspatial variability?", show_y_title = FALSE)
p_temporal <- make_question_plot("v_meantemp", "Which dimensions \n captures the true\ntemporal variability?", show_y_title = FALSE)

# ---- shared legend ----------------------------------------------------------------
legend_plot <- ggplot(dim_plot_df, aes(x = metric, y = value, fill = center_type_name, pattern = center_type_name)) +
  geom_violin_pattern(color = "black",
                       pattern_fill = "black",
                       pattern_colour = "black",
                       pattern_density = 0.15,
                       pattern_spacing = 0.03,
                       pattern_angle = 45,
                       pattern_size = 0.3) +
  scale_fill_manual(values = center_type_colors) +
  scale_pattern_manual(values = center_type_patterns) +
  labs(fill = "Centering approach", pattern = "Centering approach") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        legend.key.spacing.x = unit(1.5, "cm"),
        legend.text = element_text(margin = margin(l = 6)),
        legend.key.width = unit(1, "cm"),
        legend.title = element_blank())
# cowplot::get_legend() searches for a grob named "guide-box"; ggplot2 >= 3.5
# splits this into guide-box-top/bottom/left/right (even when only one is
# populated), so get_legend() now matches several and warns, returning
# just the first -- which may not be the one you want. Since
# legend.position = "bottom" above, grab guide-box-bottom explicitly instead.
shared_legend <- cowplot::get_plot_component(legend_plot, "guide-box-bottom", return_all = TRUE)
if (is.list(shared_legend) && length(shared_legend) == 1) shared_legend <- shared_legend[[1]]

# ---- assemble: movement first, then spatial, then temporal -----------------------
p_row <- cowplot::plot_grid(p_movement, p_spatial, p_temporal,
                            nrow = 1, align = "h", axis = "tb")

p_dim <- cowplot::plot_grid(
  p_row, shared_legend,
  ncol = 1, rel_heights = c(1, 0.1)
)

title_grob <- cowplot::ggdraw() +
  cowplot::draw_label("Cosine of the angle between estimated and true components, by ecological question and SVD dimension",
                      fontface = "bold", size = 13)

p_dim_full <- cowplot::plot_grid(title_grob, p_dim, ncol = 1, rel_heights = c(0.1, 1))

ggsave("images/sim/simulation_metrics_comparison2.png", p_dim_full, width = 12, height = 6, bg = "white")


# =================================================================================
# ---- one representative simulated dataset, used by SM Figures S1.1 and S1.3 ----
# =================================================================================
# A single illustrative replicate (rather than the full 200-replicate
# distribution used above), fit under each of the three centering models --
# used below for the SM monthly-maps figure (S1.1) and the SM
# eofs_c_full-style panel (S1.3).

# ---- coastline (for map figures) -----------------------------------------------
library(sf)
library(rnaturalearth)

coast <- ne_countries(scale = "medium", returnclass = "sf")
bbox_pad <- 0.5
coast <- suppressWarnings(sf::st_crop(coast,
  xmin = min(loc_x$x) - bbox_pad, xmax = max(loc_x$x) + bbox_pad,
  ymin = min(loc_x$y) - bbox_pad, ymax = max(loc_x$y) + bbox_pad))
# NOTE: assumes loc_x$x/$y are longitude/latitude (WGS84); adjust ne_countries()
# scale or swap in your own coastline shapefile if a different CRS is used.

seed_illustration <- 1
sim1 <- simulate_once(seed_illustration)

svd_by_model <- lapply(1:3, function(center_i) {
  rc <- center_data(sim1$Zs, center_type = center_i)
  svd(rc$data_matrix_c)
})
names(svd_by_model) <- c("Spatial centering", "Temporal centering", "AMMI model")

# ---- monthly maps of the raw simulated data (one year) -- SM Figure S1.1 -------
# Shows what a single simulated dataset actually looks like, month by month,
# analogous to Figure 1 (real data) in the manuscript.

mean_by_month <- sapply(1:12, function(mo) rowMeans(sim1$Zs[, month_t == mo, drop = FALSE]))
month_labels <- month.abb

monthly_df <- purrr::map_dfr(1:12, function(k) {
  data.frame(loc_x, value = mean_by_month[, k], month = month_labels[k])
})
monthly_df$month <- factor(monthly_df$month, levels = month_labels)

p_monthly <- ggplot(monthly_df, aes(x = x, y = y, color = value)) +
  geom_point(size = 0.8) +
  geom_sf(data = coast, inherit.aes = FALSE, fill = "grey85", color = "black", linewidth = 0.3) +
  facet_wrap(~ month, ncol = 4) +
  scale_color_viridis_c(
    option = "viridis",
    trans = "log10",
    breaks = c(0.03, 0.1, 0.3, 1, 3),
    labels = scales::label_number(accuracy = 0.01),
    name = NULL
  ) +
  coord_sf(xlim = range(loc_x$x), ylim = range(loc_x$y), expand = FALSE) +
  labs(x = NULL, y = NULL, color = "value") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_blank(),
        axis.ticks = element_blank(), legend.title = element_blank())

ggsave("images/sim/simulation_illustration_monthly_maps.png", p_monthly,
       width = 11, height = 9, bg = "white")

cat("Monthly maps saved.\n")

# =================================================================================
# ---- eofs_c_full-style figure, built from the simulated dataset -- SM Fig S1.3 -
# =================================================================================
# Reproduces the real-data eofs_c_full.png layout (maps of the first 3 spatial
# dimensions on top, temporal loadings below, one column per centering model)
# using the illustrative simulated dataset (sim1, svd_by_model already computed
# above). Each model uses its own color scale (matching the real figure, where
# AMMI/Spatial/Temporal panels each have their own legend range).

n_dim_show <- 3
switch_sign_sim <- F
which_centering_to_switch_sim <- c(1, 3)   # same convention as the real script

model_names_ordered <- c("AMMI model", "Spatial centering", "Temporal centering")
svd_by_model_ordered <- svd_by_model[c("AMMI model", "Spatial centering", "Temporal centering")]

uv_plot_list_sim <- list()

for (model_name in model_names_ordered) {

  res_svd_i <- svd_by_model_ordered[[model_name]]
  n_dim_i <- n_dim_show

  if (switch_sign_sim & which(model_names_ordered == model_name) %in% which_centering_to_switch_sim) {
    res_svd_i$u <- -res_svd_i$u
    res_svd_i$v <- -res_svd_i$v
  }

  ## ---- spatial maps: Dim.1..n_dim_i -----------------------------------------
  u_df <- purrr::map_dfr(1:n_dim_i, function(d) {
    data.frame(loc_x, value = res_svd_i$u[, d], dim = paste0("Dim.", d))
  })
  u_df$dim <- factor(u_df$dim, levels = paste0("Dim.", 1:n_dim_i))
  max_abs_u <- max(abs(u_df$value))

  u_plot <- ggplot(u_df, aes(x = x, y = y, color = value, fill = value)) +
    geom_point(size = 1) +
    geom_sf(data = coast, inherit.aes = FALSE, fill = "grey85", color = "black", linewidth = 0.3) +
    facet_wrap(~ dim, nrow = 1) +
    scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                          midpoint = 0, limits = c(-max_abs_u, max_abs_u)) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-max_abs_u, max_abs_u)) +
    coord_sf(xlim = range(loc_x$x), ylim = range(loc_x$y), expand = FALSE) +
    ggtitle(label = model_name, subtitle = "") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkgray"),
          legend.title = element_blank(),
          axis.text = element_text(size = 7),
          axis.title = element_blank())

  ## ---- temporal loadings: Dim.1..n_dim_i -------------------------------------
  v_df <- purrr::map_dfr(1:n_dim_i, function(d) {
    data.frame(date = dates_t, v = res_svd_i$v[, d], dim = paste0("Dim.", d))
  })
  v_df$dim <- factor(v_df$dim, levels = paste0("Dim.", 1:n_dim_i))
  jan1_dates <- dates_t[format(dates_t, "%m-%d") == "01-01"]

  v_plot <- ggplot(v_df) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2) +
    geom_vline(xintercept = jan1_dates, linetype = "dashed", color = "skyblue") +
    geom_line(aes(x = date, y = v, group = dim), linewidth = 1, color = "steelblue") +
    facet_wrap(~ dim, ncol = 1) +
    ylim(-0.4, 0.4) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "pt")) +
    coord_cartesian(clip = "off")

  uv_plot_list_sim[[model_name]] <- cowplot::plot_grid(u_plot, v_plot, ncol = 1, rel_heights = c(0.35, 0.65))
}

gray_line_sim <- ggdraw() +
  draw_line(x = c(0.5, 0.5), y = c(0, 1), color = "darkgray", linetype = "dashed", linewidth = 0.5)

uv_plot_full_sim <- cowplot::plot_grid(
  uv_plot_list_sim[["AMMI model"]], gray_line_sim,
  uv_plot_list_sim[["Spatial centering"]], gray_line_sim,
  uv_plot_list_sim[["Temporal centering"]],
  nrow = 1, rel_widths = c(1, 0.05, 1, 0.05, 1)
)

ggsave("images/sim/simulation_eofs_c_full.png", uv_plot_full_sim, width = 6*3, height = 9, bg = "white")

cat("Simulation eofs_c_full-style figure saved.\n")
