#################################
## Empirical Orthogonal Functions
#################################
## Requires: 00_init.R, 01_load_data.R (S_x_df, loc_x, time_step_df,
##           gridpolygon_sf, nT, nS, center_data(), plot_factors(),
##           plot_loadings())
## Reproduces manuscript Figures 2 (right panel), 4, 5, 6 and Tables 1-2,
## and Supplementary Material Figures S1.4, S1.5 (six-dimension real-data
## SVD, AMMI + spatial centering only) and S1.6 (AMMI residuals diagnostic).
## Also writes res/*.RData, consumed by 04_simulation_metrics_figures.R.

## Transform biomass data frame into wide matrix
toto <- S_x_df %>% 
  # mutate(S_x = log(S_x)) |> 
  filter(!is.na(S_x) | !is.na(x) | !is.na(y) | !is.na(cell)) %>%
  dplyr::select(cell,date,S_x) %>% 
  arrange(cell) %>% 
  pivot_wider(names_from = date,values_from = S_x) %>% 
  as.matrix() %>% t()
Zs <- t(toto[-1,])

# Plot spatial and temporal term ------------------------------------------

## Spatial centering of the spatio-temporal dataset
## --> compute spatial effect
spat_average <- apply(Zs, 1, mean,na.rm=T) # compute spatial average pattern
Zs_c <- Zs - outer(spat_average, rep(1, nT))

spatial_pattern_df = data.frame(loc_x,S_x=spat_average - mean(Zs)) %>% 
  inner_join(gridpolygon_sf,by = c("cell"="layer")) %>% 
  st_as_sf

max_abs_S_x <- max(abs(spatial_pattern_df$S_x))

spatial_pattern_plot = ggplot()+
  geom_sf(data=spatial_pattern_df,aes(col=S_x,fill=S_x))+
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
  geom_sf(data=mapBase)+
  coord_sf(xlim = c(grid_xmin,grid_xmax),
           ylim = c(grid_ymin,grid_ymax),
           expand = FALSE)+
  theme_bw()+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5,face = "bold"))+
  xlab("")+ylab("")+
  ggtitle(label = "Mean spatial pattern",
          subtitle = expression("(" * mu + bold(alpha) * " )"))
plot(spatial_pattern_plot)

## Temporal centering of the spatio-temporal dataset
## --> compute temporal effect
temp_average <- apply(Zs_c, 2, mean) # compute temporal average pattern
Zs_c <- t(t(Zs_c) - outer(temp_average,rep(1, nS))) # center on mean temporal trend

temp_average_plot <- apply(Zs, 2, mean) # compute temporal average pattern
temporal_trend_df = data.frame(time_step_df,S_x=temp_average_plot)
temporal_trend_plot = ggplot()+
  geom_line(data=temporal_trend_df,aes(x = date,y=S_x))+
  theme_classic()+
  theme(aspect.ratio = 1/4,
        plot.title = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5,face = "bold"))+
  xlab("")+ylab("")+ylim(0,NA)+
  ggtitle(label = "Mean temporal trend",
          subtitle = expression("(" * mu + bold(beta) * " )"))
plot(temporal_trend_plot)

spatial_temporal_trend_plot <- cowplot::plot_grid(spatial_pattern_plot,temporal_trend_plot,
                                                  ncol = 1,rel_heights = c(1,0.75))

ggsave(filename = paste0("images/",spp_to_plot,"/spatial_temporal_trend_plot.png"),
       width = 6,height = 6,bg = "white")

# Perform svd on the different centered matrices -----------------------------------------
n_dim <- 3 # number of dimensions to plot
n_rank <- min(dim(Zs))
switch_sign <- T # switch sign of the EOFs matrices
which_centering_to_switch <- c(1,3) # which centering is to be switched ?
uv_plot_list <- list()
per_var_plot_list <- list()
type_center_df <- data.frame(center_type = c(1,2,3),
                             center_type_name = c("Spatial centering","Temporal centering","AMMI model"))
type_center_df$center_type_name <- factor(type_center_df$center_type_name,levels = c("AMMI model","Spatial centering","Temporal centering"))
residual_family <- "laplace" # investigate the type of distribution for residuals : "normal", "student", "laplace", "skewed-normal"

# colorblind-safe color per SVD dimension (Okabe-Ito palette), used for all time series plots below
dim_colors <- c("Dim.1" = "steelblue", "Dim.2" = "steelblue", "Dim.3" = "steelblue",
                "Dim.4" = "steelblue", "Dim.5" = "steelblue", "Dim.6" = "steelblue")

n_dim <- 3
n_dim_plot <- c(3,1,3)
n_dim_filter <- c(6,1,6)

for(center_i in 1:3){
  
  ## Center the data
  #-----------------
  res_center_data <- center_data(Zs,center_type = center_i)
  Zs_c  <- res_center_data$data_matrix_c
  
  ## Perform svd
  #--------------
  res_svd_i <- svd(Zs_c)
  ncp <- FactoMineR::estim_ncp(Zs_c)$ncp
  
  L <- res_svd_i$d * res_svd_i$d
  pct <- L / sum(L)
  # plot(pct[1:15])

  ## Switch sign if necessary
  if(switch_sign & center_i %in% which_centering_to_switch){
    res_svd_i$u <- -res_svd_i$u
    res_svd_i$v <- -res_svd_i$v
  }
  
  ## Plot factors and loadings
  #---------------------------
  u_plot <- plot_factors(support = loc_x,res_svd = res_svd_i,n_dim = n_dim_plot[center_i])+
    scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0)+
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0)+
    ggtitle(label = type_center_df$center_type_name[center_i],subtitle = "")+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "darkgray"))
  v_plot <- plot_loadings(time_step_df,res_svd_i,n_dim = n_dim_plot[center_i]) # plot loadings
  v_plot <- v_plot + scale_color_manual(values = dim_colors) # colored by dimension

  ## Manual correction so that pictures have the same alignment -------------------------------
  if(center_i == 2){
    v_df <- data.frame(time_step_df,Dim=res_svd_i$v[,1:n_dim_plot[center_i]])
    v_df <- rename(v_df,Dim.1 = Dim)
    v_df <- v_df %>% 
      pivot_longer(cols = num_range("Dim.",1:n_dim)) %>% 
      rename(v = value,dim = name)
    v_df <- rbind(v_df,
      c(NA,NA,"Dim.2",NA),
      c(NA,NA,"Dim.3",NA))
    v_df$t <- as.integer(v_df$t)
    v_df$v <- as.numeric(v_df$v)
    v_plot = ggplot(v_df)+
      geom_hline(yintercept=0,linetype="dashed", color = "black", size = 0.2)+
      geom_line(size=1,aes(x=date,y=v,col=dim,group=dim))+
      scale_color_manual(values = dim_colors)+
      theme_bw()+
      theme()+
      xlab("")+ylab("")+
      theme(legend.title = element_blank(),
            aspect.ratio = c(1/4),
            legend.position = "none",
            plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "pt"))+
      coord_cartesian(clip = "off")+
      facet_wrap(.~dim,ncol = 1)+
      geom_vline(xintercept = v_df$date[which(str_detect(v_df$date,"01-01"))],
                 linetype="dashed",color = "skyblue")  
  }
  #----------------------------------------------------------------------------------------------

  v_plot <- v_plot + ylim(-0.4,0.4)
  uv_plot <- cowplot::plot_grid(u_plot,v_plot,ncol = 1,rel_heights = c(0.35,0.65)) ## Plot and save factors and loadings
  ggsave(paste0("images/",spp_to_plot,"/eofs_c",center_i,".png"),width = 9*3,height = 9,bg = "white")
  uv_plot_list[[center_i]] <- uv_plot
  
  ## Plot additional dimensions (SM Figures S1.4 - AMMI - and S1.5 - spatial
  ## centering; the SM has no six-dimension figure for temporal centering,
  ## so this block is skipped for center_i == 2)
  if(center_i %in% c(1,3)){
    plot_fact <- plot_factors(support = loc_x,res_svd = res_svd_i,n_dim = 6)+
      scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0)+
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0)
    plot_load <- plot_loadings(time_step_df,res_svd_i,n_dim = 6) +
      scale_color_manual(values = dim_colors)
    plot_fact_load <- cowplot::plot_grid(plot_fact,plot_load,nrow = 1)
    ggsave(paste0("images/",spp_to_plot,"/additionnal_dim/eofs_c",center_i,"_additional_dim.png"),width = 9*1.25,bg = "white")
  }

  ## Compute sum of squares
  #------------------------
  ss_vec <- c()
  total_ss <- sum((Zs - mean(Zs))^2)
  # spatial_ss <- ncol(Zs) * sum((rowMeans(Zs) - mean(Zs))^2)
  spatial_ss <- ncol(Zs) * sum((res_center_data$row_mean)^2)
  # temporal_ss <- nrow(Zs) * sum((colMeans(Zs) - mean(Zs))^2)
  temporal_ss <- nrow(Zs) * sum((res_center_data$col_mean)^2)
  interac_ss <- (res_svd_i$d)^2
  
  ss_df <- data.frame(ss = c(spatial_ss,temporal_ss,interac_ss),
                      dim = c("Spatial","Temporal",paste0("SVD.",1:length(interac_ss))),
                      center_type = center_i,
                      SVD_dim_in = c(NA,NA,rep("EOF",n_dim_filter[center_i]),rep("Residuals",min(dim(res_center_data$data_matrix_c)) - n_dim_filter[center_i])))
  if(center_i == 1) ss_df_full <- ss_df
  if(center_i != 1) ss_df_full <- rbind(ss_df,ss_df_full)
  
  # ## Estimate the number of dimensions to filter
  # filt_dim_est <- FactoMineR::estim_ncp(Zs_c,scale = F)
  # print(paste0("Nb of dimensions estim_ncp: ",filt_dim_est$ncp," | Kaiser Guttman: ",length(which(L > mean(L)))))
  # plot(filt_dim_est$criterion)
  ## SM Figure S1.6 is the AMMI-model residuals diagnostic only (no
  ## corresponding figure for spatial/temporal centering in the SM).
  if(analyse_residuals & center_i == 3){

    ## Look at residuals
    #-------------------
    ncp <- n_dim_filter[center_i]
    residuals_mat <- res_svd_i$u[,(ncp+1):n_rank] %*% diag(res_svd_i$d[(ncp+1):n_rank]) %*% t(res_svd_i$v[,(ncp+1):n_rank])
    residuals_df <- residuals_mat |> as.vector() |> as.data.frame()
    residuals_df <- residuals_df |>
      mutate(index = 1:nrow(residuals_df))
    colnames(residuals_df)[1] <- "residuals"
    
    x_resid <- residuals_df$residuals
    # x_resid <- residuals_df$residuals[sample(x = 1:nrow(residuals_df), size = 1000)]
    sd_res  <- sd(x_resid)
  
    # --- Student-t: estimate df from excess kurtosis (method of moments) --------
    kurt_res <- mean(((x_resid - mean(x_resid)) / sd_res)^4) - 3
    df_t     <- if (kurt_res > 0) max(4 + 6 / kurt_res, 4.01) else 30
    scale_t  <- sd_res * sqrt((df_t - 2) / df_t)
  
    # --- Skew-t: MLE fit (handles strong skew + heavy tails together) ----------
    # st_fit <- sn::selm(x_resid ~ 1, family = "ST",
    #   opt.method = "nlminb")   # try an alternative optimizer    dp_st_raw <- coef(st_fit, param.type = "DP")
    # st_fit <- sn::selm(x_resid ~ 1, family = "ST")
    # dp_st <- as.numeric(dp_st_raw)
    # names(dp_st) <- if (!is.null(names(dp_st_raw))) names(dp_st_raw) else colnames(dp_st_raw)
    # if (is.null(names(dp_st))) names(dp_st) <- c("xi", "omega", "alpha", "nu")
  
    families <- list(
      normal = list(
        label   = "Normal",
        qtheo   = function(p) as.numeric(qnorm(p, mean = 0, sd = sd_res)),
        dtheo   = function(x) dnorm(x, mean = 0, sd = sd_res)
      ),
      student = list(
        label   = sprintf("Student-t"),
        qtheo   = function(p) as.numeric(scale_t * qt(p, df = df_t)),
        dtheo   = function(x) dt(x / scale_t, df = df_t) / scale_t
      ),
      laplace = list(
        label   = "Laplace",
        qtheo   = function(p) as.numeric(qlaplace(p, scale = sd_res / sqrt(2))),
        dtheo   = function(x) dlaplace(x, scale = sd_res / sqrt(2))
      ) # ,
      # skewt = list(
      #   label = sprintf("Skew-t\n(\u03b1=%.1f, \u03bd=%.1f)", dp_st["alpha"], dp_st["nu"]),
      #   qtheo = function(p) as.numeric(sn::qst(p, dp = dp_st)),
      #   dtheo = function(x) sn::dst(x, dp = dp_st)
      # )
    )
  
    ## One row per distribution: label | QQ-plot | histogram -------------------
    n_fam <- length(families)
    clip_lim <- c(-3, 3)   # force displayed range for both plots

    dir.create(paste0("images/", spp_to_plot), showWarnings = FALSE, recursive = TRUE)
    png(paste0("images/", spp_to_plot, "/residuals_qq_hist_c", center_i, ".png"),
        width = 9 / 1.5, height = 4.5 / 1.5 * n_fam, units = "in", res = 300)

    layout(matrix(1:(n_fam * 3), nrow = n_fam, byrow = TRUE), widths = c(0.15, 0.425, 0.425))

    for (fam in families) {

      ## row label panel (rotated 90 degrees)
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

      ## Histogram + fitted density, clipped to [-3, 3]
      x_resid_clip <- x_resid[x_resid >= clip_lim[1] & x_resid <= clip_lim[2]]
      grid_x <- seq(clip_lim[1], clip_lim[2], length.out = 1000)
      dens_max <- max(
        hist(x_resid_clip, plot = FALSE, breaks = seq(clip_lim[1], clip_lim[2], length.out = 30))$density,
        density(x_resid_clip)$y,
        fam$dtheo(grid_x)
      )

      par(mar = c(4, 4, 2, 1), pty = "s")
      hist(x_resid_clip, freq = FALSE, ylim = c(0, dens_max * 1.05), xlim = clip_lim,
          breaks = seq(clip_lim[1], clip_lim[2], length.out = 30),
          main = "",
          xlab = "Residuals")
      lines(density(x_resid_clip), col = "blue", lwd = 3)
      curve(fam$dtheo(x), add = TRUE, col = "red", lwd = 2, from = clip_lim[1], to = clip_lim[2])
    }

    dev.off()
  
  }
}

## Redefine factor for sum of squares dataframe
ss_df_full$dim <- factor(ss_df_full$dim,
  levels = unique(ss_df_full$dim))
ss_df_full$center_type_name <- NA
ss_df_full$center_type_name[which(ss_df_full$center_type == 1)] <- "Spatial centering"
ss_df_full$center_type_name[which(ss_df_full$center_type == 2)] <- "Temporal centering"
ss_df_full$center_type_name[which(ss_df_full$center_type == 3)] <- "AMMI model"
ss_df_full$center_type_name <- factor(ss_df_full$center_type_name,
               levels = c("AMMI model","Spatial centering","Temporal centering"))

prop_var_plot <- ss_df_full |> 
  mutate(dim = as.character(dim)) |> 
  mutate(dim_num = as.numeric(str_remove(dim,"SVD."))) |> 
  mutate(dim = ifelse(str_detect(dim,"SVD"),"SVD",dim)) |> 
  mutate(dim = ifelse(dim == "SVD",SVD_dim_in,dim)) |> 
  mutate(dim = ifelse(dim == "EOF",paste0("Interaction Dim.",dim_num),dim)) |> 
  group_by(dim,center_type_name) |> 
  dplyr::summarise(ss = sum(ss))

prop_var_plot$dim <- factor(prop_var_plot$dim,levels = c("Spatial",
  "Temporal",
  paste0("Interaction Dim.",1:max(n_dim_filter)),
  "Residuals"))

prop_var_plot <- prop_var_plot %>%
  mutate(pattern_group = case_when(
    dim == "Spatial"   ~ "stripe",
    dim == "Temporal"  ~ "crosshatch",
    dim == "Residuals" ~ "none",
    TRUE               ~ "circle"   # all Interaction Dim.1-6
  ))

skyblue_palette <- rev(c(
  "#E3F6FF", "#CFEFFF", "#9BD7F0", "#87CEEB", "#5FA8D3", "#2E8BC0"
))

pattern_lookup <- c(
  "Spatial"          = "stripe",
  "Temporal"         = "crosshatch",
  "Interaction Dim.1" = "circle",
  "Interaction Dim.2" = "circle",
  "Interaction Dim.3" = "circle",
  "Interaction Dim.4" = "circle",
  "Interaction Dim.5" = "circle",
  "Interaction Dim.6" = "circle",
  "Residuals"        = "none"
)

lm_var_expl_plot <- ggplot(data = prop_var_plot,
                           aes(x = factor(center_type_name), y = ss / total_ss,
                               fill = dim, pattern = dim)) +   # <- pattern mapped to dim, not pattern_group
  geom_col_pattern(position = "stack",
                   color = "black",
                   pattern_fill = "black",
                   pattern_colour = "black",
                   pattern_density = 0.15,
                   pattern_spacing = 0.03,
                   pattern_angle = 45,
                   pattern_size = 0.3) +
  scale_pattern_manual(values = pattern_lookup) +   # no guide = "none" here
  scale_fill_manual(values = c("orange", "mediumseagreen", skyblue_palette, "lightgrey")) +
  labs(x = "Center Type", y = "Proportion of variance", fill = "Dimension", pattern = "Dimension") +
  theme_bw() +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(size = 10.5)) +
  ggtitle("Proportion of variance explained \n by each term of the linear models")

## Latex Table of variance explained by the AMMI model
table_perc_var  <- prop_var_plot |> 
  filter(center_type_name == "AMMI model") |> 
  dplyr::select(-center_type_name)

table_perc_var |> 
  mutate(perc_var = ss/sum(table_perc_var$ss))
## Then transform to latex with chatgpt

## Plot variance explained by each dimension of the SVD
vline_df <- data.frame(dim = n_dim_filter + 1,
  center_type_name = c("AMMI model","Temporal centering","Spatial centering"))
svd_var_expl_plot <- ss_df_full |> 
  filter(str_detect(dim,"SVD.")) |> 
  dplyr::mutate(dim = as.numeric(str_remove(dim,"SVD."))) |> 
  filter(dim <= 20) |> 
  ggplot()+
  geom_line(aes(x = dim, y = ss / total_ss))+
  geom_point(aes(x = dim, y = ss / total_ss),size=0.75)+
  facet_wrap(.~ center_type_name,scales = "free_y")+
  theme_bw()+
  theme(aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5,face = "bold"),
    legend.title = element_blank()
  )+
  ggtitle("Proportion of variance explained by the SVD dimensions \n conducted on the matrix R (interaction + residuals)")+
  ylab("Proportion of variance")+
  xlab("Dimensions")+
  geom_vline(data=vline_df,aes(xintercept = dim),color = "skyblue")

lm_var_expl_plot2 <- cowplot::plot_grid(lm_var_expl_plot,NULL,nrow = 1,rel_widths = c(0.7,0.3))
prop_var_expl_plot <- cowplot::plot_grid(lm_var_expl_plot2,svd_var_expl_plot,ncol=1)

ggsave(paste0("images/",spp_to_plot,"/prop_var_expl.png"),width = 11,height = 11,bg = "white")

## Plot all loadings and factors
gray_line <- ggdraw() + # vertical line
  draw_line(
    x = c(0.5,0.5),
    y = c(0,1),
    color = "darkgray", linetype = "dashed",
    size = 0.5
  )
uv_plot_full <- cowplot::plot_grid(uv_plot_list[[3]],gray_line,
                   uv_plot_list[[1]],gray_line,
                   uv_plot_list[[2]],nrow = 1,
                   rel_widths = c(1,0.05,1,0.05,1))
ggsave(paste0("images/",spp_to_plot,"/eofs_c_full.png"),width = 6*3,height = 9,bg = "white")

# ## Plot percentage of variance captured by each dimension
# perc_var_df_full$dim <- factor(perc_var_df_full$dim,
#                                levels = unique(perc_var_df_full$dim))
# perc_var_df_full$center_type_name <- NA
# perc_var_df_full$center_type_name[which(perc_var_df_full$center_type == 1)] <- "Spatial centering"
# perc_var_df_full$center_type_name[which(perc_var_df_full$center_type == 2)] <- "Temporal centering"
# perc_var_df_full$center_type_name[which(perc_var_df_full$center_type == 3)] <- "AMMI model"
# perc_var_df_full$center_type_name <- factor(perc_var_df_full$center_type_name,
#                                             levels = c("AMMI model","Spatial centering","Temporal centering"))
# per_var_full_plot <- ggplot(perc_var_df_full, 
#                             aes(x = dim, y = perc_var, fill = dim))+
#   geom_bar(stat = "identity",col = "black")+
#   facet_wrap(~center_type_name, nrow = 1)+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         legend.position = "none",aspect.ratio = 1)+
#   labs(x = "Dimension", y = "Percentage of variance explained")+xlab("")+
#   scale_fill_manual(values = c(  "#E69F00", "#56B4E9", # Nice distinct colors  
#                                  "#F8766D", "#00BA38", "#619CFF", # Your requested colors (3-5)
#                                  "grey" )) # Last color must be grey

# ggsave(paste0("images/perc_var_full_plot.png"),width = 9,height = 4,bg = "white")
# ---------------------------------------------------------------------------
# Save intermediate objects consumed by the simulation scripts
# (04_simulation_metrics_figures.R and 05_simulation_recovery_study.R draw
# their spatial support, temporal support, mean spatial pattern, mean
# temporal trend and a real leading spatial EOF -- used as the "movement"
# pattern -- from this real-data fit, so the simulation is calibrated on the
# actual sole dataset rather than an arbitrary synthetic shape).
#
# NOTE: the original draft of this block (commented out in the legacy
# workflow.R) called save(data = loc_x, file = ...), which is invalid --
# save() has no `data` argument, so that line silently saved nothing under
# the intended object name. Fixed below to save the objects themselves.
# `res_svd_i` here is the SVD from the last loop iteration above
# (center_i = 3, i.e. the AMMI / double-centering model), matching what the
# legacy comment intended.
# ---------------------------------------------------------------------------
save(loc_x, file = "res/loc_x.RData")
save(time_step_df, file = "res/time_step_df.RData")
save(spatial_pattern_df, file = "res/spatial_pattern_df.RData")
save(temporal_trend_df, file = "res/temporal_trend_df.RData")

movement_df <- as.data.frame(res_svd_i$u[, 1:2])
save(movement_df, file = "res/movement_df.RData")
