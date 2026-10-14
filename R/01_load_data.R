######################
## Load and shape data
######################
## Requires: 00_init.R (grid_*, resol_grid, mapBase)
## Produces: S_x_df, time_step_df, loc_x, nS, nT, gridpolygon_sf
##           + images/monthly_map_plot.png

load("data/S_x_df_Solea_solea_standardise.RData")

# Keep the needed columns, build a date, and convert the standardised
# distribution S_x into abundance by multiplying by the total abundance Ab
S_x_df <- S_x_df %>%
  dplyr::select(x, y, month, S_x, year, Year_Month, Ab) %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-")),
         S_x  = S_x * Ab)

# Time step data frame (one row per month, t = 1..nT in chronological order)
time_step_df <- data.frame(date = sort(unique(S_x_df$date)))
time_step_df$t <- seq_len(nrow(time_step_df))
nT <- nrow(time_step_df)

# Spatial domain: one row per grid cell, with a cell index
loc_x <- S_x_df %>%
  dplyr::select(x, y) %>%
  group_by(x, y) %>%
  slice(1) %>%
  ungroup() %>%
  as.data.frame()
nS <- nrow(loc_x)
loc_x$cell <- seq_len(nS)

S_x_df <- inner_join(S_x_df, loc_x, by = c("x", "y"))

## Grid polygons for plotting -------------------------------------------------
# Regular grid covering the study area, at resolution resol_grid
grid_limit <- ext(c(grid_xmin, grid_xmax, grid_ymin, grid_ymax))
grid <- rast(grid_limit)
res(grid) <- resol_grid
crs(grid) <- "+proj=longlat +datum=WGS84"
gridpolygon <- as.polygons(grid)
gridpolygon$layer <- seq_len(nrow(gridpolygon))
gridpolygon_sf <- st_as_sf(gridpolygon)

loc_x_sf <- st_as_sf(loc_x,
                     coords = c("x", "y"),
                     crs = "+proj=longlat +datum=WGS84")
loc_x_sf <- loc_x_sf[lengths(st_intersects(loc_x_sf, gridpolygon_sf)) > 0, ]

gridpolygon_sf <- st_join(gridpolygon_sf, loc_x_sf) %>%
  filter(!is.na(cell)) %>%
  arrange(layer, cell) %>%
  filter(!duplicated(paste(layer, cell))) %>%
  dplyr::select(-layer) %>%
  rename(layer = cell)

## Monthly mean distribution map -----------------------------------------------
monthly_map_plot <- S_x_df |>
  group_by(x, y, month) |>
  dplyr::summarise(S_x = mean(S_x), .groups = "drop") |>
  ggplot(aes(x = x, y = y, fill = S_x)) +
  geom_raster() +
  facet_wrap(~ month, ncol = 4) +
  scale_fill_viridis_c(
    option = "viridis",
    trans = "log10",
    breaks = c(0.03, 0.1, 0.3, 1, 3),
    labels = scales::label_number(accuracy = 0.01),
    name = NULL
  ) +
  geom_sf(data = mapBase, inherit.aes = FALSE) +
  coord_sf(xlim = c(grid_xmin, grid_xmax),
           ylim = c(grid_ymin, grid_ymax),
           expand = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

ggsave("images/monthly_map_plot.png", monthly_map_plot,
       height = 9, width = 9)
