######################
## Load and shape data
######################
## Requires: 00_init.R (grid_*, mapBase)
## Produces: S_x_df, time_step_df, loc_x, nS, nT, gridpolygon_sf
##           + images/monthly_map_plot.png

load("data/S_x_df_Solea_solea_standardise.RData")
S_x_df = S_x_df %>%
  dplyr::select(x,y,month,S_x,
                year,Year_Month,Ab) %>%
  mutate(date = as.Date(paste(year,month,"01",sep = "-"))) %>%
  mutate(S_x = S_x * Ab)

# Time step data frame
time_step_df = data.frame(date = unique(S_x_df$date),
                          t=1:length(unique(S_x_df$date)))
nT <- nrow(time_step_df)

# Spatial domain
loc_x = S_x_df %>%
  dplyr::select(x,y) %>%
  group_by(x,y) %>%
  slice(1) %>%
  as.data.frame
nS <- nrow(loc_x)

loc_x$cell = 1:nrow(loc_x)

S_x_df = inner_join(S_x_df,loc_x)

## Create grid polygon for plotting
# Coordinates limit
grid_limit <- ext(c(grid_xmin,grid_xmax,grid_ymin,grid_ymax))

# Build grid
grid <- rast(grid_limit)
res(grid) <- resol_grid
crs(grid) <- "+proj=longlat +datum=WGS84"
gridpolygon <- as.polygons(grid)
gridpolygon$layer <- c(1:length(gridpolygon$layer))
gridpolygon_sf <- st_as_sf(gridpolygon)

# Cross with gridpolygon
loc_x_sf <- st_as_sf(loc_x,
                     coords = c("x","y"),
                     crs="+proj=longlat +datum=WGS84")
loc_x_sf <- loc_x_sf[st_intersects(loc_x_sf,gridpolygon_sf) %>% lengths > 0,]
gridpolygon_sf_2 <- st_join(gridpolygon_sf,loc_x_sf) %>%
  group_by(layer,cell) %>%
  slice(1) %>%
  filter(!is.na(cell))

gridpolygon_sf = gridpolygon_sf_2 %>%
  dplyr::select(-layer) %>%
  rename(layer = cell)

monthly_map_plot <- S_x_df |>
  group_by(x,y,month) |>
  dplyr::summarise(S_x = mean(S_x)) |>
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
  coord_quickmap() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90,vjust = 0.5)
  )+
  geom_sf(data=mapBase, inherit.aes = FALSE)+
  coord_sf(xlim = c(grid_xmin,grid_xmax),
           ylim = c(grid_ymin,grid_ymax),
           expand = FALSE)

ggsave("images/monthly_map_plot.png",height = 9,width = 9)
