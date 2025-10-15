# =========================================================
# ORTHOGRAPHIC GLOBE — Africa centered + African capitals
# Fast + cached:
#   • risk points on land (st_join)        → cache/risk_*.rds
#   • projected coastlines (already XY)    → cache/coast_*.rds
# Output: globe_cs_ortho_africa.png
# =========================================================

# ------------ Packages ------------
pkgs <- c("sf","rnaturalearth","rnaturalearthdata","ggplot2",
          "dplyr","tibble","scales","viridisLite","ggforce")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, require, character.only = TRUE))

# ------------ Parameters ------------
lon0 <- 20     # Africa-centered longitude
lat0 <-  7     # Africa-centered latitude
step_deg <- 1.2         # risk grid step (bigger = faster)
coast_scale  <- "small" # coastline detail: "small" (110m) or "medium"
simplify_tol <- 0.20    # degrees for st_simplify (0–0.5 typical)
stride       <- 1       # keep every Nth vertex after simplify
out_file     <- "globe_cs_ortho_africa.png"
bbox_vals <- c(xmin=-170, ymin=-58, xmax=190, ymax=80)
force_refresh_risk  <- FALSE
force_refresh_coast <- FALSE

# ------------ Cache paths ------------
dir.create("cache", showWarnings = FALSE)
risk_cache <- file.path(
  "cache",
  sprintf("risk_lonlat_land_step%.2f_%d_%d_%d_%d.rds",
          step_deg, bbox_vals["xmin"], bbox_vals["xmax"],
          bbox_vals["ymin"], bbox_vals["ymax"])
)
coast_cache <- file.path(
  "cache",
  sprintf("coast_ortho_lon%g_lat%g_scale%s_tol%.2f_stride%d.rds",
          lon0, lat0, coast_scale, simplify_tol, stride)
)

# ------------ Helpers (no PROJ needed) ------------
deg2rad <- function(x) x*pi/180
lonlat_to_ortho <- function(lon, lat, lon0, lat0){
  lam  <- deg2rad(lon);  phi  <- deg2rad(lat)
  lam0 <- deg2rad(lon0); phi0 <- deg2rad(lat0)
  dlam <- lam - lam0
  x <- cos(phi) * sin(dlam)
  y <- cos(phi0)*sin(phi) - sin(phi0)*cos(phi)*cos(dlam)
  cosC <- sin(phi0)*sin(phi) + cos(phi0)*cos(phi)*cos(dlam)
  data.frame(x=x, y=y, visible = (cosC >= 0))
}

# ------------ Land polygons (quick) ------------
land <- rnaturalearth::ne_countries(scale="medium", returnclass="sf") |>
  dplyr::select(geometry)
bbx <- st_as_sfc(st_bbox(bbox_vals, crs = st_crs(land)))
old_s2 <- sf_use_s2(); sf_use_s2(FALSE)
land <- st_make_valid(land)
land <- suppressMessages(st_intersection(land, bbx))
sf_use_s2(old_s2)

# =========================================================
# 1) RISK SURFACE on LAND — cached
# =========================================================
if (file.exists(risk_cache) && !force_refresh_risk) {
  risk <- readRDS(risk_cache)
  message("Loaded cached risk points: ", risk_cache)
} else {
  lon <- seq(bbox_vals["xmin"], bbox_vals["xmax"], by = step_deg)
  lat <- seq(bbox_vals["ymin"], bbox_vals["ymax"], by = step_deg)
  grid <- expand.grid(lon=lon, lat=lat)

  # vivid tropics: latitude weight + gentle texture
  lat_w <- with(grid, cos(pmin(abs(lat), 75)*pi/180))^1.3
  tex   <- with(grid,
                sin((lon+25)*pi/35)*cos((lat-8)*pi/22) +
                  0.5*sin((lon-90)*pi/18)*0.4*cos((lat+6)*pi/15))
  grid$val <- scales::rescale(0.65*lat_w + 0.35*scales::rescale(tex))

  # mask to land once
  old_s2 <- sf_use_s2(); sf_use_s2(FALSE)
  grid_sf <- st_as_sf(grid, coords = c("lon","lat"), crs = 4326)
  grid_sf$val <- grid$val
  message("Computing land mask (st_join)… (will cache)")
  grid_on_land <- suppressMessages(st_join(grid_sf, land, join = st_within, left = FALSE))
  sf_use_s2(old_s2)

  coords <- st_coordinates(grid_on_land) |> as.data.frame()
  risk <- tibble::tibble(lon = coords$X, lat = coords$Y, val = grid_on_land$val)
  saveRDS(risk, risk_cache)
  message("Saved: ", risk_cache)
}

# =========================================================
# 2) FAST COASTLINES — projected & cached
# =========================================================
if (file.exists(coast_cache) && !force_refresh_coast) {
  land_v <- readRDS(coast_cache)
  message("Loaded cached projected coastlines: ", coast_cache)
} else {
  message("Building fast coastlines… (will cache)")
  coast <- rnaturalearth::ne_coastline(scale = coast_scale, returnclass = "sf")

  old_s2 <- sf_use_s2(); sf_use_s2(FALSE)
  coast <- suppressMessages(st_intersection(coast, st_as_sfc(st_bbox(bbox_vals, crs = st_crs(coast)))))
  if (simplify_tol > 0) {
    coast <- st_simplify(coast, dTolerance = simplify_tol, preserveTopology = TRUE)
  }
  sf_use_s2(old_s2)

  lc <- st_coordinates(st_cast(coast, "MULTILINESTRING"))
  if (stride > 1L && nrow(lc) > 0L) lc <- lc[seq(1, nrow(lc), by = stride), , drop = FALSE]
  coast_df <- as.data.frame(lc); names(coast_df)[1:2] <- c("lon","lat")

  P <- lonlat_to_ortho(coast_df$lon, coast_df$lat, lon0, lat0)
  coast_df$x <- P$x; coast_df$y <- P$y
  coast_df$visible <- P$visible & (coast_df$x^2 + coast_df$y^2 <= 1 + 1e-9)

  grp_cols <- intersect(c("L1","L2","L3","L4"), colnames(coast_df))
  coast_df$grp <- if (length(grp_cols)) interaction(coast_df[, grp_cols], drop = TRUE) else 1

  land_v <- coast_df[coast_df$visible, c("x","y","grp"), drop = FALSE]
  saveRDS(land_v, coast_cache)
  message("Saved projected coastlines cache: ", coast_cache)
}

# =========================================================
# 3) Project risk & build AFRICAN CAPITALS as sites
# =========================================================
# --- African capitals robustly (no st_union) ---
message("Fetching African capitals (robust)…")

# 1) Download populated places
pp <- rnaturalearth::ne_download(
  scale = 110, type = "populated_places",
  category = "cultural", returnclass = "sf"
)

# 2) Filter to national capitals if field exists, else fall back to FEATURECLA
if ("ADM0CAP" %in% names(pp)) {
  pp <- dplyr::filter(pp, ADM0CAP == 1)
} else if ("FEATURECLA" %in% names(pp)) {
  pp <- dplyr::filter(pp, grepl("Admin-0 capital", as.character(FEATURECLA)))
}

# 3) Get Africa polygons
af <- rnaturalearth::ne_countries(
  scale = "small", continent = "Africa", returnclass = "sf"
) |> dplyr::select(geometry)

# 4) Do the spatial filter via st_join with s2 OFF (no union needed)
old_s2 <- sf::sf_use_s2(); sf::sf_use_s2(FALSE)
pp  <- sf::st_make_valid(pp)
af  <- sf::st_make_valid(af)

cap_af <- suppressMessages(sf::st_join(pp, af, join = sf::st_within, left = FALSE))
sf::sf_use_s2(old_s2)

# 5) Project visible capitals to orthographic
coords_cap <- sf::st_coordinates(cap_af)
cap_df <- tibble::tibble(
  name = dplyr::coalesce(cap_af$NAME, cap_af$name),
  lon  = coords_cap[,1],
  lat  = coords_cap[,2]
)

Pcap <- lonlat_to_ortho(cap_df$lon, cap_df$lat, lon0, lat0)
sites <- cap_df[Pcap$visible & (Pcap$x^2 + Pcap$y^2 <= 1 + 1e-9), ]
sites$x <- Pcap$x[Pcap$visible & (Pcap$x^2 + Pcap$y^2 <= 1 + 1e-9)]
sites$y <- Pcap$y[Pcap$visible & (Pcap$y^2 + Pcap$y^2 <= 1 + 1e-9)]

# =========================================================
# 4) CS-style background (outside the globe) — robust hex mesh
# =========================================================
rings <- data.frame(r = c(1.05, 1.12, 1.20))
ang <- seq(0, 2*pi, length.out = 24)[-1]
rad_df <- data.frame(
  x  = cos(ang)*1.00, y  = sin(ang)*1.00,
  x2 = cos(ang)*1.28, y2 = sin(ang)*1.28
)

hex_spacing <- 0.20
xs <- seq(-1.5, 1.5, by = hex_spacing)
ys <- seq(-1.5, 1.5, by = hex_spacing*sqrt(3)/2)
hex_centers <- do.call(rbind, lapply(seq_along(ys), function(i){
  xoff <- if (i %% 2 == 0) 0.5*hex_spacing else 0
  data.frame(x0 = xs + xoff, y0 = ys[i])
}))
hex_centers <- hex_centers[(hex_centers$x0^2 + hex_centers$y0^2) > 1.0^2 &
                             abs(hex_centers$x0) <= 1.35 & abs(hex_centers$y0) <= 1.35, , drop = FALSE]

build_hex_polys <- function(centers, r){
  if (!nrow(centers)) return(data.frame(x=numeric(0), y=numeric(0), grp=integer(0)))
  ang <- seq(0, 2*pi, length.out = 7)[-7]
  do.call(rbind, lapply(seq_len(nrow(centers)), function(i){
    data.frame(
      x   = centers$x0[i] + r * cos(ang),
      y   = centers$y0[i] + r * sin(ang),
      grp = i
    )
  }))
}
hex_r   <- hex_spacing * 0.35
hex_poly <- build_hex_polys(hex_centers, hex_r)

# ================================
# Plot & export (dark planet, light outside)
# ================================

# Light outside, dark planet
outside_bg  <- "#E4F0FB"   # light blue outside the globe
planet_hex  <- "#0A2342"   # dark navy *inside* the globe
grid_col    <- scales::alpha("#2E6FA5", 0.40)
radial_col  <- scales::alpha("#2E6FA5", 0.38)
horizon_col <- "#1E3A5F"

# Start radial lines just *outside* the horizon so none draw over the planet
ang <- seq(0, 2*pi, length.out = 24)[-1]
rad_df <- data.frame(
  x  = cos(ang)*1.02, y  = sin(ang)*1.02,
  x2 = cos(ang)*1.28, y2 = sin(ang)*1.28
)

p <- ggplot() +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = outside_bg, colour = NA),
    panel.background = element_rect(fill = outside_bg, colour = NA),
    plot.margin      = margin(20,20,20,20)
  ) +
  coord_equal(xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35), expand = FALSE, clip = "off") +

  # --- OUTSIDE graphics first ---
  {if (exists("hex_poly") && nrow(hex_poly))
    geom_polygon(data = hex_poly, aes(x, y, group = grp),
                 fill = NA, colour = grid_col, linewidth = 0.25)
  } +
  geom_segment(data = rad_df, aes(x=x, y=y, xend=x2, yend=y2),
               colour = radial_col, linewidth = 0.6, lineend = "round") +
  ggforce::geom_circle(data = data.frame(r=c(1.05,1.12,1.20)),
                       aes(x0=0, y0=0, r=r),
                       colour = grid_col, linewidth = 0.7, fill = NA) +

  # --- PLANET fill (this ensures the interior is dark navy) ---
  ggforce::geom_circle(aes(x0=0, y0=0, r=1.0),
                       fill = planet_hex, colour = NA) +

  # --- INSIDE-the-planet layers on top of the dark fill ---
  geom_point(data = risk_v, aes(x, y, colour = val),
             size = 1.4, stroke = 0, alpha = 0.95) +
  scale_colour_gradientn(colours = viridisLite::plasma(11), limits = c(0,1), guide = "none") +
  geom_path(data = land_v, aes(x, y, group = grp),
            colour = "#3B5A7A", linewidth = 0.28, alpha = 0.95) +
  ggforce::geom_circle(aes(x0=0, y0=0, r=1.0),
                       colour = horizon_col, linewidth = 0.9, fill = NA) +
  geom_point(data = sites, aes(x, y),
             shape = 21, size = 3.0, stroke = 1.0,
             fill = "white", colour = "#0b7285")

ggsave(out_file, p, width = 8, height = 8, dpi = 450, bg = outside_bg)
message("Saved: ", out_file)
