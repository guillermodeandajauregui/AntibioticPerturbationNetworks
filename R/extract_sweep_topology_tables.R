# extract_sweep_topology_tables.R
# Produce dos tablotas CSV con TODAS las simulaciones:
#   1) topology_pre_all_sims.csv
#   2) topology_post_all_sims.csv

library(tidyverse)
library(igraph)

# ---- helpers -----------------------------------------------------------------

largest_component_subgraph <- function(g) {
  comp <- components(g)
  giant_id <- which.max(comp$csize)
  induced_subgraph(g, vids = V(g)[comp$membership == giant_id])
}

safe_mean_distance <- function(g) {
  if (vcount(g) <= 1) return(NA_real_)
  
  g_gc <- largest_component_subgraph(g)
  
  if (vcount(g_gc) <= 1) return(NA_real_)
  
  mean_distance(g_gc, directed = FALSE)
}

graph_metrics <- function(g) {
  
  comp <- components(g)
  gc_size <- max(comp$csize)
  
  tibble(
    n_nodes = vcount(g),
    n_edges = ecount(g),
    mean_degree = mean(degree(g)),
    density = edge_density(g),
    clustering = transitivity(g, "globalundirected"),
    mean_distance_gc = safe_mean_distance(g),
    n_components = comp$no,
    giant_component_size = gc_size,
    giant_component_frac = gc_size / vcount(g)
  )
}

read_one_sim <- function(path) {
  
  x <- readRDS(path)
  
  list(
    pre = graph_metrics(x$g_pre),
    post = graph_metrics(x$g_post)
  )
}

# ---- main --------------------------------------------------------------------

extract_sweep_tables <- function(run_dir) {
  
  index_path <- file.path(run_dir, "index.rds")
  
  index_tbl <- readRDS(index_path)
  
  index_ok <- index_tbl %>%
    mutate(raw_exists = file.exists(raw_path)) %>%
    filter(raw_exists)
  
  meta_cols <- names(index_ok)
  
  pre_tbl <- map_dfr(
    1:nrow(index_ok),
    function(i) {
      
      row <- index_ok[i,]
      
      tabs <- read_one_sim(row$raw_path)
      
      bind_cols(
        row %>% select(all_of(meta_cols)),
        tabs$pre
      )
      
    }
  )
  
  post_tbl <- map_dfr(
    1:nrow(index_ok),
    function(i) {
      
      row <- index_ok[i,]
      
      tabs <- read_one_sim(row$raw_path)
      
      bind_cols(
        row %>% select(all_of(meta_cols)),
        tabs$post
      )
      
    }
  )
  
  write_csv(
    pre_tbl,
    file.path(run_dir, "topology_pre_all_sims.csv")
  )
  
  write_csv(
    post_tbl,
    file.path(run_dir, "topology_post_all_sims.csv")
  )
  
}

# ---- example -----------------------------------------------------------------

extract_sweep_tables("outputs/sweep_er/er_2026-03-04")

extract_sweep_tables("outputs/sweep_ws/ws_2026-03-05/")

extract_sweep_tables("outputs/sweep_ba/ba_2026-03-05/")