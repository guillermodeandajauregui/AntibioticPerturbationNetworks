# extract_sweep_real_tables.R
# Produce dos tablotas CSV con TODAS las simulaciones:
#   1) topology_pre_all_sims.csv
#   2) topology_post_all_sims.csv

library(tidyverse)
library(igraph)

# ---- helpers -----------------------------------------------------------------

largest_component_subgraph <- function(g) {
  comp <- components(g)
  
  if (length(comp$csize) == 0) {
    return(g)
  }
  
  giant_id <- which.max(comp$csize)
  
  induced_subgraph(g, vids = V(g)[comp$membership == giant_id])
}

safe_mean_distance <- function(g) {
  if (vcount(g) <= 1) {
    return(NA_real_)
  }
  
  g_gc <- largest_component_subgraph(g)
  
  if (vcount(g_gc) <= 1) {
    return(NA_real_)
  }
  
  mean_distance(g_gc, directed = FALSE, unconnected = FALSE)
}

graph_metrics <- function(g) {
  comp <- components(g)
  gc_size <- if (length(comp$csize) == 0) 0 else max(comp$csize)
  
  tibble(
    n_nodes = vcount(g),
    n_edges = ecount(g),
    mean_degree = if (vcount(g) > 0) mean(degree(g)) else NA_real_,
    density = if (vcount(g) > 1) edge_density(g, loops = FALSE) else NA_real_,
    clustering = transitivity(g, type = "globalundirected"),
    mean_distance_gc = safe_mean_distance(g),
    n_components = comp$no,
    giant_component_size = gc_size,
    giant_component_frac = if (vcount(g) > 0) gc_size / vcount(g) else NA_real_
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

extract_sweep_real_tables <- function(run_dir) {
  
  index_path <-
    file.path(run_dir, "index.rds")
  
  if (!file.exists(index_path)) {
    stop("No encontré index.rds en: ", run_dir)
  }
  
  index_tbl <-
    readRDS(index_path)
  
  if (!"raw_path" %in% names(index_tbl)) {
    stop("El index.rds no tiene columna raw_path.")
  }
  
  index_ok <-
    index_tbl %>%
    mutate(
      status = if ("status" %in% names(.)) status else "ok",
      raw_exists = file.exists(raw_path)
    ) %>%
    filter(status == "ok", raw_exists)
  
  if (nrow(index_ok) == 0) {
    stop("No hay simulaciones válidas con status == 'ok' y raw_path existente.")
  }
  
  meta_cols <-
    names(index_ok)
  
  pre_tbl <-
    purrr::map_dfr(
      seq_len(nrow(index_ok)),
      function(i) {
        row <-
          index_ok[i, ]
        
        tabs <-
          read_one_sim(row$raw_path[[1]])
        
        bind_cols(
          row %>% select(all_of(meta_cols)),
          tabs$pre
        )
      }
    )
  
  post_tbl <-
    purrr::map_dfr(
      seq_len(nrow(index_ok)),
      function(i) {
        row <-
          index_ok[i, ]
        
        tabs <-
          read_one_sim(row$raw_path[[1]])
        
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
  
  invisible(
    list(
      pre = pre_tbl,
      post = post_tbl
    )
  )
}

# ---- example -----------------------------------------------------------------

extract_sweep_real_tables("outputs/sweep_real/real_2026-03-10")