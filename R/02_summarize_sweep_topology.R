# summarize_sweep_topology.R
# Resume un directorio de output de sweep y devuelve tablas agregadas
# separadas para PRE y POST, respetando la topología (ER / WS / BA).

# ---- setup -------------------------------------------------------------------
library(tidyverse)
library(igraph)

# ---- helpers -----------------------------------------------------------------
detect_network_type <- function(index_tbl) {
  if ("p_edge" %in% names(index_tbl)) {
    return("ER")
  }
  
  if (all(c("k", "beta") %in% names(index_tbl))) {
    return("WS")
  }
  
  if ("m" %in% names(index_tbl)) {
    return("BA")
  }
  
  stop("No pude detectar el tipo de red a partir de las columnas del index.")
}

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
  
  # Para evitar Inf en grafos desconectados, medimos la distancia media
  # sobre el componente gigante.
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

read_one_sim_metrics <- function(raw_path) {
  sim <- readRDS(raw_path)
  
  # Asume que cada archivo tiene g_pre y g_post, como sale de compact_sim_output()
  pre <- graph_metrics(sim$g_pre) %>%
    rename_with(~ paste0(.x, "_pre"))
  
  post <- graph_metrics(sim$g_post) %>%
    rename_with(~ paste0(.x, "_post"))
  
  bind_cols(pre, post)
}

get_group_cols <- function(index_tbl) {
  base_cols <- c("n", "k0_frac", "k0", "p", "AB", "tx")
  
  topo_cols <- case_when(
    "p_edge" %in% names(index_tbl) ~ list(c("p_edge")),
    all(c("k", "beta") %in% names(index_tbl)) ~ list(c("k", "beta")),
    "m" %in% names(index_tbl) ~ list(c("m")),
    TRUE ~ list(character(0))
  )
  
  c(base_cols, topo_cols[[1]]) %>%
    keep(~ .x %in% names(index_tbl))
}

# ---- main worker -------------------------------------------------------------
summarize_sweep_dir <- function(run_dir,
                                save_csv = TRUE,
                                out_dir = run_dir) {
  
  index_path <- file.path(run_dir, "index.rds")
  
  if (!file.exists(index_path)) {
    stop("No encontré index.rds en: ", run_dir)
  }
  
  index_tbl <- readRDS(index_path)
  
  if (!"raw_path" %in% names(index_tbl)) {
    stop("El index.rds no tiene columna raw_path.")
  }
  
  index_ok <- index_tbl %>%
    mutate(
      status = if ("status" %in% names(.)) status else "ok",
      raw_exists = file.exists(raw_path)
    ) %>%
    filter(status == "ok", raw_exists)
  
  if (nrow(index_ok) == 0) {
    stop("No hay simulaciones válidas con status == 'ok' y raw_path existente.")
  }
  
  network_type <- detect_network_type(index_ok)
  group_cols <- get_group_cols(index_ok)
  
  meta_cols <- c(
    group_cols,
    intersect(
      c("sim_id", "graph_id", "graph_rep", "graph_seed_id", "seed_sim"),
      names(index_ok)
    )
  )
  
  per_sim_metrics <- index_ok %>%
    mutate(.row_id = row_number()) %>%
    split(.$.row_id) %>%
    purrr::map_dfr(
      function(x) {
        met <- read_one_sim_metrics(x$raw_path[[1]])
        bind_cols(
          x %>%
            select(all_of(meta_cols)),
          met
        )
      }
    )
  
  pre_tbl <- per_sim_metrics %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_sims = n(),
      across(
        ends_with("_pre"),
        ~ mean(.x, na.rm = TRUE),
        .names = "{.col}"
      ),
      .groups = "drop"
    ) %>%
    arrange(across(all_of(group_cols)))
  
  post_tbl <- per_sim_metrics %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_sims = n(),
      across(
        ends_with("_post"),
        ~ mean(.x, na.rm = TRUE),
        .names = "{.col}"
      ),
      .groups = "drop"
    ) %>%
    arrange(across(all_of(group_cols)))
  
  if (isTRUE(save_csv)) {
    write_csv(
      pre_tbl,
      file.path(out_dir, paste0("table_", tolower(network_type), "_pre.csv"))
    )
    
    write_csv(
      post_tbl,
      file.path(out_dir, paste0("table_", tolower(network_type), "_post.csv"))
    )
  }
  
  list(
    network_type = network_type,
    grouping_vars = group_cols,
    per_sim = per_sim_metrics,
    pre = pre_tbl,
    post = post_tbl
  )
}

# ---- optional wrapper: varios runs -------------------------------------------
summarize_many_sweeps <- function(run_dirs,
                                  save_csv = TRUE) {
  purrr::map(run_dirs, ~ summarize_sweep_dir(.x, save_csv = save_csv))
}

# ---- example -----------------------------------------------------------------
# ER
er_res <-
  summarize_sweep_dir("outputs/sweep_er/er_2026-03-04")
# 
# er_res$pre
# er_res$post
#
# WS
ws_res <-
  summarize_sweep_dir("outputs/sweep_ws/ws_2026-03-05/")
#
# BA
ba_res <-
  summarize_sweep_dir("outputs/sweep_ba/ba_2026-03-05/")
#
# Varios de una:
# all_res <-
#   summarize_many_sweeps(
#     c(
#       "outputs/sweep_er/er_2026-03-04",
#       "outputs/sweep_ws/ws_2026-03-04",
#       "outputs/sweep_ba/ba_2026-03-04"
#     )
#   )