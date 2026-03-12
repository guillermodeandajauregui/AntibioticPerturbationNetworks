# extract_sweep_ws_topology_delta.R
# Build only the pre/post differential topology table for WS sweep.

# ---- setup -------------------------------------------------------------------
library(tidyverse)
library(igraph)

run_dir <-
  "outputs/sweep_ws/ws_2026-03-05/"   # ajusta si hace falta

index_path <-
  file.path(run_dir, "index.rds")

delta_out_file <-
  file.path(run_dir, "topology_delta_all_sims.csv")

# ---- topology metrics --------------------------------------------------------
compute_topology_metrics <-
  function(g) {
    comps <-
      components(g)
    
    gc_id <-
      which.max(comps$csize)
    
    gc_vertices <-
      V(g)[comps$membership == gc_id]
    
    g_gc <-
      induced_subgraph(g, gc_vertices)
    
    tibble(
      n_nodes = vcount(g),
      n_edges = ecount(g),
      density = edge_density(g, loops = FALSE),
      clustering = transitivity(g, type = "globalundirected"),
      n_components = comps$no,
      giant_component_size = max(comps$csize),
      giant_component_frac = max(comps$csize) / vcount(g),
      mean_distance_gc = if (vcount(g_gc) > 1) {
        mean_distance(g_gc, directed = FALSE, unconnected = FALSE)
      } else {
        NA_real_
      }
    )
  }

# ---- load index --------------------------------------------------------------
index_tbl <-
  readRDS(index_path) |>
  filter(status == "ok")

# ---- extract per sim ---------------------------------------------------------
extract_one_sim_delta <-
  function(raw_path) {
    x <-
      readRDS(raw_path)
    
    pre <-
      compute_topology_metrics(x$g_pre)
    
    post <-
      compute_topology_metrics(x$g_post)
    
    bind_cols(
      pre |>
        rename_with(~ paste0(.x, "_pre")),
      post |>
        rename_with(~ paste0(.x, "_post"))
    ) |>
      mutate(
        delta_n_nodes = n_nodes_post - n_nodes_pre,
        delta_n_edges = n_edges_post - n_edges_pre,
        delta_density = density_post - density_pre,
        delta_clustering = clustering_post - clustering_pre,
        delta_n_components = n_components_post - n_components_pre,
        delta_giant_component_size = giant_component_size_post - giant_component_size_pre,
        delta_giant_component_frac = giant_component_frac_post - giant_component_frac_pre,
        delta_mean_distance_gc = mean_distance_gc_post - mean_distance_gc_pre
      )
  }

# ---- columns to drop from index ----------------------------------------------
inutiles <-
  c(
    "run_id", "sim_id", "graph_seed_id", "cond_id",
    "seed_sim", "graph_rep", "graph_id", "raw_path",
    "error_msg", "status"
  )

index_keep <-
  index_tbl |>
  select(-any_of(inutiles))

# ---- build delta table -------------------------------------------------------
topology_delta_all_sims <-
  purrr::map_dfr(
    seq_len(nrow(index_tbl)),
    function(i) {
      bind_cols(
        index_keep[i, ],
        extract_one_sim_delta(index_tbl$raw_path[[i]])
      )
    }
  )

# ---- write output ------------------------------------------------------------
readr::write_csv(topology_delta_all_sims, delta_out_file)