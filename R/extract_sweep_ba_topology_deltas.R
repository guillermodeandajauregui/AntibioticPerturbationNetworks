library(tidyverse)
library(igraph)

# ---- paths -------------------------------------------------------------------
run_dir <-
  "outputs/sweep_ba/ba_2026-03-05/"   # ajusta si hace falta

index_path <-
  file.path(run_dir, "index.rds")

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

# ---- extract pre/post --------------------------------------------------------
extract_one_sim_topology <-
  function(raw_path) {
    x <-
      readRDS(raw_path)
    
    list(
      pre  = compute_topology_metrics(x$g_pre),
      post = compute_topology_metrics(x$g_post)
    )
  }

topo_list <-
  purrr::map(index_tbl$raw_path, extract_one_sim_topology)

topology_pre_ba <-
  purrr::map_dfr(
    seq_along(topo_list),
    function(i) {
      bind_cols(
        index_tbl[i, ],
        tibble(snapshot = "pre"),
        topo_list[[i]]$pre
      )
    }
  )

topology_post_ba <-
  purrr::map_dfr(
    seq_along(topo_list),
    function(i) {
      bind_cols(
        index_tbl[i, ],
        tibble(snapshot = "post"),
        topo_list[[i]]$post
      )
    }
  )

# ---- differential table ------------------------------------------------------
id_vars <-
  c(
    "run_id", "sim_id", "graph_id",
    "n", "m", "graph_rep", "graph_seed_id",
    "cond_id", "k0_frac", "k0",
    "p", "AB", "tx", "T_used",
    "seed_sim", "raw_path", "status", "error_msg"
  )

metric_vars <-
  c(
    "n_nodes", "n_edges", "density", "clustering",
    "n_components", "giant_component_size",
    "giant_component_frac", "mean_distance_gc"
  )

topology_delta_ba <-
  topology_pre_ba |>
  select(all_of(id_vars), all_of(metric_vars)) |>
  rename_with(~ paste0(.x, "_pre"), all_of(metric_vars)) |>
  left_join(
    topology_post_ba |>
      select(all_of(id_vars), all_of(metric_vars)) |>
      rename_with(~ paste0(.x, "_post"), all_of(metric_vars)),
    by = id_vars
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

# objetos finales en memoria:
# topology_pre_ba
# topology_post_ba
# topology_delta_ba |> 
#   vroom::vroom_write(file = paste0(run_dir, "topology_delta_ba.txt"))
