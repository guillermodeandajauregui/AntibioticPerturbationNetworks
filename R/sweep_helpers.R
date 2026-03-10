# R/sweep_helpers.R
# Helpers for parameter sweeps: dirs, ids, seeds, grids, and compact output extraction.

library(tidyverse)
library(igraph)

# ---- filesystem --------------------------------------------------------------
assert_dir_exists <- function(path) {
  if (!dir.exists(path)) stop("Directory does not exist: ", path)
  invisible(TRUE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = FALSE)
  invisible(TRUE)
}

# ---- ids ---------------------------------------------------------------------
make_graph_id <- function(n, p_edge, graph_seed_id) {
  sprintf("n%04d_p%0.3f_g%04d", n, p_edge, graph_seed_id)
}

make_sim_id <- function(graph_id, cond_id, sim_rep) {
  sprintf("%s_c%02d_s%02d", graph_id, cond_id, sim_rep)
}

# ---- seeds -------------------------------------------------------------------
# Graph seed is the global integer id 1..N_graph_total (requested).
# Sim seed is deterministic from (graph_seed_id, cond_id, sim_rep).
seed_for_sim <- function(graph_seed_id, cond_id, sim_rep) {
  as.integer(graph_seed_id) * 100000L + as.integer(cond_id) * 100L + as.integer(sim_rep)
}

# ---- parameter transforms -----------------------------------------------------
k_from_frac <- function(n, k0_frac) {
  max(0L, min(as.integer(round(n * k0_frac)), as.integer(n)))
}

# ---- grids -------------------------------------------------------------------
make_cond_grid <- function(k0_frac_values, p_values, AB_values) {
  tidyr::crossing(
    k0_frac = k0_frac_values,
    p       = p_values,     # NOTE: this is the Poisson rate parameter used by your engine
    AB      = AB_values
  ) %>%
    mutate(cond_id = row_number())
}

n_graph_total <- function(n_values, p_edge_values, n_graph) {
  length(n_values) * length(p_edge_values) * n_graph
}

# ---- compact extraction -------------------------------------------------------
counts_from_graphs <- function(graphs) {
  tibble(t = names(graphs)) %>%
    mutate(
      t = readr::parse_number(t),
      counts = map(graphs, function(g) {
        r <- igraph::vertex_attr(g, "resistant")
        r <- as.logical(r)
        
        tibble(
          resistant = c(FALSE, TRUE),
          n = c(sum(!r, na.rm = TRUE), sum(r, na.rm = TRUE))
        )
      })
    ) %>%
    unnest(counts)
}

compact_sim_output <- function(graphs, meta = list()) {
  stopifnot(is.list(graphs), length(graphs) >= 1)
  
  list(
    meta      = meta,
    g_pre     = graphs[[1]],
    g_post    = graphs[[length(graphs)]],
    counts_ts = counts_from_graphs(graphs)
  )
}

##########################################################

make_cond_grid_full <- function(k0_values,
                                seed_mode_values,
                                p_values,
                                AB_values,
                                tx_values,
                                T_used_values) {
  tidyr::crossing(
    k0 = k0_values,
    seed_mode = seed_mode_values,
    p = p_values,
    AB = AB_values,
    tx = tx_values,
    T_used = T_used_values
  ) |>
    dplyr::mutate(
      cond_id = dplyr::row_number()
    ) |>
    dplyr::select(
      cond_id,
      k0,
      seed_mode,
      p,
      AB,
      tx,
      T_used
    )
}