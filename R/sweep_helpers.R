# R/sweep_helpers.R
# Helpers for parameter sweeps (file paths, ids, seeds, grids)

library(tidyverse)

# ---- filesystem --------------------------------------------------------------
assert_dir_exists <- function(path) {
  if (!dir.exists(path)) stop("Directory does not exist: ", path)
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
# Graph seeds are global 1..N_graph_total (you requested this explicitly).
# Sim seeds are deterministic from (graph_seed_id, cond_id, sim_rep).
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
    p       = p_values,
    AB      = AB_values
  ) %>%
    mutate(cond_id = row_number())
}

n_graph_total <- function(n_values, p_edge_values, n_graph) {
  length(n_values) * length(p_edge_values) * n_graph
}