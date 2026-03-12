# run_sweep_real_raw.R
# Empirical-graph sweep; save compact sim outputs per sim; write an index table.

# ---- setup -------------------------------------------------------------------
library(tidyverse)
library(igraph)

source("R/sweep_helpers.R")

source("R/sim_resistance_helpers.R")
source("R/graph_prep.R")
source("R/sim_resistance_steps.R")
source("R/simulate_resistance_network.R")

output_dir <-
  "outputs/sweep_real"  # must exist

run_id <-
  paste0("real_", Sys.Date())

run_dir <-
  file.path(output_dir, run_id)

dir_raw <-
  file.path(run_dir, "raw")

assert_dir_exists(output_dir)
ensure_dir(run_dir)
ensure_dir(dir_raw)

# ---- input graph -------------------------------------------------------------
graph_path <-
  "data/enterotype_graph.graphml"

g0 <-
  read_graph(graph_path, format = "graphml") |>
  as.undirected(mode = "collapse") |>
  simplify(remove.multiple = TRUE, remove.loops = TRUE) |>
  decorate_resistance_graph(name_prefix = "s")

n <-
  vcount(g0)

graph_seed_id <-
  1L

# ---- parameters --------------------------------------------------------------
k0_values <-
  c(1L, 5L, 10L)

seed_mode_values <-
  c("random", "high_degree", "dispersed")

p_values <-
  c(0.001, 0.005, 0.01, 0.05, 0.1)

AB_values <-
  20L

tx_values <-
  0L

T_used_values <-
  as.integer(AB_values + tx_values + 1L)

n_sim <-
  10L

cond_grid <-
  make_cond_grid_full(
    k0_values = k0_values,
    seed_mode_values = seed_mode_values,
    p_values = p_values,
    AB_values = AB_values,
    tx_values = tx_values,
    T_used_values = T_used_values
  )

# ---- run ---------------------------------------------------------------------
index_rows <- list()
row_i <- 0L

for (cond_i in seq_len(nrow(cond_grid))) {
  
  k0 <-
    cond_grid$k0[[cond_i]]
  
  seed_mode <-
    cond_grid$seed_mode[[cond_i]]
  
  p <-
    cond_grid$p[[cond_i]]
  
  AB <-
    cond_grid$AB[[cond_i]]
  
  tx <-
    cond_grid$tx[[cond_i]]
  
  T_used <-
    cond_grid$T_used[[cond_i]]
  
  cond_id <-
    cond_grid$cond_id[[cond_i]]
  
  for (sim_rep in seq_len(n_sim)) {
    
    sim_id <-
      paste0("c", cond_id, "_r", sim_rep)
    
    sim_seed <-
      seed_for_sim(graph_seed_id, cond_id, sim_rep)
    
    raw_path <-
      file.path(dir_raw, paste0("sim_", sim_id, ".rds"))
    
    status <-
      "ok"
    
    error_msg <-
      NA_character_
    
    sim_out <-
      tryCatch(
        {
          g_init <-
            seed_resistance(
              g = g0,
              k = k0,
              seed = sim_seed,
              seed_mode = seed_mode
            )
          
          graphs <-
            simulate_resistance_network(
              g = g_init,
              p = p,
              T = T_used,
              AB = AB,
              tx = tx,
              verbose = FALSE
            )
          
          compact_sim_output(
            graphs = graphs,
            meta = list(
              run_id = run_id,
              sim_id = sim_id,
              n = n,
              graph_seed_id = graph_seed_id,
              cond_id = cond_id,
              k0 = k0,
              seed_mode = seed_mode,
              p = p,
              AB = AB,
              tx = tx,
              T_used = T_used,
              seed_sim = sim_seed
            )
          )
        },
        error = function(e) {
          status <<- "error"
          error_msg <<- conditionMessage(e)
          NULL
        }
      )
    
    if (identical(status, "ok")) {
      saveRDS(sim_out, raw_path)
    }
    
    row_i <- row_i + 1L
    
    index_rows[[row_i]] <-
      tibble(
        run_id = run_id,
        sim_id = sim_id,
        n = n,
        graph_seed_id = graph_seed_id,
        cond_id = cond_id,
        k0 = k0,
        seed_mode = seed_mode,
        p = p,
        AB = AB,
        tx = tx,
        T_used = T_used,
        seed_sim = sim_seed,
        raw_path = raw_path,
        status = status,
        error_msg = error_msg
      )
  }
}

index_tbl <-
  bind_rows(index_rows)

index_path <-
  file.path(run_dir, "index.rds")

saveRDS(index_tbl, index_path)

message("Done.")
message("Nodes: ", n)
message("Conditions: ", nrow(cond_grid))
message("Simulations: ", nrow(cond_grid) * n_sim)
message("Index saved at: ", index_path)