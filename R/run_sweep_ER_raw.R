# 01_run_sweep_ER_raw.R
# ER-only sweep; save compact sim outputs per sim; write an index table.
#
# Saves per sim:
# - g_pre (pre-antibiotic snapshot; first graph in history)
# - g_post (post-antibiotic snapshot; last graph in history)
# - counts_ts (n resistant / non-resistant over time)
#
# Does NOT save base graph separately.
# Graph generation is reproducible with global seed id 1..N_graph_total.

# ---- setup -------------------------------------------------------------------
library(tidyverse)
library(igraph)

source("R/sweep_helpers.R")

source("R/sim_resistance_helpers.R")
source("R/graph_prep.R")
source("R/sim_resistance_steps.R")
source("R/simulate_resistance_network.R")

output_dir <-
  "outputs/sweep_er"  # must exist

run_id <-
  paste0("er_", Sys.Date())

run_dir <-
  file.path(output_dir, run_id)

dir_raw <-
  file.path(run_dir, "raw")

assert_dir_exists(output_dir)
ensure_dir(run_dir)
ensure_dir(dir_raw)

# ---- parameters --------------------------------------------------------------
n_values <-
  c(100, 300, 1000)

p_edge_values <-
  c(0.02, 0.05, 0.10)

k0_frac_values <-
  c(0.05, 0.15)

# This is the Poisson rate parameter "p" used by derive_params()
p_values <-
  c(0.02, 0.08)

# Antibiotic administered once, perfectly effective (tx = 0 => kill at t = AB)
AB_values <-
  c(20, 60)

tx <-
  0L

n_graph <-
  10

n_sim <-
  10

cond_grid <-
  make_cond_grid(
    k0_frac_values = k0_frac_values,
    p_values       = p_values,
    AB_values      = AB_values
  )

N_graph_total <-
  n_graph_total(n_values, p_edge_values, n_graph)

# ---- run ---------------------------------------------------------------------
index_rows <- list()
row_i <- 0L
graph_seed_id <- 0L

for (n in n_values) {
  for (p_edge in p_edge_values) {
    for (graph_rep in seq_len(n_graph)) {
      
      graph_seed_id <- graph_seed_id + 1L
      if (graph_seed_id > N_graph_total) stop("graph_seed_id exceeded N_graph_total")
      
      graph_id <-
        make_graph_id(n, p_edge, graph_seed_id)
      
      # ER graph reproducible with global seed id 1..N
      set.seed(graph_seed_id)
      g0 <-
        sample_gnp(
          n        = n,
          p        = p_edge,
          directed = FALSE,
          loops    = FALSE
        ) %>%
        decorate_resistance_graph(name_prefix = "s")
      
      for (cond_i in seq_len(nrow(cond_grid))) {
        
        k0_frac <- cond_grid$k0_frac[[cond_i]]
        p       <- cond_grid$p[[cond_i]]
        AB      <- cond_grid$AB[[cond_i]]
        cond_id <- cond_grid$cond_id[[cond_i]]
        
        # kill at t_kill = AB + tx; keep 1 step after for post snapshot
        T_used <- as.integer(AB + tx + 1L)
        
        k0 <- k_from_frac(n, k0_frac)
        
        for (sim_rep in seq_len(n_sim)) {
          
          sim_id <-
            make_sim_id(graph_id, cond_id, sim_rep)
          
          sim_seed <-
            seed_for_sim(graph_seed_id, cond_id, sim_rep)
          
          raw_path <-
            file.path(dir_raw, paste0("sim_", sim_id, ".rds"))
          
          status <- "ok"
          error_msg <- NA_character_
          
          sim_out <-
            tryCatch(
              {
                g_init <-
                  seed_resistance(g0, k = k0, seed = sim_seed)
                
                graphs <-
                  simulate_resistance_network(
                    g       = g_init,
                    p       = p,
                    T       = T_used,
                    AB      = AB,
                    tx      = tx,
                    verbose = FALSE
                  )
                
                compact_sim_output(
                  graphs = graphs,
                  meta = list(
                    run_id        = run_id,
                    sim_id        = sim_id,
                    graph_id      = graph_id,
                    n             = n,
                    p_edge        = p_edge,
                    graph_rep     = graph_rep,
                    graph_seed_id = graph_seed_id,
                    cond_id       = cond_id,
                    k0_frac       = k0_frac,
                    k0            = k0,
                    p             = p,
                    AB            = AB,
                    tx            = tx,
                    T_used        = T_used,
                    seed_sim      = sim_seed
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
              run_id        = run_id,
              sim_id        = sim_id,
              graph_id      = graph_id,
              n             = n,
              p_edge        = p_edge,
              graph_rep     = graph_rep,
              graph_seed_id = graph_seed_id,
              cond_id       = cond_id,
              k0_frac       = k0_frac,
              k0            = k0,
              p             = p,
              AB            = AB,
              tx            = tx,
              T_used        = T_used,
              seed_sim      = sim_seed,
              raw_path      = raw_path,
              status        = status,
              error_msg     = error_msg
            )
        }
      }
    }
  }
}

index_tbl <- bind_rows(index_rows)

index_path <-
  file.path(run_dir, "index.rds")

saveRDS(index_tbl, index_path)

message("Done.")
message("Total graphs generated: ", graph_seed_id, " (expected ", N_graph_total, ")")
message("Index saved at: ", index_path)