# =============================================================================
# FILE 3: R/simulate_resistance_network.R
# Orchestrator: main time loop, calls step functions, returns snapshots.
# =============================================================================

simulate_resistance_network <- function(
    g,
    p = 0.1,
    T = 50,
    AB = 20,
    tx = 5,
    verbose = FALSE
) {
  # (0) Normalize/validate vertex attributes
  g <- ensure_vertex_attrs(g)
  
  # (0b) Derived parameters
  dp <- derive_params(p = p, AB = AB, tx = tx)
  q <- dp$q
  t_kill <- dp$t_kill
  
  # (0c) Allocate snapshot container
  graphs <- init_snapshots(T)
  
  # (1) Main time loop: apply steps, then snapshot
  for (t in seq_len(T)) {
    
    # (1A) Spread phase (only before kill time)
    if (t < t_kill) {
      g <- step_spread(g, q)
    }
    
    # (1B) Antibiotic introduction hook (informational)
    step_antibiotic_event(t, AB, verbose)
    
    # (1C) Kill phase (exactly at kill time)
    if (t == t_kill) {
      g <- step_kill(g)
      if (verbose) message("Susceptible species removed at t = ", t_kill)
    }
    
    # (1D) Snapshot at end of time step
    graphs[[t]] <- g
  }
  
  graphs
}