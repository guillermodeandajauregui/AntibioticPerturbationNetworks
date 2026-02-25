# =============================================================================
# FILE 2: R/sim_resistance_steps.R
# Step functions: clone creation, spread mechanics, antibiotic hook, kill.
# =============================================================================

add_plus_clone <- function(g, vid) {
  # Create resistant clone "s+" for vertex vid, copying all its neighbors.
  # Idempotent: if "s+" already exists, does nothing.
  
  s_name <- V(g)$name[vid]
  splus  <- paste0(s_name, "+")
  
  if (splus %in% V(g)$name) return(g)
  
  nbr_ids   <- neighbors(g, vid)
  nbr_names <- V(g)$name[nbr_ids]
  
  g <- add_vertices(g, 1, attr = list(name = splus, resistant = TRUE))
  
  if (length(nbr_names) > 0) {
    g <- add_edges(g, rbind(rep(splus, length(nbr_names)), nbr_names))
  }
  
  g
}

get_spread_candidates <- function(g) {
  # Susceptible vertices that are 1-hop neighbors of resistant ones.
  
  res_ids <- which(V(g)$resistant)
  if (length(res_ids) == 0) return(integer(0))
  
  neigh_ids <- unique(unlist(neighborhood(g, order = 1, nodes = res_ids)))
  
  cand_ids <- setdiff(neigh_ids, res_ids)
  if (length(cand_ids) == 0) return(integer(0))
  
  cand_ids <- cand_ids[!V(g)$resistant[cand_ids]]
  if (length(cand_ids) == 0) return(integer(0))
  
  cand_ids
}

sample_acquisitions <- function(cand_ids, q) {
  # Bernoulli trial per candidate: acquire if runif(1) < q
  if (length(cand_ids) == 0) return(integer(0))
  cand_ids[runif(length(cand_ids)) < q]
}

apply_acquisitions <- function(g, acquired) {
  # For each acquired susceptible, create its "s+" clone.
  if (length(acquired) == 0) return(g)
  
  for (vid in acquired) {
    g <- add_plus_clone(g, vid)
  }
  
  g
}

step_spread <- function(g, q) {
  # One spread step: candidates -> sample -> apply
  cand_ids <- get_spread_candidates(g)
  acquired <- sample_acquisitions(cand_ids, q)
  apply_acquisitions(g, acquired)
}

step_antibiotic_event <- function(t, AB, verbose) {
  # Informational hook only (no mutation).
  if (verbose && t == AB) message("AB introduced at t = ", AB)
  invisible(NULL)
}

step_kill <- function(g) {
  # Remove all susceptible nodes.
  sus_ids <- which(!V(g)$resistant)
  if (length(sus_ids) == 0) return(g)
  delete_vertices(g, sus_ids)
}