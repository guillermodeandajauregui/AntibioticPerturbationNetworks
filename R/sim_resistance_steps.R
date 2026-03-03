# =============================================================================
# FILE 2: R/sim_resistance_steps.R
# Step functions: clone creation, spread mechanics, antibiotic hook, kill.
# =============================================================================

add_plus_clone <- function(g, vid) {
  # Create a resistant clone of vertex `vid`, copying all its neighbors.
  # The clone inherits the same `name` and is distinguished by `resistant = TRUE`.
  
  s_name <- V(g)$name[vid]
  nbr_ids <- neighbors(g, vid)
  
  # Add new vertex (clone)
  g <- add_vertices(g, 1, attr = list(name = s_name, resistant = TRUE))
  clone_id <- vcount(g)  # the newly added vertex index
  
  # Copy edges by vertex indices (NOT by names)
  if (length(nbr_ids) > 0) {
    edge_vec <- as.vector(rbind(rep(clone_id, length(nbr_ids)), as.integer(nbr_ids)))
    g <- add_edges(g, edge_vec)
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