# =============================================================================
# FILE: R/graph_prep.R
# Graph preparation: decorate graphs + seed initial resistance states.
# =============================================================================

decorate_resistance_graph <- function(g, name_prefix = "s") {
  stopifnot(is_igraph(g))
  
  # ---- name ---------------------------------------------------------------
  if (is.null(V(g)$name)) {
    V(g)$name <- paste0(name_prefix, seq_len(vcount(g)))
  } else {
    V(g)$name <- as.character(V(g)$name)
  }
  
  if (anyNA(V(g)$name)) stop("V(g)$name contains NA.")
  if (anyDuplicated(V(g)$name)) stop("V(g)$name must be unique.")
  
  # ---- resistant ----------------------------------------------------------
  if (is.null(V(g)$resistant)) {
    V(g)$resistant <- rep(FALSE, vcount(g))
  } else {
    V(g)$resistant <- as.logical(V(g)$resistant)
    if (length(V(g)$resistant) != vcount(g)) stop("V(g)$resistant length mismatch.")
    if (anyNA(V(g)$resistant)) stop("V(g)$resistant contains NA.")
  }
  
  # ---- cloned -------------------------------------------------------------
  if (is.null(V(g)$cloned)) {
    V(g)$cloned <- rep(FALSE, vcount(g))
  } else {
    V(g)$cloned <- as.logical(V(g)$cloned)
    if (length(V(g)$cloned) != vcount(g)) stop("V(g)$cloned length mismatch.")
    if (anyNA(V(g)$cloned)) stop("V(g)$cloned contains NA.")
  }
  
  g
}

pick_seed_nodes_random <- function(g, k0) {
  stopifnot(is_igraph(g))
  
  n <- vcount(g)
  if (n == 0L) return(integer(0))
  
  k0 <- as.integer(k0)
  if (is.na(k0) || k0 < 0L) stop("k0 must be a non-negative integer.")
  
  size <- min(k0, n)
  if (size == 0L) return(integer(0))
  
  sample.int(n = n, size = size, replace = FALSE)
}

pick_seed_nodes_high_degree <- function(g, k0) {
  stopifnot(is_igraph(g))
  
  n <- vcount(g)
  if (n == 0L) return(integer(0))
  
  k0 <- as.integer(k0)
  if (is.na(k0) || k0 < 0L) stop("k0 must be a non-negative integer.")
  
  size <- min(k0, n)
  if (size == 0L) return(integer(0))
  
  deg <- degree(g)
  
  ord <- order(deg, decreasing = TRUE)
  ord[seq_len(size)]
}

pick_seed_nodes_dispersed <- function(g, k0) {
  stopifnot(is_igraph(g))
  
  n <- vcount(g)
  if (n == 0L) return(integer(0))
  
  k0 <- as.integer(k0)
  if (is.na(k0) || k0 < 0L) stop("k0 must be a non-negative integer.")
  
  size <- min(k0, n)
  if (size == 0L) return(integer(0))
  if (size == 1L) return(sample.int(n, 1))
  
  dmat <- distances(g, v = V(g), to = V(g), mode = "all")
  deg <- degree(g)
  
  # primer nodo: hub
  selected <- integer(size)
  selected[1] <- which.max(deg)
  
  remaining <- setdiff(seq_len(n), selected[1])
  
  for (i in 2:size) {
    min_dist_to_selected <-
      apply(dmat[remaining, selected[seq_len(i - 1)], drop = FALSE], 1, min)
    
    # prioriza nodos más lejos; desempata por mayor grado
    best_dist <- max(min_dist_to_selected)
    candidates <- remaining[min_dist_to_selected == best_dist]
    
    if (length(candidates) > 1L) {
      cand_deg <- deg[candidates]
      candidates <- candidates[cand_deg == max(cand_deg)]
    }
    
    if (length(candidates) > 1L) {
      next_node <- sample(candidates, 1)
    } else {
      next_node <- candidates
    }
    
    selected[i] <- next_node
    remaining <- setdiff(remaining, next_node)
  }
  
  selected
}

seed_resistance <- function(g,
                            k = 3,
                            seed = NULL,
                            replace = FALSE,
                            seed_mode = "random") {
  stopifnot(is_igraph(g))
  
  if (is.null(V(g)$resistant)) {
    stop("V(g)$resistant is missing. Call decorate_resistance_graph(g) first.")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  seed_mode <- match.arg(
    seed_mode,
    choices = c("random", "high_degree", "dispersed")
  )
  
  if (isTRUE(replace) && seed_mode != "random") {
    warning("replace is only used for seed_mode = 'random'. Ignoring replace.")
  }
  
  n <- vcount(g)
  if (n == 0L) return(g)
  
  k <- as.integer(k)
  if (is.na(k) || k < 0L) stop("k must be a non-negative integer.")
  
  idx <-
    switch(
      seed_mode,
      random = {
        size <- min(k, n)
        if (size == 0L) integer(0) else sample.int(n = n, size = size, replace = replace)
      },
      high_degree = pick_seed_nodes_high_degree(g, k0 = k),
      dispersed = pick_seed_nodes_dispersed(g, k0 = k)
    )
  
  if (length(idx) > 0L) {
    V(g)$resistant[idx] <- TRUE
  }
  
  g
}