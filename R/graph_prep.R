# =============================================================================
# FILE: R/graph_prep.R
# Graph preparation: decorate graphs + seed initial resistance states.
# =============================================================================

decorate_resistance_graph <- function(g, name_prefix = "s") {
  # Ensure the graph has the required vertex attributes for the resistance model.
  #
  # Guarantees:
  #   - V(g)$name is a unique character vector
  #   - V(g)$resistant is a logical vector (all FALSE by default if missing)
  #
  # This function does NOT change topology; it only sets/validates attributes.
  
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
    # be strict: coerce to logical and validate length
    V(g)$resistant <- as.logical(V(g)$resistant)
    if (length(V(g)$resistant) != vcount(g)) stop("V(g)$resistant length mismatch.")
    if (anyNA(V(g)$resistant)) stop("V(g)$resistant contains NA.")
  }
  
  g
}

seed_resistance <- function(g, k = 3, seed = NULL, replace = FALSE) {
  # Seed initial resistant nodes by sampling k vertices uniformly at random.
  #
  # Behavior:
  #   - Sets exactly (up to) k nodes to resistant == TRUE.
  #   - If k > vcount(g), it seeds all nodes.
  #   - If replace = FALSE (default), it samples without replacement.
  #
  # Notes:
  #   - This function assumes the graph already has V(g)$resistant.
  #     If not, call decorate_resistance_graph(g) first.
  
  stopifnot(is_igraph(g))
  
  if (is.null(V(g)$resistant)) {
    stop("V(g)$resistant is missing. Call decorate_resistance_graph(g) first.")
  }
  
  n <- vcount(g)
  if (n == 0) return(g)
  
  k <- as.integer(k)
  if (is.na(k) || k < 0) stop("k must be a non-negative integer.")
  
  if (!is.null(seed)) set.seed(seed)
  
  size <- min(k, n)
  
  # Sample vertex ids (1..n). If size == 0, nothing to do.
  if (size == 0) return(g)
  
  idx <- sample.int(n = n, size = size, replace = replace)
  
  V(g)$resistant[idx] <- TRUE
  g
}