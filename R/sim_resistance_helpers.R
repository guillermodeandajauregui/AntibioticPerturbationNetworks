# =============================================================================
# FILE 1: R/sim_resistance_helpers.R
# Small utilities: input normalization, derived params, snapshot allocation.
# =============================================================================

ensure_vertex_attrs <- function(g) {
  stopifnot(is_igraph(g))
  
  if (is.null(V(g)$name)) {
    V(g)$name <- as.character(seq_len(vcount(g)))
  }
  
  if (is.null(V(g)$resistant)) {
    V(g)$resistant <- rep(FALSE, vcount(g))
  }
  
  #if (anyNA(V(g)$name)) stop("V(g)$name contains NA.")
  #if (anyDuplicated(V(g)$name)) stop("V(g)$name must be unique.")
  if (length(V(g)$resistant) != vcount(g)) stop("V(g)$resistant length mismatch.")
  
  g
}

derive_params <- function(p, AB, tx) {
  q <- 1 - exp(-p)     # P(Poisson(p) >= 1)
  t_kill <- AB + tx
  list(q = q, t_kill = t_kill)
}

init_snapshots <- function(T) {
  setNames(vector("list", T), paste0("t", seq_len(T)))
}