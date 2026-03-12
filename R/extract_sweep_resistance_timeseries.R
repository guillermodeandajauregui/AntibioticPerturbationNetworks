# extract_sweep_resistance_timeseries.R
# Read per-simulation counts_ts objects from a sweep run, attach metadata,
# stack all simulations into a single table, and write one CSV.
#
# Expected inputs in run_dir:
# - sim_index.csv
# - raw/ ... per-simulation .rds files referenced by sim_index$counts_ts_file
#
# Output:
# - tables/resistance_timeseries_all_sims.csv

# ---- setup -------------------------------------------------------------------
library(tidyverse)
library(vroom)

# ---- user params --------------------------------------------------------------
run_dir <-
  "outputs/sweep_er/er_2026-03-04"

index_file <-
  file.path(run_dir, "sim_index.csv")

dir_tables <-
  file.path(run_dir, "tables")

output_file <-
  file.path(dir_tables, "resistance_timeseries_all_sims.csv")

# ---- checks ------------------------------------------------------------------
if (!dir.exists(run_dir)) {
  stop("run_dir does not exist: ", run_dir)
}

if (!file.exists(index_file)) {
  stop("sim_index.csv not found at: ", index_file)
}

if (!dir.exists(dir_tables)) {
  stop("tables directory does not exist: ", dir_tables)
}

# ---- helpers -----------------------------------------------------------------
read_counts_ts_one <- function(index_row) {
  
  counts_path <-
    file.path(run_dir, index_row$counts_ts_file)
  
  if (!file.exists(counts_path)) {
    warning("Missing counts_ts file: ", counts_path)
    return(NULL)
  }
  
  x <-
    readRDS(counts_path)
  
  if (!is.data.frame(x)) {
    warning("counts_ts is not a data.frame/tibble: ", counts_path)
    return(NULL)
  }
  
  x <-
    as_tibble(x)
  
  # Standardize expected columns if possible
  if (!"t" %in% colnames(x)) {
    warning("counts_ts missing column `t`: ", counts_path)
    return(NULL)
  }
  
  # Common expected columns from the simulator
  possible_count_cols <-
    c("n_resistant", "n_sensitive", "n_non_resistant",
      "frac_resistant", "frac_sensitive", "frac_non_resistant")
  
  # Attach metadata from index
  x |>
    mutate(
      sim_id = index_row$sim_id,
      graph_seed = dplyr::coalesce(index_row$graph_seed, NA),
      sim_seed = dplyr::coalesce(index_row$sim_seed, NA)
    ) |>
    bind_cols(
      index_row |>
        select(-counts_ts_file)
    ) |>
    relocate(
      sim_id, t
    )
}

# ---- read index ---------------------------------------------------------------
sim_index <-
  vroom::vroom(
    file = index_file,
    show_col_types = FALSE,
    progress = FALSE
  ) |>
  as_tibble()

required_cols <-
  c("sim_id", "counts_ts_file")

missing_required <-
  setdiff(required_cols, colnames(sim_index))

if (length(missing_required) > 0) {
  stop(
    "sim_index.csv is missing required columns: ",
    paste(missing_required, collapse = ", ")
  )
}

# ---- extract -----------------------------------------------------------------
resistance_ts_all <-
  purrr::map_dfr(
    .x = seq_len(nrow(sim_index)),
    .f = function(i) {
      read_counts_ts_one(sim_index[i, ])
    }
  )

if (nrow(resistance_ts_all) == 0) {
  stop("No time series could be read.")
}

# ---- optional standardization ------------------------------------------------
# If the simulator used n_non_resistant instead of n_sensitive, normalize name.
if ("n_non_resistant" %in% colnames(resistance_ts_all) &&
    !"n_sensitive" %in% colnames(resistance_ts_all)) {
  resistance_ts_all <-
    resistance_ts_all |>
    rename(n_sensitive = n_non_resistant)
}

if ("frac_non_resistant" %in% colnames(resistance_ts_all) &&
    !"frac_sensitive" %in% colnames(resistance_ts_all)) {
  resistance_ts_all <-
    resistance_ts_all |>
    rename(frac_sensitive = frac_non_resistant)
}

# Recompute fractions if absent and counts are available
if (!"frac_resistant" %in% colnames(resistance_ts_all) &&
    all(c("n_resistant", "n_sensitive") %in% colnames(resistance_ts_all))) {
  resistance_ts_all <-
    resistance_ts_all |>
    mutate(
      n_total = n_resistant + n_sensitive,
      frac_resistant = if_else(n_total > 0, n_resistant / n_total, NA_real_),
      frac_sensitive = if_else(n_total > 0, n_sensitive / n_total, NA_real_)
    )
}

# ---- arrange -----------------------------------------------------------------
id_first <-
  c(
    "sim_id", "t",
    "n", "p_edge", "m", "k", "beta",
    "k0_frac", "k0", "seed_mode",
    "p", "AB", "tx", "T_used",
    "graph_seed", "sim_seed"
  )

id_first <-
  id_first[id_first %in% colnames(resistance_ts_all)]

other_cols <-
  setdiff(colnames(resistance_ts_all), id_first)

resistance_ts_all <-
  resistance_ts_all |>
  select(all_of(id_first), all_of(other_cols)) |>
  arrange(sim_id, t)

# ---- write -------------------------------------------------------------------
vroom::vroom_write(
  x = resistance_ts_all,
  file = output_file,
  delim = ","
)

message("Wrote: ", output_file)
message("Rows: ", nrow(resistance_ts_all))
message("Cols: ", ncol(resistance_ts_all))