# extract_resistance_timeseries_real.R
# Extract counts_ts from compact sim outputs produced by run_sweep_real_raw.R

# ---- setup -------------------------------------------------------------------
library(tidyverse)
library(vroom)

#run_dir <- "outputs/sweep_real/real_2026-03-10" #borked run
run_dir <- "outputs/sweep_real/real_2026-03-12/"

index_path <-
  file.path(run_dir, "index.rds")

dir_tables <-
  file.path(run_dir, "tables")

output_file <-
  file.path(dir_tables, "resistance_timeseries_all_sims.csv")

# ---- checks ------------------------------------------------------------------
if (!dir.exists(run_dir)) {
  stop("run_dir does not exist: ", run_dir)
}

if (!file.exists(index_path)) {
  stop("index.rds not found: ", index_path)
}

if (!dir.exists(dir_tables)) {
  dir.create(dir_tables, recursive = TRUE)
}

# ---- read index --------------------------------------------------------------
index_tbl <-
  readRDS(index_path) |>
  as_tibble()

index_ok <-
  index_tbl |>
  filter(status == "ok")

if (nrow(index_ok) == 0) {
  stop("No successful simulations found in index.rds")
}

# ---- helper ------------------------------------------------------------------
read_one_counts_ts <- function(row) {
  
  f <-
    row$raw_path
  
  if (!file.exists(f)) {
    warning("Missing raw file: ", f)
    return(NULL)
  }
  
  sim_out <-
    readRDS(f)
  
  counts_ts <-
    sim_out$counts_ts |>
    as_tibble()
  
  bind_cols(
    row |>
      select(
        run_id, sim_id,
        n, graph_seed_id,
        cond_id, k0, seed_mode,
        p, AB, tx, T_used, seed_sim
      ) |>
      slice(rep(1, nrow(counts_ts))),
    counts_ts
  )
}

# ---- extract -----------------------------------------------------------------
resistance_ts_all <-
  map_dfr(
    seq_len(nrow(index_ok)),
    ~read_one_counts_ts(index_ok[.x, ])
  ) |>
  arrange(sim_id, t)

# ---- write -------------------------------------------------------------------
vroom::vroom_write(
  resistance_ts_all,
  output_file
)

message("Wrote: ", output_file)
message("Rows: ", nrow(resistance_ts_all))
message("Cols: ", ncol(resistance_ts_all))