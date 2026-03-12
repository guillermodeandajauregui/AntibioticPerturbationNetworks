# plot_all_100nodes_ab60_resistance_timeseries_v2.R
# One combined figure for all 100-node simulations with AB = 60
# Line = network structure
# Facets = transmission probability p and initial resistant fraction k0_frac

library(tidyverse)
library(vroom)
library(glue)

# ---- setup -------------------------------------------------------------------
output_dir <-
  "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- helper ------------------------------------------------------------------
# Variables:
# - n_input: tamaño total de la red
# - time: paso de tiempo
# - resistant_state: TRUE si esa fila corresponde a nodos resistentes
# - n_count: número de nodos resistentes en esa fila
# - pc_resistant: porcentaje de nodos resistentes
# - k0_frac: fracción inicial de nodos resistentes
# - p: probabilidad de transmisión de resistencia por paso
# - AB: tiempo de administración de antibiótico
# - structure_label: etiqueta de la estructura de red
#     * ER: incluye p_edge
#     * WS: incluye k y beta
#     * BA: incluye m

prep_ts <- function(path, network_type, n_count_col) {
  
  inutiles <-
    c("run_id", "graph_id", "sim_id", "graph_seed_id", "cond_id", "seed_sim", "graph_rep")
  
  x <-
    vroom::vroom(
      path,
      show_col_types = FALSE
    ) |>
    select(-any_of(inutiles)) |>
    rename(
      time = t,
      resistant_state = resistant,
      n_input = `n...4`,
      n_count = all_of(n_count_col)
    ) |>
    filter(resistant_state) |>
    mutate(
      pc_resistant = 100 * n_count / n_input,
      network_type = network_type
    )
  
  if (network_type == "ER") {
    x |>
      group_by(
        time,
        n_input,
        p_edge,
        k0_frac,
        p,
        AB
      ) |>
      summarise(
        mean_resistant = mean(pc_resistant, na.rm = TRUE),
        sd_resistant = sd(pc_resistant, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        sd_resistant = replace_na(sd_resistant, 0),
        structure_label = glue("ER (p_edge = {p_edge})")
      )
  } else if (network_type == "WS") {
    x |>
      group_by(
        time,
        n_input,
        k,
        beta,
        k0_frac,
        p,
        AB
      ) |>
      summarise(
        mean_resistant = mean(pc_resistant, na.rm = TRUE),
        sd_resistant = sd(pc_resistant, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        sd_resistant = replace_na(sd_resistant, 0),
        structure_label = glue("WS (k = {k}, beta = {beta})")
      )
  } else if (network_type == "BA") {
    x |>
      group_by(
        time,
        n_input,
        m,
        k0_frac,
        p,
        AB
      ) |>
      summarise(
        mean_resistant = mean(pc_resistant, na.rm = TRUE),
        sd_resistant = sd(pc_resistant, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        sd_resistant = replace_na(sd_resistant, 0),
        structure_label = glue("BA (m = {m})")
      )
  }
}

# ---- load --------------------------------------------------------------------
er_sum <-
  prep_ts(
    path = "outputs/sweep_er/er_2026-03-04/tables/resistance_timeseries_all_sims.csv",
    network_type = "ER",
    n_count_col = "n...18"
  )

ws_sum <-
  prep_ts(
    path = "outputs/sweep_ws/ws_2026-03-05/tables/resistance_timeseries_all_sims.csv",
    network_type = "WS",
    n_count_col = "n...19"
  )

ba_sum <-
  prep_ts(
    path = "outputs/sweep_ba/ba_2026-03-05/tables/resistance_timeseries_all_sims.csv",
    network_type = "BA",
    n_count_col = "n...18"
  )

all_sum <-
  bind_rows(er_sum, ws_sum, ba_sum)

# ---- filter target slice -----------------------------------------------------
plot_df <-
  all_sum |>
  filter(
    n_input == 100,
    AB == 60
  ) |>
  mutate(
    p = factor(p),
    k0_frac = factor(k0_frac),
    structure_label = factor(structure_label)
  )

# ---- plot --------------------------------------------------------------------
p_fig <-
  ggplot(
    plot_df,
    aes(
      x = time,
      y = mean_resistant,
      color = structure_label,
      group = structure_label
    )
  ) +
  geom_line(linewidth = 1) +
  facet_grid(
    k0_frac ~ p,
    labeller = label_both,
    scales = "free_x"
  ) +
  labs(
    x = "Time step",
    y = "% resistant population",
    color = "Network structure",
    title = "Resistance dynamics across network structures",
    subtitle = "n = 100, AB = 60"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey70")
  )

p_fig

# ---- save --------------------------------------------------------------------
ggsave(
  filename = file.path(output_dir, "all_100nodes_ab60_resistance_timeseries_v2.pdf"),
  plot = p_fig,
  width = 15,
  height = 10,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(output_dir, "all_100nodes_ab60_resistance_timeseries_v2.png"),
  plot = p_fig,
  width = 15,
  height = 10,
  units = "in",
  dpi = 96
)