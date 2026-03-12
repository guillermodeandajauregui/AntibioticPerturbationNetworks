# 04_plot_ws_resistance_timeseries.R
# Generate one figure per Watts-Strogatz topology setting
# and save each as PDF and PNG in results/

library(tidyverse)
library(vroom)
library(glue)

# ---- setup -------------------------------------------------------------------
input_file <-
  "outputs/sweep_ws/ws_2026-03-05/tables/resistance_timeseries_all_sims.csv"

output_dir <-
  "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- load data ----------------------------------------------------------------
# Variables:
# - k: número de vecinos iniciales en la red Watts-Strogatz
# - beta: probabilidad de rewiring en la red Watts-Strogatz
# - k0_frac: fracción inicial de nodos resistentes
# - p: probabilidad de transmisión de resistencia por paso de tiempo
# - AB: tiempo en que se administra el antibiótico
# - n_input: tamaño total de la red (número total de nodos)
# - n_count: número de nodos en el estado indicado en esa fila
# - pc_resistant: porcentaje de población resistente

s_ws <-
  vroom::vroom(
    input_file,
    show_col_types = FALSE
  )

inutiles <-
  c("run_id", "graph_id", "sim_id", "graph_seed_id", "cond_id", "seed_sim", "graph_rep")

ws_ts <-
  s_ws |>
  select(-any_of(inutiles)) |>
  rename(
    time = t,
    resistant_state = resistant,
    n_input = `n...4`,
    n_count = `n...19`
  ) |>
  filter(resistant_state) |>
  mutate(
    pc_resistant = 100 * n_count / n_input
  )

ws_sum <-
  ws_ts |>
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
    n_reps = n(),
    .groups = "drop"
  ) |>
  mutate(
    sd_resistant = replace_na(sd_resistant, 0)
  )

# ---- combinations to plot ----------------------------------------------------
plot_specs <-
  ws_sum |>
  distinct(k, n_input) |>
  arrange(k, n_input)

# ---- plot + save -------------------------------------------------------------
for (i in seq_len(nrow(plot_specs))) {
  
  k_fix <-
    plot_specs$k[[i]]
  
  n_fix <-
    plot_specs$n_input[[i]]
  
  plot_df <-
    ws_sum |>
    filter(
      k == k_fix,
      n_input == n_fix
    ) |>
    mutate(
      k0_frac = factor(k0_frac),
      beta = factor(beta),
      p = factor(p),
      AB = factor(AB)
    )
  
  p_fig <-
    ggplot(
      plot_df,
      aes(
        x = time,
        y = mean_resistant,
        color = k0_frac,
        fill = k0_frac,
        group = k0_frac
      )
    ) +
    geom_ribbon(
      aes(
        ymin = pmax(mean_resistant - sd_resistant, 0),
        ymax = pmin(mean_resistant + sd_resistant, 100)
      ),
      alpha = 0.18,
      color = NA
    ) +
    geom_line(linewidth = 1) +
    facet_grid(
      beta ~ p + AB,
      labeller = label_both,
      scales = "free_x"
    ) +
    labs(
      x = "Time step",
      y = "% resistant population",
      color = "k0_frac",
      fill = "k0_frac"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70")
    ) +
    ggtitle(
      label = "Watts-Strogatz graph",
      subtitle = glue("k = {k_fix}, n = {n_fix}")
    )
  
  file_stub <-
    glue("ws_resistance_timeseries_k_{k_fix}_n_{n_fix}")
  
  ggsave(
    filename = file.path(output_dir, glue("{file_stub}.pdf")),
    plot = p_fig,
    width = 15,
    height = 10,
    units = "in",
    device = cairo_pdf
  )
  
  ggsave(
    filename = file.path(output_dir, glue("{file_stub}.png")),
    plot = p_fig,
    width = 15,
    height = 10,
    units = "in",
    dpi = 96
  )
}