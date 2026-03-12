# 06_plot_er_resistance_timeseries.R
# Generate one figure per Erdos-Renyi topology setting
# and save each as PDF and PNG in results/

library(tidyverse)
library(vroom)
library(glue)

# ---- setup -------------------------------------------------------------------
input_file <-
  "outputs/sweep_er/er_2026-03-04/tables/resistance_timeseries_all_sims.csv"

output_dir <-
  "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- load data ----------------------------------------------------------------
# Variables:
# - p_edge: probabilidad de conexión entre pares de nodos en la red Erdos-Renyi
# - k0_frac: fracción inicial de nodos resistentes
# - p: probabilidad de transmisión de resistencia por paso de tiempo
# - AB: tiempo en que se administra el antibiótico
# - n_input: tamaño total de la red (número total de nodos)
# - n_count: número de nodos en el estado indicado en esa fila
# - pc_resistant: porcentaje de población resistente

s_er <-
  vroom::vroom(
    input_file,
    show_col_types = FALSE
  )

inutiles <-
  c("run_id", "graph_id", "sim_id", "graph_seed_id", "cond_id", "seed_sim", "graph_rep")

er_ts <-
  s_er |>
  select(-any_of(inutiles)) |>
  rename(
    time = t,
    resistant_state = resistant,
    n_input = `n...4`,
    n_count = `n...18`
  ) |>
  filter(resistant_state) |>
  mutate(
    pc_resistant = 100 * n_count / n_input
  )

er_sum <-
  er_ts |>
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
    n_reps = n(),
    .groups = "drop"
  ) |>
  mutate(
    sd_resistant = replace_na(sd_resistant, 0)
  )

# ---- combinations to plot ----------------------------------------------------
plot_specs <-
  er_sum |>
  distinct(p_edge, n_input) |>
  arrange(p_edge, n_input)

# ---- plot + save -------------------------------------------------------------
for (i in seq_len(nrow(plot_specs))) {
  
  p_edge_fix <-
    plot_specs$p_edge[[i]]
  
  n_fix <-
    plot_specs$n_input[[i]]
  
  plot_df <-
    er_sum |>
    filter(
      p_edge == p_edge_fix,
      n_input == n_fix
    ) |>
    mutate(
      k0_frac = factor(k0_frac),
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
      . ~ p + AB,
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
      label = "Erdos-Renyi graph",
      subtitle = glue("p_edge = {p_edge_fix}, n = {n_fix}")
    )
  
  file_stub <-
    glue("er_resistance_timeseries_pedge_{p_edge_fix}_n_{n_fix}")
  
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