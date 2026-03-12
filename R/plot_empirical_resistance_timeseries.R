# plot_empirical_resistance_timeseries.R
# Generate one figure per empirical-network setting
# and save each as PDF and PNG in results/

library(tidyverse)
library(vroom)
library(glue)

# ---- setup -------------------------------------------------------------------
input_file <-
  "outputs/sweep_real/real_2026-03-12/tables/resistance_timeseries_all_sims.csv"

output_dir <-
  "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- load data ----------------------------------------------------------------
# Variables:
# - k0: número inicial de nodos resistentes
# - seed_mode: modo de siembra inicial de resistencia
# - p: probabilidad de transmisión de resistencia por paso de tiempo
# - AB: tiempo en que se administra el antibiótico
# - tx: duración o ventana de tratamiento con antibiótico
# - T_used: número total de pasos simulados
# - n_input: tamaño total de la red (número total de nodos)
# - n_count: número de nodos en el estado indicado en esa fila
# - pc_resistant: porcentaje de población resistente

s_real <-
  vroom::vroom(
    input_file,
    show_col_types = FALSE
  )

inutiles <-
  c("run_id", "sim_id", "graph_seed_id", "cond_id", "seed_sim")

real_ts <-
  s_real |>
  select(-any_of(inutiles)) |>
  rename(
    n_input = `n...3`,
    time = t,
    resistant_state = resistant,
    n_count = `n...15`
  ) |>
  filter(resistant_state) |>
  mutate(
    pc_resistant = 100 * n_count / n_input
  )

real_sum <-
  real_ts |>
  group_by(
    time,
    n_input,
    k0,
    seed_mode,
    p,
    AB,
    tx,
    T_used
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

plot_specs <-
  real_sum |>
  distinct(k0, seed_mode, n_input, tx, T_used) |>
  arrange(k0, seed_mode, n_input, tx, T_used)

for (i in seq_len(nrow(plot_specs))) {
  
  k0_fix <-
    plot_specs$k0[[i]]
  
  seed_mode_fix <-
    plot_specs$seed_mode[[i]]
  
  n_fix <-
    plot_specs$n_input[[i]]
  
  tx_fix <-
    plot_specs$tx[[i]]
  
  T_used_fix <-
    plot_specs$T_used[[i]]
  
  plot_df <-
    real_sum |>
    filter(
      k0 == k0_fix,
      seed_mode == seed_mode_fix,
      n_input == n_fix,
      tx == tx_fix,
      T_used == T_used_fix
    ) |>
    mutate(
      p = factor(p),
      AB = factor(AB)
    )
  
  p_fig <-
    ggplot(
      plot_df,
      aes(
        x = time,
        y = mean_resistant,
        group = 1
      )
    ) +
    geom_ribbon(
      aes(
        ymin = pmax(mean_resistant - sd_resistant, 0),
        ymax = pmin(mean_resistant + sd_resistant, 100)
      ),
      alpha = 0.18,
      fill = "#4B9CD3",
      color = NA
    ) +
    geom_line(
      linewidth = 1,
      color = "#1F4E79"
    ) +
    facet_grid(
      . ~ p + AB,
      labeller = label_both,
      scales = "free_x"
    ) +
    labs(
      x = "Time step",
      y = "% resistant population"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70")
    ) +
    ggtitle(
      label = "Empirical graph",
      subtitle = glue(
        "k0 = {k0_fix}, seed_mode = {seed_mode_fix}, n = {n_fix}, tx = {tx_fix}, T_used = {T_used_fix}"
      )
    )
  
  seed_mode_safe <-
    seed_mode_fix |>
    str_replace_all("[^[:alnum:]]+", "_") |>
    str_to_lower()
  
  file_stub <-
    glue(
      "empirical_resistance_timeseries_k0_{k0_fix}_seed_{seed_mode_safe}_n_{n_fix}_tx_{tx_fix}_T_{T_used_fix}"
    )
  
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
