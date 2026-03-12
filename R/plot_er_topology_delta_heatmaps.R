# plot_er_topology_delta_heatmaps.R
#failed script. use full version
# ---- setup -------------------------------------------------------------------
library(tidyverse)

input_file <-
  "outputs/sweep_er/er_2026-03-04/topology_delta_all_sims.csv"

output_dir <-
  "results"


x <-
  vroom::vroom(input_file, show_col_types = FALSE)

# ---- summarise ---------------------------------------------------------------
df_sum <-
  x |>
  filter(n == 100) |>
  mutate(
    pct_delta_n_edges =
      100 * (n_edges_post - n_edges_pre) / n_edges_pre
  ) |>
  group_by(p_edge, k0_frac, p, AB) |>
  summarise(
    mean_pct_delta_n_edges = mean(pct_delta_n_edges, na.rm = TRUE),
    sd_pct_delta_n_edges = sd(pct_delta_n_edges, na.rm = TRUE),
    
    mean_delta_density = mean(delta_density, na.rm = TRUE),
    sd_delta_density = sd(delta_density, na.rm = TRUE),
    
    mean_delta_clustering = mean(delta_clustering, na.rm = TRUE),
    sd_delta_clustering = sd(delta_clustering, na.rm = TRUE),
    
    mean_delta_n_components = mean(delta_n_components, na.rm = TRUE),
    sd_delta_n_components = sd(delta_n_components, na.rm = TRUE),
    
    mean_delta_giant_component_size = mean(delta_giant_component_size, na.rm = TRUE),
    sd_delta_giant_component_size = sd(delta_giant_component_size, na.rm = TRUE),
    
    mean_delta_giant_component_frac = mean(delta_giant_component_frac, na.rm = TRUE),
    sd_delta_giant_component_frac = sd(delta_giant_component_frac, na.rm = TRUE),
    
    mean_delta_mean_distance_gc = mean(delta_mean_distance_gc, na.rm = TRUE),
    sd_delta_mean_distance_gc = sd(delta_mean_distance_gc, na.rm = TRUE),
    .groups = "drop"
  ) |>
  distinct() |>
  arrange(p_edge) |>
  mutate(
    er_panel = LETTERS[seq_len(n())]
  ) |>
  mutate(
    er_panel = factor(er_panel, levels = unique(er_panel)),
    AB = factor(AB),
    k0_frac = factor(k0_frac),
    p = factor(p)
  )

panel_key <-
  df_sum |>
  distinct(er_panel, p_edge) |>
  mutate(
    key_txt = paste0(er_panel, " = p_edge ", p_edge)
  ) |>
  pull(key_txt) |>
  paste(collapse = "   |   ")

# ---- plot specs --------------------------------------------------------------
plot_specs <-
  tribble(
    ~mean_col,                         ~sd_col,                         ~file_stub,                    ~fill_name,                         ~title_text,
    "mean_pct_delta_n_edges",          "sd_pct_delta_n_edges",          "er_pct_delta_n_edges",        "Mean %Δ\nedge count",              "Erdős-Rényi networks: mean percent change in number of edges",
    "mean_delta_density",              "sd_delta_density",              "er_delta_density",            "Mean Δ\ndensity",                  "Erdős-Rényi networks: mean change in density",
    "mean_delta_clustering",           "sd_delta_clustering",           "er_delta_clustering",         "Mean Δ\nclustering",               "Erdős-Rényi networks: mean change in clustering coefficient",
    "mean_delta_n_components",         "sd_delta_n_components",         "er_delta_n_components",       "Mean Δ\ncomponents",               "Erdős-Rényi networks: mean change in number of components",
    "mean_delta_giant_component_size", "sd_delta_giant_component_size", "er_delta_gc_size",            "Mean Δ\ngiant comp. size",         "Erdős-Rényi networks: mean change in giant component size",
    "mean_delta_giant_component_frac", "sd_delta_giant_component_frac", "er_delta_gc_frac",            "Mean Δ\ngiant comp. fraction",     "Erdős-Rényi networks: mean change in giant component fraction",
    "mean_delta_mean_distance_gc",     "sd_delta_mean_distance_gc",     "er_delta_mean_distance_gc",   "Mean Δ\npath length",              "Erdős-Rényi networks: mean change in path length of the giant component"
  )

# ---- helper ------------------------------------------------------------------
make_delta_heatmap <-
  function(data,
           mean_col,
           sd_col,
           fill_name,
           title_text,
           panel_caption_text) {
    
    plot_data <-
      data |>
      transmute(
        er_panel,
        AB,
        k0_frac,
        p,
        mean_value = .data[[mean_col]],
        sd_value = .data[[sd_col]],
        label = paste0(
          sprintf("%.2f", .data[[mean_col]]),
          "\n±",
          sprintf("%.2f", .data[[sd_col]])
        )
      )
    
    lim <-
      max(abs(plot_data$mean_value), na.rm = TRUE)
    
    if (!is.finite(lim) || lim == 0) {
      lim <- 1
    }
    
    ggplot(plot_data, aes(x = k0_frac, y = p, fill = mean_value)) +
      geom_tile(color = "white", linewidth = 0.7) +
      geom_text(aes(label = label), size = 3.3, lineheight = 0.95) +
      facet_grid(
        rows = vars(er_panel),
        cols = vars(AB),
        labeller = labeller(
          er_panel = label_value,
          AB = function(x) paste("Antibiotic time AB =", x)
        )
      ) +
      scale_fill_gradient2(
        low = "#b2182b",
        mid = "white",
        high = "#2166ac",
        midpoint = 0,
        limits = c(-lim, lim),
        name = fill_name
      ) +
      labs(
        title = title_text,
        subtitle = "Numbers inside tiles show mean ± SD across simulation replicates",
        caption = panel_caption_text,
        x = "Initial resistant fraction (k0 / N)",
        y = "Transmission probability (p)"
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        plot.caption = element_text(
          hjust = 0,
          size = 10,
          margin = margin(t = 12)
        )
      )
  }

# ---- export ------------------------------------------------------------------
for (i in seq_len(nrow(plot_specs))) {
  
  p_plot <-
    make_delta_heatmap(
      data = df_sum,
      mean_col = plot_specs$mean_col[[i]],
      sd_col = plot_specs$sd_col[[i]],
      fill_name = plot_specs$fill_name[[i]],
      title_text = plot_specs$title_text[[i]],
      panel_caption_text = paste0("Panel key: ", panel_key)
    )
  
  ggsave(
    filename = file.path(output_dir, paste0(plot_specs$file_stub[[i]], ".pdf")),
    plot = p_plot,
    width = 15,
    height = 10,
    units = "in"
  )
  
  ggsave(
    filename = file.path(output_dir, paste0(plot_specs$file_stub[[i]], ".png")),
    plot = p_plot,
    width = 15,
    height = 10,
    units = "in",
    dpi = 96
  )
}