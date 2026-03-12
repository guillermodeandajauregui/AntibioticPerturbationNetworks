# plot_ba_topology_delta_heatmaps.R

# ---- setup -------------------------------------------------------------------
library(tidyverse)

input_file <-
  "outputs/sweep_ba/ba_2026-03-05/topology_delta_ba.txt"

output_dir <-
  "results"


x <-
  vroom::vroom(input_file)

# ---- summarise ---------------------------------------------------------------
df_sum <-
  x |>
  filter(n == 100) |>
  mutate(
    pct_delta_n_edges =
      100 * (n_edges_post - n_edges_pre) / n_edges_pre
  ) |>
  group_by(m, k0_frac, p, AB) |>
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
  mutate(
    m = factor(m),
    AB = factor(AB),
    k0_frac = factor(k0_frac),
    p = factor(p)
  )

# ---- plot specs --------------------------------------------------------------
plot_specs <-
  tribble(
    ~mean_col,                      ~sd_col,                      ~file_stub,                  ~fill_name,                         ~title_text,
    "mean_pct_delta_n_edges",       "sd_pct_delta_n_edges",       "ba_pct_delta_n_edges",      "Mean %Δ\nedge count",              "Barabási-Albert networks: mean percent change in number of edges",
    "mean_delta_density",           "sd_delta_density",           "ba_delta_density",          "Mean Δ\ndensity",                  "Barabási-Albert networks: mean change in density",
    "mean_delta_clustering",        "sd_delta_clustering",        "ba_delta_clustering",       "Mean Δ\nclustering",               "Barabási-Albert networks: mean change in clustering coefficient",
    "mean_delta_n_components",      "sd_delta_n_components",      "ba_delta_n_components",     "Mean Δ\ncomponents",               "Barabási-Albert networks: mean change in number of components",
    "mean_delta_giant_component_size", "sd_delta_giant_component_size", "ba_delta_gc_size",   "Mean Δ\ngiant comp. size",         "Barabási-Albert networks: mean change in giant component size",
    "mean_delta_giant_component_frac", "sd_delta_giant_component_frac", "ba_delta_gc_frac",   "Mean Δ\ngiant comp. fraction",     "Barabási-Albert networks: mean change in giant component fraction",
    "mean_delta_mean_distance_gc",  "sd_delta_mean_distance_gc",  "ba_delta_mean_distance_gc", "Mean Δ\npath length",              "Barabási-Albert networks: mean change in path length of the giant component"
  )

# ---- helper ------------------------------------------------------------------
make_delta_heatmap <-
  function(data,
           mean_col,
           sd_col,
           fill_name,
           title_text) {
    
    plot_data <-
      data |>
      transmute(
        m = m,
        AB = AB,
        k0_frac = k0_frac,
        p = p,
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
      geom_text(aes(label = label), size = 3.8, lineheight = 0.95) +
      facet_grid(
        rows = vars(m),
        cols = vars(AB),
        labeller = labeller(
          m = function(x) paste("BA attachment m =", x),
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
        subtitle = "Tiles and labels show mean ± SD across simulation replicates",
        x = "Initial resistant fraction (k0 / N)",
        y = "Transmission probability (p)"
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")
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
      title_text = plot_specs$title_text[[i]]
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