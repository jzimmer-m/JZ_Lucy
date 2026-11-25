suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(tibble)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_file_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/"))
} else {
  getwd()
}

repo_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/")
auc_path <- file.path(repo_dir, "outputs", "bl251_auc_summary.csv")
protein_path <- file.path(repo_dir, "data", "processed", "bl251_invitro_protein_expression.csv")
comparison_path <- file.path(repo_dir, "outputs", "bl251_comparison_values_40ng.csv")
color_path <- file.path(repo_dir, "config", "bl251_lnp_colors.csv")
figure_dir <- file.path(repo_dir, "outputs", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(
  theme_classic(base_size = 10, base_family = "sans") +
    theme(
      axis.line = element_line(color = "grey35", linewidth = 0.3),
      axis.ticks = element_line(color = "grey35", linewidth = 0.3),
      axis.text = element_text(color = "black", size = 9),
      axis.text.x = element_text(size = 8),
      axis.title = element_text(color = "black", face = "bold", size = 12),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      legend.key.height = grid::unit(0.14, "in"),
      legend.key.width = grid::unit(0.18, "in"),
      panel.border = element_rect(color = "grey65", fill = NA, linewidth = 0.35),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      panel.spacing.x = grid::unit(0.48, "in"),
      panel.spacing.y = grid::unit(0.24, "in"),
      plot.margin = margin(6, 6, 6, 6),
      strip.background = element_rect(fill = "grey96", color = "grey65", linewidth = 0.35),
      strip.text = element_text(color = "black", face = "bold", size = 11)
  )
)

log_breaks <- function(limits) {
  limits <- limits[is.finite(limits) & limits > 0]
  if (length(limits) == 0L) {
    return(numeric())
  }

  first_power <- floor(log10(min(limits)))
  last_power <- floor(log10(max(limits)))
  breaks <- as.vector(outer(c(1, 2, 5), 10 ^ seq.int(first_power, last_power)))
  sort(breaks[breaks >= min(limits) & breaks <= max(limits)])
}

x_log_breaks <- function(limits) {
  limits <- limits[is.finite(limits) & limits > 0]
  if (length(limits) == 0L) {
    return(numeric())
  }

  first_power <- ceiling(log10(min(limits)))
  last_power <- floor(log10(max(limits)))
  if (first_power > last_power) {
    return(numeric())
  }

  breaks <- as.vector(outer(c(1, 5), 10 ^ seq.int(first_power, last_power)))
  sort(breaks[breaks <= max(limits)])
}

label_round_ticks <- function(x) {
  vapply(
    x,
    function(value) {
      if (!is.finite(value) || value <= 0) {
        return("")
      }

      format(value, scientific = FALSE, trim = TRUE)
    },
    character(1)
  )
}

protein_log_breaks <- function(limits) {
  limits <- limits[is.finite(limits) & limits > 0]
  if (length(limits) == 0L) {
    return(numeric())
  }
  breaks <- pretty(limits, n = 6)
  breaks[breaks > 0 & breaks >= min(limits) & breaks <= max(limits)]
}

axis_value_label <- function(metric, dataset) {
  label <- c(
    "in vitro Scar-seq" = "In vitro Scar-seq cyto abundance AUC",
    "in vivo Scar-seq" = "In vivo Scar-seq cyto abundance AUC",
    "in vitro protein expression" = "In vitro protein expression at 40 ng, 24 h (Alexa 488 mean per well)"
  )[dataset]

  if (!is.na(label)) {
    return(unname(label))
  }

  paste(metric, dataset)
}

figure_specs <- tribble(
  ~x_dataset, ~y_dataset, ~stem,
  "in vitro Scar-seq", "in vitro protein expression", "bl251_auc_invitro_scarseq_vs_invitro_protein",
  "in vivo Scar-seq", "in vitro protein expression", "bl251_auc_invivo_scarseq_vs_invitro_protein"
)

make_pairs <- function(values, x_dataset, y_dataset) {
  x_values <- values %>%
    filter(.data$dataset == x_dataset) %>%
    transmute(
      lnp = .data$lnp,
      x_context = .data$context,
      x_metric = .data$metric,
      x_value = .data$value
    )

  y_values <- values %>%
    filter(.data$dataset == y_dataset) %>%
    transmute(
      lnp = .data$lnp,
      y_context = .data$context,
      y_metric = .data$metric,
      y_value = .data$value
    )

  inner_join(x_values, y_values, by = "lnp", relationship = "many-to-many") %>%
    filter(is.finite(.data$x_value), is.finite(.data$y_value), .data$x_value > 0, .data$y_value > 0) %>%
    mutate(
      x_context = factor(.data$x_context, levels = sort(unique(.data$x_context))),
      y_context = factor(.data$y_context, levels = sort(unique(.data$y_context)))
    )
}

make_stats <- function(pairs) {
  pairs %>%
    group_by(.data$x_context, .data$y_context) %>%
    summarize(
      rho = cor(.data$x_value, .data$y_value, method = "spearman"),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(label = sprintf("rho = %.2f\nn = %d", .data$rho, .data$n))
}

plot_pairs <- function(pairs, stats, colors, x_dataset, y_dataset) {
  x_metric <- pairs$x_metric[[1]]
  y_metric <- pairs$y_metric[[1]]
  x_range <- range(log10(pairs$x_value))
  y_range <- range(log10(pairs$y_value))
  stats$x_label <- 10 ^ (x_range[[1]] + 0.02 * diff(x_range))
  stats$y_label <- 10 ^ (y_range[[2]] - 0.02 * diff(y_range))

  ggplot(pairs, aes(x = .data$x_value, y = .data$y_value)) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "black",
      linewidth = 0.45,
      na.rm = TRUE
    ) +
    geom_point(
      aes(fill = .data$lnp),
      shape = 21,
      color = "white",
      stroke = 0.25,
      size = 2.4,
      alpha = 0.95
    ) +
    geom_text(
      data = stats,
      aes(x = .data$x_label, y = .data$y_label, label = .data$label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      lineheight = 0.92,
      size = 2.6
    ) +
    facet_grid(y_context ~ x_context) +
    scale_x_log10(breaks = x_log_breaks, labels = label_round_ticks) +
    scale_y_log10(
      breaks = if (y_dataset == "in vitro protein expression") protein_log_breaks else log_breaks,
      labels = label_round_ticks
    ) +
    scale_fill_manual(values = colors, breaks = names(colors), drop = FALSE) +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))) +
    labs(
      x = axis_value_label(x_metric, x_dataset),
      y = axis_value_label(y_metric, y_dataset)
    )
}

scar_values <- read_csv(auc_path, show_col_types = FALSE) %>%
  filter(.data$dataset != "in vitro protein expression") %>%
  transmute(
    metric = .data$metric,
    dataset = .data$dataset,
    context = .data$context,
    lnp = .data$lnp,
    measure = "time_auc",
    value = .data$auc
  )
write_csv(
  scar_values %>% transmute(metric, dataset, context, lnp, auc = value),
  auc_path
)
protein_values <- read_csv(protein_path, show_col_types = FALSE) %>%
  filter(.data$dose_ng == 40, is.finite(.data$protein_expression)) %>%
  mutate(lnp = if_else(.data$lnp %in% c(
    "86/56 11", "86/56 11 member pool", "86/56 11-plex", "86/56 20", "86/56 20-plex"
  ), "86/56", .data$lnp)) %>%
  group_by(.data$dataset, .data$context, .data$lnp) %>%
  summarize(value = mean(.data$protein_expression), .groups = "drop") %>%
  mutate(metric = "In vitro protein expression at 40 ng", measure = "single_dose_mean") %>%
  select(metric, dataset, context, lnp, measure, value)
if (nrow(protein_values) == 0L) {
  stop("No 40 ng protein-expression measurements found")
}
values <- bind_rows(scar_values, protein_values)
write_csv(values, comparison_path)
colors <- read_csv(color_path, show_col_types = FALSE)
color_values <- setNames(colors$color, colors$lnp)

for (i in seq_len(nrow(figure_specs))) {
  spec <- figure_specs[i, ]
  pairs <- make_pairs(values, spec$x_dataset, spec$y_dataset)
  stats <- make_stats(pairs)
  p <- plot_pairs(pairs, stats, color_values, spec$x_dataset, spec$y_dataset)

  n_x <- n_distinct(pairs$x_context)
  n_y <- n_distinct(pairs$y_context)
  width <- 2.60 * n_x + 2.0
  height <- 1.80 * n_y + 1.6

  ggsave(
    filename = file.path(figure_dir, paste0(spec$stem, ".pdf")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  ggsave(
    filename = file.path(figure_dir, paste0(spec$stem, ".png")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}
