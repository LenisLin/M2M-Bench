theme_m2m <- function(base_family = "Arial", base_size = 10.5) {
  ggplot2::theme_minimal(base_family = base_family, base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#EEF2F6", linewidth = 0.16),
      panel.spacing = grid::unit(8, "mm"),
      panel.spacing.y = grid::unit(8.5, "mm"),
      strip.background = ggplot2::element_rect(fill = "#F2F5F8", color = "#E2E8EF", linewidth = 0.32),
      strip.text = ggplot2::element_text(
        face = "bold",
        color = "#1F2D3A",
        size = base_size - 0.6,
        margin = ggplot2::margin(3.2, 3, 3.2, 3)
      ),
      strip.text.x = ggplot2::element_text(
        margin = ggplot2::margin(2.8, 3, 3.8, 3)
      ),
      strip.text.y = ggplot2::element_text(
        margin = ggplot2::margin(2.8, 4, 2.8, 3)
      ),
      strip.switch.pad.grid = grid::unit(1.5, "mm"),
      axis.title.x = ggplot2::element_text(
        face = "bold",
        color = "#1F2D3A",
        size = base_size - 0.2,
        margin = ggplot2::margin(t = 7)
      ),
      axis.title.y = ggplot2::element_text(
        face = "bold",
        color = "#1F2D3A",
        size = base_size - 0.2,
        margin = ggplot2::margin(r = 7)
      ),
      axis.text = ggplot2::element_text(color = "#42576D", size = base_size - 1.25),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 0.7, color = "#111111"),
      plot.tag.position = c(0, 1),
      plot.margin = ggplot2::margin(9, 9, 7, 9),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", color = "#1F2D3A", size = base_size - 0.8),
      legend.text = ggplot2::element_text(color = "#3B5068", size = base_size - 1.1),
      legend.key = ggplot2::element_blank(),
      legend.key.size = grid::unit(3.8, "mm"),
      legend.box.spacing = grid::unit(1, "mm"),
      legend.spacing.x = grid::unit(2.1, "mm"),
      legend.margin = ggplot2::margin(1, 0, 0, 0)
    )
}
