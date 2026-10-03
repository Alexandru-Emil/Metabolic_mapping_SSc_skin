# Cell-level Figure 7A-B plots adapted from
# data_preprocessing/imc/workflows/Interaction/Expr_hasneigh_plots_violin_nostatistics.Rmd.
# The deposited scores and saved neighbor groups are used without aggregation.
methi_neighbor_score_plot <- function(d, score) {
  groups <- c("Without Methi_Fib neighbour", "With Methi_Fib neighbour")
  if (!setequal(unique(d$neighbor_group), groups)) stop("Both Methi_Fib neighbor groups are required")
  if (!score %in% c("Glycolysis_z_c", "TCA_z_c")) stop("Unexpected metabolic score")
  d$neighbor_group <- factor(d$neighbor_group, levels = groups)
  baseline <- median(d[[score]][d$neighbor_group == groups[1]])
  ggplot(d, aes(x = neighbor_group, y = .data[[score]], fill = neighbor_group)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_point(position = position_jitter(width = 0.05, seed = 1), size = 0.1, alpha = 0.1) +
    geom_hline(yintercept = baseline, linetype = "dashed", color = "blue") +
    scale_fill_manual(values = setNames(c("#B2DF8A", "#FB9A99"), groups), guide = "none") +
    scale_x_discrete(labels = c("-", "+")) +
    labs(title = if(score == "Glycolysis_z_c") "Glycolysis score" else "TCA/OXPHOS score",
      x = expression(Met^hi * "_Fib as neighbours"), y = "Expression (z-rescaled)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black"))
}
methi_neighbor_violins <- function(d, population) {
  plots <- lapply(c("Glycolysis_z_c", "TCA_z_c"), function(score) methi_neighbor_score_plot(d, score))
  patchwork::wrap_plots(plots, nrow = 1) +
    patchwork::plot_annotation(title = paste(population, "neighboring Methi_Fib"),
      theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5)))
}
