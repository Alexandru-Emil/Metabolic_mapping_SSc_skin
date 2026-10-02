# Plot expressions and functions retained from the project analysis scripts.
# Only source-data bindings are replaced by the panel adapters.
original_modules <- list()
original_modules[["six"]] <- list(
  setup = expression(
fib_levels <- c("Met_hi_Fib", "Other_Fib"),
fib_colors <- c(Met_hi_Fib = "#F8766D", Other_Fib = "#00BFC4"),
ec_levels <- c("Met_hi_EC", "Other_EC"),
ec_colors <- c(Met_hi_EC = "#E69F00", Other_EC = "#56B4E9"),
score_palette <- c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C"),
score_limits <- c(-1, 1),
score_breaks <- seq(-1, 1, by = 0.5),
heatmap_palette <- viridisLite::viridis(100),
metabolic_markers <- c("GLUT1", "HK1", "PFKL1_PFKM", "PKM2", "LDHA", "GAPDH", "CS", "OGDH", "ATP5A", "SDHA", "CD98", "ACAC", 
    "CPT1a", "G6PD", "pAMPK", "p_mTOR", "NOX4", "HIF1a", "PGC1a", "mtTFA", "pNRF2"),
glycolysis_markers <- c("GLUT1", "HK1", "PFKL1_PFKM", "PKM2", "LDHA", "GAPDH"),
tca_oxphos_markers <- c("CS", "OGDH", "ATP5A", "SDHA"),
protein_markers <- c("aSMA", "Cdh11", "FAP"),
save_plot_set <- function(plot, basename, width, height) {
    png_file <- file.path(output_dir, paste0(basename, ".png"))
    pdf_file <- file.path(output_dir, paste0(basename, ".pdf"))
    svg_file <- file.path(output_dir, paste0(basename, ".svg"))
    ggsave(png_file, plot = plot, width = width, height = height, dpi = 450, bg = "white")
    grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE, bg = "white")
    print(plot)
    grDevices::dev.off()
    grDevices::svg(svg_file, width = width, height = height, bg = "white")
    print(plot)
    grDevices::dev.off()
    invisible(c(png_file, pdf_file, svg_file))
},
draw_grob <- function(grob) {
    grid::grid.newpage()
    grid::grid.draw(grob)
},
save_grob_set <- function(grob, basename, width, height) {
    png_file <- file.path(output_dir, paste0(basename, ".png"))
    pdf_file <- file.path(output_dir, paste0(basename, ".pdf"))
    svg_file <- file.path(output_dir, paste0(basename, ".svg"))
    grDevices::png(png_file, width = width, height = height, units = "in", res = 450, bg = "white")
    draw_grob(grob)
    grDevices::dev.off()
    grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE, bg = "white")
    draw_grob(grob)
    grDevices::dev.off()
    grDevices::svg(svg_file, width = width, height = height, bg = "white")
    draw_grob(grob)
    grDevices::dev.off()
    invisible(c(png_file, pdf_file, svg_file))
},
save_plot_aliases <- function(plot, basenames, width, height) {
    invisible(lapply(basenames, function(current_basename) {
        save_plot_set(plot, current_basename, width, height)
    }))
},
save_grob_aliases <- function(grob, basenames, width, height) {
    invisible(lapply(basenames, function(current_basename) {
        save_grob_set(grob, current_basename, width, height)
    }))
},
format_p <- function(p_value) {
    dplyr::case_when(is.na(p_value) ~ "p = NA", p_value < 0.001 ~ "p < 0.001", TRUE ~ paste0("p = ", formatC(p_value, digits = 3, 
        format = "f")))
},
paired_test <- function(data, feature_col, value_col, group_col, group_1, group_2) {
    test_data <- data %>% dplyr::transmute(donor = as.character(donor), feature = as.character(.data[[feature_col]]), group = as.character(.data[[group_col]]), 
        value = as.numeric(.data[[value_col]]))
    dplyr::bind_rows(lapply(unique(test_data$feature), function(current_feature) {
        wide <- test_data %>% dplyr::filter(feature == current_feature, group %in% c(group_1, group_2)) %>% dplyr::select(donor, 
            group, value) %>% tidyr::pivot_wider(names_from = group, values_from = value) %>% dplyr::filter(is.finite(.data[[group_1]]), 
            is.finite(.data[[group_2]]))
        test_result <- if (nrow(wide) >= 2) {
            suppressWarnings(stats::wilcox.test(wide[[group_1]], wide[[group_2]], paired = TRUE, exact = FALSE))
        }
        else {
            NULL
        }
        delta <- wide[[group_1]] - wide[[group_2]]
        data.frame(feature = current_feature, group_1 = group_1, group_2 = group_2, n_paired_donors = nrow(wide), mean_delta_group_1_minus_group_2 = if (length(delta)) 
            mean(delta)
        else NA_real_, median_delta_group_1_minus_group_2 = if (length(delta)) 
            median(delta)
        else NA_real_, p_value = if (is.null(test_result)) 
            NA_real_
        else test_result$p.value, stringsAsFactors = FALSE)
    }))
},
base_umap_theme <- theme_classic(base_size = 6.8) + theme(axis.title = element_text(size = 6.5), axis.text = element_text(size = 5.7, 
    color = "grey20"), axis.ticks = element_line(linewidth = 0.3), axis.line = element_line(linewidth = 0.35), legend.title = element_text(size = 6), 
    legend.text = element_text(size = 6), legend.key.height = grid::unit(0.2, "cm"), legend.key.width = grid::unit(0.22, 
        "cm"), legend.spacing.y = grid::unit(0.02, "cm"), legend.margin = margin(0, 0, 0, 1), plot.title = element_text(size = 7.4, 
        face = "bold", hjust = 0.5, lineheight = 0.95, margin = margin(b = 1.5)), plot.subtitle = element_blank(), plot.margin = margin(2, 
        2, 2, 2)),
make_score_umap <- function(data, score_col, title) {
    ggplot(data, aes(x = UMAP1, y = UMAP2, color = .data[[score_col]])) + geom_point(size = 0.52, alpha = 0.92, stroke = 0) + 
        scale_color_gradientn(colors = score_palette, limits = score_limits, breaks = score_breaks, oob = scales::squish, 
            name = "Score", guide = guide_colorbar(barwidth = grid::unit(0.13, "cm"), barheight = grid::unit(0.86, "cm"), 
                title.position = "top", title.hjust = 0.5)) + labs(title = title, x = "UMAP1", y = "UMAP2") + base_umap_theme + 
        theme(legend.position = "right", legend.title = element_text(size = 5.8, hjust = 0.5), legend.text = element_text(size = 5.4), 
            legend.margin = margin(0, 0, 0, 0), plot.margin = margin(2, 1, 2, 1))
},
make_fib_score_box <- function(current_pathway) {
    panel_data <- fib_score_sample %>% dplyr::filter(as.character(pathway) == current_pathway)
    panel_stat <- fib_score_stat_positions %>% dplyr::filter(as.character(pathway) == current_pathway)
    ggplot(panel_data, aes(x = metfiblabel, y = mean_score, fill = metfiblabel)) + geom_boxplot(width = 0.52, outlier.shape = NA, 
        alpha = 0.7, color = "grey30", linewidth = 0.3) + geom_point(position = position_jitter(width = 0.045, height = 0), 
        size = 0.78, color = "black", alpha = 0.9) + geom_text(data = panel_stat, aes(x = 1.5, y = label_y, label = label), 
        inherit.aes = FALSE, size = 1.65) + scale_fill_manual(values = fib_colors, drop = FALSE) + labs(title = current_pathway, 
        x = NULL, y = "Mean IMC protein score") + theme_bw(base_size = 6.4) + theme(axis.text.x = element_text(size = 6, 
        angle = 20, hjust = 1), axis.text.y = element_text(size = 5.6), axis.title.y = element_text(size = 6.2, margin = margin(r = 1.5)), 
        axis.ticks = element_line(linewidth = 0.3), panel.border = element_rect(fill = NA, color = "grey25", linewidth = 0.35), 
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25), panel.grid.minor = element_line(color = "grey95", 
            linewidth = 0.2), legend.position = "none", plot.title = element_text(size = 7, face = "bold", hjust = 0.5, margin = margin(b = 0.5)), 
        plot.margin = margin(2, 3, 3, 8))
}
  ),
  expressions = list(
"fib_levels" = quote(c("Met_hi_Fib", "Other_Fib")),
"fib_colors" = quote(c(Met_hi_Fib = "#F8766D", Other_Fib = "#00BFC4")),
"ec_levels" = quote(c("Met_hi_EC", "Other_EC")),
"ec_colors" = quote(c(Met_hi_EC = "#E69F00", Other_EC = "#56B4E9")),
"score_palette" = quote(c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C")),
"score_limits" = quote(c(-1, 1)),
"score_breaks" = quote(seq(-1, 1, by = 0.5)),
"heatmap_palette" = quote(viridisLite::viridis(100)),
"metabolic_markers" = quote(c("GLUT1", "HK1", "PFKL1_PFKM", "PKM2", "LDHA", "GAPDH", "CS", "OGDH", "ATP5A", "SDHA", "CD98", "ACAC", "CPT1a", "G6PD", 
    "pAMPK", "p_mTOR", "NOX4", "HIF1a", "PGC1a", "mtTFA", "pNRF2")),
"glycolysis_markers" = quote(c("GLUT1", "HK1", "PFKL1_PFKM", "PKM2", "LDHA", "GAPDH")),
"tca_oxphos_markers" = quote(c("CS", "OGDH", "ATP5A", "SDHA")),
"protein_markers" = quote(c("aSMA", "Cdh11", "FAP")),
"save_plot_set" = quote(function(plot, basename, width, height) {
    png_file <- file.path(output_dir, paste0(basename, ".png"))
    pdf_file <- file.path(output_dir, paste0(basename, ".pdf"))
    svg_file <- file.path(output_dir, paste0(basename, ".svg"))
    ggsave(png_file, plot = plot, width = width, height = height, dpi = 450, bg = "white")
    grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE, bg = "white")
    print(plot)
    grDevices::dev.off()
    grDevices::svg(svg_file, width = width, height = height, bg = "white")
    print(plot)
    grDevices::dev.off()
    invisible(c(png_file, pdf_file, svg_file))
}),
"draw_grob" = quote(function(grob) {
    grid::grid.newpage()
    grid::grid.draw(grob)
}),
"save_grob_set" = quote(function(grob, basename, width, height) {
    png_file <- file.path(output_dir, paste0(basename, ".png"))
    pdf_file <- file.path(output_dir, paste0(basename, ".pdf"))
    svg_file <- file.path(output_dir, paste0(basename, ".svg"))
    grDevices::png(png_file, width = width, height = height, units = "in", res = 450, bg = "white")
    draw_grob(grob)
    grDevices::dev.off()
    grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE, bg = "white")
    draw_grob(grob)
    grDevices::dev.off()
    grDevices::svg(svg_file, width = width, height = height, bg = "white")
    draw_grob(grob)
    grDevices::dev.off()
    invisible(c(png_file, pdf_file, svg_file))
}),
"save_plot_aliases" = quote(function(plot, basenames, width, height) {
    invisible(lapply(basenames, function(current_basename) {
        save_plot_set(plot, current_basename, width, height)
    }))
}),
"save_grob_aliases" = quote(function(grob, basenames, width, height) {
    invisible(lapply(basenames, function(current_basename) {
        save_grob_set(grob, current_basename, width, height)
    }))
}),
"format_p" = quote(function(p_value) {
    dplyr::case_when(is.na(p_value) ~ "p = NA", p_value < 0.001 ~ "p < 0.001", TRUE ~ paste0("p = ", formatC(p_value, digits = 3, 
        format = "f")))
}),
"paired_test" = quote(function(data, feature_col, value_col, group_col, group_1, group_2) {
    test_data <- data %>% dplyr::transmute(donor = as.character(donor), feature = as.character(.data[[feature_col]]), group = as.character(.data[[group_col]]), 
        value = as.numeric(.data[[value_col]]))
    dplyr::bind_rows(lapply(unique(test_data$feature), function(current_feature) {
        wide <- test_data %>% dplyr::filter(feature == current_feature, group %in% c(group_1, group_2)) %>% dplyr::select(donor, 
            group, value) %>% tidyr::pivot_wider(names_from = group, values_from = value) %>% dplyr::filter(is.finite(.data[[group_1]]), 
            is.finite(.data[[group_2]]))
        test_result <- if (nrow(wide) >= 2) {
            suppressWarnings(stats::wilcox.test(wide[[group_1]], wide[[group_2]], paired = TRUE, exact = FALSE))
        }
        else {
            NULL
        }
        delta <- wide[[group_1]] - wide[[group_2]]
        data.frame(feature = current_feature, group_1 = group_1, group_2 = group_2, n_paired_donors = nrow(wide), mean_delta_group_1_minus_group_2 = if (length(delta)) 
            mean(delta)
        else NA_real_, median_delta_group_1_minus_group_2 = if (length(delta)) 
            median(delta)
        else NA_real_, p_value = if (is.null(test_result)) 
            NA_real_
        else test_result$p.value, stringsAsFactors = FALSE)
    }))
}),
"base_umap_theme" = quote(theme_classic(base_size = 6.8) + theme(axis.title = element_text(size = 6.5), axis.text = element_text(size = 5.7, color = "grey20"), 
    axis.ticks = element_line(linewidth = 0.3), axis.line = element_line(linewidth = 0.35), legend.title = element_text(size = 6), 
    legend.text = element_text(size = 6), legend.key.height = grid::unit(0.2, "cm"), legend.key.width = grid::unit(0.22, 
        "cm"), legend.spacing.y = grid::unit(0.02, "cm"), legend.margin = margin(0, 0, 0, 1), plot.title = element_text(size = 7.4, 
        face = "bold", hjust = 0.5, lineheight = 0.95, margin = margin(b = 1.5)), plot.subtitle = element_blank(), plot.margin = margin(2, 
        2, 2, 2))),
"p_fib_umap" = quote(ggplot(fib_umap, aes(x = UMAP1, y = UMAP2, color = metfiblabel)) + geom_point(size = 0.55, alpha = 0.9, stroke = 0) + scale_color_manual(values = fib_colors, 
    drop = FALSE) + labs(title = "Metabolic phenotypes of fibroblasts (FibMET)", x = "UMAP1", y = "UMAP2", color = NULL) + 
    base_umap_theme + guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1), keyheight = grid::unit(0.28, 
    "cm"))) + theme(legend.position = "right", plot.title = element_text(size = 5.9, face = "bold", hjust = 0.5, lineheight = 0.95, 
    margin = margin(b = 1)))),
"make_score_umap" = quote(function(data, score_col, title) {
    ggplot(data, aes(x = UMAP1, y = UMAP2, color = .data[[score_col]])) + geom_point(size = 0.52, alpha = 0.92, stroke = 0) + 
        scale_color_gradientn(colors = score_palette, limits = score_limits, breaks = score_breaks, oob = scales::squish, 
            name = "Score", guide = guide_colorbar(barwidth = grid::unit(0.13, "cm"), barheight = grid::unit(0.86, "cm"), 
                title.position = "top", title.hjust = 0.5)) + labs(title = title, x = "UMAP1", y = "UMAP2") + base_umap_theme + 
        theme(legend.position = "right", legend.title = element_text(size = 5.8, hjust = 0.5), legend.text = element_text(size = 5.4), 
            legend.margin = margin(0, 0, 0, 0), plot.margin = margin(2, 1, 2, 1))
}),
"p_fib_glycolysis" = quote(make_score_umap(fib_umap_scores, "glycolysis_score", "Glycolysis score")),
"p_fib_tca" = quote(make_score_umap(fib_umap_scores, "tca_oxphos_score", "TCA/OXPHOS score")),
"p_fib_scores" = quote((p_fib_glycolysis + p_fib_tca + patchwork::plot_layout(guides = "keep") + patchwork::plot_annotation(title = "Metabolic scores", 
    theme = theme(plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5, margin = margin(b = 1)))))),
"fib_protein_stat_positions" = quote(fib_protein_sample %>% dplyr::group_by(marker) %>% dplyr::summarise(ymax = max(mean_expression_z, na.rm = TRUE), spread = diff(range(mean_expression_z, 
    na.rm = TRUE)), .groups = "drop") %>% dplyr::mutate(spread = ifelse(spread > 0, spread, 0.2), label_y = ymax + 0.17 * 
    spread) %>% dplyr::left_join(fib_protein_stats %>% dplyr::select(marker, p_value, n_paired_donors), by = "marker") %>% 
    dplyr::mutate(label = format_p(p_value))),
"p_fib_protein" = quote(ggplot(fib_protein_sample, aes(x = metfiblabel, y = mean_expression_z, fill = metfiblabel)) + geom_violin(width = 0.64, alpha = 0.76, 
    color = "grey25", linewidth = 0.28, trim = TRUE) + geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "grey20", 
    linewidth = 0.28) + geom_point(aes(group = donor), position = position_jitter(width = 0.045, height = 0), size = 0.82, 
    color = "black", alpha = 0.9) + geom_text(data = fib_protein_stat_positions, aes(x = 1.5, y = label_y, label = label), 
    inherit.aes = FALSE, size = 1.55) + scale_fill_manual(values = fib_colors, drop = FALSE) + facet_wrap(~marker, nrow = 1, 
    scales = "free_y") + labs(title = "Protein expression", x = NULL, y = "Mean expression z score") + theme_bw(base_size = 6.2) + 
    theme(axis.text.x = element_text(size = 5.5, angle = 42, hjust = 1), axis.text.y = element_text(size = 5.2), axis.title.y = element_text(size = 6), 
        axis.ticks = element_line(linewidth = 0.3), panel.border = element_rect(fill = NA, color = "grey25", linewidth = 0.35), 
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25), panel.grid.minor = element_line(color = "grey95", 
            linewidth = 0.2), strip.background = element_blank(), strip.text = element_text(size = 6.6, face = "bold", margin = margin(b = 0.5)), 
        legend.position = "none", plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5, margin = margin(b = 1)), 
        plot.margin = margin(2, 2, 3, 2))),
"fib_heatmap_annotation_colors" = quote(list(mols = (grDevices::colorRampPalette(c("#F7F0F5", "#CC79A7")))(100), metfiblabel = fib_colors)),
"fib_heatmap_object" = quote(pheatmap::pheatmap(fib_heatmap_matrix, color = heatmap_palette, breaks = seq(-1, 0.5, length.out = length(heatmap_palette) + 
    1), cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = fib_heatmap_annotation, annotation_colors = fib_heatmap_annotation_colors, 
    show_colnames = FALSE, show_rownames = TRUE, border_color = NA, legend = TRUE, legend_breaks = c(-1, -0.5, 0, 0.5), legend_labels = c("-1.0", 
        "-0.5", "0.0", "0.5"), annotation_legend = TRUE, annotation_names_col = TRUE, fontsize = 6, fontsize_row = 5.2, fontsize_col = 5.2, 
    cellwidth = 42, cellheight = 5.7, treeheight_row = 0, treeheight_col = 0, silent = TRUE)),
"fib_score_stat_positions" = quote(fib_score_sample %>% dplyr::group_by(pathway) %>% dplyr::summarise(ymax = max(mean_score, na.rm = TRUE), spread = diff(range(mean_score, 
    na.rm = TRUE)), .groups = "drop") %>% dplyr::mutate(spread = ifelse(spread > 0, spread, 0.2), label_y = ymax + 0.16 * 
    spread) %>% dplyr::left_join(fib_score_stats %>% dplyr::select(pathway, p_value, n_paired_donors), by = "pathway") %>% 
    dplyr::mutate(label = format_p(p_value))),
"make_fib_score_box" = quote(function(current_pathway) {
    panel_data <- fib_score_sample %>% dplyr::filter(as.character(pathway) == current_pathway)
    panel_stat <- fib_score_stat_positions %>% dplyr::filter(as.character(pathway) == current_pathway)
    ggplot(panel_data, aes(x = metfiblabel, y = mean_score, fill = metfiblabel)) + geom_boxplot(width = 0.52, outlier.shape = NA, 
        alpha = 0.7, color = "grey30", linewidth = 0.3) + geom_point(position = position_jitter(width = 0.045, height = 0), 
        size = 0.78, color = "black", alpha = 0.9) + geom_text(data = panel_stat, aes(x = 1.5, y = label_y, label = label), 
        inherit.aes = FALSE, size = 1.65) + scale_fill_manual(values = fib_colors, drop = FALSE) + labs(title = current_pathway, 
        x = NULL, y = "Mean IMC protein score") + theme_bw(base_size = 6.4) + theme(axis.text.x = element_text(size = 6, 
        angle = 20, hjust = 1), axis.text.y = element_text(size = 5.6), axis.title.y = element_text(size = 6.2, margin = margin(r = 1.5)), 
        axis.ticks = element_line(linewidth = 0.3), panel.border = element_rect(fill = NA, color = "grey25", linewidth = 0.35), 
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25), panel.grid.minor = element_line(color = "grey95", 
            linewidth = 0.2), legend.position = "none", plot.title = element_text(size = 7, face = "bold", hjust = 0.5, margin = margin(b = 0.5)), 
        plot.margin = margin(2, 3, 3, 8))
}),
"p_fib_sample_glycolysis" = quote(make_fib_score_box("Glycolysis")),
"p_fib_sample_tca" = quote(make_fib_score_box("TCA/OXPHOS")),
"p_fib_sample_scores" = quote((p_fib_sample_glycolysis/p_fib_sample_tca + patchwork::plot_layout(heights = c(1, 1)))),
"p_ec_umap" = quote(ggplot(ec_umap, aes(x = UMAP1, y = UMAP2, color = metEClabel)) + geom_point(size = 0.5, alpha = 0.92, stroke = 0) + scale_color_manual(values = ec_colors, 
    drop = FALSE) + labs(title = "Metabolic phenotypes of\nendothelial cells (End-Met)", x = "UMAP_EC_1", y = "UMAP_EC_2", 
    color = NULL) + base_umap_theme + guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1), keyheight = grid::unit(0.28, 
    "cm"))) + theme(legend.position = "right", plot.title = element_text(size = 5.9, face = "bold", hjust = 0.5, lineheight = 0.95, 
    margin = margin(b = 1)))),
"ec_heatmap_object" = quote(pheatmap::pheatmap(ec_heatmap_matrix, color = heatmap_palette, breaks = seq(-1, 1, length.out = length(heatmap_palette) + 
    1), cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = ec_heatmap_annotation, annotation_colors = ec_heatmap_annotation_colors, 
    show_colnames = FALSE, show_rownames = TRUE, border_color = NA, legend = TRUE, legend_breaks = c(-1, -0.5, 0, 0.5, 1), 
    legend_labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"), annotation_legend = TRUE, annotation_names_col = TRUE, fontsize = 6, 
    fontsize_row = 5.2, fontsize_col = 5.2, cellwidth = 42, cellheight = 4.1, treeheight_row = 0, treeheight_col = 0, silent = TRUE))
  ),
  plots = list(
"p_fib_umap" = quote(ggplot(fib_umap, aes(x = UMAP1, y = UMAP2, color = metfiblabel)) + geom_point(size = 0.55, alpha = 0.9, stroke = 0) + scale_color_manual(values = fib_colors, 
    drop = FALSE) + labs(title = "Metabolic phenotypes of fibroblasts (FibMET)", x = "UMAP1", y = "UMAP2", color = NULL) + 
    base_umap_theme + guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1), keyheight = grid::unit(0.28, 
    "cm"))) + theme(legend.position = "right", plot.title = element_text(size = 5.9, face = "bold", hjust = 0.5, lineheight = 0.95, 
    margin = margin(b = 1)))),
"p_fib_protein" = quote(ggplot(fib_protein_sample, aes(x = metfiblabel, y = mean_expression_z, fill = metfiblabel)) + geom_violin(width = 0.64, alpha = 0.76, 
    color = "grey25", linewidth = 0.28, trim = TRUE) + geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "grey20", 
    linewidth = 0.28) + geom_point(aes(group = donor), position = position_jitter(width = 0.045, height = 0), size = 0.82, 
    color = "black", alpha = 0.9) + geom_text(data = fib_protein_stat_positions, aes(x = 1.5, y = label_y, label = label), 
    inherit.aes = FALSE, size = 1.55) + scale_fill_manual(values = fib_colors, drop = FALSE) + facet_wrap(~marker, nrow = 1, 
    scales = "free_y") + labs(title = "Protein expression", x = NULL, y = "Mean expression z score") + theme_bw(base_size = 6.2) + 
    theme(axis.text.x = element_text(size = 5.5, angle = 42, hjust = 1), axis.text.y = element_text(size = 5.2), axis.title.y = element_text(size = 6), 
        axis.ticks = element_line(linewidth = 0.3), panel.border = element_rect(fill = NA, color = "grey25", linewidth = 0.35), 
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25), panel.grid.minor = element_line(color = "grey95", 
            linewidth = 0.2), strip.background = element_blank(), strip.text = element_text(size = 6.6, face = "bold", margin = margin(b = 0.5)), 
        legend.position = "none", plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5, margin = margin(b = 1)), 
        plot.margin = margin(2, 2, 3, 2))),
"p_ec_umap" = quote(ggplot(ec_umap, aes(x = UMAP1, y = UMAP2, color = metEClabel)) + geom_point(size = 0.5, alpha = 0.92, stroke = 0) + scale_color_manual(values = ec_colors, 
    drop = FALSE) + labs(title = "Metabolic phenotypes of\nendothelial cells (End-Met)", x = "UMAP_EC_1", y = "UMAP_EC_2", 
    color = NULL) + base_umap_theme + guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1), keyheight = grid::unit(0.28, 
    "cm"))) + theme(legend.position = "right", plot.title = element_text(size = 5.9, face = "bold", hjust = 0.5, lineheight = 0.95, 
    margin = margin(b = 1))))
  )
)
original_modules[["signal"]] <- list(
  setup = expression(
`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || is.na(x[1])) 
        y
    else x
},
canonical_marker <- function(x) {
    replacements <- c(Glut1 = "GLUT1", NRF2_p = "pNRF2", TFAM = "mtTFA", PDGFRA = "PDGFRa", S6_p = "pS6", pRibosom_S6 = "pS6", 
        `PGC-1a` = "PGC1a", `E-cadherin` = "E-cadherin", E_cadherin = "E-cadherin", Collagen = "COL1A1", PFKL1 = "PFKL/PFKM", 
        PFKL1_PFKM = "PFKL/PFKM", `Ki-67` = "Ki67")
    out <- x
    matched <- match(x, names(replacements))
    out[!is.na(matched)] <- unname(replacements[matched[!is.na(matched)]])
    out
},
get_sample_id <- function(spe) {
    cd <- as.data.frame(colData(spe))
    candidates <- c("sample_id", "patient_id", "roi")
    selected <- candidates[candidates %in% names(cd)][1]
    if (is.na(selected)) {
        stop("No sample, patient, or ROI identifier found in colData.")
    }
    as.character(cd[[selected]])
},
fit_marker_snr <- function(transformed, raw, sample_id, marker, cohort, min_cells_per_component = 20L) {
    keep <- is.finite(transformed) & is.finite(raw) & !is.na(sample_id)
    transformed <- transformed[keep]
    raw <- raw[keep]
    sample_id <- sample_id[keep]
    fit <- tryCatch(suppressWarnings(Mclust(transformed, G = 2, verbose = FALSE)), error = function(e) NULL)
    if (is.null(fit) || length(unique(fit$classification)) != 2) {
        global <- tibble(cohort = cohort, marker = marker, marker_canonical = canonical_marker(marker), n_cells = length(raw), 
            model_name = NA_character_, signal_component = NA_integer_, background_component = NA_integer_, n_signal = NA_integer_, 
            n_background = NA_integer_, fraction_signal = NA_real_, signal_mean = NA_real_, background_mean = NA_real_, snr = NA_real_, 
            log2_signal = NA_real_, log2_snr = NA_real_, fit_status = "mixture_fit_failed")
        return(list(global = global, sample = tibble()))
    }
    component <- as.integer(fit$classification)
    component_means <- tapply(raw, component, mean, na.rm = TRUE)
    signal_component <- as.integer(names(which.max(component_means)))
    background_component <- as.integer(names(which.min(component_means)))
    is_signal <- component == signal_component
    is_background <- component == background_component
    signal_mean <- mean(raw[is_signal], na.rm = TRUE)
    background_mean <- mean(raw[is_background], na.rm = TRUE)
    snr <- if (is.finite(background_mean) && background_mean > 0) {
        signal_mean/background_mean
    }
    else {
        NA_real_
    }
    global <- tibble(cohort = cohort, marker = marker, marker_canonical = canonical_marker(marker), n_cells = length(raw), 
        model_name = fit$modelName %||% NA_character_, signal_component = signal_component, background_component = background_component, 
        n_signal = sum(is_signal), n_background = sum(is_background), fraction_signal = mean(is_signal), signal_mean = signal_mean, 
        background_mean = background_mean, snr = snr, log2_signal = ifelse(signal_mean > 0, log2(signal_mean), NA_real_), 
        log2_snr = ifelse(snr > 0, log2(snr), NA_real_), fit_status = "ok")
    sample_levels <- unique(sample_id)
    sample_stats <- lapply(sample_levels, function(current_sample) {
        idx <- sample_id == current_sample
        cur_signal <- raw[idx & is_signal]
        cur_background <- raw[idx & is_background]
        cur_signal_mean <- if (length(cur_signal)) 
            mean(cur_signal, na.rm = TRUE)
        else NA_real_
        cur_background_mean <- if (length(cur_background)) 
            mean(cur_background, na.rm = TRUE)
        else NA_real_
        cur_snr <- if (is.finite(cur_signal_mean) && is.finite(cur_background_mean) && cur_background_mean > 0) {
            cur_signal_mean/cur_background_mean
        }
        else {
            NA_real_
        }
        tibble(cohort = cohort, sample_id = current_sample, marker = marker, marker_canonical = canonical_marker(marker), 
            n_cells = sum(idx), n_signal = length(cur_signal), n_background = length(cur_background), fraction_signal = length(cur_signal)/sum(idx), 
            signal_mean = cur_signal_mean, background_mean = cur_background_mean, snr = cur_snr, log2_signal = ifelse(cur_signal_mean > 
                0, log2(cur_signal_mean), NA_real_), log2_snr = ifelse(cur_snr > 0, log2(cur_snr), NA_real_), sufficient_components = length(cur_signal) >= 
                min_cells_per_component && length(cur_background) >= min_cells_per_component)
    })
    list(global = global, sample = bind_rows(sample_stats))
},
analyze_cohort <- function(spe, cohort) {
    stopifnot("counts" %in% assayNames(spe), "asinh" %in% assayNames(spe))
    markers <- setdiff(rownames(spe), technical_markers)
    sample_id <- get_sample_id(spe)
    counts_mat <- assay(spe, "counts")
    transformed_mat <- assay(spe, "asinh")
    results <- lapply(markers, function(marker) {
        fit_marker_snr(transformed = as.numeric(transformed_mat[marker, ]), raw = as.numeric(counts_mat[marker, ]), sample_id = sample_id, 
            marker = marker, cohort = cohort)
    })
    list(global = bind_rows(lapply(results, `[[`, "global")), sample = bind_rows(lapply(results, `[[`, "sample")))
},
is_displayed_marker <- function(cohort, marker_canonical) {
    !(cohort == "Second cohort" & marker_canonical == "TCF21")
},
palette <- c(`First cohort` = "#246A9A", `Second cohort` = "#D97706"),
neutral_dark <- "#33383D",
theme_qc <- function(base_size = 11.5) {
    theme_bw(base_size = base_size) + theme(plot.title = element_text(face = "bold", size = base_size + 2, margin = margin(b = 4)), 
        plot.subtitle = element_text(color = "#4A5157", size = base_size - 0.2, margin = margin(b = 7)), plot.caption = element_text(color = "#596168", 
            hjust = 0, size = base_size - 1.2, margin = margin(t = 4)), strip.background = element_rect(fill = "#F3F4F5", 
            color = "#B8BEC3"), strip.text = element_text(face = "bold", size = base_size), panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = "#E7E9EB", linewidth = 0.3), axis.title = element_text(size = base_size), 
        axis.text = element_text(size = base_size - 0.7), legend.title = element_text(size = base_size), legend.text = element_text(size = base_size - 
            0.5), legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 2, 0), plot.title.position = "plot", 
        plot.margin = margin(6, 8, 6, 6))
},
theme_manuscript_panel <- function(base_size = 8.2) {
    theme_bw(base_size = base_size) + theme(plot.title = element_blank(), plot.subtitle = element_blank(), strip.background = element_rect(fill = "#F3F4F5", 
        color = "#AEB4B9", linewidth = 0.35), strip.text = element_text(face = "bold", size = 8.2), panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = "#E5E8EA", linewidth = 0.28), panel.border = element_rect(color = "#44484C", 
            linewidth = 0.45), axis.title = element_text(size = 8.2), axis.text = element_text(size = 7.3), axis.ticks = element_line(linewidth = 0.35), 
        axis.ticks.length = grid::unit(1.4, "mm"), legend.position = "none", plot.margin = margin(4, 4, 4, 5, unit = "pt"))
},
save_plot <- function(plot, stem, width, height, dpi = 400) {
    ggsave(file.path(plot_dir, paste0(stem, ".pdf")), plot, width = width, height = height)
    ggsave(file.path(plot_dir, paste0(stem, ".svg")), plot, width = width, height = height, device = svglite::svglite)
    ggsave(file.path(plot_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
},
make_heatmap <- function(data, cohort_name) {
    ggplot(data %>% filter(cohort == cohort_name), aes(x = sample_short, y = marker_canonical, fill = log2_snr)) + geom_tile(color = "white", 
        linewidth = 0.15) + scale_fill_gradient(low = "white", high = palette[[cohort_name]], limits = heat_limits, na.value = "#E8EAEC") + 
        labs(title = paste(cohort_name, "sample-level SNR"), subtitle = "Common, unclipped color scale across cohorts", x = NULL, 
            y = NULL, fill = expression(log[2] * " SNR")) + theme_qc(11.5) + theme(legend.position = "right", axis.text.x = element_text(angle = 45, 
        hjust = 1, vjust = 1, size = 8.5), axis.text.y = element_text(size = 8.9), panel.grid = element_blank()) + guides(fill = guide_colorbar(barheight = grid::unit(35, 
        "mm"), barwidth = grid::unit(5, "mm")))
}
  ),
  expressions = list(
"%||%" = quote(function(x, y) {
    if (is.null(x) || length(x) == 0 || is.na(x[1])) 
        y
    else x
}),
"canonical_marker" = quote(function(x) {
    replacements <- c(Glut1 = "GLUT1", NRF2_p = "pNRF2", TFAM = "mtTFA", PDGFRA = "PDGFRa", S6_p = "pS6", pRibosom_S6 = "pS6", 
        `PGC-1a` = "PGC1a", `E-cadherin` = "E-cadherin", E_cadherin = "E-cadherin", Collagen = "COL1A1", PFKL1 = "PFKL/PFKM", 
        PFKL1_PFKM = "PFKL/PFKM", `Ki-67` = "Ki67")
    out <- x
    matched <- match(x, names(replacements))
    out[!is.na(matched)] <- unname(replacements[matched[!is.na(matched)]])
    out
}),
"get_sample_id" = quote(function(spe) {
    cd <- as.data.frame(colData(spe))
    candidates <- c("sample_id", "patient_id", "roi")
    selected <- candidates[candidates %in% names(cd)][1]
    if (is.na(selected)) {
        stop("No sample, patient, or ROI identifier found in colData.")
    }
    as.character(cd[[selected]])
}),
"fit_marker_snr" = quote(function(transformed, raw, sample_id, marker, cohort, min_cells_per_component = 20L) {
    keep <- is.finite(transformed) & is.finite(raw) & !is.na(sample_id)
    transformed <- transformed[keep]
    raw <- raw[keep]
    sample_id <- sample_id[keep]
    fit <- tryCatch(suppressWarnings(Mclust(transformed, G = 2, verbose = FALSE)), error = function(e) NULL)
    if (is.null(fit) || length(unique(fit$classification)) != 2) {
        global <- tibble(cohort = cohort, marker = marker, marker_canonical = canonical_marker(marker), n_cells = length(raw), 
            model_name = NA_character_, signal_component = NA_integer_, background_component = NA_integer_, n_signal = NA_integer_, 
            n_background = NA_integer_, fraction_signal = NA_real_, signal_mean = NA_real_, background_mean = NA_real_, snr = NA_real_, 
            log2_signal = NA_real_, log2_snr = NA_real_, fit_status = "mixture_fit_failed")
        return(list(global = global, sample = tibble()))
    }
    component <- as.integer(fit$classification)
    component_means <- tapply(raw, component, mean, na.rm = TRUE)
    signal_component <- as.integer(names(which.max(component_means)))
    background_component <- as.integer(names(which.min(component_means)))
    is_signal <- component == signal_component
    is_background <- component == background_component
    signal_mean <- mean(raw[is_signal], na.rm = TRUE)
    background_mean <- mean(raw[is_background], na.rm = TRUE)
    snr <- if (is.finite(background_mean) && background_mean > 0) {
        signal_mean/background_mean
    }
    else {
        NA_real_
    }
    global <- tibble(cohort = cohort, marker = marker, marker_canonical = canonical_marker(marker), n_cells = length(raw), 
        model_name = fit$modelName %||% NA_character_, signal_component = signal_component, background_component = background_component, 
        n_signal = sum(is_signal), n_background = sum(is_background), fraction_signal = mean(is_signal), signal_mean = signal_mean, 
        background_mean = background_mean, snr = snr, log2_signal = ifelse(signal_mean > 0, log2(signal_mean), NA_real_), 
        log2_snr = ifelse(snr > 0, log2(snr), NA_real_), fit_status = "ok")
    sample_levels <- unique(sample_id)
    sample_stats <- lapply(sample_levels, function(current_sample) {
        idx <- sample_id == current_sample
        cur_signal <- raw[idx & is_signal]
        cur_background <- raw[idx & is_background]
        cur_signal_mean <- if (length(cur_signal)) 
            mean(cur_signal, na.rm = TRUE)
        else NA_real_
        cur_background_mean <- if (length(cur_background)) 
            mean(cur_background, na.rm = TRUE)
        else NA_real_
        cur_snr <- if (is.finite(cur_signal_mean) && is.finite(cur_background_mean) && cur_background_mean > 0) {
            cur_signal_mean/cur_background_mean
        }
        else {
            NA_real_
        }
        tibble(cohort = cohort, sample_id = current_sample, marker = marker, marker_canonical = canonical_marker(marker), 
            n_cells = sum(idx), n_signal = length(cur_signal), n_background = length(cur_background), fraction_signal = length(cur_signal)/sum(idx), 
            signal_mean = cur_signal_mean, background_mean = cur_background_mean, snr = cur_snr, log2_signal = ifelse(cur_signal_mean > 
                0, log2(cur_signal_mean), NA_real_), log2_snr = ifelse(cur_snr > 0, log2(cur_snr), NA_real_), sufficient_components = length(cur_signal) >= 
                min_cells_per_component && length(cur_background) >= min_cells_per_component)
    })
    list(global = global, sample = bind_rows(sample_stats))
}),
"analyze_cohort" = quote(function(spe, cohort) {
    stopifnot("counts" %in% assayNames(spe), "asinh" %in% assayNames(spe))
    markers <- setdiff(rownames(spe), technical_markers)
    sample_id <- get_sample_id(spe)
    counts_mat <- assay(spe, "counts")
    transformed_mat <- assay(spe, "asinh")
    results <- lapply(markers, function(marker) {
        fit_marker_snr(transformed = as.numeric(transformed_mat[marker, ]), raw = as.numeric(counts_mat[marker, ]), sample_id = sample_id, 
            marker = marker, cohort = cohort)
    })
    list(global = bind_rows(lapply(results, `[[`, "global")), sample = bind_rows(lapply(results, `[[`, "sample")))
}),
"is_displayed_marker" = quote(function(cohort, marker_canonical) {
    !(cohort == "Second cohort" & marker_canonical == "TCF21")
}),
"palette" = quote(c(`First cohort` = "#246A9A", `Second cohort` = "#D97706")),
"neutral_dark" = quote("#33383D"),
"theme_qc" = quote(function(base_size = 11.5) {
    theme_bw(base_size = base_size) + theme(plot.title = element_text(face = "bold", size = base_size + 2, margin = margin(b = 4)), 
        plot.subtitle = element_text(color = "#4A5157", size = base_size - 0.2, margin = margin(b = 7)), plot.caption = element_text(color = "#596168", 
            hjust = 0, size = base_size - 1.2, margin = margin(t = 4)), strip.background = element_rect(fill = "#F3F4F5", 
            color = "#B8BEC3"), strip.text = element_text(face = "bold", size = base_size), panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = "#E7E9EB", linewidth = 0.3), axis.title = element_text(size = base_size), 
        axis.text = element_text(size = base_size - 0.7), legend.title = element_text(size = base_size), legend.text = element_text(size = base_size - 
            0.5), legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 2, 0), plot.title.position = "plot", 
        plot.margin = margin(6, 8, 6, 6))
}),
"theme_manuscript_panel" = quote(function(base_size = 8.2) {
    theme_bw(base_size = base_size) + theme(plot.title = element_blank(), plot.subtitle = element_blank(), strip.background = element_rect(fill = "#F3F4F5", 
        color = "#AEB4B9", linewidth = 0.35), strip.text = element_text(face = "bold", size = 8.2), panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = "#E5E8EA", linewidth = 0.28), panel.border = element_rect(color = "#44484C", 
            linewidth = 0.45), axis.title = element_text(size = 8.2), axis.text = element_text(size = 7.3), axis.ticks = element_line(linewidth = 0.35), 
        axis.ticks.length = grid::unit(1.4, "mm"), legend.position = "none", plot.margin = margin(4, 4, 4, 5, unit = "pt"))
}),
"save_plot" = quote(function(plot, stem, width, height, dpi = 400) {
    ggsave(file.path(plot_dir, paste0(stem, ".pdf")), plot, width = width, height = height)
    ggsave(file.path(plot_dir, paste0(stem, ".svg")), plot, width = width, height = height, device = svglite::svglite)
    ggsave(file.path(plot_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
}),
"manuscript_signal_plot" = quote(ggplot(manuscript_signal_data, aes(x = log2_signal, y = log2_snr)) + geom_point(aes(fill = cohort), shape = 21, size = 1.55, 
    stroke = 0.38, color = neutral_dark) + scale_fill_manual(values = palette, guide = "none") + ggrepel::geom_text_repel(data = manuscript_signal_labels, 
    aes(label = marker_canonical), size = 1.95, seed = 220224, max.overlaps = Inf, min.segment.length = 0, box.padding = 0.1, 
    point.padding = 0.07, segment.color = "#9EA5AA", segment.size = 0.2, force = 2.8, force_pull = 0.18, max.iter = 1e+05, 
    max.time = 20) + facet_wrap(vars(cohort), ncol = 1) + scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) + 
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) + labs(x = expression(log[2] * " mean signal intensity"), 
    y = expression(log[2] * " signal-to-noise ratio")) + theme_manuscript_panel(8.2) + theme(axis.text = element_text(size = 7.4), 
    strip.text = element_text(size = 8.4), panel.spacing.y = grid::unit(2.2, "mm"), plot.margin = margin(6, 3, 4, 7, unit = "pt"))),
"manuscript_background_plot" = quote(ggplot(dumbbell_data, aes(y = marker_canonical)) + geom_segment(aes(x = log2_background, xend = log2_signal, yend = marker_canonical), 
    color = "#B9BEC2", linewidth = 0.62, lineend = "round") + geom_point(aes(x = log2_background), shape = 21, fill = "white", 
    color = "#70777D", size = 1.7, stroke = 0.38) + geom_point(aes(x = log2_signal, fill = cohort), shape = 21, color = neutral_dark, 
    size = 1.88, stroke = 0.38) + facet_grid(. ~ cohort) + scale_fill_manual(values = palette, guide = "none") + scale_x_continuous(breaks = c(-6, 
    0, 6), expand = expansion(mult = c(0.04, 0.05))) + labs(x = expression(log[2] * " mean component intensity"), y = NULL) + 
    theme_manuscript_panel(8.2) + theme(axis.text.x = element_text(size = 6.8), axis.text.y = element_text(size = 6.4), strip.text = element_text(size = 7.8), 
    panel.spacing.x = grid::unit(2.2, "mm"), plot.margin = margin(6, 4, 4, 3, unit = "pt"))),
"make_heatmap" = quote(function(data, cohort_name) {
    ggplot(data %>% filter(cohort == cohort_name), aes(x = sample_short, y = marker_canonical, fill = log2_snr)) + geom_tile(color = "white", 
        linewidth = 0.15) + scale_fill_gradient(low = "white", high = palette[[cohort_name]], limits = heat_limits, na.value = "#E8EAEC") + 
        labs(title = paste(cohort_name, "sample-level SNR"), subtitle = "Common, unclipped color scale across cohorts", x = NULL, 
            y = NULL, fill = expression(log[2] * " SNR")) + theme_qc(11.5) + theme(legend.position = "right", axis.text.x = element_text(angle = 45, 
        hjust = 1, vjust = 1, size = 8.5), axis.text.y = element_text(size = 8.9), panel.grid = element_blank()) + guides(fill = guide_colorbar(barheight = grid::unit(35, 
        "mm"), barwidth = grid::unit(5, "mm")))
})
  ),
  plots = list(
"reference_base" = quote(ggplot(plot_global_snr %>% filter(fit_status == "ok", is.finite(log2_signal), is.finite(log2_snr)), aes(x = log2_signal, 
    y = log2_snr)) + geom_point(aes(fill = cohort), shape = 21, size = 3, stroke = 0.55, color = neutral_dark) + scale_fill_manual(values = palette, 
    guide = "none")),
"interval_plot" = quote(ggplot(interval_data, aes(x = median_log2_snr, y = marker_canonical, xmin = q25_log2_snr, xmax = q75_log2_snr, color = cohort)) + 
    geom_vline(xintercept = 0, color = neutral_mid, linewidth = 0.45) + geom_errorbarh(position = position_dodge(width = 0.55), 
    height = 0, linewidth = 0.95) + geom_point(position = position_dodge(width = 0.55), size = 2.8) + scale_color_manual(values = palette) + 
    labs(title = "Per-sample SNR consistency", subtitle = "Point: median; line: interquartile range across samples", x = expression(log[2] * 
        " signal-to-noise ratio"), y = NULL, color = NULL) + theme_qc(11.5) + theme(axis.text.y = element_text(size = 10), 
    legend.position = "top")),
"dumbbell_plot" = quote(ggplot(dumbbell_data, aes(y = marker_canonical)) + geom_segment(aes(x = log2_background, xend = log2_signal, yend = marker_canonical), 
    color = neutral_light, linewidth = 1.4, lineend = "round") + geom_point(aes(x = log2_background), shape = 21, fill = "white", 
    color = neutral_mid, size = 2.8) + geom_point(aes(x = log2_signal, fill = cohort), shape = 21, color = neutral_dark, 
    size = 3.2, stroke = 0.55) + facet_wrap(vars(cohort), nrow = 1) + scale_fill_manual(values = palette, guide = "none") + 
    labs(title = "Empirical signal and background", subtitle = "Open: lower-intensity component; filled: higher-intensity component", 
        x = expression(log[2] * " mean untransformed single-cell intensity"), y = NULL) + theme_qc(11.5) + theme(axis.text.y = element_text(size = 10))),
"manuscript_signal_plot" = quote(ggplot(manuscript_signal_data, aes(x = log2_signal, y = log2_snr)) + geom_point(aes(fill = cohort), shape = 21, size = 1.55, 
    stroke = 0.38, color = neutral_dark) + scale_fill_manual(values = palette, guide = "none") + ggrepel::geom_text_repel(data = manuscript_signal_labels, 
    aes(label = marker_canonical), size = 1.95, seed = 220224, max.overlaps = Inf, min.segment.length = 0, box.padding = 0.1, 
    point.padding = 0.07, segment.color = "#9EA5AA", segment.size = 0.2, force = 2.8, force_pull = 0.18, max.iter = 1e+05, 
    max.time = 20) + facet_wrap(vars(cohort), ncol = 1) + scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) + 
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) + labs(x = expression(log[2] * " mean signal intensity"), 
    y = expression(log[2] * " signal-to-noise ratio")) + theme_manuscript_panel(8.2) + theme(axis.text = element_text(size = 7.4), 
    strip.text = element_text(size = 8.4), panel.spacing.y = grid::unit(2.2, "mm"), plot.margin = margin(6, 3, 4, 7, unit = "pt"))),
"manuscript_background_plot" = quote(ggplot(dumbbell_data, aes(y = marker_canonical)) + geom_segment(aes(x = log2_background, xend = log2_signal, yend = marker_canonical), 
    color = "#B9BEC2", linewidth = 0.62, lineend = "round") + geom_point(aes(x = log2_background), shape = 21, fill = "white", 
    color = "#70777D", size = 1.7, stroke = 0.38) + geom_point(aes(x = log2_signal, fill = cohort), shape = 21, color = neutral_dark, 
    size = 1.88, stroke = 0.38) + facet_grid(. ~ cohort) + scale_fill_manual(values = palette, guide = "none") + scale_x_continuous(breaks = c(-6, 
    0, 6), expand = expansion(mult = c(0.04, 0.05))) + labs(x = expression(log[2] * " mean component intensity"), y = NULL) + 
    theme_manuscript_panel(8.2) + theme(axis.text.x = element_text(size = 6.8), axis.text.y = element_text(size = 6.4), strip.text = element_text(size = 7.8), 
    panel.spacing.x = grid::unit(2.2, "mm"), plot.margin = margin(6, 4, 4, 3, unit = "pt"))),
"concordance_plot" = quote(ggplot(shared_global, aes(x = `log2_snr__First cohort`, y = `log2_snr__Second cohort`)) + geom_abline(intercept = 0, slope = 1, 
    color = neutral_mid, linetype = 2, linewidth = 0.5) + geom_point(shape = 21, fill = "#D7A12E", color = neutral_dark, 
    size = 3, stroke = 0.55) + ggrepel::geom_text_repel(aes(label = marker_canonical), size = 2.9, seed = 220224, max.overlaps = Inf, 
    min.segment.length = 0, box.padding = 0.25, point.padding = 0.16, segment.color = "#A3A9AE", segment.size = 0.28) + annotate("text", 
    x = -Inf, y = Inf, hjust = -0.08, vjust = 1.1, label = sprintf("Spearman rho = %.2f; n = %d", cross_cohort_spearman, 
        nrow(shared_global)), size = 3.2) + coord_equal() + labs(title = "Cross-cohort SNR concordance", x = expression("First cohort " * 
    log[2] * " SNR"), y = expression("Second cohort " * log[2] * " SNR")) + theme_qc(11.5) + theme(legend.position = "none"))
  )
)
original_modules[["density"]] <- list(
  setup = expression(
`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || is.na(x[1])) 
        y
    else x
},
bootstrap_spearman_ci <- function(x, y, iterations = 5000L, seed = 1L) {
    set.seed(seed)
    n <- length(x)
    boot <- replicate(iterations, {
        idx <- sample.int(n, n, replace = TRUE)
        if (length(unique(x[idx])) < 2 || length(unique(y[idx])) < 2) {
            return(NA_real_)
        }
        suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
    })
    unname(quantile(boot, c(0.025, 0.975), na.rm = TRUE))
},
spearman_row <- function(data, density_col, score_col, min_met_hi_cells, scope) {
    current <- data %>% filter(n_met_hi_fib >= min_met_hi_cells)
    if (scope == "SSc only") {
        current <- current %>% filter(Progression_skin != "healthy")
    }
    test <- suppressWarnings(cor.test(current[[density_col]], current[[score_col]], method = "spearman", exact = FALSE))
    ci <- bootstrap_spearman_ci(current[[density_col]], current[[score_col]], seed = 1000 + min_met_hi_cells + match(density_col, 
        c("mean_delaunay_neighbors_80um", "mean_distance_5nn_um")) * 10 + match(score_col, c("mean_glycolysis_score", "mean_tca_oxphos_score")) * 
        100)
    data.frame(scope = scope, minimum_met_hi_fib_cells = min_met_hi_cells, density_metric = density_col, metabolic_score = score_col, 
        n_donors = nrow(current), rho = unname(test$estimate), ci_low = ci[1], ci_high = ci[2], p_value = test$p.value)
},
group_test <- function(data, value_col, analysis_name) {
    formula <- reformulate("Progression_skin", response = value_col)
    test <- kruskal.test(formula, data = data)
    data.frame(analysis = analysis_name, metric = value_col, n_healthy = sum(data$Progression_skin == "healthy"), n_progressive = sum(data$Progression_skin == 
        "progressive"), n_stable = sum(data$Progression_skin == "stable"), statistic = unname(test$statistic), df = unname(test$parameter), 
        p_value = test$p.value)
},
group_labels <- c(healthy = "Healthy", progressive = "Progressive SSc", stable = "Stable SSc"),
group_colors <- c(healthy = "#6CC3F4", progressive = "#C03830", stable = "#F6BF93"),
ink <- "#30363B",
grid_color <- "#E5E8EA",
theme_reviewer <- function(base_size = 11) {
    theme_classic(base_size = base_size) + theme(plot.title = element_text(face = "bold", size = rel(1.18)), plot.subtitle = element_text(color = "#4B535A", 
        size = rel(0.92)), plot.caption = element_text(color = "#596168", hjust = 0, size = rel(0.78)), strip.background = element_rect(fill = "#F3F4F5", 
        color = "#AEB4B9"), strip.text = element_text(face = "bold", size = rel(0.93)), panel.grid.major = element_line(color = grid_color, 
        linewidth = 0.3), panel.grid.minor = element_blank(), axis.line = element_line(color = ink, linewidth = 0.4), axis.title = element_text(face = "bold"), 
        legend.position = "top", legend.justification = "left", legend.box = "horizontal", plot.margin = margin(9, 15, 9, 
            9))
},
theme_manuscript_density <- function(base_size = 8.4) {
    theme_bw(base_size = base_size) + theme(plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(), 
        strip.background = element_rect(fill = "#F3F4F5", color = "#AEB4B9", linewidth = 0.35), strip.text = element_text(face = "bold", 
            size = 8.2, lineheight = 0.95), panel.grid.minor = element_blank(), panel.grid.major = element_line(color = grid_color, 
            linewidth = 0.28), panel.border = element_rect(color = ink, linewidth = 0.45), axis.title = element_text(size = 8.6, 
            face = "bold"), axis.text = element_text(size = 7.4, color = ink), axis.ticks = element_line(linewidth = 0.35), 
        axis.ticks.length = grid::unit(1.3, "mm"), legend.position = "top", legend.justification = "center", legend.box = "horizontal", 
        legend.box.just = "center", legend.title = element_text(size = 7.5), legend.text = element_text(size = 7.1), legend.key.height = grid::unit(3.4, 
            "mm"), legend.key.width = grid::unit(4, "mm"), legend.spacing.x = grid::unit(1.3, "mm"), legend.margin = margin(0, 
            0, 2, 0), plot.margin = margin(4, 5, 4, 5, unit = "pt"))
},
save_plot <- function(plot, stem, width, height, dpi = 450) {
    ggsave(file.path(plot_dir, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
    ggsave(file.path(plot_dir, paste0(stem, ".svg")), plot, width = width, height = height, device = svglite::svglite)
    ggsave(file.path(plot_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
}
  ),
  expressions = list(
"%||%" = quote(function(x, y) {
    if (is.null(x) || length(x) == 0 || is.na(x[1])) 
        y
    else x
}),
"bootstrap_spearman_ci" = quote(function(x, y, iterations = 5000L, seed = 1L) {
    set.seed(seed)
    n <- length(x)
    boot <- replicate(iterations, {
        idx <- sample.int(n, n, replace = TRUE)
        if (length(unique(x[idx])) < 2 || length(unique(y[idx])) < 2) {
            return(NA_real_)
        }
        suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
    })
    unname(quantile(boot, c(0.025, 0.975), na.rm = TRUE))
}),
"spearman_row" = quote(function(data, density_col, score_col, min_met_hi_cells, scope) {
    current <- data %>% filter(n_met_hi_fib >= min_met_hi_cells)
    if (scope == "SSc only") {
        current <- current %>% filter(Progression_skin != "healthy")
    }
    test <- suppressWarnings(cor.test(current[[density_col]], current[[score_col]], method = "spearman", exact = FALSE))
    ci <- bootstrap_spearman_ci(current[[density_col]], current[[score_col]], seed = 1000 + min_met_hi_cells + match(density_col, 
        c("mean_delaunay_neighbors_80um", "mean_distance_5nn_um")) * 10 + match(score_col, c("mean_glycolysis_score", "mean_tca_oxphos_score")) * 
        100)
    data.frame(scope = scope, minimum_met_hi_fib_cells = min_met_hi_cells, density_metric = density_col, metabolic_score = score_col, 
        n_donors = nrow(current), rho = unname(test$estimate), ci_low = ci[1], ci_high = ci[2], p_value = test$p.value)
}),
"group_test" = quote(function(data, value_col, analysis_name) {
    formula <- reformulate("Progression_skin", response = value_col)
    test <- kruskal.test(formula, data = data)
    data.frame(analysis = analysis_name, metric = value_col, n_healthy = sum(data$Progression_skin == "healthy"), n_progressive = sum(data$Progression_skin == 
        "progressive"), n_stable = sum(data$Progression_skin == "stable"), statistic = unname(test$statistic), df = unname(test$parameter), 
        p_value = test$p.value)
}),
"group_labels" = quote(c(healthy = "Healthy", progressive = "Progressive SSc", stable = "Stable SSc")),
"group_colors" = quote(c(healthy = "#6CC3F4", progressive = "#C03830", stable = "#F6BF93")),
"ink" = quote("#30363B"),
"grid_color" = quote("#E5E8EA"),
"theme_reviewer" = quote(function(base_size = 11) {
    theme_classic(base_size = base_size) + theme(plot.title = element_text(face = "bold", size = rel(1.18)), plot.subtitle = element_text(color = "#4B535A", 
        size = rel(0.92)), plot.caption = element_text(color = "#596168", hjust = 0, size = rel(0.78)), strip.background = element_rect(fill = "#F3F4F5", 
        color = "#AEB4B9"), strip.text = element_text(face = "bold", size = rel(0.93)), panel.grid.major = element_line(color = grid_color, 
        linewidth = 0.3), panel.grid.minor = element_blank(), axis.line = element_line(color = ink, linewidth = 0.4), axis.title = element_text(face = "bold"), 
        legend.position = "top", legend.justification = "left", legend.box = "horizontal", plot.margin = margin(9, 15, 9, 
            9))
}),
"theme_manuscript_density" = quote(function(base_size = 8.4) {
    theme_bw(base_size = base_size) + theme(plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(), 
        strip.background = element_rect(fill = "#F3F4F5", color = "#AEB4B9", linewidth = 0.35), strip.text = element_text(face = "bold", 
            size = 8.2, lineheight = 0.95), panel.grid.minor = element_blank(), panel.grid.major = element_line(color = grid_color, 
            linewidth = 0.28), panel.border = element_rect(color = ink, linewidth = 0.45), axis.title = element_text(size = 8.6, 
            face = "bold"), axis.text = element_text(size = 7.4, color = ink), axis.ticks = element_line(linewidth = 0.35), 
        axis.ticks.length = grid::unit(1.3, "mm"), legend.position = "top", legend.justification = "center", legend.box = "horizontal", 
        legend.box.just = "center", legend.title = element_text(size = 7.5), legend.text = element_text(size = 7.1), legend.key.height = grid::unit(3.4, 
            "mm"), legend.key.width = grid::unit(4, "mm"), legend.spacing.x = grid::unit(1.3, "mm"), legend.margin = margin(0, 
            0, 2, 0), plot.margin = margin(4, 5, 4, 5, unit = "pt"))
}),
"save_plot" = quote(function(plot, stem, width, height, dpi = 450) {
    ggsave(file.path(plot_dir, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
    ggsave(file.path(plot_dir, paste0(stem, ".svg")), plot, width = width, height = height, device = svglite::svglite)
    ggsave(file.path(plot_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
})
  ),
  plots = list(
"correlation_plot" = quote(ggplot(correlation_plot_data, aes(x = density_value, y = score_value)) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, 
    linewidth = 0.48, linetype = 2, color = "#777F86") + geom_point(aes(fill = Progression_skin, size = n_met_hi_fib), shape = 21, 
    color = ink, stroke = 0.38, alpha = 0.95) + ggrepel::geom_text_repel(aes(label = Donor), size = 2.25, seed = 20260729, 
    max.overlaps = Inf, min.segment.length = 0, box.padding = 0.16, point.padding = 0.11, segment.color = "#A3A9AE", segment.size = 0.2, 
    force = 2.2, force_pull = 0.25, max.iter = 1e+05, max.time = 20) + geom_label(data = correlation_annotations, aes(x = annotation_x, 
    y = annotation_y, label = label, hjust = annotation_hjust, vjust = annotation_vjust), inherit.aes = FALSE, size = 2.3, 
    lineheight = 0.95, fill = "white", color = ink, linewidth = 0.18, label.padding = grid::unit(0.12, "lines"), label.r = grid::unit(0.08, 
        "lines")) + facet_wrap(vars(panel_label), ncol = 2, scales = "free") + scale_fill_manual(values = group_colors, breaks = names(group_labels), 
    labels = unname(group_labels), name = NULL) + scale_size_continuous(range = c(1.45, 3.3), trans = "sqrt", breaks = c(5, 
    20, 100, 500, 1000), name = "Met-hi Fib cells") + guides(size = guide_legend(order = 1, title.position = "left", title.theme = element_text(size = 7.5, 
    vjust = 0.5)), fill = guide_legend(order = 2, title.position = "left", override.aes = list(size = 2.2))) + labs(x = "Local-density metric", 
    y = "Mean metabolic score") + theme_manuscript_density(8.4) + theme(strip.text = element_text(size = 8), panel.spacing = grid::unit(2.6, 
    "mm"), legend.position = "top", legend.justification = "center", legend.box.just = "center", plot.margin = margin(3, 
    4, 3, 10, unit = "pt"))),
"cell_count_plot" = quote(ggplot(count_plot_data, aes(x = Progression_skin, y = value, fill = Progression_skin)) + geom_boxplot(width = 0.56, outlier.shape = NA, 
    alpha = 0.42, color = ink, linewidth = 0.45) + geom_point(position = position_jitter(width = 0.1, height = 0, seed = 20260729), 
    shape = 21, size = 2.35, color = ink, stroke = 0.38) + geom_label(data = count_annotations, aes(x = -Inf, y = Inf, label = label), 
    inherit.aes = FALSE, hjust = -0.04, vjust = 1.08, size = 2.45, fill = "white", color = ink, linewidth = 0.18, label.padding = grid::unit(0.12, 
        "lines"), label.r = grid::unit(0.08, "lines")) + facet_wrap(vars(metric_label), nrow = 1, scales = "free_y") + scale_fill_manual(values = group_colors, 
    guide = "none") + scale_x_discrete(labels = c(healthy = "Healthy", progressive = "Progressive\nSSc", stable = "Stable\nSSc")) + 
    scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0.05, 0.14))) + labs(x = NULL, y = NULL) + 
    theme_manuscript_density(8.4) + theme(axis.text.x = element_text(size = 7.6, angle = 0, hjust = 0.5, vjust = 1, lineheight = 0.92), 
    strip.text = element_text(size = 8.3), panel.spacing = grid::unit(4.2, "mm"), plot.margin = margin(3, 4, 4, 5, unit = "pt"))),
"local_density_plot" = quote(ggplot(local_density_plot_data, aes(x = Progression_skin, y = value, fill = Progression_skin)) + geom_boxplot(width = 0.56, 
    outlier.shape = NA, alpha = 0.35, color = ink, linewidth = 0.45) + geom_point(aes(size = n_met_hi_fib), position = position_jitter(width = 0.1, 
    height = 0, seed = 20260729), shape = 21, color = ink, stroke = 0.45) + geom_label(data = local_density_annotations, 
    aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.05, vjust = 1.1, size = 3, fill = "white", color = ink, 
    linewidth = 0.2) + facet_wrap(vars(metric_label), nrow = 1, scales = "free_y") + scale_fill_manual(values = group_colors, 
    guide = "none") + scale_x_discrete(labels = group_labels) + scale_size_continuous(range = c(2.3, 4.8), trans = "sqrt", 
    breaks = c(5, 20, 100, 500, 1000), name = "Met-hi Fib cells") + labs(title = "Local cellularity around Met-hi Fib by clinical group", 
    subtitle = "Each point is one donor; point size indicates contributing Met-hi Fib cells", x = NULL, y = NULL) + theme_reviewer(10.8) + 
    theme(axis.text.x = element_text(angle = 20, hjust = 1), panel.spacing = unit(1.2, "lines"))),
"threshold_plot" = quote(ggplot(threshold_plot_data, aes(x = minimum_label, y = comparison, fill = rho)) + geom_tile(color = "white", linewidth = 0.7) + 
    geom_text(aes(label = cell_label), size = 3, lineheight = 0.95, color = ink) + scale_fill_gradient2(low = "#3F78A8", 
    mid = "white", high = "#D97706", midpoint = 0, limits = c(-1, 1), name = "Spearman \317\201") + labs(title = "Sensitivity to minimum Met-hi Fib cell counts", 
    subtitle = "Each column repeats the donor-level analysis after excluding sparsely represented donor means", x = "Minimum Met-hi Fib cells per donor", 
    y = NULL, caption = "Raw two-sided Spearman p-values are displayed; Holm-adjusted values are available in the source table.") + 
    theme_reviewer(10.5) + theme(legend.position = "right", panel.grid = element_blank(), axis.text.x = element_text(angle = 25, 
    hjust = 1)))
  )
)
original_modules[["ec"]] <- list(
  setup = expression(
ec_lv2_levels <- c("ACKR1_EC", "ACTA2_EC", "cycling_EC", "EPC", "HEY1_EC", "lymphatic_EC"),
score_levels <- c("CoreMatrisome score", "Collagens score", "ECM Glycoproteins score", "Proteoglycans score"),
score_labels <- c(`CoreMatrisome score` = "Core\nmatrisome", `Collagens score` = "Collagens", `ECM Glycoproteins score` = "ECM\nglycoproteins", 
    `Proteoglycans score` = "Proteoglycans"),
ec_colors <- c(ACKR1_EC = "#4C78A8", ACTA2_EC = "#E45756", cycling_EC = "#72B7B2", EPC = "#F58518", HEY1_EC = "#54A24B", 
    lymphatic_EC = "#B279A2"),
to_dgC <- function(x) {
    if (inherits(x, "dgCMatrix")) 
        return(x)
    x <- as(x, "CsparseMatrix")
    as(x, "dgCMatrix")
},
rescale01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) {
        return(rep(0.5, length(x)))
    }
    (x - rng[1])/(rng[2] - rng[1])
},
format_p <- function(p_value) {
    dplyr::case_when(is.na(p_value) ~ "p = NA", p_value < 0.001 ~ "p < 0.001", TRUE ~ paste0("p = ", formatC(p_value, format = "f", 
        digits = 3)))
},
save_plot_set <- function(plot, basename, width, height) {
    png_file <- file.path(output_dir, paste0(basename, ".png"))
    pdf_file <- file.path(output_dir, paste0(basename, ".pdf"))
    svg_file <- file.path(output_dir, paste0(basename, ".svg"))
    ggplot2::ggsave(png_file, plot, width = width, height = height, dpi = 300, bg = "white")
    grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE, bg = "white")
    print(plot)
    grDevices::dev.off()
    grDevices::svg(svg_file, width = width, height = height, bg = "white")
    print(plot)
    grDevices::dev.off()
    invisible(c(png_file, pdf_file, svg_file))
},
make_ecm_dotplot <- function(summary_df, scale_by_score = FALSE) {
    plot_df <- summary_df %>% dplyr::mutate(score_name = factor(as.character(score_name), levels = rev(score_levels)), lv2_anno_new = factor(as.character(lv2_anno_new), 
        levels = ec_lv2_levels))
    if (isTRUE(scale_by_score)) {
        plot_df <- plot_df %>% dplyr::group_by(score_name) %>% dplyr::mutate(fill_score = rescale01(mean_score), size_score = fill_score) %>% 
            dplyr::ungroup()
        legend_title <- "Score-scaled mean"
        plot_title <- "EC ECM scores across Xenium EC annotations, scaled"
    }
    else {
        plot_df <- plot_df %>% dplyr::mutate(fill_score = mean_score, size_score = rescale01(mean_score))
        legend_title <- "Mean AUCell score"
        plot_title <- "EC ECM scores across Xenium EC annotations"
    }
    ggplot(plot_df, aes(x = lv2_anno_new, y = score_name)) + geom_point(aes(fill = fill_score, size = size_score), shape = 21, 
        color = "grey35", stroke = 0.35) + scale_y_discrete(labels = score_labels) + scale_fill_gradient(low = "white", high = "#D7301F", 
        name = legend_title) + scale_size(range = c(4.5, 11.5), guide = "none") + labs(title = plot_title, x = NULL, y = NULL) + 
        theme_bw(base_size = 14) + theme(axis.text.x = element_text(size = 12, face = "bold", angle = 35, hjust = 1), axis.text.y = element_text(size = 12.5, 
        lineheight = 0.95), legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 10.5), 
        panel.grid.major = element_line(color = "grey88"), panel.grid.minor = element_blank(), plot.title = element_text(size = 15.5, 
            hjust = 0.5, margin = margin(b = 6)), plot.margin = margin(5.5, 6, 5.5, 5.5))
}
  ),
  expressions = list(
"ec_lv2_levels" = quote(c("ACKR1_EC", "ACTA2_EC", "cycling_EC", "EPC", "HEY1_EC", "lymphatic_EC")),
"score_levels" = quote(c("CoreMatrisome score", "Collagens score", "ECM Glycoproteins score", "Proteoglycans score")),
"score_labels" = quote(c(`CoreMatrisome score` = "Core\nmatrisome", `Collagens score` = "Collagens", `ECM Glycoproteins score` = "ECM\nglycoproteins", 
    `Proteoglycans score` = "Proteoglycans")),
"ec_colors" = quote(c(ACKR1_EC = "#4C78A8", ACTA2_EC = "#E45756", cycling_EC = "#72B7B2", EPC = "#F58518", HEY1_EC = "#54A24B", lymphatic_EC = "#B279A2")),
"to_dgC" = quote(function(x) {
    if (inherits(x, "dgCMatrix")) 
        return(x)
    x <- as(x, "CsparseMatrix")
    as(x, "dgCMatrix")
}),
"rescale01" = quote(function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) {
        return(rep(0.5, length(x)))
    }
    (x - rng[1])/(rng[2] - rng[1])
}),
"format_p" = quote(function(p_value) {
    dplyr::case_when(is.na(p_value) ~ "p = NA", p_value < 0.001 ~ "p < 0.001", TRUE ~ paste0("p = ", formatC(p_value, format = "f", 
        digits = 3)))
}),
"save_plot_set" = quote(function(plot, basename, width, height) {
    png_file <- file.path(output_dir, paste0(basename, ".png"))
    pdf_file <- file.path(output_dir, paste0(basename, ".pdf"))
    svg_file <- file.path(output_dir, paste0(basename, ".svg"))
    ggplot2::ggsave(png_file, plot, width = width, height = height, dpi = 300, bg = "white")
    grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE, bg = "white")
    print(plot)
    grDevices::dev.off()
    grDevices::svg(svg_file, width = width, height = height, bg = "white")
    print(plot)
    grDevices::dev.off()
    invisible(c(png_file, pdf_file, svg_file))
}),
"make_ecm_dotplot" = quote(function(summary_df, scale_by_score = FALSE) {
    plot_df <- summary_df %>% dplyr::mutate(score_name = factor(as.character(score_name), levels = rev(score_levels)), lv2_anno_new = factor(as.character(lv2_anno_new), 
        levels = ec_lv2_levels))
    if (isTRUE(scale_by_score)) {
        plot_df <- plot_df %>% dplyr::group_by(score_name) %>% dplyr::mutate(fill_score = rescale01(mean_score), size_score = fill_score) %>% 
            dplyr::ungroup()
        legend_title <- "Score-scaled mean"
        plot_title <- "EC ECM scores across Xenium EC annotations, scaled"
    }
    else {
        plot_df <- plot_df %>% dplyr::mutate(fill_score = mean_score, size_score = rescale01(mean_score))
        legend_title <- "Mean AUCell score"
        plot_title <- "EC ECM scores across Xenium EC annotations"
    }
    ggplot(plot_df, aes(x = lv2_anno_new, y = score_name)) + geom_point(aes(fill = fill_score, size = size_score), shape = 21, 
        color = "grey35", stroke = 0.35) + scale_y_discrete(labels = score_labels) + scale_fill_gradient(low = "white", high = "#D7301F", 
        name = legend_title) + scale_size(range = c(4.5, 11.5), guide = "none") + labs(title = plot_title, x = NULL, y = NULL) + 
        theme_bw(base_size = 14) + theme(axis.text.x = element_text(size = 12, face = "bold", angle = 35, hjust = 1), axis.text.y = element_text(size = 12.5, 
        lineheight = 0.95), legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 10.5), 
        panel.grid.major = element_line(color = "grey88"), panel.grid.minor = element_blank(), plot.title = element_text(size = 15.5, 
            hjust = 0.5, margin = margin(b = 6)), plot.margin = margin(5.5, 6, 5.5, 5.5))
})
  ),
  plots = list(
"p_heatmap_raw" = quote(ggplot(heatmap_raw_df, aes(x = lv2_anno_new, y = score_name, fill = mean_score)) + geom_tile(color = "white", linewidth = 0.7) + 
    scale_y_discrete(labels = score_labels) + scale_fill_gradient(low = "white", high = "#D7301F", name = "Mean AUCell") + 
    labs(title = "Mean EC ECM AUCell scores", x = NULL, y = NULL) + theme_bw(base_size = 14) + theme(axis.text.x = element_text(size = 12, 
    face = "bold", angle = 35, hjust = 1), axis.text.y = element_text(size = 12.5, lineheight = 0.95), panel.grid = element_blank(), 
    plot.title = element_text(size = 15.5, hjust = 0.5))),
"p_heatmap_scaled" = quote(ggplot(heatmap_scaled_df, aes(x = lv2_anno_new, y = score_name, fill = mean_z)) + geom_tile(color = "white", linewidth = 0.7) + 
    scale_y_discrete(labels = score_labels) + scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#D7301F", midpoint = 0, 
    name = "Row z") + labs(title = "Mean EC ECM AUCell scores, row-scaled", x = NULL, y = NULL) + theme_bw(base_size = 14) + 
    theme(axis.text.x = element_text(size = 12, face = "bold", angle = 35, hjust = 1), axis.text.y = element_text(size = 12.5, 
        lineheight = 0.95), panel.grid = element_blank(), plot.title = element_text(size = 15.5, hjust = 0.5))),
"p_violin" = quote(score_long %>% dplyr::mutate(score_name = factor(as.character(score_name), levels = score_levels), lv2_anno_new = factor(as.character(lv2_anno_new), 
    levels = ec_lv2_levels)) %>% ggplot(aes(x = lv2_anno_new, y = score, fill = lv2_anno_new)) + geom_violin(width = 0.86, 
    alpha = 0.72, color = "grey30", trim = TRUE) + geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "grey15", 
    linewidth = 0.4) + facet_wrap(~score_name, ncol = 2, scales = "free_y", labeller = labeller(score_name = score_labels)) + 
    scale_fill_manual(values = ec_colors) + labs(title = "Single-cell EC ECM score distributions", x = NULL, y = "AUCell score", 
    fill = NULL) + theme_bw(base_size = 14) + theme(axis.text.x = element_text(size = 11, face = "bold", angle = 35, hjust = 1), 
    strip.text = element_text(size = 12.5, face = "bold"), legend.position = "none", plot.title = element_text(size = 16, 
        hjust = 0.5))),
"p_density" = quote(score_long %>% dplyr::mutate(score_name = factor(as.character(score_name), levels = score_levels)) %>% ggplot(aes(x = score, 
    color = lv2_anno_new, fill = lv2_anno_new)) + geom_density(alpha = 0.13, linewidth = 0.75, adjust = 1.05) + facet_wrap(~score_name, 
    ncol = 2, scales = "free", labeller = labeller(score_name = score_labels)) + scale_color_manual(values = ec_colors) + 
    scale_fill_manual(values = ec_colors) + labs(title = "EC ECM score density by Xenium annotation", x = "AUCell score", 
    y = "Density", color = NULL, fill = NULL) + theme_bw(base_size = 14) + theme(strip.text = element_text(size = 12.5, face = "bold"), 
    legend.position = "bottom", plot.title = element_text(size = 16, hjust = 0.5)))
  )
)
original_modules[["fib"]] <- list(
  setup = expression(
lv3_preferred_order <- c("papillary_Fib", "PI16_Fib", "CXCL12_Fib", "COL8A1_Fib", "CCL19_Fib", "COMP_Fib", "COCH_Fib", "NGFR_Fib"),
metfib_levels <- c("Other_Fib", "Met_hi_Fib"),
metfib_colors <- c(Other_Fib = "#4C78A8", Met_hi_Fib = "#E45756"),
score_levels <- c("CoreMatrisome score", "Collagens score", "ECM Glycoproteins score", "Proteoglycans score"),
score_labels <- c(`CoreMatrisome score` = "Core\nmatrisome", `Collagens score` = "Collagens", `ECM Glycoproteins score` = "ECM\nglycoproteins", 
    `Proteoglycans score` = "Proteoglycans"),
to_dgC <- function(x) {
    if (inherits(x, "dgCMatrix")) 
        return(x)
    x <- as(x, "CsparseMatrix")
    as(x, "dgCMatrix")
},
rescale01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) {
        return(rep(0.5, length(x)))
    }
    (x - rng[1])/(rng[2] - rng[1])
},
format_p <- function(p_value) {
    dplyr::case_when(is.na(p_value) ~ "p = NA", p_value < 0.001 ~ "p < 0.001", TRUE ~ paste0("p = ", formatC(p_value, format = "f", 
        digits = 3)))
},
save_plot_set <- function(plot, basename, width, height) {
    png_file <- file.path(output_dir, paste0(basename, ".png"))
    pdf_file <- file.path(output_dir, paste0(basename, ".pdf"))
    svg_file <- file.path(output_dir, paste0(basename, ".svg"))
    ggplot2::ggsave(png_file, plot, width = width, height = height, dpi = 300, bg = "white")
    grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE, bg = "white")
    print(plot)
    grDevices::dev.off()
    grDevices::svg(svg_file, width = width, height = height, bg = "white")
    print(plot)
    grDevices::dev.off()
    invisible(c(png_file, pdf_file, svg_file))
},
make_split_dotplot <- function(summary_df, scale_by_score = FALSE) {
    plot_df <- summary_df %>% dplyr::mutate(score_name = factor(as.character(score_name), levels = rev(score_levels)), lv3_anno = factor(as.character(lv3_anno), 
        levels = lv3_levels), metfiblabel = factor(as.character(metfiblabel), levels = metfib_levels))
    if (isTRUE(scale_by_score)) {
        plot_df <- plot_df %>% dplyr::group_by(score_name) %>% dplyr::mutate(fill_score = rescale01(mean_score), size_score = fill_score) %>% 
            dplyr::ungroup()
        legend_title <- "Score-scaled mean"
        plot_title <- "Fibroblast ECM scores by Xenium label and IMC metfiblabel, scaled"
    }
    else {
        plot_df <- plot_df %>% dplyr::mutate(fill_score = mean_score, size_score = rescale01(mean_score))
        legend_title <- "Mean AUCell score"
        plot_title <- "Fibroblast ECM scores by Xenium label and IMC metfiblabel"
    }
    ggplot(plot_df, aes(x = lv3_anno, y = score_name)) + geom_point(aes(fill = fill_score, size = size_score, color = metfiblabel), 
        shape = 21, stroke = 0.85, position = position_dodge(width = 0.68)) + scale_y_discrete(labels = score_labels) + scale_fill_gradient(low = "white", 
        high = "#D7301F", name = legend_title) + scale_color_manual(values = metfib_colors, name = NULL) + scale_size(range = c(3.8, 
        10.5), guide = "none") + labs(title = plot_title, x = NULL, y = NULL) + theme_bw(base_size = 14) + theme(axis.text.x = element_text(size = 11.5, 
        face = "bold", angle = 35, hjust = 1), axis.text.y = element_text(size = 12.5, lineheight = 0.95), legend.position = "right", 
        legend.title = element_text(size = 12), legend.text = element_text(size = 10.5), panel.grid.major = element_line(color = "grey88"), 
        panel.grid.minor = element_blank(), plot.title = element_text(size = 15.5, hjust = 0.5, margin = margin(b = 6)), 
        plot.margin = margin(5.5, 6, 5.5, 5.5))
},
make_presentation_mean_dot_matrix <- function(summary_df, scale_by_score = FALSE) {
    plot_df <- summary_df %>% dplyr::mutate(score_name = factor(as.character(score_name), levels = score_levels), lv3_anno = factor(as.character(lv3_anno), 
        levels = rev(lv3_levels)), metfiblabel = factor(as.character(metfiblabel), levels = metfib_levels))
    if (isTRUE(scale_by_score)) {
        plot_df <- plot_df %>% dplyr::group_by(score_name) %>% dplyr::mutate(fill_score = rescale01(mean_score), size_score = fill_score) %>% 
            dplyr::ungroup()
        legend_title <- "Scaled mean\nwithin score"
        plot_title <- "ECM score presentation matrix, scaled"
    }
    else {
        plot_df <- plot_df %>% dplyr::mutate(fill_score = mean_score, size_score = mean_score)
        legend_title <- "Mean\nAUCell score"
        plot_title <- "ECM score presentation matrix"
    }
    ggplot(plot_df, aes(x = metfiblabel, y = lv3_anno)) + geom_point(aes(fill = fill_score, size = size_score), shape = 21, 
        color = "grey25", stroke = 0.35) + facet_wrap(~score_name, ncol = 2, labeller = labeller(score_name = score_labels)) + 
        scale_x_discrete(labels = c(Other_Fib = "Other", Met_hi_Fib = "Met_hi")) + scale_fill_gradient(low = "white", high = "#B30000", 
        name = legend_title) + scale_size(range = c(5, 13), guide = "none") + labs(title = plot_title, x = NULL, y = NULL) + 
        theme_bw(base_size = 16) + theme(panel.grid = element_blank(), axis.text.x = element_text(size = 14, face = "bold"), 
        axis.text.y = element_text(size = 13.5, face = "bold"), strip.background = element_rect(fill = "grey95", color = "grey35"), 
        strip.text = element_text(size = 13.5, face = "bold"), legend.title = element_text(size = 12.5), legend.text = element_text(size = 11.5), 
        legend.position = "right", plot.title = element_text(size = 17, hjust = 0.5), plot.margin = margin(7, 7, 7, 7))
},
make_presentation_delta_dot_matrix <- function(delta_df, scale_by_score = FALSE) {
    plot_df <- delta_df %>% dplyr::mutate(score_name = factor(as.character(score_name), levels = score_levels), lv3_anno = factor(as.character(lv3_anno), 
        levels = rev(lv3_levels)))
    if (isTRUE(scale_by_score)) {
        plot_df <- plot_df %>% dplyr::group_by(score_name) %>% dplyr::mutate(max_abs_delta = max(abs(mean_delta_met_hi_vs_other), 
            na.rm = TRUE), delta_fill = dplyr::if_else(is.finite(max_abs_delta) & max_abs_delta > 0, mean_delta_met_hi_vs_other/max_abs_delta, 
            0), size_score = abs(delta_fill)) %>% dplyr::ungroup()
        legend_title <- "Scaled delta\nwithin score"
        plot_title <- "ECM score delta presentation matrix, scaled"
        fill_limits <- c(-1, 1)
    }
    else {
        plot_df <- plot_df %>% dplyr::mutate(delta_fill = mean_delta_met_hi_vs_other, size_score = abs(mean_delta_met_hi_vs_other))
        legend_title <- "Mean delta\nMet_hi - Other"
        plot_title <- "ECM score delta presentation matrix"
        max_abs <- max(abs(plot_df$delta_fill), na.rm = TRUE)
        fill_limits <- c(-max_abs, max_abs)
    }
    ggplot(plot_df, aes(x = score_name, y = lv3_anno)) + geom_point(aes(fill = delta_fill, size = size_score), shape = 21, 
        color = "grey25", stroke = 0.35) + scale_x_discrete(labels = score_labels) + scale_fill_gradient2(low = "#3B6FB6", 
        mid = "white", high = "#B30000", midpoint = 0, limits = fill_limits, name = legend_title) + scale_size(range = c(4.5, 
        13), guide = "none") + labs(title = plot_title, x = NULL, y = NULL) + theme_bw(base_size = 16) + theme(panel.grid = element_blank(), 
        axis.text.x = element_text(size = 12.5, face = "bold", lineheight = 0.95), axis.text.y = element_text(size = 13.5, 
            face = "bold"), legend.title = element_text(size = 12.5), legend.text = element_text(size = 11.5), legend.position = "right", 
        plot.title = element_text(size = 17, hjust = 0.5), plot.margin = margin(7, 7, 7, 7))
}
  ),
  expressions = list(
"lv3_preferred_order" = quote(c("papillary_Fib", "PI16_Fib", "CXCL12_Fib", "COL8A1_Fib", "CCL19_Fib", "COMP_Fib", "COCH_Fib", "NGFR_Fib")),
"metfib_levels" = quote(c("Other_Fib", "Met_hi_Fib")),
"metfib_colors" = quote(c(Other_Fib = "#4C78A8", Met_hi_Fib = "#E45756")),
"score_levels" = quote(c("CoreMatrisome score", "Collagens score", "ECM Glycoproteins score", "Proteoglycans score")),
"score_labels" = quote(c(`CoreMatrisome score` = "Core\nmatrisome", `Collagens score` = "Collagens", `ECM Glycoproteins score` = "ECM\nglycoproteins", 
    `Proteoglycans score` = "Proteoglycans")),
"to_dgC" = quote(function(x) {
    if (inherits(x, "dgCMatrix")) 
        return(x)
    x <- as(x, "CsparseMatrix")
    as(x, "dgCMatrix")
}),
"rescale01" = quote(function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) {
        return(rep(0.5, length(x)))
    }
    (x - rng[1])/(rng[2] - rng[1])
}),
"format_p" = quote(function(p_value) {
    dplyr::case_when(is.na(p_value) ~ "p = NA", p_value < 0.001 ~ "p < 0.001", TRUE ~ paste0("p = ", formatC(p_value, format = "f", 
        digits = 3)))
}),
"save_plot_set" = quote(function(plot, basename, width, height) {
    png_file <- file.path(output_dir, paste0(basename, ".png"))
    pdf_file <- file.path(output_dir, paste0(basename, ".pdf"))
    svg_file <- file.path(output_dir, paste0(basename, ".svg"))
    ggplot2::ggsave(png_file, plot, width = width, height = height, dpi = 300, bg = "white")
    grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE, bg = "white")
    print(plot)
    grDevices::dev.off()
    grDevices::svg(svg_file, width = width, height = height, bg = "white")
    print(plot)
    grDevices::dev.off()
    invisible(c(png_file, pdf_file, svg_file))
}),
"make_split_dotplot" = quote(function(summary_df, scale_by_score = FALSE) {
    plot_df <- summary_df %>% dplyr::mutate(score_name = factor(as.character(score_name), levels = rev(score_levels)), lv3_anno = factor(as.character(lv3_anno), 
        levels = lv3_levels), metfiblabel = factor(as.character(metfiblabel), levels = metfib_levels))
    if (isTRUE(scale_by_score)) {
        plot_df <- plot_df %>% dplyr::group_by(score_name) %>% dplyr::mutate(fill_score = rescale01(mean_score), size_score = fill_score) %>% 
            dplyr::ungroup()
        legend_title <- "Score-scaled mean"
        plot_title <- "Fibroblast ECM scores by Xenium label and IMC metfiblabel, scaled"
    }
    else {
        plot_df <- plot_df %>% dplyr::mutate(fill_score = mean_score, size_score = rescale01(mean_score))
        legend_title <- "Mean AUCell score"
        plot_title <- "Fibroblast ECM scores by Xenium label and IMC metfiblabel"
    }
    ggplot(plot_df, aes(x = lv3_anno, y = score_name)) + geom_point(aes(fill = fill_score, size = size_score, color = metfiblabel), 
        shape = 21, stroke = 0.85, position = position_dodge(width = 0.68)) + scale_y_discrete(labels = score_labels) + scale_fill_gradient(low = "white", 
        high = "#D7301F", name = legend_title) + scale_color_manual(values = metfib_colors, name = NULL) + scale_size(range = c(3.8, 
        10.5), guide = "none") + labs(title = plot_title, x = NULL, y = NULL) + theme_bw(base_size = 14) + theme(axis.text.x = element_text(size = 11.5, 
        face = "bold", angle = 35, hjust = 1), axis.text.y = element_text(size = 12.5, lineheight = 0.95), legend.position = "right", 
        legend.title = element_text(size = 12), legend.text = element_text(size = 10.5), panel.grid.major = element_line(color = "grey88"), 
        panel.grid.minor = element_blank(), plot.title = element_text(size = 15.5, hjust = 0.5, margin = margin(b = 6)), 
        plot.margin = margin(5.5, 6, 5.5, 5.5))
}),
"make_presentation_mean_dot_matrix" = quote(function(summary_df, scale_by_score = FALSE) {
    plot_df <- summary_df %>% dplyr::mutate(score_name = factor(as.character(score_name), levels = score_levels), lv3_anno = factor(as.character(lv3_anno), 
        levels = rev(lv3_levels)), metfiblabel = factor(as.character(metfiblabel), levels = metfib_levels))
    if (isTRUE(scale_by_score)) {
        plot_df <- plot_df %>% dplyr::group_by(score_name) %>% dplyr::mutate(fill_score = rescale01(mean_score), size_score = fill_score) %>% 
            dplyr::ungroup()
        legend_title <- "Scaled mean\nwithin score"
        plot_title <- "ECM score presentation matrix, scaled"
    }
    else {
        plot_df <- plot_df %>% dplyr::mutate(fill_score = mean_score, size_score = mean_score)
        legend_title <- "Mean\nAUCell score"
        plot_title <- "ECM score presentation matrix"
    }
    ggplot(plot_df, aes(x = metfiblabel, y = lv3_anno)) + geom_point(aes(fill = fill_score, size = size_score), shape = 21, 
        color = "grey25", stroke = 0.35) + facet_wrap(~score_name, ncol = 2, labeller = labeller(score_name = score_labels)) + 
        scale_x_discrete(labels = c(Other_Fib = "Other", Met_hi_Fib = "Met_hi")) + scale_fill_gradient(low = "white", high = "#B30000", 
        name = legend_title) + scale_size(range = c(5, 13), guide = "none") + labs(title = plot_title, x = NULL, y = NULL) + 
        theme_bw(base_size = 16) + theme(panel.grid = element_blank(), axis.text.x = element_text(size = 14, face = "bold"), 
        axis.text.y = element_text(size = 13.5, face = "bold"), strip.background = element_rect(fill = "grey95", color = "grey35"), 
        strip.text = element_text(size = 13.5, face = "bold"), legend.title = element_text(size = 12.5), legend.text = element_text(size = 11.5), 
        legend.position = "right", plot.title = element_text(size = 17, hjust = 0.5), plot.margin = margin(7, 7, 7, 7))
}),
"make_presentation_delta_dot_matrix" = quote(function(delta_df, scale_by_score = FALSE) {
    plot_df <- delta_df %>% dplyr::mutate(score_name = factor(as.character(score_name), levels = score_levels), lv3_anno = factor(as.character(lv3_anno), 
        levels = rev(lv3_levels)))
    if (isTRUE(scale_by_score)) {
        plot_df <- plot_df %>% dplyr::group_by(score_name) %>% dplyr::mutate(max_abs_delta = max(abs(mean_delta_met_hi_vs_other), 
            na.rm = TRUE), delta_fill = dplyr::if_else(is.finite(max_abs_delta) & max_abs_delta > 0, mean_delta_met_hi_vs_other/max_abs_delta, 
            0), size_score = abs(delta_fill)) %>% dplyr::ungroup()
        legend_title <- "Scaled delta\nwithin score"
        plot_title <- "ECM score delta presentation matrix, scaled"
        fill_limits <- c(-1, 1)
    }
    else {
        plot_df <- plot_df %>% dplyr::mutate(delta_fill = mean_delta_met_hi_vs_other, size_score = abs(mean_delta_met_hi_vs_other))
        legend_title <- "Mean delta\nMet_hi - Other"
        plot_title <- "ECM score delta presentation matrix"
        max_abs <- max(abs(plot_df$delta_fill), na.rm = TRUE)
        fill_limits <- c(-max_abs, max_abs)
    }
    ggplot(plot_df, aes(x = score_name, y = lv3_anno)) + geom_point(aes(fill = delta_fill, size = size_score), shape = 21, 
        color = "grey25", stroke = 0.35) + scale_x_discrete(labels = score_labels) + scale_fill_gradient2(low = "#3B6FB6", 
        mid = "white", high = "#B30000", midpoint = 0, limits = fill_limits, name = legend_title) + scale_size(range = c(4.5, 
        13), guide = "none") + labs(title = plot_title, x = NULL, y = NULL) + theme_bw(base_size = 16) + theme(panel.grid = element_blank(), 
        axis.text.x = element_text(size = 12.5, face = "bold", lineheight = 0.95), axis.text.y = element_text(size = 13.5, 
            face = "bold"), legend.title = element_text(size = 12.5), legend.text = element_text(size = 11.5), legend.position = "right", 
        plot.title = element_text(size = 17, hjust = 0.5), plot.margin = margin(7, 7, 7, 7))
})
  ),
  plots = list(
"p_delta_heatmap" = quote(delta_summary %>% dplyr::mutate(score_name = factor(as.character(score_name), levels = rev(score_levels)), lv3_anno = factor(as.character(lv3_anno), 
    levels = lv3_levels)) %>% ggplot(aes(x = lv3_anno, y = score_name, fill = mean_delta_met_hi_vs_other)) + geom_tile(color = "white", 
    linewidth = 0.7) + scale_y_discrete(labels = score_labels) + scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#D7301F", 
    midpoint = 0, name = "Mean delta\nMet_hi - Other") + labs(title = "ECM score delta by Xenium fibroblast label", x = NULL, 
    y = NULL) + theme_bw(base_size = 14) + theme(axis.text.x = element_text(size = 11.5, face = "bold", angle = 35, hjust = 1), 
    axis.text.y = element_text(size = 12.5, lineheight = 0.95), panel.grid = element_blank(), plot.title = element_text(size = 15.5, 
        hjust = 0.5))),
"p_delta_bar" = quote(delta_summary %>% dplyr::mutate(score_name = factor(as.character(score_name), levels = score_levels), lv3_anno = factor(as.character(lv3_anno), 
    levels = lv3_levels)) %>% ggplot(aes(x = lv3_anno, y = mean_delta_met_hi_vs_other, fill = mean_delta_met_hi_vs_other > 
    0)) + geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) + geom_col(width = 0.72, color = "grey30", linewidth = 0.25) + 
    facet_wrap(~score_name, ncol = 2, scales = "free_y", labeller = labeller(score_name = score_labels)) + scale_fill_manual(values = c(`TRUE` = "#D7301F", 
    `FALSE` = "#3B6FB6"), guide = "none") + labs(title = "Mean ECM score delta, Met_hi_Fib minus Other_Fib", x = NULL, y = "Mean AUCell delta") + 
    theme_bw(base_size = 14) + theme(axis.text.x = element_text(size = 10.5, face = "bold", angle = 35, hjust = 1), strip.text = element_text(size = 12.5, 
    face = "bold"), plot.title = element_text(size = 15.5, hjust = 0.5))),
"p_violin" = quote(score_long %>% dplyr::mutate(score_name = factor(as.character(score_name), levels = score_levels), lv3_anno = factor(as.character(lv3_anno), 
    levels = lv3_levels), metfiblabel = factor(as.character(metfiblabel), levels = metfib_levels)) %>% ggplot(aes(x = lv3_anno, 
    y = score, fill = metfiblabel)) + geom_violin(position = position_dodge(width = 0.82), width = 0.78, alpha = 0.72, color = "grey30", 
    trim = TRUE) + geom_boxplot(position = position_dodge(width = 0.82), width = 0.14, outlier.shape = NA, fill = "white", 
    color = "grey15", linewidth = 0.35) + facet_wrap(~score_name, ncol = 2, scales = "free_y", labeller = labeller(score_name = score_labels)) + 
    scale_fill_manual(values = metfib_colors) + labs(title = "Single-cell fibroblast ECM scores split by IMC metfiblabel", 
    x = NULL, y = "AUCell score", fill = NULL) + theme_bw(base_size = 14) + theme(axis.text.x = element_text(size = 10.5, 
    face = "bold", angle = 35, hjust = 1), strip.text = element_text(size = 12.5, face = "bold"), legend.position = "top", 
    plot.title = element_text(size = 16, hjust = 0.5))),
"p_density" = quote(score_long %>% dplyr::mutate(score_name = factor(as.character(score_name), levels = score_levels)) %>% ggplot(aes(x = score, 
    color = metfiblabel, fill = metfiblabel)) + geom_density(alpha = 0.18, linewidth = 0.8, adjust = 1.05) + facet_grid(score_name ~ 
    lv3_anno, scales = "free", labeller = labeller(score_name = score_labels)) + scale_color_manual(values = metfib_colors) + 
    scale_fill_manual(values = metfib_colors) + labs(title = "ECM score density within each Xenium fibroblast label", x = "AUCell score", 
    y = "Density", color = NULL, fill = NULL) + theme_bw(base_size = 12) + theme(axis.text.x = element_text(size = 8, angle = 35, 
    hjust = 1), axis.text.y = element_text(size = 8), strip.text = element_text(size = 9.5, face = "bold"), legend.position = "top", 
    plot.title = element_text(size = 15.5, hjust = 0.5)))
  )
)
original_modules[["energy"]] <- list(
  setup = expression(
make_imc_protein_score <- function(markers) {
    markers <- intersect(markers, rownames(spe))
    colMeans(as.matrix(SummarizedExperiment::assay(spe, score_assay)[markers, , drop = FALSE]), na.rm = TRUE)
},
recode_group <- function(group, config) {
    group <- as.character(group)
    if (!is.null(config$group_recode)) {
        idx <- match(group, names(config$group_recode))
        matched <- !is.na(idx)
        group[matched] <- unname(config$group_recode[idx[matched]])
    }
    group
},
rescale01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) {
        return(rep(0.5, length(x)))
    }
    (x - rng[1])/(rng[2] - rng[1])
},
format_p <- function(p_value) {
    dplyr::case_when(is.na(p_value) ~ "p = NA", p_value < 0.001 ~ "p < 0.001", TRUE ~ paste0("p = ", formatC(p_value, format = "f", 
        digits = 3)))
},
save_plot_set <- function(plot, basename, width, height) {
    png_file <- file.path(output_dir, paste0(basename, ".png"))
    pdf_file <- file.path(output_dir, paste0(basename, ".pdf"))
    svg_file <- file.path(output_dir, paste0(basename, ".svg"))
    ggsave(png_file, plot, width = width, height = height, dpi = 300, bg = "white")
    grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE, bg = "white")
    print(plot)
    grDevices::dev.off()
    grDevices::svg(svg_file, width = width, height = height, bg = "white")
    print(plot)
    grDevices::dev.off()
    invisible(c(png_file, pdf_file, svg_file))
},
summarise_score_long <- function(score_long) {
    score_long %>% dplyr::group_by(pathway, group) %>% dplyr::summarise(n_cells = dplyr::n(), n_samples = dplyr::n_distinct(sample_id[!is.na(sample_id) & 
        sample_id != ""]), mean_score = mean(score, na.rm = TRUE), median_score = median(score, na.rm = TRUE), q25 = stats::quantile(score, 
        0.25, na.rm = TRUE), q75 = stats::quantile(score, 0.75, na.rm = TRUE), sd_score = stats::sd(score, na.rm = TRUE), 
        se_score = sd_score/sqrt(n_cells), ci95 = 1.96 * se_score, .groups = "drop")
},
summarise_sample_scores <- function(score_long) {
    score_long %>% dplyr::filter(!is.na(sample_id), sample_id != "") %>% dplyr::group_by(sample_id, pathway, group) %>% dplyr::summarise(n_cells = dplyr::n(), 
        mean_score = mean(score, na.rm = TRUE), median_score = median(score, na.rm = TRUE), .groups = "drop")
},
make_imc_score_long <- function(config) {
    score_wide_all %>% dplyr::mutate(group = recode_group(.data[[config$group_col]], config)) %>% dplyr::filter(!is.na(group), 
        group %in% config$group_levels) %>% tidyr::pivot_longer(cols = dplyr::all_of(pathway_levels), names_to = "pathway", 
        values_to = "score") %>% dplyr::mutate(group = factor(group, levels = config$group_levels), pathway = factor(pathway, 
        levels = pathway_levels), source = "IMC protein") %>% dplyr::select(cell, sample_id, group, pathway, score, source)
},
load_scmetab_score_long <- function(config) {
    score_file <- file.path(scmetab_output_dir, config$scmetab_cell_file)
    if (!file.exists(score_file)) {
        warning("Missing scMetabolism score file: ", score_file)
        return(data.frame())
    }
    cell_scores <- read.csv(score_file, check.names = FALSE, stringsAsFactors = FALSE)
    required_cols <- c("cell", config$group_col, "pathway", "score")
    missing_cols <- setdiff(required_cols, colnames(cell_scores))
    if (length(missing_cols) > 0) {
        stop("Missing columns in ", score_file, ": ", paste(missing_cols, collapse = ", "))
    }
    score_long <- cell_scores %>% dplyr::transmute(cell = as.character(cell), group = recode_group(.data[[config$group_col]], 
        config), pathway = dplyr::case_when(pathway == "Glycolysis / Gluconeogenesis" ~ "Glycolysis", pathway == "Citrate cycle (TCA cycle)" ~ 
        "TCA / OXPHOS", TRUE ~ as.character(pathway)), score = as.numeric(score), source = "Xenium scMetabolism") %>% dplyr::filter(!is.na(group), 
        group %in% config$group_levels, pathway %in% pathway_levels)
    joined_file <- file.path(scmetab_output_dir, config$scmetab_joined_file)
    if (file.exists(joined_file)) {
        joined <- read.csv(joined_file, check.names = FALSE, stringsAsFactors = FALSE)
        if ("cell" %in% colnames(joined)) {
            if ("sample_id" %in% colnames(joined)) {
                joined$sample_id_lookup <- as.character(joined$sample_id)
            }
            else if ("sample_id_xenium" %in% colnames(joined)) {
                joined$sample_id_lookup <- as.character(joined$sample_id_xenium)
            }
            else {
                joined$sample_id_lookup <- NA_character_
            }
            sample_lookup <- joined %>% dplyr::transmute(cell = as.character(cell), sample_id = sample_id_lookup) %>% dplyr::distinct(cell, 
                .keep_all = TRUE)
            score_long <- score_long %>% dplyr::left_join(sample_lookup, by = "cell")
        }
        else {
            score_long$sample_id <- NA_character_
        }
    }
    else {
        score_long$sample_id <- NA_character_
    }
    score_long %>% dplyr::mutate(group = factor(group, levels = config$group_levels), pathway = factor(pathway, levels = pathway_levels)) %>% 
        dplyr::select(cell, sample_id, group, pathway, score, source)
},
make_energy_dotplot <- function(summary_df, config, title, legend_title, scale_by_column = FALSE) {
    plot_df <- summary_df %>% dplyr::mutate(pathway = factor(as.character(pathway), levels = rev(pathway_levels)), group = factor(as.character(group), 
        levels = config$group_levels))
    if (isTRUE(scale_by_column)) {
        plot_df <- plot_df %>% dplyr::group_by(pathway) %>% dplyr::mutate(fill_score = rescale01(mean_score), size_score = fill_score) %>% 
            dplyr::ungroup()
        legend_title <- "Column-scaled mean"
    }
    else {
        plot_df <- plot_df %>% dplyr::mutate(fill_score = mean_score, size_score = rescale01(mean_score))
    }
    ggplot(plot_df, aes(x = group, y = pathway)) + geom_point(aes(fill = fill_score, size = size_score), shape = 21, color = "grey35", 
        stroke = 0.35) + scale_y_discrete(labels = pathway_labels) + scale_fill_gradient(low = "white", high = "#D7301F", 
        name = legend_title) + scale_size(range = c(7, 14), guide = "none") + labs(title = title, x = NULL, y = NULL) + theme_bw(base_size = 15) + 
        theme(axis.text.x = element_text(size = 14, face = "bold"), axis.text.y = element_text(size = 14, lineheight = 0.95), 
            legend.position = "right", legend.title = element_text(size = 12.5), legend.text = element_text(size = 11.5), 
            panel.grid.major = element_line(color = "grey88"), panel.grid.minor = element_blank(), plot.title = element_text(size = 16, 
                hjust = 0.5, margin = margin(b = 6)), plot.margin = margin(5.5, 6, 5.5, 5.5))
},
make_violin_plot <- function(score_long, title, y_label) {
    ggplot(score_long, aes(x = group, y = score, fill = group)) + geom_violin(width = 0.82, alpha = 0.8, color = "grey30", 
        trim = TRUE) + geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", color = "grey15") + facet_wrap(~pathway, 
        nrow = 1, scales = "free_y") + scale_fill_manual(values = group_colors) + labs(title = title, x = NULL, y = y_label, 
        fill = NULL) + theme_bw(base_size = 13) + theme(axis.text.x = element_text(size = 12, face = "bold", angle = 20, 
        hjust = 1), strip.text = element_text(size = 12, face = "bold"), legend.position = "none", plot.title = element_text(hjust = 0.5))
},
make_density_plot <- function(score_long, title, x_label) {
    ggplot(score_long, aes(x = score, fill = group, color = group)) + geom_density(alpha = 0.32, linewidth = 0.8, adjust = 1.1) + 
        facet_wrap(~pathway, nrow = 1, scales = "free") + scale_fill_manual(values = group_colors) + scale_color_manual(values = group_colors_dark) + 
        labs(title = title, x = x_label, y = "Density", fill = NULL, color = NULL) + theme_bw(base_size = 13) + theme(strip.text = element_text(size = 12, 
        face = "bold"), legend.position = "top", plot.title = element_text(hjust = 0.5))
},
make_sample_plot <- function(sample_summary, title, y_label) {
    if (nrow(sample_summary) == 0) {
        return(ggplot() + annotate("text", x = 0, y = 0, label = "No sample IDs available", size = 4) + labs(title = title) + 
            theme_void())
    }
    ggplot(sample_summary, aes(x = group, y = mean_score, color = group)) + geom_point(position = position_jitter(width = 0.08, 
        height = 0), size = 2.8, alpha = 0.85) + stat_summary(fun = mean, geom = "crossbar", width = 0.42, color = "grey20", 
        linewidth = 0.35) + facet_wrap(~pathway, nrow = 1, scales = "free_y") + scale_color_manual(values = group_colors) + 
        labs(title = title, subtitle = "Each dot is one sample/group mean", x = NULL, y = y_label, color = NULL) + theme_bw(base_size = 13) + 
        theme(axis.text.x = element_text(size = 12, face = "bold", angle = 20, hjust = 1), strip.text = element_text(size = 12, 
            face = "bold"), legend.position = "none", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
},
cor_summary <- function(df) {
    df <- df %>% dplyr::filter(is.finite(imc_protein_score), is.finite(scmetabolism_score))
    if (nrow(df) < 3 || length(unique(df$imc_protein_score)) < 2 || length(unique(df$scmetabolism_score)) < 2) {
        return(data.frame(n = nrow(df), spearman_rho = NA_real_, spearman_p = NA_real_, pearson_r = NA_real_, pearson_p = NA_real_))
    }
    spearman_test <- suppressWarnings(stats::cor.test(df$imc_protein_score, df$scmetabolism_score, method = "spearman", exact = FALSE))
    pearson_test <- suppressWarnings(stats::cor.test(df$imc_protein_score, df$scmetabolism_score, method = "pearson"))
    data.frame(n = nrow(df), spearman_rho = unname(spearman_test$estimate), spearman_p = spearman_test$p.value, pearson_r = unname(pearson_test$estimate), 
        pearson_p = pearson_test$p.value)
},
make_score_bundle <- function(score_long, config, source_id, source_label, y_label, raw_legend_title) {
    score_summary <- summarise_score_long(score_long)
    sample_summary <- summarise_sample_scores(score_long)
    p_dot_raw <- make_energy_dotplot(score_summary, config, paste0(config$label, " ", source_label, " scores"), raw_legend_title, 
        scale_by_column = FALSE)
    p_dot_scaled <- make_energy_dotplot(score_summary, config, paste0(config$label, " ", source_label, " scores - scaled"), 
        "Column-scaled mean", scale_by_column = TRUE)
    p_violin <- make_violin_plot(score_long, paste0(config$label, " ", source_label, " score distributions"), y_label)
    p_density <- make_density_plot(score_long, paste0(config$label, " ", source_label, " score densities"), y_label)
    p_sample <- make_sample_plot(sample_summary, paste0(config$label, " ", source_label, " sample means"), paste0("Sample mean ", 
        y_label))
    basename <- paste0(config$file_prefix, "_", source_id)
    save_plot_set(p_dot_raw, paste0(basename, "_dotplot_raw_reversed_axes"), width = 5.6, height = 3.4)
    save_plot_set(p_dot_scaled, paste0(basename, "_dotplot_scaled_by_column_reversed_axes"), width = 5.6, height = 3.4)
    save_plot_set(p_violin, paste0(basename, "_violin"), width = 7.2, height = 4.4)
    save_plot_set(p_density, paste0(basename, "_density"), width = 7.2, height = 4.4)
    save_plot_set(p_sample, paste0(basename, "_sample_means"), width = 7.2, height = 4.4)
    write.csv(score_long, file.path(output_dir, paste0(basename, "_cell_scores.csv")), row.names = FALSE)
    write.csv(score_summary, file.path(output_dir, paste0(basename, "_summary.csv")), row.names = FALSE)
    write.csv(sample_summary, file.path(output_dir, paste0(basename, "_sample_summary.csv")), row.names = FALSE)
    list(score_long = score_long, summary = score_summary, sample_summary = sample_summary, plots = list(dot_raw = p_dot_raw, 
        dot_scaled = p_dot_scaled, violin = p_violin, density = p_density, sample = p_sample))
},
make_correlation_bundle <- function(config) {
    joined_file <- file.path(scmetab_output_dir, config$scmetab_joined_file)
    if (!file.exists(joined_file)) {
        warning("Missing scMetabolism joined file: ", joined_file)
        return(list(joined = data.frame(), long = data.frame(), stats = data.frame(), plots = list()))
    }
    scmetab_cell <- read.csv(joined_file, check.names = FALSE, stringsAsFactors = FALSE)
    required_cols <- c("KEY", config$group_col, "xenium_glycolysis_score", "xenium_tca_score")
    missing_cols <- setdiff(required_cols, colnames(scmetab_cell))
    if (length(missing_cols) > 0) {
        stop("Missing columns in ", joined_file, ": ", paste(missing_cols, collapse = ", "))
    }
    scmetab_cell$xenium_cell <- if ("cell" %in% colnames(scmetab_cell)) {
        as.character(scmetab_cell$cell)
    }
    else {
        NA_character_
    }
    scmetab_cell$scmetab_sample_id <- if ("sample_id" %in% colnames(scmetab_cell)) {
        as.character(scmetab_cell$sample_id)
    }
    else if ("sample_id_xenium" %in% colnames(scmetab_cell)) {
        as.character(scmetab_cell$sample_id_xenium)
    }
    else {
        NA_character_
    }
    scmetab_keep <- scmetab_cell %>% dplyr::transmute(KEY = as.character(KEY), xenium_cell, scmetab_sample_id, scmetab_group = recode_group(.data[[config$group_col]], 
        config), scmetabolism_glycolysis_score = xenium_glycolysis_score, scmetabolism_tca_score = xenium_tca_score)
    imc_protein_keep <- score_wide_all %>% dplyr::mutate(KEY = as.character(cell), group = recode_group(.data[[config$group_col]], 
        config)) %>% dplyr::select(KEY, imc_cell = cell, imc_sample_id = sample_id, group, imc_protein_glycolysis_score = Glycolysis, 
        imc_protein_tca_oxphos_score = `TCA / OXPHOS`)
    joined <- imc_protein_keep %>% dplyr::inner_join(scmetab_keep, by = "KEY") %>% dplyr::filter(group %in% config$group_levels, 
        scmetab_group %in% config$group_levels, group == scmetab_group) %>% dplyr::mutate(group = factor(group, levels = config$group_levels))
    cor_long <- dplyr::bind_rows(joined %>% dplyr::transmute(KEY, imc_cell, xenium_cell, imc_sample_id, scmetab_sample_id, 
        group, comparison = "Glycolysis", imc_protein_score = imc_protein_glycolysis_score, scmetabolism_score = scmetabolism_glycolysis_score), 
        joined %>% dplyr::transmute(KEY, imc_cell, xenium_cell, imc_sample_id, scmetab_sample_id, group, comparison = "TCA / OXPHOS", 
            imc_protein_score = imc_protein_tca_oxphos_score, scmetabolism_score = scmetabolism_tca_score)) %>% dplyr::filter(is.finite(imc_protein_score), 
        is.finite(scmetabolism_score)) %>% dplyr::mutate(group = factor(as.character(group), levels = config$group_levels), 
        comparison = factor(comparison, levels = pathway_levels))
    cor_stats_by_group <- cor_long %>% dplyr::group_by(comparison, group) %>% dplyr::group_modify(~cor_summary(.x)) %>% dplyr::ungroup() %>% 
        dplyr::mutate(level = "group")
    cor_stats_all <- cor_long %>% dplyr::group_by(comparison) %>% dplyr::group_modify(~cor_summary(.x)) %>% dplyr::ungroup() %>% 
        dplyr::mutate(group = "All matched cells", level = "all")
    cor_stats <- dplyr::bind_rows(cor_stats_by_group, cor_stats_all) %>% dplyr::select(level, comparison, group, dplyr::everything())
    cor_plot_stats <- cor_stats %>% dplyr::filter(level == "group") %>% dplyr::mutate(comparison = factor(as.character(comparison), 
        levels = pathway_levels), group = factor(as.character(group), levels = config$group_levels), label = paste0("Spearman rho = ", 
        formatC(spearman_rho, format = "f", digits = 2), "\n", format_p(spearman_p), "\nn = ", n))
    if (nrow(cor_long) > 0) {
        p_points <- ggplot(cor_long, aes(x = imc_protein_score, y = scmetabolism_score, color = group)) + geom_point(alpha = 0.45, 
            size = 1.25) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.65) + geom_label(data = cor_plot_stats, 
            aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.04, vjust = 1.08, size = 3.1, label.size = 0.2, 
            fill = "white", alpha = 0.9) + facet_grid(group ~ comparison, scales = "free") + scale_color_manual(values = group_colors) + 
            labs(title = paste0(config$label, " IMC protein score vs Xenium scMetabolism"), x = "IMC protein score", y = "Xenium scMetabolism AUCell score", 
                color = NULL) + theme_bw(base_size = 12) + theme(legend.position = "none", strip.text = element_text(size = 11.5, 
            face = "bold"), plot.title = element_text(hjust = 0.5))
        p_bins <- ggplot(cor_long, aes(x = imc_protein_score, y = scmetabolism_score)) + geom_bin2d(bins = 32) + geom_smooth(method = "lm", 
            formula = y ~ x, se = FALSE, color = "grey15", linewidth = 0.65) + geom_label(data = cor_plot_stats, aes(x = -Inf, 
            y = Inf, label = label), inherit.aes = FALSE, hjust = -0.04, vjust = 1.08, size = 3.1, label.size = 0.2, fill = "white", 
            alpha = 0.9) + facet_grid(group ~ comparison, scales = "free") + scale_fill_gradient(low = "white", high = "#D7301F", 
            name = "Matched cells") + labs(title = paste0(config$label, " binned IMC protein vs scMetabolism correlation"), 
            x = "IMC protein score", y = "Xenium scMetabolism AUCell score") + theme_bw(base_size = 12) + theme(strip.text = element_text(size = 11.5, 
            face = "bold"), plot.title = element_text(hjust = 0.5))
    }
    else {
        p_points <- ggplot() + annotate("text", x = 0, y = 0, label = "No matched cells available", size = 4) + labs(title = paste0(config$label, 
            " IMC protein score vs Xenium scMetabolism")) + theme_void()
        p_bins <- ggplot() + annotate("text", x = 0, y = 0, label = "No matched cells available", size = 4) + labs(title = paste0(config$label, 
            " binned IMC protein vs scMetabolism correlation")) + theme_void()
    }
    basename <- paste0(config$file_prefix, "_imc_protein_vs_scmetabolism")
    save_plot_set(p_points, paste0(basename, "_correlation_points"), width = 8.2, height = 5.2)
    save_plot_set(p_bins, paste0(basename, "_correlation_bins"), width = 8.2, height = 5.2)
    write.csv(joined, file.path(output_dir, paste0(basename, "_cell_joined.csv")), row.names = FALSE)
    write.csv(cor_long, file.path(output_dir, paste0(basename, "_correlation_long.csv")), row.names = FALSE)
    write.csv(cor_stats, file.path(output_dir, paste0(basename, "_correlation_stats.csv")), row.names = FALSE)
    list(joined = joined, long = cor_long, stats = cor_stats, plots = list(points = p_points, bins = p_bins))
},
build_analysis_result <- function(config) {
    imc_long <- make_imc_score_long(config)
    scmetab_long <- load_scmetab_score_long(config)
    list(imc = make_score_bundle(imc_long, config, source_id = "imc_protein", source_label = "IMC protein", y_label = "Mean zrescaled IMC protein score", 
        raw_legend_title = "Mean IMC protein score"), scmetabolism = make_score_bundle(scmetab_long, config, source_id = "scmetabolism", 
        source_label = "Xenium scMetabolism", y_label = "scMetabolism AUCell score", raw_legend_title = "Mean AUCell score"), 
        correlation = make_correlation_bundle(config))
}
  ),
  expressions = list(
"make_imc_protein_score" = quote(function(markers) {
    markers <- intersect(markers, rownames(spe))
    colMeans(as.matrix(SummarizedExperiment::assay(spe, score_assay)[markers, , drop = FALSE]), na.rm = TRUE)
}),
"recode_group" = quote(function(group, config) {
    group <- as.character(group)
    if (!is.null(config$group_recode)) {
        idx <- match(group, names(config$group_recode))
        matched <- !is.na(idx)
        group[matched] <- unname(config$group_recode[idx[matched]])
    }
    group
}),
"rescale01" = quote(function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) {
        return(rep(0.5, length(x)))
    }
    (x - rng[1])/(rng[2] - rng[1])
}),
"format_p" = quote(function(p_value) {
    dplyr::case_when(is.na(p_value) ~ "p = NA", p_value < 0.001 ~ "p < 0.001", TRUE ~ paste0("p = ", formatC(p_value, format = "f", 
        digits = 3)))
}),
"save_plot_set" = quote(function(plot, basename, width, height) {
    png_file <- file.path(output_dir, paste0(basename, ".png"))
    pdf_file <- file.path(output_dir, paste0(basename, ".pdf"))
    svg_file <- file.path(output_dir, paste0(basename, ".svg"))
    ggsave(png_file, plot, width = width, height = height, dpi = 300, bg = "white")
    grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE, bg = "white")
    print(plot)
    grDevices::dev.off()
    grDevices::svg(svg_file, width = width, height = height, bg = "white")
    print(plot)
    grDevices::dev.off()
    invisible(c(png_file, pdf_file, svg_file))
}),
"summarise_score_long" = quote(function(score_long) {
    score_long %>% dplyr::group_by(pathway, group) %>% dplyr::summarise(n_cells = dplyr::n(), n_samples = dplyr::n_distinct(sample_id[!is.na(sample_id) & 
        sample_id != ""]), mean_score = mean(score, na.rm = TRUE), median_score = median(score, na.rm = TRUE), q25 = stats::quantile(score, 
        0.25, na.rm = TRUE), q75 = stats::quantile(score, 0.75, na.rm = TRUE), sd_score = stats::sd(score, na.rm = TRUE), 
        se_score = sd_score/sqrt(n_cells), ci95 = 1.96 * se_score, .groups = "drop")
}),
"summarise_sample_scores" = quote(function(score_long) {
    score_long %>% dplyr::filter(!is.na(sample_id), sample_id != "") %>% dplyr::group_by(sample_id, pathway, group) %>% dplyr::summarise(n_cells = dplyr::n(), 
        mean_score = mean(score, na.rm = TRUE), median_score = median(score, na.rm = TRUE), .groups = "drop")
}),
"make_imc_score_long" = quote(function(config) {
    score_wide_all %>% dplyr::mutate(group = recode_group(.data[[config$group_col]], config)) %>% dplyr::filter(!is.na(group), 
        group %in% config$group_levels) %>% tidyr::pivot_longer(cols = dplyr::all_of(pathway_levels), names_to = "pathway", 
        values_to = "score") %>% dplyr::mutate(group = factor(group, levels = config$group_levels), pathway = factor(pathway, 
        levels = pathway_levels), source = "IMC protein") %>% dplyr::select(cell, sample_id, group, pathway, score, source)
}),
"load_scmetab_score_long" = quote(function(config) {
    score_file <- file.path(scmetab_output_dir, config$scmetab_cell_file)
    if (!file.exists(score_file)) {
        warning("Missing scMetabolism score file: ", score_file)
        return(data.frame())
    }
    cell_scores <- read.csv(score_file, check.names = FALSE, stringsAsFactors = FALSE)
    required_cols <- c("cell", config$group_col, "pathway", "score")
    missing_cols <- setdiff(required_cols, colnames(cell_scores))
    if (length(missing_cols) > 0) {
        stop("Missing columns in ", score_file, ": ", paste(missing_cols, collapse = ", "))
    }
    score_long <- cell_scores %>% dplyr::transmute(cell = as.character(cell), group = recode_group(.data[[config$group_col]], 
        config), pathway = dplyr::case_when(pathway == "Glycolysis / Gluconeogenesis" ~ "Glycolysis", pathway == "Citrate cycle (TCA cycle)" ~ 
        "TCA / OXPHOS", TRUE ~ as.character(pathway)), score = as.numeric(score), source = "Xenium scMetabolism") %>% dplyr::filter(!is.na(group), 
        group %in% config$group_levels, pathway %in% pathway_levels)
    joined_file <- file.path(scmetab_output_dir, config$scmetab_joined_file)
    if (file.exists(joined_file)) {
        joined <- read.csv(joined_file, check.names = FALSE, stringsAsFactors = FALSE)
        if ("cell" %in% colnames(joined)) {
            if ("sample_id" %in% colnames(joined)) {
                joined$sample_id_lookup <- as.character(joined$sample_id)
            }
            else if ("sample_id_xenium" %in% colnames(joined)) {
                joined$sample_id_lookup <- as.character(joined$sample_id_xenium)
            }
            else {
                joined$sample_id_lookup <- NA_character_
            }
            sample_lookup <- joined %>% dplyr::transmute(cell = as.character(cell), sample_id = sample_id_lookup) %>% dplyr::distinct(cell, 
                .keep_all = TRUE)
            score_long <- score_long %>% dplyr::left_join(sample_lookup, by = "cell")
        }
        else {
            score_long$sample_id <- NA_character_
        }
    }
    else {
        score_long$sample_id <- NA_character_
    }
    score_long %>% dplyr::mutate(group = factor(group, levels = config$group_levels), pathway = factor(pathway, levels = pathway_levels)) %>% 
        dplyr::select(cell, sample_id, group, pathway, score, source)
}),
"make_energy_dotplot" = quote(function(summary_df, config, title, legend_title, scale_by_column = FALSE) {
    plot_df <- summary_df %>% dplyr::mutate(pathway = factor(as.character(pathway), levels = rev(pathway_levels)), group = factor(as.character(group), 
        levels = config$group_levels))
    if (isTRUE(scale_by_column)) {
        plot_df <- plot_df %>% dplyr::group_by(pathway) %>% dplyr::mutate(fill_score = rescale01(mean_score), size_score = fill_score) %>% 
            dplyr::ungroup()
        legend_title <- "Column-scaled mean"
    }
    else {
        plot_df <- plot_df %>% dplyr::mutate(fill_score = mean_score, size_score = rescale01(mean_score))
    }
    ggplot(plot_df, aes(x = group, y = pathway)) + geom_point(aes(fill = fill_score, size = size_score), shape = 21, color = "grey35", 
        stroke = 0.35) + scale_y_discrete(labels = pathway_labels) + scale_fill_gradient(low = "white", high = "#D7301F", 
        name = legend_title) + scale_size(range = c(7, 14), guide = "none") + labs(title = title, x = NULL, y = NULL) + theme_bw(base_size = 15) + 
        theme(axis.text.x = element_text(size = 14, face = "bold"), axis.text.y = element_text(size = 14, lineheight = 0.95), 
            legend.position = "right", legend.title = element_text(size = 12.5), legend.text = element_text(size = 11.5), 
            panel.grid.major = element_line(color = "grey88"), panel.grid.minor = element_blank(), plot.title = element_text(size = 16, 
                hjust = 0.5, margin = margin(b = 6)), plot.margin = margin(5.5, 6, 5.5, 5.5))
}),
"make_violin_plot" = quote(function(score_long, title, y_label) {
    ggplot(score_long, aes(x = group, y = score, fill = group)) + geom_violin(width = 0.82, alpha = 0.8, color = "grey30", 
        trim = TRUE) + geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", color = "grey15") + facet_wrap(~pathway, 
        nrow = 1, scales = "free_y") + scale_fill_manual(values = group_colors) + labs(title = title, x = NULL, y = y_label, 
        fill = NULL) + theme_bw(base_size = 13) + theme(axis.text.x = element_text(size = 12, face = "bold", angle = 20, 
        hjust = 1), strip.text = element_text(size = 12, face = "bold"), legend.position = "none", plot.title = element_text(hjust = 0.5))
}),
"make_density_plot" = quote(function(score_long, title, x_label) {
    ggplot(score_long, aes(x = score, fill = group, color = group)) + geom_density(alpha = 0.32, linewidth = 0.8, adjust = 1.1) + 
        facet_wrap(~pathway, nrow = 1, scales = "free") + scale_fill_manual(values = group_colors) + scale_color_manual(values = group_colors_dark) + 
        labs(title = title, x = x_label, y = "Density", fill = NULL, color = NULL) + theme_bw(base_size = 13) + theme(strip.text = element_text(size = 12, 
        face = "bold"), legend.position = "top", plot.title = element_text(hjust = 0.5))
}),
"make_sample_plot" = quote(function(sample_summary, title, y_label) {
    if (nrow(sample_summary) == 0) {
        return(ggplot() + annotate("text", x = 0, y = 0, label = "No sample IDs available", size = 4) + labs(title = title) + 
            theme_void())
    }
    ggplot(sample_summary, aes(x = group, y = mean_score, color = group)) + geom_point(position = position_jitter(width = 0.08, 
        height = 0), size = 2.8, alpha = 0.85) + stat_summary(fun = mean, geom = "crossbar", width = 0.42, color = "grey20", 
        linewidth = 0.35) + facet_wrap(~pathway, nrow = 1, scales = "free_y") + scale_color_manual(values = group_colors) + 
        labs(title = title, subtitle = "Each dot is one sample/group mean", x = NULL, y = y_label, color = NULL) + theme_bw(base_size = 13) + 
        theme(axis.text.x = element_text(size = 12, face = "bold", angle = 20, hjust = 1), strip.text = element_text(size = 12, 
            face = "bold"), legend.position = "none", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
}),
"cor_summary" = quote(function(df) {
    df <- df %>% dplyr::filter(is.finite(imc_protein_score), is.finite(scmetabolism_score))
    if (nrow(df) < 3 || length(unique(df$imc_protein_score)) < 2 || length(unique(df$scmetabolism_score)) < 2) {
        return(data.frame(n = nrow(df), spearman_rho = NA_real_, spearman_p = NA_real_, pearson_r = NA_real_, pearson_p = NA_real_))
    }
    spearman_test <- suppressWarnings(stats::cor.test(df$imc_protein_score, df$scmetabolism_score, method = "spearman", exact = FALSE))
    pearson_test <- suppressWarnings(stats::cor.test(df$imc_protein_score, df$scmetabolism_score, method = "pearson"))
    data.frame(n = nrow(df), spearman_rho = unname(spearman_test$estimate), spearman_p = spearman_test$p.value, pearson_r = unname(pearson_test$estimate), 
        pearson_p = pearson_test$p.value)
}),
"make_score_bundle" = quote(function(score_long, config, source_id, source_label, y_label, raw_legend_title) {
    score_summary <- summarise_score_long(score_long)
    sample_summary <- summarise_sample_scores(score_long)
    p_dot_raw <- make_energy_dotplot(score_summary, config, paste0(config$label, " ", source_label, " scores"), raw_legend_title, 
        scale_by_column = FALSE)
    p_dot_scaled <- make_energy_dotplot(score_summary, config, paste0(config$label, " ", source_label, " scores - scaled"), 
        "Column-scaled mean", scale_by_column = TRUE)
    p_violin <- make_violin_plot(score_long, paste0(config$label, " ", source_label, " score distributions"), y_label)
    p_density <- make_density_plot(score_long, paste0(config$label, " ", source_label, " score densities"), y_label)
    p_sample <- make_sample_plot(sample_summary, paste0(config$label, " ", source_label, " sample means"), paste0("Sample mean ", 
        y_label))
    basename <- paste0(config$file_prefix, "_", source_id)
    save_plot_set(p_dot_raw, paste0(basename, "_dotplot_raw_reversed_axes"), width = 5.6, height = 3.4)
    save_plot_set(p_dot_scaled, paste0(basename, "_dotplot_scaled_by_column_reversed_axes"), width = 5.6, height = 3.4)
    save_plot_set(p_violin, paste0(basename, "_violin"), width = 7.2, height = 4.4)
    save_plot_set(p_density, paste0(basename, "_density"), width = 7.2, height = 4.4)
    save_plot_set(p_sample, paste0(basename, "_sample_means"), width = 7.2, height = 4.4)
    write.csv(score_long, file.path(output_dir, paste0(basename, "_cell_scores.csv")), row.names = FALSE)
    write.csv(score_summary, file.path(output_dir, paste0(basename, "_summary.csv")), row.names = FALSE)
    write.csv(sample_summary, file.path(output_dir, paste0(basename, "_sample_summary.csv")), row.names = FALSE)
    list(score_long = score_long, summary = score_summary, sample_summary = sample_summary, plots = list(dot_raw = p_dot_raw, 
        dot_scaled = p_dot_scaled, violin = p_violin, density = p_density, sample = p_sample))
}),
"make_correlation_bundle" = quote(function(config) {
    joined_file <- file.path(scmetab_output_dir, config$scmetab_joined_file)
    if (!file.exists(joined_file)) {
        warning("Missing scMetabolism joined file: ", joined_file)
        return(list(joined = data.frame(), long = data.frame(), stats = data.frame(), plots = list()))
    }
    scmetab_cell <- read.csv(joined_file, check.names = FALSE, stringsAsFactors = FALSE)
    required_cols <- c("KEY", config$group_col, "xenium_glycolysis_score", "xenium_tca_score")
    missing_cols <- setdiff(required_cols, colnames(scmetab_cell))
    if (length(missing_cols) > 0) {
        stop("Missing columns in ", joined_file, ": ", paste(missing_cols, collapse = ", "))
    }
    scmetab_cell$xenium_cell <- if ("cell" %in% colnames(scmetab_cell)) {
        as.character(scmetab_cell$cell)
    }
    else {
        NA_character_
    }
    scmetab_cell$scmetab_sample_id <- if ("sample_id" %in% colnames(scmetab_cell)) {
        as.character(scmetab_cell$sample_id)
    }
    else if ("sample_id_xenium" %in% colnames(scmetab_cell)) {
        as.character(scmetab_cell$sample_id_xenium)
    }
    else {
        NA_character_
    }
    scmetab_keep <- scmetab_cell %>% dplyr::transmute(KEY = as.character(KEY), xenium_cell, scmetab_sample_id, scmetab_group = recode_group(.data[[config$group_col]], 
        config), scmetabolism_glycolysis_score = xenium_glycolysis_score, scmetabolism_tca_score = xenium_tca_score)
    imc_protein_keep <- score_wide_all %>% dplyr::mutate(KEY = as.character(cell), group = recode_group(.data[[config$group_col]], 
        config)) %>% dplyr::select(KEY, imc_cell = cell, imc_sample_id = sample_id, group, imc_protein_glycolysis_score = Glycolysis, 
        imc_protein_tca_oxphos_score = `TCA / OXPHOS`)
    joined <- imc_protein_keep %>% dplyr::inner_join(scmetab_keep, by = "KEY") %>% dplyr::filter(group %in% config$group_levels, 
        scmetab_group %in% config$group_levels, group == scmetab_group) %>% dplyr::mutate(group = factor(group, levels = config$group_levels))
    cor_long <- dplyr::bind_rows(joined %>% dplyr::transmute(KEY, imc_cell, xenium_cell, imc_sample_id, scmetab_sample_id, 
        group, comparison = "Glycolysis", imc_protein_score = imc_protein_glycolysis_score, scmetabolism_score = scmetabolism_glycolysis_score), 
        joined %>% dplyr::transmute(KEY, imc_cell, xenium_cell, imc_sample_id, scmetab_sample_id, group, comparison = "TCA / OXPHOS", 
            imc_protein_score = imc_protein_tca_oxphos_score, scmetabolism_score = scmetabolism_tca_score)) %>% dplyr::filter(is.finite(imc_protein_score), 
        is.finite(scmetabolism_score)) %>% dplyr::mutate(group = factor(as.character(group), levels = config$group_levels), 
        comparison = factor(comparison, levels = pathway_levels))
    cor_stats_by_group <- cor_long %>% dplyr::group_by(comparison, group) %>% dplyr::group_modify(~cor_summary(.x)) %>% dplyr::ungroup() %>% 
        dplyr::mutate(level = "group")
    cor_stats_all <- cor_long %>% dplyr::group_by(comparison) %>% dplyr::group_modify(~cor_summary(.x)) %>% dplyr::ungroup() %>% 
        dplyr::mutate(group = "All matched cells", level = "all")
    cor_stats <- dplyr::bind_rows(cor_stats_by_group, cor_stats_all) %>% dplyr::select(level, comparison, group, dplyr::everything())
    cor_plot_stats <- cor_stats %>% dplyr::filter(level == "group") %>% dplyr::mutate(comparison = factor(as.character(comparison), 
        levels = pathway_levels), group = factor(as.character(group), levels = config$group_levels), label = paste0("Spearman rho = ", 
        formatC(spearman_rho, format = "f", digits = 2), "\n", format_p(spearman_p), "\nn = ", n))
    if (nrow(cor_long) > 0) {
        p_points <- ggplot(cor_long, aes(x = imc_protein_score, y = scmetabolism_score, color = group)) + geom_point(alpha = 0.45, 
            size = 1.25) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.65) + geom_label(data = cor_plot_stats, 
            aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.04, vjust = 1.08, size = 3.1, label.size = 0.2, 
            fill = "white", alpha = 0.9) + facet_grid(group ~ comparison, scales = "free") + scale_color_manual(values = group_colors) + 
            labs(title = paste0(config$label, " IMC protein score vs Xenium scMetabolism"), x = "IMC protein score", y = "Xenium scMetabolism AUCell score", 
                color = NULL) + theme_bw(base_size = 12) + theme(legend.position = "none", strip.text = element_text(size = 11.5, 
            face = "bold"), plot.title = element_text(hjust = 0.5))
        p_bins <- ggplot(cor_long, aes(x = imc_protein_score, y = scmetabolism_score)) + geom_bin2d(bins = 32) + geom_smooth(method = "lm", 
            formula = y ~ x, se = FALSE, color = "grey15", linewidth = 0.65) + geom_label(data = cor_plot_stats, aes(x = -Inf, 
            y = Inf, label = label), inherit.aes = FALSE, hjust = -0.04, vjust = 1.08, size = 3.1, label.size = 0.2, fill = "white", 
            alpha = 0.9) + facet_grid(group ~ comparison, scales = "free") + scale_fill_gradient(low = "white", high = "#D7301F", 
            name = "Matched cells") + labs(title = paste0(config$label, " binned IMC protein vs scMetabolism correlation"), 
            x = "IMC protein score", y = "Xenium scMetabolism AUCell score") + theme_bw(base_size = 12) + theme(strip.text = element_text(size = 11.5, 
            face = "bold"), plot.title = element_text(hjust = 0.5))
    }
    else {
        p_points <- ggplot() + annotate("text", x = 0, y = 0, label = "No matched cells available", size = 4) + labs(title = paste0(config$label, 
            " IMC protein score vs Xenium scMetabolism")) + theme_void()
        p_bins <- ggplot() + annotate("text", x = 0, y = 0, label = "No matched cells available", size = 4) + labs(title = paste0(config$label, 
            " binned IMC protein vs scMetabolism correlation")) + theme_void()
    }
    basename <- paste0(config$file_prefix, "_imc_protein_vs_scmetabolism")
    save_plot_set(p_points, paste0(basename, "_correlation_points"), width = 8.2, height = 5.2)
    save_plot_set(p_bins, paste0(basename, "_correlation_bins"), width = 8.2, height = 5.2)
    write.csv(joined, file.path(output_dir, paste0(basename, "_cell_joined.csv")), row.names = FALSE)
    write.csv(cor_long, file.path(output_dir, paste0(basename, "_correlation_long.csv")), row.names = FALSE)
    write.csv(cor_stats, file.path(output_dir, paste0(basename, "_correlation_stats.csv")), row.names = FALSE)
    list(joined = joined, long = cor_long, stats = cor_stats, plots = list(points = p_points, bins = p_bins))
}),
"build_analysis_result" = quote(function(config) {
    imc_long <- make_imc_score_long(config)
    scmetab_long <- load_scmetab_score_long(config)
    list(imc = make_score_bundle(imc_long, config, source_id = "imc_protein", source_label = "IMC protein", y_label = "Mean zrescaled IMC protein score", 
        raw_legend_title = "Mean IMC protein score"), scmetabolism = make_score_bundle(scmetab_long, config, source_id = "scmetabolism", 
        source_label = "Xenium scMetabolism", y_label = "scMetabolism AUCell score", raw_legend_title = "Mean AUCell score"), 
        correlation = make_correlation_bundle(config))
})
  ),
  plots = list(

  )
)
original_modules[["acta"]] <- list(
  setup = expression(
group_levels <- c("Other_EC", "ACTA2_EC"),
pathway_lookup <- c(`Glycolysis / Gluconeogenesis` = "Glycolysis", `Citrate cycle (TCA cycle)` = "TCA / OXPHOS"),
pathway_levels <- unname(pathway_lookup),
pathway_labels <- c(Glycolysis = "Glycolysis", `TCA / OXPHOS` = "TCA /\nOXPHOS"),
group_colors <- c(Other_EC = "#4C78A8", ACTA2_EC = "#E45756"),
to_dgC <- function(x) {
    if (inherits(x, "dgCMatrix")) 
        return(x)
    x <- as(x, "CsparseMatrix")
    as(x, "dgCMatrix")
},
rescale01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) {
        return(rep(0.5, length(x)))
    }
    (x - rng[1])/(rng[2] - rng[1])
},
format_p <- function(p_value) {
    dplyr::case_when(is.na(p_value) ~ "p = NA", p_value < 0.001 ~ "p < 0.001", TRUE ~ paste0("p = ", formatC(p_value, format = "f", 
        digits = 3)))
},
save_plot_set <- function(plot, basename, width, height) {
    png_file <- file.path(output_dir, paste0(basename, ".png"))
    pdf_file <- file.path(output_dir, paste0(basename, ".pdf"))
    svg_file <- file.path(output_dir, paste0(basename, ".svg"))
    ggplot2::ggsave(png_file, plot, width = width, height = height, dpi = 300, bg = "white")
    grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE, bg = "white")
    print(plot)
    grDevices::dev.off()
    grDevices::svg(svg_file, width = width, height = height, bg = "white")
    print(plot)
    grDevices::dev.off()
    invisible(c(png_file, pdf_file, svg_file))
},
make_dotplot <- function(summary_df, scale_by_column = FALSE) {
    plot_df <- summary_df %>% dplyr::mutate(pathway = factor(as.character(pathway), levels = rev(pathway_levels)), acta2_ec_group = factor(as.character(acta2_ec_group), 
        levels = group_levels))
    if (isTRUE(scale_by_column)) {
        plot_df <- plot_df %>% dplyr::group_by(pathway) %>% dplyr::mutate(fill_score = rescale01(mean_score), size_score = fill_score) %>% 
            dplyr::ungroup()
        legend_title <- "Column-scaled mean"
        plot_title <- "ACTA2_EC scMetabolism scores, scaled"
    }
    else {
        plot_df <- plot_df %>% dplyr::mutate(fill_score = mean_score, size_score = rescale01(mean_score))
        legend_title <- "Mean AUCell score"
        plot_title <- "ACTA2_EC scMetabolism scores"
    }
    ggplot(plot_df, aes(x = acta2_ec_group, y = pathway)) + geom_point(aes(fill = fill_score, size = size_score), shape = 21, 
        color = "grey35", stroke = 0.35) + scale_y_discrete(labels = pathway_labels) + scale_fill_gradient(low = "white", 
        high = "#D7301F", name = legend_title) + scale_size(range = c(7, 14), guide = "none") + labs(title = plot_title, 
        x = NULL, y = NULL) + theme_bw(base_size = 15) + theme(axis.text.x = element_text(size = 14, face = "bold"), axis.text.y = element_text(size = 14, 
        lineheight = 0.95), legend.position = "right", legend.title = element_text(size = 12.5), legend.text = element_text(size = 11.5), 
        panel.grid.major = element_line(color = "grey88"), panel.grid.minor = element_blank(), plot.title = element_text(size = 16, 
            hjust = 0.5, margin = margin(b = 6)), plot.margin = margin(5.5, 6, 5.5, 5.5))
}
  ),
  expressions = list(
"group_levels" = quote(c("Other_EC", "ACTA2_EC")),
"pathway_lookup" = quote(c(`Glycolysis / Gluconeogenesis` = "Glycolysis", `Citrate cycle (TCA cycle)` = "TCA / OXPHOS")),
"pathway_levels" = quote(unname(pathway_lookup)),
"pathway_labels" = quote(c(Glycolysis = "Glycolysis", `TCA / OXPHOS` = "TCA /\nOXPHOS")),
"group_colors" = quote(c(Other_EC = "#4C78A8", ACTA2_EC = "#E45756")),
"to_dgC" = quote(function(x) {
    if (inherits(x, "dgCMatrix")) 
        return(x)
    x <- as(x, "CsparseMatrix")
    as(x, "dgCMatrix")
}),
"rescale01" = quote(function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) {
        return(rep(0.5, length(x)))
    }
    (x - rng[1])/(rng[2] - rng[1])
}),
"format_p" = quote(function(p_value) {
    dplyr::case_when(is.na(p_value) ~ "p = NA", p_value < 0.001 ~ "p < 0.001", TRUE ~ paste0("p = ", formatC(p_value, format = "f", 
        digits = 3)))
}),
"save_plot_set" = quote(function(plot, basename, width, height) {
    png_file <- file.path(output_dir, paste0(basename, ".png"))
    pdf_file <- file.path(output_dir, paste0(basename, ".pdf"))
    svg_file <- file.path(output_dir, paste0(basename, ".svg"))
    ggplot2::ggsave(png_file, plot, width = width, height = height, dpi = 300, bg = "white")
    grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE, bg = "white")
    print(plot)
    grDevices::dev.off()
    grDevices::svg(svg_file, width = width, height = height, bg = "white")
    print(plot)
    grDevices::dev.off()
    invisible(c(png_file, pdf_file, svg_file))
}),
"make_dotplot" = quote(function(summary_df, scale_by_column = FALSE) {
    plot_df <- summary_df %>% dplyr::mutate(pathway = factor(as.character(pathway), levels = rev(pathway_levels)), acta2_ec_group = factor(as.character(acta2_ec_group), 
        levels = group_levels))
    if (isTRUE(scale_by_column)) {
        plot_df <- plot_df %>% dplyr::group_by(pathway) %>% dplyr::mutate(fill_score = rescale01(mean_score), size_score = fill_score) %>% 
            dplyr::ungroup()
        legend_title <- "Column-scaled mean"
        plot_title <- "ACTA2_EC scMetabolism scores, scaled"
    }
    else {
        plot_df <- plot_df %>% dplyr::mutate(fill_score = mean_score, size_score = rescale01(mean_score))
        legend_title <- "Mean AUCell score"
        plot_title <- "ACTA2_EC scMetabolism scores"
    }
    ggplot(plot_df, aes(x = acta2_ec_group, y = pathway)) + geom_point(aes(fill = fill_score, size = size_score), shape = 21, 
        color = "grey35", stroke = 0.35) + scale_y_discrete(labels = pathway_labels) + scale_fill_gradient(low = "white", 
        high = "#D7301F", name = legend_title) + scale_size(range = c(7, 14), guide = "none") + labs(title = plot_title, 
        x = NULL, y = NULL) + theme_bw(base_size = 15) + theme(axis.text.x = element_text(size = 14, face = "bold"), axis.text.y = element_text(size = 14, 
        lineheight = 0.95), legend.position = "right", legend.title = element_text(size = 12.5), legend.text = element_text(size = 11.5), 
        panel.grid.major = element_line(color = "grey88"), panel.grid.minor = element_blank(), plot.title = element_text(size = 16, 
            hjust = 0.5, margin = margin(b = 6)), plot.margin = margin(5.5, 6, 5.5, 5.5))
})
  ),
  plots = list(
"p_density" = quote(ggplot(score_long, aes(x = score, fill = acta2_ec_group, color = acta2_ec_group)) + geom_density(alpha = 0.28, linewidth = 0.8, 
    adjust = 1.05) + geom_vline(data = score_summary, aes(xintercept = median_score, color = acta2_ec_group), linetype = "dashed", 
    linewidth = 0.55, show.legend = FALSE) + facet_wrap(~pathway, nrow = 1, scales = "free_x") + scale_fill_manual(values = group_colors) + 
    scale_color_manual(values = group_colors) + labs(title = "Single-cell score distributions", x = "AUCell score", y = "Density", 
    fill = NULL, color = NULL) + theme_bw(base_size = 14) + theme(strip.text = element_text(size = 13.5, face = "bold"), 
    axis.text = element_text(size = 12), legend.position = "top", plot.title = element_text(size = 16, hjust = 0.5))),
"p_ecdf" = quote(ggplot(score_long, aes(x = score, color = acta2_ec_group)) + stat_ecdf(linewidth = 0.9) + facet_wrap(~pathway, nrow = 1, 
    scales = "free_x") + scale_color_manual(values = group_colors) + labs(title = "Single-cell empirical cumulative distributions", 
    x = "AUCell score", y = "Fraction of cells", color = NULL) + theme_bw(base_size = 14) + theme(strip.text = element_text(size = 13.5, 
    face = "bold"), axis.text = element_text(size = 12), legend.position = "top", plot.title = element_text(size = 16, hjust = 0.5))),
"p_cell_jitter" = quote(score_long %>% dplyr::mutate(pathway = factor(as.character(pathway), levels = pathway_levels), acta2_ec_group = factor(as.character(acta2_ec_group), 
    levels = group_levels)) %>% ggplot(aes(x = acta2_ec_group, y = score, color = acta2_ec_group)) + geom_jitter(width = 0.16, 
    height = 0, alpha = 0.34, size = 1.15) + geom_boxplot(width = 0.22, outlier.shape = NA, fill = "white", color = "grey20", 
    linewidth = 0.55) + facet_wrap(~pathway, nrow = 1, scales = "free_y") + scale_color_manual(values = group_colors) + labs(title = "Single-cell scores", 
    x = NULL, y = "AUCell score", color = NULL) + theme_bw(base_size = 14) + theme(axis.text.x = element_text(size = 12.5, 
    face = "bold", angle = 15, hjust = 1), strip.text = element_text(size = 13.5, face = "bold"), legend.position = "none", 
    plot.title = element_text(size = 16, hjust = 0.5)))
  )
)
original_modules[["gsea"]] <- list(
  setup = expression(
group_levels <- c("OXPHOS", "Glycolysis / carbon catabolism", "ECM / fibrosis"),
extract_exact_terms <- function(contrast_key, treatment, terms, groups) {
    contrast_data <- gsea[gsea$Contrast_key == contrast_key, , drop = FALSE]
    matches <- match(tolower(terms), tolower(contrast_data$Term))
    if (anyNA(matches)) {
        stop("Missing terms for ", treatment, ": ", paste(terms[is.na(matches)], collapse = ", "))
    }
    selected <- contrast_data[matches, , drop = FALSE]
    if (anyDuplicated(tolower(selected$Term))) {
        stop("Duplicate requested terms found for ", treatment, ".")
    }
    selected$Treatment <- treatment
    selected$Category <- groups
    selected$Treatment_order <- seq_along(terms)
    selected
},
extract_union_terms <- function(contrast_key, treatment) {
    selected <- extract_exact_terms(contrast_key = contrast_key, treatment = treatment, terms = union_terms, groups = union_annotation$Category)
    selected$Union_order <- union_annotation$Union_order
    selected
},
prepare_plot_data <- function(data, term_order) {
    data$Treatment <- factor(data$Treatment, levels = c("2DG + OX", "LDHA + OX"))
    data$Category <- factor(data$Category, levels = group_levels)
    term_labels <- setNames(vapply(term_order, function(term) paste(strwrap(term, width = 39), collapse = "\n"), character(1)), 
        term_order)
    data$Term_label <- factor(unname(term_labels[data$Term]), levels = rev(unname(term_labels[term_order])))
    data$minus_log10_FDR <- -log10(pmax(data$FDR, .Machine$double.xmin))
    data$FDR_status <- ifelse(data$FDR <= 0.25, "FDR <= 0.25", "FDR > 0.25")
    data
},
common_theme <- function(base_size = 9) {
    theme_bw(base_size = base_size, base_family = "Arial") + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_line(color = "#D9DEE2", linewidth = 0.28), panel.border = element_rect(color = "#4D555A", 
            linewidth = 0.45), strip.background = element_rect(fill = "#EEF1F3", color = "#9FAAB1", linewidth = 0.4), strip.text.y.left = element_text(face = "bold", 
            color = "#202326", angle = 0, size = 7.4), axis.text = element_text(color = "#202326"), axis.text.y = element_text(size = 7.7), 
        axis.title = element_text(color = "#202326", size = 8.5), plot.title = element_text(face = "bold", size = 10.5, color = "#202326"), 
        plot.subtitle = element_text(size = 7.8, color = "#4F5960"), plot.title.position = "plot", legend.position = "bottom", 
        legend.box = "vertical", legend.direction = "horizontal", legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7), legend.key.height = grid::unit(3.5, "mm"), legend.key.width = grid::unit(5, 
            "mm"), legend.spacing.x = grid::unit(2, "mm"), legend.margin = margin(t = -3, unit = "pt"), plot.margin = margin(7, 
            8, 6, 7, unit = "pt"))
},
plot_scales <- function(title, subtitle, scale_data, include_treatment_legend = FALSE) {
    significant_fdr <- scale_data$minus_log10_FDR[scale_data$FDR <= 0.25]
    if (!length(significant_fdr)) {
        stop("No terms at FDR <= 0.25 are available for the red color scale.")
    }
    fdr_limits <- range(significant_fdr, na.rm = TRUE)
    if (diff(fdr_limits) == 0) {
        fdr_limits <- fdr_limits + c(-0.05, 0.05)
    }
    fdr_breaks <- seq(fdr_limits[1], fdr_limits[2], length.out = 3)
    has_nonsignificant <- any(scale_data$FDR > 0.25)
    guides_list <- list(fill = guide_colorbar(order = 1, title.position = "top", barwidth = grid::unit(25, "mm"), barheight = grid::unit(2.8, 
        "mm")), size = guide_legend(order = 2, title.position = "top", nrow = 1, override.aes = list(shape = 21, fill = "#F4B3B0", 
        color = "#202326")))
    if (has_nonsignificant) {
        guides_list$colour <- guide_legend(order = 3, title.position = "top", nrow = 1, override.aes = list(shape = 21, fill = "#B8BDC3", 
            colour = "#202326", size = 3.5))
    }
    if (include_treatment_legend) {
        guides_list$shape <- guide_legend(order = 4, title.position = "top", nrow = 1, override.aes = list(size = 3.5, fill = "#C83E3A"))
    }
    list(geom_vline(xintercept = 0, color = "#4D555A", linewidth = 0.35, linetype = "dashed"), scale_x_continuous(name = "Normalized enrichment score (NES)", 
        limits = x_limits, breaks = x_breaks, expand = expansion(mult = c(0.03, 0.03))), scale_size_continuous(name = "GO gene-set size", 
        range = c(2.6, 6), limits = size_limits, breaks = size_breaks), scale_fill_gradient(name = expression(-log[10] * 
        "(FDR)"), low = "#FADBD8", high = "#B2182B", limits = fdr_limits, breaks = fdr_breaks, labels = scales::label_number(accuracy = 0.1)), 
        if (has_nonsignificant) {
            scale_color_manual(name = "FDR status", values = c(`FDR > 0.25` = "#202326"), drop = TRUE)
        }, facet_grid(rows = vars(Category), scales = "free_y", space = "free_y", switch = "y", drop = TRUE, labeller = as_labeller(c(OXPHOS = "OXPHOS", 
            `Glycolysis / carbon catabolism` = "Glycolysis /\ncarbon catabolism", `ECM / fibrosis` = "ECM /\nfibrosis"))), 
        labs(title = title, subtitle = subtitle, y = NULL), do.call(guides, guides_list), coord_cartesian(clip = "off"), 
        common_theme())
},
single_treatment_plot <- function(data, title) {
    significant <- data[data$FDR <= 0.25, , drop = FALSE]
    nonsignificant <- data[data$FDR > 0.25, , drop = FALSE]
    plot <- ggplot(data, aes(x = NES, y = Term_label, size = Set_size)) + geom_point(data = significant, aes(fill = minus_log10_FDR), 
        shape = 21, color = "#202326", stroke = 0.45, alpha = 0.96)
    if (nrow(nonsignificant)) {
        plot <- plot + geom_point(data = nonsignificant, aes(color = FDR_status), shape = 21, fill = "#B8BDC3", stroke = 0.45, 
            alpha = 0.96)
    }
    plot + plot_scales(title = title, subtitle = paste0("edgeR exact test; ComBat-seq PCS-adjusted counts\n", "red scale uses this plot's FDR range; gray = FDR > 0.25"), 
        scale_data = data)
},
save_plot <- function(plot, stem, width_mm, height_mm) {
    pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
    png_file <- file.path(output_dir, paste0(stem, ".png"))
    ggsave(pdf_file, plot = plot, width = width_mm, height = height_mm, units = "mm", device = grDevices::cairo_pdf, bg = "white")
    ggsave(png_file, plot = plot, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
    c(pdf_file, png_file)
}
  ),
  expressions = list(
"group_levels" = quote(c("OXPHOS", "Glycolysis / carbon catabolism", "ECM / fibrosis")),
"extract_exact_terms" = quote(function(contrast_key, treatment, terms, groups) {
    contrast_data <- gsea[gsea$Contrast_key == contrast_key, , drop = FALSE]
    matches <- match(tolower(terms), tolower(contrast_data$Term))
    if (anyNA(matches)) {
        stop("Missing terms for ", treatment, ": ", paste(terms[is.na(matches)], collapse = ", "))
    }
    selected <- contrast_data[matches, , drop = FALSE]
    if (anyDuplicated(tolower(selected$Term))) {
        stop("Duplicate requested terms found for ", treatment, ".")
    }
    selected$Treatment <- treatment
    selected$Category <- groups
    selected$Treatment_order <- seq_along(terms)
    selected
}),
"extract_union_terms" = quote(function(contrast_key, treatment) {
    selected <- extract_exact_terms(contrast_key = contrast_key, treatment = treatment, terms = union_terms, groups = union_annotation$Category)
    selected$Union_order <- union_annotation$Union_order
    selected
}),
"prepare_plot_data" = quote(function(data, term_order) {
    data$Treatment <- factor(data$Treatment, levels = c("2DG + OX", "LDHA + OX"))
    data$Category <- factor(data$Category, levels = group_levels)
    term_labels <- setNames(vapply(term_order, function(term) paste(strwrap(term, width = 39), collapse = "\n"), character(1)), 
        term_order)
    data$Term_label <- factor(unname(term_labels[data$Term]), levels = rev(unname(term_labels[term_order])))
    data$minus_log10_FDR <- -log10(pmax(data$FDR, .Machine$double.xmin))
    data$FDR_status <- ifelse(data$FDR <= 0.25, "FDR <= 0.25", "FDR > 0.25")
    data
}),
"common_theme" = quote(function(base_size = 9) {
    theme_bw(base_size = base_size, base_family = "Arial") + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_line(color = "#D9DEE2", linewidth = 0.28), panel.border = element_rect(color = "#4D555A", 
            linewidth = 0.45), strip.background = element_rect(fill = "#EEF1F3", color = "#9FAAB1", linewidth = 0.4), strip.text.y.left = element_text(face = "bold", 
            color = "#202326", angle = 0, size = 7.4), axis.text = element_text(color = "#202326"), axis.text.y = element_text(size = 7.7), 
        axis.title = element_text(color = "#202326", size = 8.5), plot.title = element_text(face = "bold", size = 10.5, color = "#202326"), 
        plot.subtitle = element_text(size = 7.8, color = "#4F5960"), plot.title.position = "plot", legend.position = "bottom", 
        legend.box = "vertical", legend.direction = "horizontal", legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7), legend.key.height = grid::unit(3.5, "mm"), legend.key.width = grid::unit(5, 
            "mm"), legend.spacing.x = grid::unit(2, "mm"), legend.margin = margin(t = -3, unit = "pt"), plot.margin = margin(7, 
            8, 6, 7, unit = "pt"))
}),
"plot_scales" = quote(function(title, subtitle, scale_data, include_treatment_legend = FALSE) {
    significant_fdr <- scale_data$minus_log10_FDR[scale_data$FDR <= 0.25]
    if (!length(significant_fdr)) {
        stop("No terms at FDR <= 0.25 are available for the red color scale.")
    }
    fdr_limits <- range(significant_fdr, na.rm = TRUE)
    if (diff(fdr_limits) == 0) {
        fdr_limits <- fdr_limits + c(-0.05, 0.05)
    }
    fdr_breaks <- seq(fdr_limits[1], fdr_limits[2], length.out = 3)
    has_nonsignificant <- any(scale_data$FDR > 0.25)
    guides_list <- list(fill = guide_colorbar(order = 1, title.position = "top", barwidth = grid::unit(25, "mm"), barheight = grid::unit(2.8, 
        "mm")), size = guide_legend(order = 2, title.position = "top", nrow = 1, override.aes = list(shape = 21, fill = "#F4B3B0", 
        color = "#202326")))
    if (has_nonsignificant) {
        guides_list$colour <- guide_legend(order = 3, title.position = "top", nrow = 1, override.aes = list(shape = 21, fill = "#B8BDC3", 
            colour = "#202326", size = 3.5))
    }
    if (include_treatment_legend) {
        guides_list$shape <- guide_legend(order = 4, title.position = "top", nrow = 1, override.aes = list(size = 3.5, fill = "#C83E3A"))
    }
    list(geom_vline(xintercept = 0, color = "#4D555A", linewidth = 0.35, linetype = "dashed"), scale_x_continuous(name = "Normalized enrichment score (NES)", 
        limits = x_limits, breaks = x_breaks, expand = expansion(mult = c(0.03, 0.03))), scale_size_continuous(name = "GO gene-set size", 
        range = c(2.6, 6), limits = size_limits, breaks = size_breaks), scale_fill_gradient(name = expression(-log[10] * 
        "(FDR)"), low = "#FADBD8", high = "#B2182B", limits = fdr_limits, breaks = fdr_breaks, labels = scales::label_number(accuracy = 0.1)), 
        if (has_nonsignificant) {
            scale_color_manual(name = "FDR status", values = c(`FDR > 0.25` = "#202326"), drop = TRUE)
        }, facet_grid(rows = vars(Category), scales = "free_y", space = "free_y", switch = "y", drop = TRUE, labeller = as_labeller(c(OXPHOS = "OXPHOS", 
            `Glycolysis / carbon catabolism` = "Glycolysis /\ncarbon catabolism", `ECM / fibrosis` = "ECM /\nfibrosis"))), 
        labs(title = title, subtitle = subtitle, y = NULL), do.call(guides, guides_list), coord_cartesian(clip = "off"), 
        common_theme())
}),
"single_treatment_plot" = quote(function(data, title) {
    significant <- data[data$FDR <= 0.25, , drop = FALSE]
    nonsignificant <- data[data$FDR > 0.25, , drop = FALSE]
    plot <- ggplot(data, aes(x = NES, y = Term_label, size = Set_size)) + geom_point(data = significant, aes(fill = minus_log10_FDR), 
        shape = 21, color = "#202326", stroke = 0.45, alpha = 0.96)
    if (nrow(nonsignificant)) {
        plot <- plot + geom_point(data = nonsignificant, aes(color = FDR_status), shape = 21, fill = "#B8BDC3", stroke = 0.45, 
            alpha = 0.96)
    }
    plot + plot_scales(title = title, subtitle = paste0("edgeR exact test; ComBat-seq PCS-adjusted counts\n", "red scale uses this plot's FDR range; gray = FDR > 0.25"), 
        scale_data = data)
}),
"save_plot" = quote(function(plot, stem, width_mm, height_mm) {
    pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
    png_file <- file.path(output_dir, paste0(stem, ".png"))
    ggsave(pdf_file, plot = plot, width = width_mm, height = height_mm, units = "mm", device = grDevices::cairo_pdf, bg = "white")
    ggsave(png_file, plot = plot, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
    c(pdf_file, png_file)
})
  ),
  plots = list(
"figure_combined" = quote(ggplot(plot_combined, aes(x = NES, y = Term_label, size = Set_size, shape = Treatment, group = Treatment)) + geom_point(data = plot_combined[plot_combined$Treatment == 
    "2DG + OX" & plot_combined$FDR <= 0.25, , drop = FALSE], aes(fill = minus_log10_FDR), position = position_nudge(y = 0.13), 
    color = "#202326", stroke = 0.45, alpha = 0.96) + geom_point(data = plot_combined[plot_combined$Treatment == "LDHA + OX" & 
    plot_combined$FDR <= 0.25, , drop = FALSE], aes(fill = minus_log10_FDR), position = position_nudge(y = -0.13), color = "#202326", 
    stroke = 0.45, alpha = 0.96) + geom_point(data = plot_combined[plot_combined$Treatment == "2DG + OX" & plot_combined$FDR > 
    0.25, , drop = FALSE], aes(color = FDR_status), position = position_nudge(y = 0.13), fill = "#B8BDC3", stroke = 0.45, 
    alpha = 0.96) + geom_point(data = plot_combined[plot_combined$Treatment == "LDHA + OX" & plot_combined$FDR > 0.25, , 
    drop = FALSE], aes(color = FDR_status), position = position_nudge(y = -0.13), fill = "#B8BDC3", stroke = 0.45, alpha = 0.96) + 
    scale_shape_manual(name = "Treatment", values = c(`2DG + OX` = 21, `LDHA + OX` = 24)) + plot_scales(title = "2DG + OX and LDHA + OX vs control", 
    subtitle = paste0("edgeR exact test; union of both term lists\nred scale uses this ", "plot's FDR range; gray = FDR > 0.25"), 
    scale_data = plot_combined, include_treatment_legend = TRUE))
  )
)
