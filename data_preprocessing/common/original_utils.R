# Scientific helper functions from the original IMC utils.R; installer omitted.
umap_nn_fx <- function(sem, umapreduc, n_neighbors = 100) {
    umapd <- data.table(sem@reductions[[umapreduc]]@cell.embeddings, keep.rownames = TRUE)
    setnames(umapd, c(names(umapd)[1:3]), c("cell_ID", "UMAP_1", "UMAP_2"))
    nn_umap <- RANN::nn2(umapd[, .(UMAP_1, UMAP_2)], k = n_neighbors + 1)$nn.idx
    nn_umap <- data.table::melt(cbind(umapd[, .(cell_ID)], data.table(nn_umap)), id.vars = c("cell_ID", 
        "V1"))
    colnames(nn_umap) <- c("cell_ID1", "cell_ID1_idx", "neighbor", "cell_ID2_idx")
    nn_umap <- merge(nn_umap, nn_umap[, .(cell_ID2 = cell_ID1, cell_ID2_idx = cell_ID1_idx)][, unique(.SD)], 
        by = "cell_ID2_idx")
    wumap <- Matrix::sparseMatrix(i = c(unique(nn_umap$cell_ID1_idx), nn_umap$cell_ID2_idx), j = c(unique(nn_umap$cell_ID1_idx), 
        nn_umap$cell_ID1_idx), x = 1)
    rownames(wumap) <- colnames(wumap) <- nn_umap[order(nn_umap$cell_ID1_idx), unique(cell_ID1)]
    mumap <- Matrix::sparseMatrix(i = 1:ncol(wumap), j = 1:ncol(wumap), x = 1/Matrix::colSums(wumap))
    dimnames(mumap) <- dimnames(wumap)
    smoother <- mumap %*% wumap
    smoother <- smoother[, colnames(sem)]
    return(smoother)
}
extract_and_clean_genes_fx <- function(data, gene_col = "gene", cluster_col = "cluster", top_n = 10) {
    genes_per_cluster <- tapply(data[[gene_col]], list(data[[cluster_col]]), function(i) paste0(i[1:top_n], 
        collapse = ","))
    genes_unlisted <- unname(unlist(strsplit(unlist(genes_per_cluster), ",")))
    genes_cleaned <- genes_unlisted[genes_unlisted != "NA"]
    genes_cleaned <- genes_cleaned[!is.na(genes_cleaned)]
    return(genes_cleaned)
}
subset_and_count_cells_fx <- function(spe, cluster_col) {
    unique_clusters <- unique(colData(spe)[[cluster_col]])
    cell_counts <- numeric(length(unique_clusters))
    names(cell_counts) <- unique_clusters
    for (cluster in unique_clusters) {
        subset1_spe <- spe[, colData(spe)[[cluster_col]] == cluster]
        cell_counts[as.character(cluster)] <- ncol(subset1_spe)
    }
    return(cell_counts)
}
sort_and_cumulate_fx <- function(cell_counts, matching_values_sort) {
    sorted_cell_counts <- cell_counts[match(matching_values_sort, names(cell_counts))]
    cumulative_counts <- sorted_cell_counts
    for (i in 2:length(sorted_cell_counts)) {
        cumulative_counts[i] <- cumulative_counts[i] + cumulative_counts[i - 1]
    }
    return(cumulative_counts)
}
create_beautiful_radarchart_fx <- function(data, color = "#00AFBB", vlabels = colnames(data), vlcex = 0.7, 
    caxislabels = NULL, title = NULL, ...) {
    radarchart(data, axistype = 1, pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1, 
        cglcol = "grey", cglty = 1, cglwd = 0.8, axislabcol = "grey", vlcex = vlcex, vlabels = vlabels, 
        caxislabels = caxislabels, title = title, ...)
}
output_dir_svg_plots_fx <- function(ts) {
    output_dir_base <- getwd()
    if (!is.null(knitr::current_input())) {
        md_filename <- tools::file_path_sans_ext(basename(knitr::current_input()))
    }
    else {
        md_filename <- "interactive_session"
    }
    output_dir <- file.path(output_dir_base, paste0(md_filename, "_svg_files_", ts))
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    return(output_dir)
}
get_common_boxplot_components_fx <- function() {
    list(geom_point(position = position_dodge2(width = 0.5, padding = 1.5)), theme_minimal(), theme(axis.text.x = element_text(size = 14, 
        face = "bold", angle = 45, hjust = 1), axis.text.y = element_text(size = 14, face = "bold"), 
        legend.text = element_text(size = 14, face = "bold"), axis.title = element_text(size = 14, face = "bold"), 
        legend.title = element_text(size = 14, face = "bold"), plot.title = element_text(size = 16, face = "bold"), 
        strip.text = element_text(size = 14, face = "bold")))
}
get_common_UMAP_components_fx <- function(p) {
    p + guides(color = guide_legend(override.aes = list(size = 5), title = NULL)) + theme(legend.text = element_text(size = 16), 
        legend.key.size = unit(1.5, "lines"), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) + coord_fixed(ratio = 1)
}
sanitize_filename_fx <- function(filename) {
    filename <- gsub(" ", "_", filename)
    filename <- gsub("[/\\:*?\"<>|+]", "", filename)
    filename <- gsub("[^A-Za-z0-9_\\-]", "", filename)
    return(filename)
}
plot_conti_ht_fx <- function(tbl, xlab, ylab) {
    tmp <- reshape2::melt(tbl)
    p <- ggplot(tmp, aes(x = Var2, y = Var1, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", 
        high = "blue") + labs(x = xlab, y = ylab, fill = "prop") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, 
        hjust = 1))
    return(p)
}
replace_duplicate_colors_fx <- function(color_vector) {
    unique_colors <- unique(color_vector)
    if (length(unique_colors) < length(color_vector)) {
        all_unique_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(color_vector) * 2)
        remaining_colors <- setdiff(all_unique_colors, unique_colors)
        used_colors <- setNames(character(0), character(0))
        for (i in names(color_vector)) {
            if (color_vector[i] %in% used_colors) {
                color_vector[i] <- remaining_colors[1]
                remaining_colors <- remaining_colors[-1]
            }
            used_colors <- c(used_colors, color_vector[i])
        }
    }
    return(color_vector)
}
generate_spatial_plots_fx <- function(spe, cnames, color_vector_nsp, color_vector_bank) {
    plot_nsp <- plotColData(spe, x = "x", y = "y", point_size = 0.6, colour_by = cnames[1]) + scale_color_manual(values = color_vector_nsp) + 
        labs(title = paste0("Non-spatial clusters \n", celltype_used))
    plot_bank <- plotColData(spe, x = "x", y = "y", point_size = 0.6, colour_by = cnames[2]) + scale_color_manual(values = color_vector_bank) + 
        labs(title = paste0("BANKSY clusters \n", celltype_used))
    return(list(plot_nsp = plot_nsp, plot_bank = plot_bank))
}
rescale_quantiles_asinh_fx <- function(x, cofactor = 1) {
    x_transformed <- asinh(x/cofactor)
    q5 <- quantile(x_transformed, 0.05, na.rm = TRUE)
    q95 <- quantile(x_transformed, 0.95, na.rm = TRUE)
    x_rescaled <- (x_transformed - q5)/(q95 - q5)
    x_rescaled[x_rescaled < 0] <- 0
    x_rescaled[x_rescaled > 1] <- 1
    return(x_rescaled)
}
plot_violin_spatial_fx <- function(sfe, feature, sample_id = "all", size = 0.5) {
    violin <- plotColData(sfe, feature, point_fun = function(...) list())
    spatial <- plotSpatialFeature(sfe, feature, colGeometryName = "centroids", scattermore = TRUE, sample_id = sample_id, 
        size = size)
    violin + spatial + plot_layout(widths = c(1, 2))
}
apply_thresholds_fx <- function(sfe, donor, dapi_thresholds, area_thresholds) {
    threshold_DAPI <- unname(dapi_thresholds[donor])
    threshold_Area <- unname(area_thresholds[donor])
    pass_DAPI <- sfe$Mean.DAPI > threshold_DAPI
    pass_Area <- sfe$area < threshold_Area
    sfe$pass_both_DAPI_area <- pass_DAPI & pass_Area
    return(sfe)
}
apply_thresholds_to_list_fx <- function(sfe_list, dapi_thresholds, area_thresholds) {
    for (i in seq_along(sfe_list)) {
        donor <- names(sfe_list)[i]
        sfe_list[[i]] <- apply_thresholds_fx(sfe_list[[i]], donor, dapi_thresholds, area_thresholds)
    }
    return(sfe_list)
}
exclude_sfe_noncoding_genes_fx <- function(sfe, patterns = c("^Neg", "^SystemControl")) {
    if (is.null(rownames(sfe))) {
        stop("'sfe' does not have row names. Please assign row names before using this function.")
    }
    n_before <- nrow(sfe)
    combined_pattern <- paste(patterns, collapse = "|")
    rows_to_exclude <- str_detect(rownames(sfe), combined_pattern)
    n_excluded <- sum(rows_to_exclude)
    if (n_excluded > 0) {
        sfe_subset <- sfe[!rows_to_exclude, ]
        n_after <- nrow(sfe_subset)
        message(sprintf("Excluding non-coding genes:\n- Genes before exclusion: %d\n- Genes excluded: %d\n- Genes after exclusion: %d", 
            n_before, n_excluded, n_after))
        return(sfe_subset)
    }
    else {
        message(sprintf("No genes matching the specified patterns were found in 'sfe'.\n- Genes before: %d\n- Genes after: %d", 
            n_before, n_before))
        return(sfe)
    }
}
check_coldata_compatibility_fx <- function(s4_object, auto_convert = FALSE, verbose = TRUE) {
    if (!"colData" %in% slotNames(s4_object)) {
        stop("The provided object does not have a 'colData' slot.")
    }
    col_data <- colData(s4_object)
    if (is.null(col_data)) {
        stop("The 'colData' slot is NULL.")
    }
    problematic_cols <- list()
    cols_to_drop <- logical(length = ncol(col_data))
    names(cols_to_drop) <- colnames(col_data)
    for (col_name in colnames(col_data)) {
        column_data <- col_data[[col_name]]
        col_class <- class(column_data)
        if (is.list(column_data) || isS4(column_data) || any(sapply(column_data, isS4))) {
            problematic_cols[[col_name]] <- col_class
            if (auto_convert) {
                if (is.list(column_data) && all(sapply(column_data, length) == 1)) {
                  col_data[[col_name]] <- unlist(column_data)
                  if (verbose) {
                    message(sprintf("Auto-converted column '%s' from list to vector.", col_name))
                  }
                }
                else {
                  cols_to_drop[col_name] <- TRUE
                  if (verbose) {
                    message(sprintf("Column '%s' could not be auto-converted and will be dropped.", col_name))
                  }
                }
            }
        }
    }
    if (auto_convert && any(cols_to_drop)) {
        col_data <- col_data[, !cols_to_drop, drop = FALSE]
        if (verbose) {
            message("Dropped columns: ", paste(names(cols_to_drop)[cols_to_drop], collapse = ", "))
        }
    }
    if (length(problematic_cols) > 0) {
        if (verbose) {
            message("\342\232\240\357\270\217 The following columns in 'colData' may cause compatibility issues:")
            for (col in names(problematic_cols)) {
                message(sprintf("- '%s' (class: %s)", col, paste(problematic_cols[[col]], collapse = ", ")))
            }
            if (!auto_convert) {
                message("Consider converting these columns to atomic vectors or handling them before conversion.")
            }
        }
    }
    else {
        if (verbose) {
            message("\342\234\205 No compatibility issues detected in 'colData'.")
        }
    }
    return(list(problematic_cols = problematic_cols, coldata_fixed = col_data))
}
