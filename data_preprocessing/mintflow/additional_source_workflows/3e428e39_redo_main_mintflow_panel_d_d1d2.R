source("data_preprocessing/common/config.R")
#!/usr/bin/env Rscript

source(project_path("external/make_main_mintflow_4panel_report_with_per_donor_spatial_maps.R"), local = FALSE)

baseline_cols <- c(
  "Met_hi_Fib" = "#4a0b0b",
  "Other_Fib" = "grey76",
  "Met_hi_EC" = "#2b7bba",
  "Met_int_EC" = "#55a868",
  "Other" = "grey88"
)

make_two_layer_plot <- function(ec_df, tissue_df, title, point_scale = 1) {
  baseline_map <- ggplot(tissue_df, aes(x = x_centroid, y = y_centroid, color = baseline_label)) +
    geom_point(size = 0.72 * point_scale, alpha = 0.86) +
    scale_color_manual(values = baseline_cols, drop = FALSE, name = "Baseline label") +
    coord_equal() +
    scale_y_reverse() +
    labs(title = "D1  Baseline niche", x = "x", y = "y") +
    theme_paper(8) +
    theme(legend.position = "right")

  response_map <- ggplot() +
    geom_point(
      data = tissue_df %>% filter(combined_met_label != "Met_hi_Fib"),
      aes(x = x_centroid, y = y_centroid),
      color = "grey88",
      size = 0.32 * point_scale,
      alpha = 0.45
    ) +
    geom_point(
      data = tissue_df %>% filter(combined_met_label == "Met_hi_Fib"),
      aes(x = x_centroid, y = y_centroid),
      shape = 21,
      color = "grey10",
      fill = NA,
      size = 0.92 * point_scale,
      stroke = 0.34 * point_scale,
      alpha = 0.9
    ) +
    geom_point(
      data = ec_df,
      aes(x = x_centroid, y = y_centroid, color = delta_prob),
      size = 1.45 * point_scale,
      alpha = 0.98
    ) +
    scale_color_gradient2(
      low = "#2b59c3",
      mid = "grey94",
      high = "#c9332b",
      midpoint = 0,
      limits = c(-vmax, vmax),
      oob = scales::squish,
      name = "Delta Met_hi_EC-like\nprobability"
    ) +
    coord_equal() +
    scale_y_reverse() +
    labs(
      title = "D2  EC response after replacement",
      subtitle = "EC color: replacement - baseline. Black outlines: former Met_hi_Fib positions.",
      x = "x",
      y = "y"
    ) +
    theme_paper(8) +
    theme(
      legend.position = "right",
      plot.subtitle = element_text(size = 7.2, color = "grey25")
    )

  (baseline_map | response_map) +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(face = "bold", size = 10, hjust = 0))
    )
}

panel_d <- make_two_layer_plot(
  ec_df = paired_df,
  tissue_df = all_cells_df,
  title = "D  Baseline niche and EC response to Met_hi_Fib replacement",
  point_scale = 1
)

global_two_layer <- make_two_layer_plot(
  ec_df = paired_df,
  tissue_df = all_cells_df,
  title = "Baseline niche and EC response to Met_hi_Fib replacement",
  point_scale = 1
)

donor_outputs <- donor_outputs %>%
  rowwise() %>%
  mutate(.plot_done = {
    donor <- donor_id
    donor_slug <- sub("^two_layer_spatial_by_donor_", "", sub("\\.png$", "", two_layer_png))
    tissue_d <- all_cells_df %>% filter(donor_id == donor)
    ec_d <- paired_df %>% filter(donor_id == donor)
    two_layer_plot <- make_two_layer_plot(
      ec_df = ec_d,
      tissue_df = tissue_d,
      title = paste0("Baseline niche and EC response: ", donor),
      point_scale = 1.25
    )
    save_plot(two_layer_plot, paste0("two_layer_spatial_by_donor_", donor_slug), 12.6, 5.6, dpi = 500)
    TRUE
  }) %>%
  ungroup() %>%
  select(-.plot_done)

save_plot(panel_d, "panel_D_baseline_niche_and_ec_response", 12.8, 5.7, dpi = 450)
save_plot(panel_d, "panel_D_spatial_delta_map_faceted_by_original_ec_state", 12.8, 5.7, dpi = 450)
save_plot(global_two_layer, "supplemental_two_layer_spatial_niche_and_response", 12.8, 5.7, dpi = 450)

combined <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_layout(widths = c(0.76, 1.52), heights = c(0.84, 1.22))
save_plot(combined, "main_mintflow_replacement_4panel_figure", 15.2, 10.6, dpi = 450)

donor_tabs <- paste0(
  "<div class='tabs'>",
  paste0(
    "<button class='tablink' onclick=\"openTab(event, '", donor_outputs$tab_id, "')\">",
    html_escape(donor_outputs$donor_id),
    "</button>",
    collapse = ""
  ),
  "</div>",
  paste0(
    "<div id='", donor_outputs$tab_id, "' class='tabcontent'>",
    "<h3>", html_escape(donor_outputs$donor_id), "</h3>",
    "<p>Cells: ", donor_outputs$n_cells,
    "; baseline Met_hi_Fib: ", donor_outputs$n_met_hi_fib,
    "; ECs with paired MintFlow probabilities: ", donor_outputs$n_ec, ".</p>",
    "<p><a href='", donor_outputs$two_layer_pdf, "'>Baseline niche plus EC response PDF</a> | ",
    "<a href='", donor_outputs$response_pdf, "'>Standalone EC response PDF</a></p>",
    "<h4>Baseline niche and EC response side by side</h4>",
    "<img src='", donor_outputs$two_layer_png, "'>",
    "<h4>Supplemental standalone EC response with full tissue context</h4>",
    "<img src='", donor_outputs$response_png, "'>",
    "</div>",
    collapse = ""
  )
)

html_file <- file.path(out_dir, "main_mintflow_replacement_4panel_report.html")
html <- paste0(
  "<!doctype html><html><head><meta charset='utf-8'>",
  "<title>Main MintFlow Met_hi_Fib Replacement Figure</title>",
  "<style>body{font-family:Arial,sans-serif;margin:32px;color:#222;line-height:1.45} img{border:1px solid #ddd;margin-bottom:24px;max-width:100%} code{background:#f4f4f4;padding:1px 4px}.tabs{display:flex;flex-wrap:wrap;gap:6px;margin:16px 0}.tablink{border:1px solid #bbb;background:#f4f4f4;padding:8px 12px;cursor:pointer}.tablink.active{background:#222;color:#fff}.tabcontent{display:none;border-top:1px solid #ddd;padding-top:12px}</style>",
  "<script>function openTab(evt,id){var i,t,b;t=document.getElementsByClassName('tabcontent');for(i=0;i<t.length;i++){t[i].style.display='none';}b=document.getElementsByClassName('tablink');for(i=0;i<b.length;i++){b[i].className=b[i].className.replace(' active','');}document.getElementById(id).style.display='block';evt.currentTarget.className+=' active';}window.onload=function(){var b=document.getElementsByClassName('tablink');if(b.length){b[0].click();}}</script>",
  "</head><body>",
  "<h1>Main MintFlow Met_hi_Fib Replacement Figure</h1>",
  "<p>This report reorganizes the key MintFlow result from <code>mintflow_met_hi_fib_perturbation_ec_report_with_frequency_dotplots.html</code> into a manuscript-focused 4-panel layout. Panel D is split into D1 baseline niche and D2 EC response after replacement: former Met_hi_Fib positions are overlaid as black outlines, all other cells are faint gray, and ECs are colored by replacement minus baseline Met_hi_EC-like probability.</p>",
  "<p><b>Source folder:</b> <code>", html_escape(final_dir), "</code></p>",
  "<p><b>Coordinate/label source:</b> <code>", html_escape(obs_path), "</code></p>",
  "<h2>Combined Figure</h2>",
  "<h3>main_mintflow_replacement_4panel_figure.png</h3><img src='main_mintflow_replacement_4panel_figure.png'>",
  "<h2>Separate Panels</h2>",
  "<h3>panel_A_perturbation_schematic.png</h3><img src='panel_A_perturbation_schematic.png'>",
  "<h3>panel_B_replacement_ec_state_frequency.png</h3><img src='panel_B_replacement_ec_state_frequency.png'>",
  "<h3>panel_C_per_cell_probability_shift_by_original_ec_state.png</h3><img src='panel_C_per_cell_probability_shift_by_original_ec_state.png'>",
  "<h3>panel_D_baseline_niche_and_ec_response.png</h3><img src='panel_D_baseline_niche_and_ec_response.png'>",
  "<h2>Per-Donor Spatial Maps</h2>",
  "<p>Each tab shows one donor/sample separately, avoiding overlap from plotting all donors in one coordinate system.</p>",
  donor_tabs,
  "<h2>Supplemental Two-Layer Spatial Map</h2>",
  "<p>This plot first shows the baseline tissue niche and then the EC response after Met_hi_Fib replacement.</p>",
  "<h3>supplemental_two_layer_spatial_niche_and_response.png</h3><img src='supplemental_two_layer_spatial_niche_and_response.png'>",
  zoom_html,
  "</body></html>"
)
writeLines(html, html_file)

message(html_file)
message(file.path(out_dir, "main_mintflow_replacement_4panel_figure.pdf"))
message(file.path(out_dir, "panel_D_baseline_niche_and_ec_response.pdf"))
message(file.path(out_dir, "spatial_map_outputs_by_donor.csv"))
