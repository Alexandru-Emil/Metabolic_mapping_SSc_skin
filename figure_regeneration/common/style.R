suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
theme_panel <- function(base_size = 10) {
  theme_classic(base_size = base_size) + theme(
    plot.title = element_text(face = "bold", hjust = .5),
    axis.text = element_text(colour = "black"),
    legend.title = element_text(size = base_size),
    strip.background = element_blank(), strip.text = element_text(face = "bold"))
}
group_colours <- c(Healthy = "#6CC3F4", `SSc stable` = "#F6BF93", `SSc progressive` = "#C03830", SSc = "#C03830", healthy = "#6CC3F4", stable = "#F6BF93", progressive = "#C03830", Met_hi_Fib = "#F8766D", Other_Fib = "#00BFC4", Met_hi_EC = "#E69F00", Methi_EC = "#E69F00", Other_EC = "#56B4E9")
group_colours <- c(group_colours, `Stable SSc` = "#F6BF93", `Progressive SSc` = "#C03830")
cluster_colours <- c(stromal_1 = "#BD2424", stromal_2 = "#C2DAB7", stromal_3 = "#FC8615", stromal_4 = "#7DB9D1", stromal_5 = "#BE72A5", stromal_6 = "#146B9E", stromal_7 = "#F6D32F", stromal_8 = "#277D4B", CD31_1 = "#FBB4AE", CD31_2 = "#B3CDE3", CD31_3 = "#CCEBC5", CD31_4 = "#DECBE4", CD31_5 = "#FED9A6", CD45_CD68_1 = "#E41A1C", CD45_CD68_2 = "#377EB8", CD45_CD68_3 = "#4DAF4A", CD45_CD68_4 = "#984EA3")
ordered_groups <- function(x) {
  preferred <- c("Healthy", "healthy", "Stable SSc", "SSc stable", "stable", "Progressive SSc", "SSc progressive", "progressive", "SSc", "Met_hi_Fib", "Other_Fib", "Met_hi_EC", "Methi_EC", "Other_EC")
  factor(x, levels = c(intersect(preferred, unique(x)), setdiff(unique(x), preferred)))
}
colours_for <- function(x, palette = group_colours) {
  lev <- if (is.factor(x)) levels(x) else unique(as.character(x))
  col <- setNames(scales::hue_pal()(length(lev)), lev)
  known <- intersect(lev, names(palette)); col[known] <- palette[known]; col
}
display_labels <- function(x) {
  out <- x
  if (exists("state_labels", inherits = TRUE)) {
    k <- match(x, names(state_labels)); out[!is.na(k)] <- state_labels[k[!is.na(k)]]
  }
  out <- gsub("<br>", "\n", out, fixed = TRUE)
  out <- gsub("<[^>]+>", "", out)
  out <- gsub("Met_hi_", "Met (hi) ", out, fixed = TRUE)
  out <- gsub("_", " ", out, fixed = TRUE)
  out
}
