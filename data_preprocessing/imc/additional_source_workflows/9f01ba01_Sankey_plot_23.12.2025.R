
## Sankey plots:
##   (1) gates -> lv2_anno_from_sfe
##   (2) (stromal) Rphenograph_500 -> lv3_anno_from_sfe
## ============================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggsankey)
})

## ----------------------------
## USER INPUTS (edit these)
## ----------------------------
# Objects assumed already in memory:
#   spe : SpatialExperiment with transferred columns

GATES_COL <- "gates"
LV2_COL   <- "lv2_anno_from_sfe"
LV3_COL   <- "lv3_anno_from_sfe"
CLUST_COL <- "Rphenograph_500"
STROMAL_LABEL <- "Stromal"

## ----------------------------
## 1) gates -> lv2 sankey
## ----------------------------
df_lv2 <- as.data.frame(colData(spe)) %>%
  filter(!is.na(.data[[LV2_COL]])) %>%
  dplyr::select(all_of(c(GATES_COL, LV2_COL)))

sankey_lv2 <- make_long(df_lv2, !!sym(GATES_COL), !!sym(LV2_COL))

p_lv2 <- ggplot(
  sankey_lv2,
  aes(x = x, next_x = next_x,
      node = node, next_node = next_node,
      fill = factor(node),
      label = node)
) +
  geom_sankey(flow.alpha = 0.4, node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  theme_sankey(base_size = 16) +
  theme(legend.position = "none") +
  labs(x = NULL)

print(p_lv2)

## ----------------------------
## 2) stromal: Rphenograph_500 -> lv3 sankey
## ----------------------------
df_lv3 <- as.data.frame(colData(spe)) %>%
  filter(.data[[GATES_COL]] == STROMAL_LABEL,
         !is.na(.data[[LV3_COL]]),
         !is.na(.data[[CLUST_COL]])) %>%
  dplyr::select(all_of(c(CLUST_COL, LV3_COL)))

sankey_lv3 <- make_long(df_lv3, !!sym(CLUST_COL), !!sym(LV3_COL))

p_lv3 <- ggplot(
  sankey_lv3,
  aes(x = x, next_x = next_x,
      node = node, next_node = next_node,
      fill = factor(node),
      label = node)
) +
  geom_sankey(flow.alpha = 0.4, node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  theme_sankey(base_size = 16) +
  theme(legend.position = "none") +
  labs(x = NULL)

print(p_lv3)

## ----------------------------
## Optional: save plots
## ----------------------------
# ggsave("sankey_gates_to_lv2.png", p_lv2, width = 10, height = 6, dpi = 300)
# ggsave("sankey_cluster_to_lv3_stromal.png", p_lv3, width = 10, height = 6, dpi = 300)
