# Figure regeneration from deposited tables

Scripts for plotting numerical figure panels from the deposited source tables. Each `Figure_<id>/Panel_<letter>.R` reads the corresponding `Figure_<id>/Panel_<letter>_values.csv`. See the [panel index](PANEL_INDEX.md).

Use R 4.4.1 with the packages listed in [environment.lock](environment.lock). Run from the repository root:

```sh
Rscript --vanilla figure_regeneration/render_all.R DATA_DIR OUTPUT_DIR
Rscript --vanilla figure_regeneration/Figure_3/Panel_C.R DATA_DIR OUTPUT_DIR
```

`DATA_DIR` contains the figure folders and panel CSV files. Image and schematic panels require the original image assets.
