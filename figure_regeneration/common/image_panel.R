# Schematics and representative images use the original supplied panel asset.
image_panel_main <- function(figure, panel, args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)<2)stop("Usage: Rscript PANEL_SCRIPT ASSET_DIR OUTPUT_DIR")
  path <- file.path(args[1],paste0("Figure_",figure),paste0("Panel_",panel,"_image.png"))
  if(!file.exists(path))stop("This panel is an image/schematic. Supply the original asset: ",path)
  if(!requireNamespace("png",quietly=TRUE))stop("Install png before exporting image assets")
  img <- png::readPNG(path)
  out <- file.path(args[2],paste0("Figure_",figure));dir.create(out,recursive=TRUE,showWarnings=FALSE)
  grDevices::cairo_pdf(file.path(out,paste0("Panel_",panel,".pdf")),width=8,height=8*dim(img)[1]/dim(img)[2])
  tryCatch({grid::grid.newpage();grid::grid.raster(img)},finally=grDevices::dev.off())
  file.copy(path,file.path(out,paste0("Panel_",panel,".png")),overwrite=TRUE)
}
