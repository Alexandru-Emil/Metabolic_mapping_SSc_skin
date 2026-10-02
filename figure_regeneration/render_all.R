# Run from the repository root: Rscript figure_regeneration/render_all.R DATA_DIR OUTPUT_DIR
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2)stop("Usage: Rscript figure_regeneration/render_all.R DATA_DIR OUTPUT_DIR")
source("figure_regeneration/common/render_panel.R")
source("figure_regeneration/verify_source_data.R")
verification<-verify_source_data(args[1])
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
results <- list()
for(key in names(panel_registry)) {
  ids <- strsplit(key,"/",fixed=TRUE)[[1]]
  message("Rendering Figure ",ids[1]," panel ",ids[2])
  source(file.path("figure_regeneration",paste0("Figure_",ids[1]),paste0("Panel_",ids[2],".R")),local=new.env(parent=globalenv()))
  spec<-panel_registry[[key]]
  results[[key]]<-list(figure=ids[1],panel=ids[2],rows=spec$rows,table=spec$table,type=spec$type)
}
log <- do.call(rbind,lapply(results,as.data.frame))
write.csv(log,file.path(args[2],"rendered_panels.csv"),row.names=FALSE)
write.csv(verification,file.path(args[2],"verified_source_tables.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(args[2],"sessionInfo.txt"))
cat("Rendered ",length(results)," numerical panels.\n",sep="")
