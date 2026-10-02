# Run from the repository root: Rscript figure_regeneration/render_all.R DATA_DIR OUTPUT_DIR
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2)stop("Usage: Rscript figure_regeneration/render_all.R DATA_DIR OUTPUT_DIR")
source("figure_regeneration/common/panel_registry.R")
source("figure_regeneration/verify_source_data.R")
verification<-verify_source_data(args[1])
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
repo_dir<-normalizePath(getwd(),winslash="/")
data_dir<-normalizePath(args[1],winslash="/")
output_dir<-normalizePath(args[2],winslash="/")
dir.create(file.path(output_dir,"logs"),showWarnings=FALSE)
rscript<-file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe" else "Rscript")
jobs<-lapply(names(panel_registry),function(key){
  ids<-strsplit(key,"/",fixed=TRUE)[[1]];spec<-panel_registry[[key]]
  list(figure=ids[1],panel=ids[2],rows=spec$rows,table=spec$table,type=spec$type,
       script=file.path("figure_regeneration",paste0("Figure_",ids[1]),paste0("Panel_",ids[2],".R")))
})
run_named_panel<-function(job,rscript,data_dir,output_dir,repo_dir){
  setwd(repo_dir)
  log_file<-file.path(output_dir,"logs",paste0("Figure_",job$figure,"_Panel_",job$panel,".log"))
  status<-system2(rscript,c("--vanilla",shQuote(job$script),shQuote(data_dir),shQuote(output_dir)),stdout=log_file,stderr=log_file)
  if(status!=0)stop("Figure ",job$figure," panel ",job$panel," failed; see ",log_file)
  job$script<-NULL
  job
}
# Each named script starts in its own R session. This prevents graphics-device
# and font state from earlier panels affecting later plots in a full run.
workers<-as.integer(Sys.getenv("METABOLIC_PLOT_WORKERS","4"))
if(is.na(workers)||workers<1L||workers>8L)stop("METABOLIC_PLOT_WORKERS must be between 1 and 8")
cluster<-parallel::makeCluster(workers)
results<-tryCatch({
  completed<-list()
  batches<-split(jobs,(seq_along(jobs)-1L)%/%workers)
  for(batch in batches){
    completed<-c(completed,parallel::clusterApply(cluster,batch,run_named_panel,rscript=rscript,
      data_dir=data_dir,output_dir=output_dir,repo_dir=repo_dir))
    message("Rendered ",length(completed)," / ",length(jobs)," numerical panels")
  }
  completed
},finally=parallel::stopCluster(cluster))
log <- do.call(rbind,lapply(results,as.data.frame))
write.csv(log,file.path(args[2],"rendered_panels.csv"),row.names=FALSE)
write.csv(verification,file.path(args[2],"verified_source_tables.csv"),row.names=FALSE)
writeLines(c("Each panel runs in a separate R session. See Figure_<id>/Panel_<id>_sessionInfo.txt for its loaded package versions.",
             "Orchestrator session:",capture.output(sessionInfo())),file.path(args[2],"sessionInfo.txt"))
cat("Rendered ",length(results)," numerical panels.\n",sep="")
