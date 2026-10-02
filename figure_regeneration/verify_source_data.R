# Verify the deposited numerical inputs before regenerating panels.
verify_source_data <- function(data_dir) {
  if(!requireNamespace("digest",quietly=TRUE))stop("Install digest for SHA-256 verification")
  source("figure_regeneration/common/panel_registry.R",local=TRUE)
  source("figure_regeneration/common/source_table_hashes.R",local=TRUE)
  results<-list()
  for(key in names(panel_registry)) {
    spec<-panel_registry[[key]];path<-file.path(data_dir,spec$table)
    if(!file.exists(path))stop("Missing panel table: ",spec$table)
    d<-read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
    stopifnot(nrow(d)==spec$rows,identical(names(d),spec$columns))
    expected<-unname(source_table_hashes[spec$table])
    observed<-digest::digest(file=path,algo="sha256")
    if(length(expected)!=1||is.na(expected)||observed!=expected)stop("Source-table hash mismatch: ",spec$table)
    if(all(c("frequency_fraction","frequency_percent")%in%names(d)))stopifnot(all(abs(d$frequency_fraction-d$frequency_percent/100)<1e-12,na.rm=TRUE))
    if(all(c("cluster_cells","parent_cells","frequency_fraction")%in%names(d)))stopifnot(all(abs(d$frequency_fraction-d$cluster_cells/d$parent_cells)<1e-12,na.rm=TRUE))
    if(all(c("cell_count","cluster_total","fraction")%in%names(d)))stopifnot(all(abs(d$fraction-d$cell_count/d$cluster_total)<1e-12,na.rm=TRUE))
    if("UMAP1"%in%names(d))stopifnot(all(is.finite(d$UMAP1)),all(is.finite(d$UMAP2)))
    if("recorded_Gender"%in%names(d))stopifnot(all(d$recorded_Gender%in%c("F","M")))
    results[[key]]<-data.frame(panel=key,rows=nrow(d),sha256=observed)
  }
  validation<-read.csv(file.path(data_dir,"Figure_5/Panel_B_values.csv"));stopifnot(length(unique(validation$donor))==6,nrow(validation)==1155,sum(validation$metfiblabel=="Met_hi_Fib")==622)
  projection<-read.csv(file.path(data_dir,"Figure_S9/Panel_D_values.csv"));stopifnot(nrow(projection)==21812)
  collagen<-read.csv(file.path(data_dir,"Figure_S11/Panel_B_values.csv"));stopifnot(identical(sort(as.integer(table(collagen$condition))),rep(8L,3)))
  do.call(rbind,results)
}
if(sys.nframe()==0L) {
  args<-commandArgs(trailingOnly=TRUE);if(length(args)!=1)stop("Usage: Rscript figure_regeneration/verify_source_data.R DATA_DIR")
  results<-verify_source_data(args[1]);cat("PASS: ",nrow(results)," panel tables; ",sum(results$rows)," numerical rows.\n",sep="")
}
