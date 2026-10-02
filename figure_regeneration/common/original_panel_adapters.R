# These adapters supply deposited tables to the retained original plot code.
source("figure_regeneration/common/original_plot_code.R")
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(tibble)})
original_environment <- function(module) {
  env <- new.env(parent=globalenv())
  eval(original_modules[[module]]$setup,env)
  env
}
original_expression <- function(module,name,env) {
  expr <- original_modules[[module]]$expressions[[name]]
  if(is.null(expr))expr<-original_modules[[module]]$plots[[name]]
  if(is.null(expr))stop("Original expression not found: ",module,"/",name)
  eval(expr,env)
}
original_six_donor <- function(d,key) {
  e <- original_environment("six")
  e$matched_donor_levels <- paste0("Validation",1:6)
  if(key %in% c("5/B","S14/A")) {
    group <- if(key=="5/B")"metfiblabel" else "metEClabel"
    d[[group]]<-factor(d[[group]],levels=if(group=="metfiblabel")e$fib_levels else e$ec_levels)
    if(key=="5/B")e$fib_umap<-d else e$ec_umap<-d
    return(original_expression("six",if(key=="5/B")"p_fib_umap" else "p_ec_umap",e))
  }
  if(key=="5/C") {
    e$fib_umap_scores<-d
    for(n in c("p_fib_glycolysis","p_fib_tca","p_fib_scores"))e[[n]]<-original_expression("six",n,e)
    return(e$p_fib_scores)
  }
  if(key=="5/E") {
    d$metfiblabel<-factor(d$metfiblabel,levels=e$fib_levels);d$marker<-factor(d$marker,levels=e$protein_markers)
    e$fib_protein_sample<-d
    stat<-e$paired_test(d,"marker","mean_expression_z","metfiblabel","Met_hi_Fib","Other_Fib")
    e$fib_protein_stats<-stat %>% mutate(marker=factor(feature,levels=e$protein_markers))
    e$fib_protein_stat_positions<-original_expression("six","fib_protein_stat_positions",e)
    return(original_expression("six","p_fib_protein",e))
  }
  if(key=="S9/B") {
    d$metfiblabel<-factor(d$metfiblabel,levels=e$fib_levels)
    e$fib_score_sample<-d
    stat<-e$paired_test(d,"pathway","mean_score","metfiblabel","Met_hi_Fib","Other_Fib")
    e$fib_score_stats<-stat %>% mutate(pathway=feature)
    e$fib_score_stat_positions<-original_expression("six","fib_score_stat_positions",e)
    for(n in c("p_fib_sample_glycolysis","p_fib_sample_tca","p_fib_sample_scores"))e[[n]]<-original_expression("six",n,e)
    return(e$p_fib_sample_scores)
  }
  if(key %in% c("S9/A","S14/B")) {
    fib<-key=="S9/A";group<-if(fib)"metfiblabel" else "metEClabel";lev<-if(fib)e$fib_levels else e$ec_levels
    matrix<-heatmap_matrix(d,"marker",group,"mean_expression_z")[e$metabolic_markers,lev,drop=FALSE]
    counts<-d[!duplicated(d[[group]]),c(group,"n_cells")];nn<-setNames(counts$n_cells,counts[[group]])
    if(fib) {
      e$fib_heatmap_matrix<-matrix
      e$fib_heatmap_annotation<-data.frame(metfiblabel=factor(lev,levels=lev),mols=as.numeric(nn[lev]),row.names=lev)
      e$fib_heatmap_annotation_colors<-original_expression("six","fib_heatmap_annotation_colors",e)
      return(original_expression("six","fib_heatmap_object",e)$gtable)
    }
    e$ec_heatmap_matrix<-matrix
    e$ec_heatmap_annotation<-data.frame(ids=factor(lev,levels=lev),ncells=as.numeric(nn[lev]),row.names=lev)
    e$ec_heatmap_annotation_colors<-list(ncells=colorRampPalette(c("#F7F0F5","#CC79A7"))(100),ids=e$ec_colors)
    return(original_expression("six","ec_heatmap_object",e)$gtable)
  }
  stop("Unsupported six-donor panel")
}
original_signal <- function(d,key) {
  e<-original_environment("signal")
  e$manuscript_signal_data<-d[d$fit_status=="ok"&is.finite(d$log2_signal)&is.finite(d$log2_snr),]
  e$manuscript_signal_labels<-e$manuscript_signal_data
  e$dumbbell_data<-d %>% mutate(log2_background=log2(background_mean),marker_canonical=factor(marker_canonical,levels=rev(unique(marker_canonical))))
  original_expression("signal",if(key=="S2/A")"manuscript_signal_plot" else "manuscript_background_plot",e)
}
original_energy <- function(d,key) {
  e<-original_environment("energy")
  e$pathway_levels<-unique(d$pathway);e$pathway_labels<-setNames(e$pathway_levels,e$pathway_levels)
  e$group_colors<-colours_for(d$group)
  if(key%in%c("S9/C","S14/C")) {
    tab<-d %>% group_by(pathway,group) %>% summarise(mean_score=mean(score),.groups="drop")
    return(e$make_energy_dotplot(tab,list(group_levels=if(key=="S9/C")c("Other_Fib","Met_hi_Fib")else c("Other_EC","Met_hi_EC")),"Xenium-based scMetabolism scores","Mean AUCell score"))
  }
  e$make_violin_plot(d,"Xenium-based scMetabolism scores","AUCell score")
}
make_original_panel <- function(d,spec,key,data_dir) {
  if(key=="S14/H") {
    e<-original_environment("acta")
    tab<-d%>%group_by(pathway,acta2_ec_group)%>%summarise(mean_score=mean(score),.groups="drop")
    return(e$make_dotplot(tab,FALSE))
  }
  if(key=="5/G") {
    e<-original_environment("fib")
    d<-d[d$score_name=="CoreMatrisome score",]
    e$score_levels<-"CoreMatrisome score"
    e$lv3_levels<-e$lv3_preferred_order
    return(e$make_presentation_mean_dot_matrix(d,FALSE))
  }
  if(key=="5/D") {
    # IMC_heatmap_sankey_new.Rmd, selected-gene plot p2: unchanged geometry/scales.
    d$id<-factor(d$id,levels=c("Other_Fib","Met_hi_Fib"))
    d$features.plot<-factor(d$features.plot,levels=unique(d$features.plot))
    return(ggplot(d,aes(x=id,y=features.plot))+geom_point(aes(size=pct.exp,color=avg.exp.scaled))+
      scale_color_gradientn(colors=rev(RColorBrewer::brewer.pal(7,"RdBu")))+scale_size(range=c(1,8))+
      scale_y_discrete(limits=rev)+scale_x_discrete(labels=expression("Other_Fib",plain(Met)^hi*"_Fib"))+
      guides(colour=guide_colourbar(order=1),size=guide_legend(order=2))+
      labs(x=NULL,y=NULL,color="Average expression",size="Percent expressed")+theme_bw()+
      theme(axis.text.x=element_text(size=18,angle=0,hjust=.5),axis.text.y=element_text(size=18),legend.title=element_text(size=18),legend.text=element_text(size=16)))
  }
  if(key=="S14/F") {
    e<-original_environment("ec")
    e$score_long<-d[d$score_name%in%c("CoreMatrisome score","Collagens score"),]
    e$score_levels<-c("CoreMatrisome score","Collagens score")
    return(original_expression("ec","p_violin",e)+facet_wrap(~score_name,ncol=1,scales="free_y",labeller=labeller(score_name=e$score_labels)))
  }
  if(key %in% c("6/B","6/C"))return(original_rnaseq_heatmap(d,key))
  if(key %in% c("5/B","5/C","5/E","S9/A","S9/B","S14/A","S14/B"))return(original_six_donor(d,key))
  if(key %in% c("S2/A","S2/B"))return(original_signal(d,key))
  if(key %in% c("S9/C","S14/C"))return(original_energy(d,key))
  if(spec$type=="gsea") {
    e<-original_environment("gsea")
    both<-bind_rows(read.csv(file.path(data_dir,"Figure_6/Panel_D_values.csv")),read.csv(file.path(data_dir,"Figure_6/Panel_E_values.csv")))
    e$x_limits<-c(floor(min(c(0,both$NES))*10)/10-.08,ceiling(max(c(0,both$NES))*10)/10+.08)
    e$x_breaks<-pretty(e$x_limits,n=6);e$size_limits<-range(both$Set_size);e$size_breaks<-scales::breaks_pretty(n=3)(e$size_limits)
    plot_data<-e$prepare_plot_data(d,d$Term)
    return(e$single_treatment_plot(plot_data,as.character(d$Treatment[1])))
  }
  if(grepl("^S10/|^S11/|^S12/",key))return(original_functional_bar(d,key))
  if(key=="2/C" || (spec$type=="box" && "frequency_fraction"%in%names(d)))return(original_imc_box(d,spec$x,spec$y,spec$facet,frequency=key!="2/C"))
  # The original IMC notebooks use dittoDimPlot for saved embeddings and scores.
  if(spec$type %in% c("umap","score_umap") && !grepl("^5/|^S14/A",key)) {
    require_packages_panel(c("SingleCellExperiment","dittoSeq"))
    obj<-SingleCellExperiment::SingleCellExperiment(assays=list(counts=matrix(0,1,nrow(d))),colData=S4Vectors::DataFrame(d))
    colnames(obj)<-if("cell"%in%names(d))d$cell else paste0("PlotPoint",seq_len(nrow(d)))
    SingleCellExperiment::reducedDim(obj,"UMAP")<-as.matrix(d[,c("UMAP1","UMAP2")])
    if(spec$type=="umap") {
    return(dittoSeq::dittoDimPlot(obj,var=spec$colour,reduction.use="UMAP",size=.2,color.panel=colours_for(d[[spec$colour]],c(cluster_colours,group_colours)),do.label=FALSE)+labs(colour=NULL)+scale_colour_manual(values=colours_for(d[[spec$colour]],c(cluster_colours,group_colours)),labels=display_labels))
    }
    return(wrap_plots(lapply(spec$values,function(v)dittoSeq::dittoDimPlot(obj,var=v,reduction.use="UMAP",size=.2)+ggtitle(display_labels(v))),ncol=2))
  }
  NULL
}
original_rnaseq_heatmap <- function(d,key) {
  # build_heatmap in 15_all_feature_strict_deg_figure_heatmaps_edgeR_exact.R.
  require_packages_panel("pheatmap")
  d<-d[order(d$Cluster_rank),]
  mat<-heatmap_matrix(d,"Feature_ID","sample","row_z_score")
  sample_metadata<-d[!duplicated(d$sample),c("sample","condition")]
  condition<-sample_metadata$condition[match(colnames(mat),sample_metadata$sample)]
  treatment<-unique(condition[condition!="Control"])
  direction<-ifelse(rowMeans(mat[,condition!="Control",drop=FALSE])>rowMeans(mat[,condition=="Control",drop=FALSE]),"Up","Down")
  anno_col<-data.frame(Condition=factor(condition,levels=c("Control",treatment)),row.names=colnames(mat))
  anno_row<-data.frame(Direction=direction,row.names=rownames(mat))
  mat<-pmin(pmax(mat,-2.5),2.5)
  code<-if(key=="6/B")"2DG" else "LDHA"
  # The deposited Cluster_rank preserves the original displayed gene order.
  pheatmap::pheatmap(mat,color=colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(101),
    breaks=seq(-2.5,2.5,length.out=102),cluster_rows=FALSE,cluster_cols=FALSE,gaps_col=5,
    annotation_row=anno_row,annotation_col=anno_col,
    annotation_colors=list(Condition=setNames(c("#59636E",if(code=="2DG")"#C83E3A" else "#8F4E93"),c("Control",treatment)),Direction=c(Down="#2166AC",Up="#B2182B")),
    annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,
    labels_col=c(paste0("C",1:5),paste0(code,1:5)),angle_col=45,border_color=NA,
    legend_breaks=c(-2,0,2),legend_labels=c("-2","0","2"),silent=TRUE)$gtable
}
original_imc_box <- function(d,x,y,facet=NULL,frequency=TRUE) {
  # Geometry, comparisons, and components from Met_expr.Rmd / Freq_clust.Rmd.
  d[[x]]<-ordered_groups(d[[x]])
  p<-ggplot(d,aes(.data[[x]],.data[[y]],fill=.data[[x]]))+geom_boxplot(outlier.shape=NA)+
    geom_point(position=position_dodge2(width=.5,padding=1.5))+theme_minimal()+
    theme(axis.text.x=element_text(size=14,face="bold",angle=45,hjust=1),axis.text.y=element_text(size=14,face="bold"),legend.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"),legend.title=element_text(size=14,face="bold"),plot.title=element_text(size=16,face="bold"),strip.text=element_text(size=14,face="bold"))+
    scale_fill_manual(values=colours_for(d[[x]]))+labs(x=NULL,y=y,fill=NULL)
  if(!is.null(facet))p<-p+facet_wrap(reformulate(facet),scales="free_y")
  tests<-list()
  groups<-if(is.null(facet))list(d)else split(d,d[[facet]])
  for(a in groups) {
    a<-a[is.finite(a[[y]])&!is.na(a[[x]]),];a[[x]]<-droplevels(a[[x]])
    if(nlevels(a[[x]])<2)next
    if(nlevels(a[[x]])==2)res<-rstatix::wilcox_test(a,reformulate(x,y)) else res<-rstatix::dunn_test(a,reformulate(x,y),p.adjust.method="none")
    if(!"p.adj"%in%names(res))res$p.adj<-res$p
    res$label<-ifelse(res$p.adj<.001,sprintf("p = %.2e",res$p.adj),ifelse(res$p.adj<.01,sprintf("p = %.5f",res$p.adj),ifelse(res$p.adj<.1,sprintf("p = %.4f",res$p.adj),sprintf("p = %.3f",res$p.adj))))
    span<-diff(range(a[[y]]));if(span==0)span<-.1
    unit_scale<-if(y=="frequency_percent")100 else 1
    offsets<-if(frequency)c(.1,.15,.1)*unit_scale else c(.25,.5,.25)
    res$y.position<-max(a[[y]])+offsets[seq_len(nrow(res))]
    if(!is.null(facet))res[[facet]]<-a[[facet]][1]
    tests[[length(tests)+1]]<-res
  }
  if(length(tests))p<-p+ggpubr::stat_pvalue_manual(bind_rows(tests),label="label",inherit.aes=FALSE,bracket.shorten=.05)
  if(length(unique(d[[x]]))>2) {
    kw<-lapply(groups,function(a){a<-a[is.finite(a[[y]]),];v<-kruskal.test(reformulate(x,y),data=a)$p.value
      z<-data.frame(label=paste0("Kruskal-Wallis, p = ",format.pval(v,digits=2)),y=max(a[[y]])+if(frequency){if(y=="frequency_percent")25 else .25}else .75)
      if(!is.null(facet))z[[facet]]<-a[[facet]][1];z})
    p<-p+geom_text(data=bind_rows(kw),aes(x=1.5,y=y,label=label),inherit.aes=FALSE,size=3.5)
  }
  p
}
original_network <- function(d) {
  require_packages_panel(c("ggraph","tidygraph","igraph"))
  # Literal manual layout in Interaction tests plots beautiful plots.Rmd.
  nodes <- data.frame(name=c("CD31_5","CD45_CD68_1","CD45_CD68_2","CD45_CD68_4","stromal_1","stromal_2","stromal_4","stromal_8","stromal_5","CD31_1","CD31_2"),
    x=c(450,-450,-200,-100,-514.7458,-220.6153,692.6986,349.6034,250,500,250),
    y=c(-350,150,350,200,-135.97703,-408.62595,25,145.17038,-100,-200,-300))
  plots <- lapply(split(d,d$Disease),function(edges) {
    graph<-igraph::graph_from_data_frame(edges[,c("from_label","to_label","InteractType","AbsMeanSig")],directed=FALSE,vertices=nodes)
    ggraph::ggraph(graph,layout="manual",x=x,y=y)+
      ggraph::geom_edge_fan(aes(colour=InteractType,width=AbsMeanSig),show.legend=TRUE)+
      ggraph::scale_edge_width_continuous(range=c(.5,3),limits=c(.15,1))+
      ggraph::scale_edge_color_manual(values=c(Attraction="#E64B35CC",Avoidance="#4DBBD577"))+
      ggraph::geom_node_point(aes(colour=name),size=6)+
      ggraph::geom_node_text(aes(label=display_labels(name)),repel=TRUE,size=2.7)+
      scale_colour_manual(values=colours_for(nodes$name,cluster_colours),guide="none")+
      ggtitle(edges$Disease[1])+theme_void()
  })
  wrap_plots(plots,ncol=2,guides="collect")
}
original_functional_bar <- function(d,key) {
  if(grepl("^S10/",key)) {
    x<-"dose_label";y<-"normalized_viability";facet<-NULL
    colors<-c("#769BCF","#A785D3","#DE9BDD","#B3D77F","#378EAE","#6E5399")
  } else if(key=="S11/A") {
    d<-d[d$Condition%in%c("Control","OX+2DG","OX+LDHA"),]
    x<-"Condition";y<-"normalized_value";facet<-"Channel";colors<-c("#909090","#67C1E4","#D77820")
  } else if(key=="S11/B") {
    x<-"condition";y<-"normalized_collagen_percent";facet<-NULL;colors<-c("#2020F4","#EF1919","#16C82D")
  } else {
    x<-"Condition";y<-"FoldChange";facet<-NULL;colors<-c("#909090","#67C1E4","#D77820")
  }
  lev<-unique(d[[x]]);d[[x]]<-factor(d[[x]],levels=lev)
  formula<-if(is.null(facet))reformulate(x,y) else reformulate(c(x,facet),y)
  tab<-aggregate(formula,d,function(v)c(mean=mean(v),sem=sd(v)/sqrt(length(v))))
  tab$mean<-tab[[y]][,"mean"];tab$sem<-tab[[y]][,"sem"]
  p<-ggplot(tab,aes(.data[[x]],mean,fill=.data[[x]]))+
    geom_col(width=.7,colour="grey35",linewidth=.4,alpha=.7)+
    geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),width=.25,colour="grey30",linewidth=.4)+
    geom_point(data=d,aes(.data[[x]],.data[[y]],colour=.data[[x]]),position=position_jitter(width=.12,height=0,seed=1),size=1.8,inherit.aes=FALSE)+
    scale_fill_manual(values=setNames(colors,lev),guide="none")+scale_colour_manual(values=setNames(colors,lev),guide="none")+
    labs(x=NULL,y=y)+theme_panel()+theme(axis.text.x=element_text(angle=45,hjust=1))
  if(!is.null(facet))p<-p+facet_wrap(reformulate(facet),nrow=1)
  if(grepl("^S10/",key)) {
    p<-p+annotate("rect",xmin=1.5,xmax=2.5,ymin=0,ymax=145,fill=NA,colour="red",linetype="dotted")+coord_cartesian(ylim=c(0,150))
  }
  if(key=="S11/A") {
    stats<-lapply(split(d,d$Channel),function(a){
      t<-tidyr::pivot_wider(a,names_from=Condition,values_from=normalized_value)
      data.frame(Channel=a$Channel[1],group1="Control",group2=c("OX+2DG","OX+LDHA"),p=vapply(c("OX+2DG","OX+LDHA"),function(v)t.test(t$Control,t[[v]],paired=TRUE)$p.value,numeric(1)),y.position=c(105,115))
    });stats<-do.call(rbind,stats);stats$label<-formatC(stats$p,digits=3,format="f")
    p<-p+geom_hline(yintercept=100,linetype="dashed",colour="grey50")+ggpubr::stat_pvalue_manual(stats,label="label",inherit.aes=FALSE)
  }
  if(key=="S11/B") {
    stats<-rstatix::dunn_test(d,normalized_collagen_percent~condition,p.adjust.method="bonferroni")
    stats<-stats[stats$group1=="Control",];stats$label<-formatC(stats$p.adj,digits=4,format="f");stats$y.position<-c(120,145)
    p<-p+ggpubr::stat_pvalue_manual(stats,label="label",inherit.aes=FALSE)
  }
  p
}
