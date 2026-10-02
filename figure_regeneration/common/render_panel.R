source("figure_regeneration/common/style.R")
source("figure_regeneration/common/panel_registry.R")
source("figure_regeneration/common/original_panel_adapters.R")

read_panel <- function(spec, data_dir) {
  path <- file.path(data_dir, spec$table)
  if (!file.exists(path)) stop("Source table not found: ", path)
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(d), spec$columns)) stop("Unexpected column names/order: ", spec$table)
  if (nrow(d) != spec$rows) stop("Unexpected row count: ", spec$table)
  d
}
heatmap_matrix <- function(d, row, col, value) {
  if (anyDuplicated(d[c(row, col)])) stop("Duplicate heatmap keys")
  rn <- unique(as.character(d[[row]])); cn <- unique(as.character(d[[col]]))
  mat <- matrix(NA_real_, length(rn), length(cn), dimnames = list(rn, cn))
  mat[cbind(match(d[[row]], rn), match(d[[col]], cn))] <- d[[value]]
  mat
}
draw_heatmap <- function(d, spec) {
  require_packages_panel(c("pheatmap"))
  mat <- heatmap_matrix(d, spec$row, spec$col, spec$value)
  palette <- if (isTRUE(spec$diverging)) colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100) else viridisLite::viridis(100)
  limits <- range(mat, na.rm = TRUE)
  if (limits[1] == limits[2]) limits <- limits + c(-.01, .01)
  anno <- NULL; anno_colors <- NULL
  if (spec$col %in% c("group", "metfiblabel", "metEClabel")) {
    anno <- data.frame(State = colnames(mat), row.names = colnames(mat))
    anno_colors <- list(State = colours_for(colnames(mat), c(cluster_colours, group_colours)))
    if ("n_cells" %in% names(d)) {
      n <- d[!duplicated(d[[spec$col]]), c(spec$col, "n_cells")]
      anno$`N cells` <- n$n_cells[match(colnames(mat), n[[spec$col]])]
    }
  }
  pheatmap::pheatmap(mat, color = palette, breaks = seq(limits[1], limits[2], length.out = 101),
    scale = "none", cluster_rows = isTRUE(spec$cluster_rows) && nrow(mat)>1 && all(is.finite(mat)),
    cluster_cols = isTRUE(spec$cluster_cols) && ncol(mat)>1 && all(is.finite(mat)),
    annotation_col = anno, annotation_colors = anno_colors, border_color = NA,
    labels_col = display_labels(colnames(mat)), labels_row = rownames(mat),
    main = if(is.null(spec$title)) NA_character_ else spec$title, fontsize = 9, angle_col = "45", silent = TRUE)$gtable
}
require_packages_panel <- function(x) {
  absent <- x[!vapply(x, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Missing plotting dependencies: ", paste(absent, collapse = ", "))
}
point_plot <- function(d, x, y, colour, title = NULL, continuous = FALSE, size = .2) {
  if (continuous) {
    p <- ggplot(d, aes(x = .data[[x]], y = .data[[y]], colour = .data[[colour]])) +
      geom_point(size = size, alpha = .65) + scale_colour_viridis_c()
  } else {
    d[[colour]] <- factor(d[[colour]], levels = sort(unique(d[[colour]])))
    p <- ggplot(d, aes(x = .data[[x]], y = .data[[y]], colour = .data[[colour]])) +
      geom_point(size = size, alpha = .65) +
      scale_colour_manual(values = colours_for(d[[colour]], c(cluster_colours, group_colours)), labels = display_labels)
  }
  p + labs(title = title, colour = NULL, x = x, y = y) + theme_panel()
}
box_panel <- function(d, x, y, facet = NULL, title = NULL, paired = NULL, violin = FALSE, test = "none", adjustment = "none") {
  require_packages_panel(c("ggpubr", "rstatix"))
  d[[x]] <- ordered_groups(d[[x]])
  p <- ggplot(d, aes(x = .data[[x]], y = .data[[y]], fill = .data[[x]]))
  if (violin) p <- p + geom_violin(trim = TRUE, alpha = .7)
  if (!is.null(paired)) p <- p + geom_line(aes(group = .data[[paired]]), colour = "grey65", linewidth = .35)
  p <- p + geom_boxplot(outlier.shape = NA, width = if (violin) .18 else .65) +
    geom_point(position = position_dodge2(width = .5, padding = 1.5), size = 1.3) +
    scale_fill_manual(values = colours_for(d[[x]]), labels = display_labels) +
    labs(title = title, x = NULL, y = y, fill = NULL) + theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), strip.text = element_text(face = "bold"), legend.position = "none")
  if (!is.null(facet)) p <- p + facet_wrap(stats::reformulate(facet), scales = "free_y")
  panels <- if (is.null(facet)) list(d) else split(d, d[[facet]])
  annotations <- list()
  for (dd in panels) {
    dd <- dd[is.finite(dd[[y]]) & !is.na(dd[[x]]), ]; dd[[x]] <- droplevels(dd[[x]])
    ng <- nlevels(dd[[x]])
    if (test == "none" || ng < 2L) next
    formula <- stats::reformulate(x, y)
    if(test=="t") {
      result <- rstatix::t_test(dd, formula, paired=!is.null(paired))
    } else if (test == "wilcox" || ng == 2L) {
      result <- tryCatch(rstatix::wilcox_test(dd, formula, paired = !is.null(paired), exact = FALSE), error = function(e) NULL)
    } else result <- tryCatch(rstatix::dunn_test(dd, formula, p.adjust.method = adjustment), error = function(e) NULL)
    if (is.null(result)) next
    if (!"p.adj" %in% names(result)) result$p.adj <- result$p
    result$label <- paste0("p = ", format.pval(result$p.adj, digits = 2, eps = .0001))
    span <- diff(range(dd[[y]], na.rm = TRUE)); if (!is.finite(span) || span == 0) span <- .1
    result$y.position <- max(dd[[y]], na.rm = TRUE) + span * (.15 + .17 * seq_len(nrow(result)))
    if (!is.null(facet)) result[[facet]] <- as.character(dd[[facet]][1])
    annotations[[length(annotations)+1L]] <- result
  }
  if (length(annotations)) p <- p + ggpubr::stat_pvalue_manual(do.call(rbind, annotations), label = "label", hide.ns = FALSE, tip.length = .01, inherit.aes = FALSE)
  p
}
radar_panel <- function(d, value, title) {
  require_packages_panel("ggradar")
  tab <- tidyr::pivot_wider(d[,c("group","cluster",value)],names_from=cluster,values_from=all_of(value))
  minimum <- min(as.matrix(tab[,-1]),na.rm=TRUE);maximum<-max(as.matrix(tab[,-1]),na.rm=TRUE)
  return(ggradar::ggradar(as.data.frame(tab),grid.min=minimum,grid.mid=(minimum+maximum)/2,grid.max=maximum,label.gridline.min=FALSE,label.gridline.mid=FALSE,label.gridline.max=FALSE,group.line.width=1,group.point.size=3,group.colours=unname(colours_for(tab$group)),fill=TRUE,fill.alpha=.5,background.circle.colour="white",gridline.mid.colour="grey",legend.position="bottom")+ggtitle(title))
  states <- sort(unique(d$cluster)); groups <- unique(d$group)
  r <- expand.grid(cluster = states, group = groups, stringsAsFactors = FALSE)
  r$value <- d[[value]][match(paste(r$cluster,r$group), paste(d$cluster,d$group))]
  minimum <- min(r$value, na.rm = TRUE); maximum <- max(r$value, na.rm = TRUE)
  # ggplot polar radius must be nonnegative; the axis labels retain the stored scale.
  r$radius <- r$value - minimum; r$cluster <- factor(r$cluster, levels = states); r$group <- ordered_groups(r$group)
  ggplot(r, aes(cluster, radius, group = group, colour = group, fill = group)) +
    geom_polygon(alpha = .15, linewidth = .7, na.rm = TRUE) + geom_point(size = 2) + coord_polar() +
    scale_y_continuous(labels = function(x) round(x + minimum, 2), limits = c(0, maximum-minimum)) +
    scale_colour_manual(values = colours_for(r$group)) + scale_fill_manual(values = colours_for(r$group)) +
    scale_x_discrete(labels = display_labels) + labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) + theme(axis.text.y = element_text(size = 7), legend.position = "bottom")
}
gene_dot <- function(d, spec) {
  gene <- spec$row; state <- spec$col
  d[[gene]] <- factor(d[[gene]], levels = unique(d[[gene]]))
  d[[state]] <- factor(d[[state]], levels = unique(d[[state]]))
  ggplot(d, aes(x = .data[[gene]], y = .data[[state]], size = .data[[spec$size]], colour = .data[[spec$value]])) +
    geom_point() + scale_size(range=c(0,6), name = "Percent Expressed") +
    scale_colour_gradient(low = "lightgrey", high = "darkred", name = "Average Expression") +
    labs(x = "Features", y = "Identity", title = spec$title) + cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
sankey_panel <- function(d) {
  require_packages_panel("ggsankey")
  # Aggregate counts are expanded to the cell-to-state pairs consumed by make_long.
  pairs<-d[rep(seq_len(nrow(d)),d$Freq),c("Xenium","IMC")]
  long<-ggsankey::make_long(pairs,Xenium,IMC)
  if(all(c("papillary_Fib","PI16_Fib")%in%pairs$Xenium)) {
    # The Figure 5F palette orders RNA annotations before the two IMC states.
    palette_order<-c("CCL19_Fib","COCH_Fib","COL8A1_Fib","COMP_Fib","CXCL12_Fib","NGFR_Fib","papillary_Fib","PI16_Fib","Met_hi_Fib","Other_Fib")
    return(ggplot(long,aes(x=x,next_x=next_x,node=node,next_node=next_node,fill=node,label=node))+
      ggsankey::geom_sankey(flow.alpha=.3,node.color="grey30")+
      ggsankey::geom_sankey_text(aes(hjust=ifelse(x=="Xenium",1.6,-.6)),size=3,color="black")+
      scale_fill_manual(values=setNames(scales::hue_pal()(10),palette_order))+
      scale_x_discrete(expand=expansion(add=c(.55,.55)))+theme_void()+theme(legend.position="none"))
  }
  ggplot(long,aes(x=x,next_x=next_x,node=node,next_node=next_node,fill=node,label=node))+
    ggsankey::geom_sankey(flow.alpha=.7,node.color="grey30")+
    ggsankey::geom_sankey_label(size=3,color="black")+theme_minimal()+labs(x=NULL,y=NULL)+
    theme(axis.text.x=element_text(size=10,face="bold"),legend.position="none")
}
network_panel <- function(d) {
  nodes <- unique(c(d$from_label, d$to_label)); angle <- seq(0,2*pi,length.out = length(nodes)+1)[-length(nodes)-1]
  xy <- data.frame(state = nodes, x = cos(angle), y = sin(angle))
  d$x <- xy$x[match(d$from_label, xy$state)]; d$y <- xy$y[match(d$from_label, xy$state)]
  d$xend <- xy$x[match(d$to_label, xy$state)]; d$yend <- xy$y[match(d$to_label, xy$state)]
  xy <- merge(xy, data.frame(Disease = unique(d$Disease)), by = NULL)
  ggplot() + geom_segment(data=d, aes(x,y,xend=xend,yend=yend,linewidth=AbsMeanSig, colour=MeanSig),alpha=.7) +
    geom_point(data=xy,aes(x,y,fill=state),shape=21,size=5) +
    geom_text(data=xy,aes(x*1.17,y*1.17,label=display_labels(state)),size=2.4) +
    facet_wrap(~Disease) + scale_colour_gradient2(low="#2166AC",mid="grey90",high="#B2182B") +
    scale_fill_manual(values=colours_for(nodes,cluster_colours)) + coord_equal(clip="off") + theme_void() + theme(legend.position="bottom")
}
ridge_panel <- function(d, spec) {
  # Density expression from Dist_to_cell_types_density_heatmap.Rmd.
  require_packages_panel(c("ggsci","ggtext"))
  d$cluster<-factor(d$cluster,levels=sort(unique(d$cluster)))
  ggplot(d,aes(distance_um,after_stat(density),colour=disease,fill=disease))+
    geom_density(alpha=.05,linewidth=1,bounds=c(0,2000))+
    scale_x_continuous(breaks=seq(0,160,40))+
    scale_colour_manual(values=setNames(rev(ggsci::pal_npg("nrc")(2)),c("Healthy","SSc")))+
    scale_fill_manual(values=setNames(rev(ggsci::pal_npg("nrc")(2)),c("Healthy","SSc")))+
    coord_cartesian(expand=FALSE,xlim=c(0,120),ylim=c(0,.04))+
    facet_wrap(~cluster,ncol=4,labeller=as_labeller(display_labels))+
    geom_vline(xintercept=40,linewidth=1)+theme_bw()+labs(x="Distance (µm)",y="Density",title=spec$title)+
    theme(strip.text=element_text(size=16,face="bold",lineheight=.8),legend.key.size=grid::unit(.8,"cm"),
      panel.spacing=grid::unit(.5,"cm"),plot.margin=margin(10,10,10,10),legend.title=element_blank(),
      legend.text=element_text(size=16),axis.text=element_text(size=10),axis.title=element_text(size=18))
}
niche_panel <- function(d,spec) {
  # Proportion heatmap from Dist_to_cell_types_density_heatmap.Rmd.
  require_packages_panel(c("ComplexHeatmap","circlize","RColorBrewer"))
  plots<-lapply(split(d,d$disease),function(a){
    a<-a[order(a$cluster),];mat<-heatmap_matrix(a,"cluster","niche","fraction")
    mat<-mat[,intersect(c("Distant","Within Niche"),colnames(mat)),drop=FALSE]
    h<-ComplexHeatmap::Heatmap(mat,name="matrix",col=circlize::colorRamp2(seq(0,1,length.out=100),colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdYlBu")))(100)),
      na_col="lightgrey",cluster_rows=FALSE,cluster_columns=FALSE,row_labels=display_labels(rownames(mat)),column_title=a$disease[1],
      rect_gp=grid::gpar(col="grey",lwd=.1),row_title_gp=grid::gpar(fontsize=10,fontface="bold"),column_title_gp=grid::gpar(fontsize=10,fontface="bold"))
    grid::grid.grabExpr(ComplexHeatmap::draw(h,heatmap_legend_side="right"),width=3,height=4)
  })
  wrap_plots(lapply(plots,function(p)wrap_elements(full=p)),ncol=2)
}
correlation_panel <- function(d, spec) {
  if ("current_mRSS" %in% names(d)) {
    rho <- cor(d$current_mRSS,d$frequency_percent,method="spearman")
    pval <- suppressWarnings(cor.test(d$current_mRSS,d$frequency_percent,method="spearman",exact=FALSE)$p.value)
    return(ggplot(d,aes(current_mRSS,frequency_percent)) + geom_smooth(method="lm",formula=y~x,colour="grey40",fill="grey80") +
      geom_point(size=2,colour="#C03830") + annotate("text",x=Inf,y=Inf,hjust=1.05,vjust=1.3,label=sprintf("rho = %.2f; p = %.3g",rho,pval),size=3) +
      labs(x="Total mRSS",y="Frequency (%)",title=spec$title) + theme_panel())
  }
  if ("clinical_var" %in% names(d)) {
    ggplot(d,aes(clinical_var,label,colour=estimate,size=-log10(pmax(p.value,1e-300)))) + geom_point() +
      scale_colour_gradient2(low="#2166AC",mid="white",high="#B2182B",limits=c(-1,1)) +
      scale_y_discrete(labels=display_labels) + labs(x=NULL,y=NULL,size="-log10(p)",colour="Spearman rho") + theme_panel() + theme(axis.text.x=element_text(angle=45,hjust=1))
  } else {
    d$interaction <- paste(display_labels(d$base_cluster),display_labels(d$neighbor_cluster),sep=" / ")
    ggplot(d,aes(rho,reorder(interaction,rho),size=-log10(pmax(p_spearman_raw,1e-300)),colour=rho)) +
      geom_vline(xintercept=0,colour="grey80") + geom_point() + scale_colour_gradient2(low="#2166AC",mid="white",high="#B2182B",limits=c(-1,1)) +
      labs(x="Spearman rho",y=NULL,title=spec$title,size="-log10(p)") + theme_panel()
  }
}
make_panel <- function(d, spec) {
  type <- spec$type
  if (type == "heatmap") return(draw_heatmap(d,spec))
  if (type == "heatmap_split") {
    pieces <- lapply(split(d,d$disease),function(a)draw_heatmap(a,spec))
    return(wrap_plots(lapply(pieces,function(x)wrap_elements(full=x))))
  }
  if (type == "umap") return(point_plot(d,"UMAP1","UMAP2",spec$colour,spec$title,size=spec$point_size))
  if (type == "mds") return(point_plot(d,"x","y","group",spec$title,size=3))
  if (type == "score_umap") {
    plots <- lapply(spec$values,function(v)point_plot(d,"UMAP1","UMAP2",v,display_labels(v),TRUE,size=.2))
    if (isTRUE(spec$validation_scale)) plots <- lapply(plots,function(p)p+scale_colour_gradientn(colours=c("#2C7BB6","#ABD9E9","#FFFFBF","#FDAE61","#D7191C"),limits=c(-1,1),oob=scales::squish))
    return(wrap_plots(plots,ncol=2))
  }
  if (type == "projection") return(wrap_plots(lapply(split(d,d$cohort),function(a)point_plot(a,"UMAP1","UMAP2","plotted_group",unique(a$cohort),size=.2)),ncol=2))
  if (type == "box") return(box_panel(d,spec$x,spec$y,spec$facet,spec$title,spec$paired,spec$violin,spec$test,spec$adjustment))
  if (type == "score_boxes") {
    cols <- intersect(c("Glycolysis","TCA_OXPHOS"),names(d))
    a <- do.call(rbind,lapply(cols,function(c)data.frame(donor=d$donor,group=d$disease,outcome=c,value=d[[c]])))
    return(box_panel(a,"group","value","outcome",spec$title,test="wilcox"))
  }
  if (type == "radar") return(wrap_plots(radar_panel(d,"median_frequency_fraction","Frequency"),radar_panel(d,"median_standardized_frequency","Standardized frequency"),ncol=2))
  if (type == "gene_dot") return(gene_dot(d,spec))
  if (type == "sankey") return(sankey_panel(d))
  if (type == "score_dot") {
    if ("score_name" %in% names(d)) {
      # The manuscript panel shows the ECM score; remaining scores stay in the source table.
      d <- d[d$score_name == spec$score_name, ]
      return(ggplot(d,aes(metfiblabel,lv3_anno,colour=mean_score))+geom_point(size=6)+scale_colour_gradient(low="grey95",high="#B30000")+labs(x=NULL,y=NULL,colour="Mean score")+theme_panel()+theme(axis.text.x=element_text(angle=45,hjust=1)))
    }
    a <- aggregate(d$score,list(pathway=d$pathway,group=d$group),mean); names(a)[3] <- "score"
    return(ggplot(a,aes(group,pathway,colour=score))+geom_point(size=6)+scale_colour_gradient(low="grey95",high="#B30000")+labs(x=NULL,y=NULL,colour="Mean AUCell score")+theme_panel())
  }
  if (type == "gsea") return(ggplot(d,aes(NES,reorder(Term,NES),size=Core_gene_count,colour=FDR))+geom_point()+scale_colour_gradient(low="#B2182B",high="#2166AC")+labs(y=NULL,title=spec$title)+theme_panel())
  if (type == "ridge") return(ridge_panel(d,spec))
  if (type == "niche") return(niche_panel(d,spec))
  if (type == "network") return(original_network(d))
  if (type == "interaction_dot") return(ggplot(d,aes(to_label,from_label,size=AbsMeanSig,colour=InteractType))+geom_point()+facet_wrap(~Disease)+scale_colour_manual(values=c(Attraction="#E64B35CC",Avoidance="#4DBBD577",NoInteraction="grey90"))+scale_size_continuous(range=c(0,7),limits=c(.05,1))+scale_x_discrete(labels=display_labels)+scale_y_discrete(labels=display_labels)+theme_minimal()+theme(axis.text.x=element_text(angle=90,hjust=1)))
  if (type == "difference") return(ggplot(d,aes(neigh_var,cluster_var,size=-log10(pmax(p,1e-300)),colour=difference))+geom_point()+scale_colour_gradient2(low="#2166AC",mid="white",high="#B2182B")+scale_y_discrete(labels=display_labels)+labs(x=NULL,y=NULL,colour="Difference (%)",size="-log10(p)")+theme_panel()+theme(axis.text.x=element_text(angle=90,hjust=1)))
  if (type == "correlation") return(correlation_panel(d,spec))
  if (type == "qc_signal") return(ggplot(d,aes(marker,cohort,fill=.data[[spec$value]]))+geom_tile()+scale_fill_viridis_c()+labs(x=NULL,y=NULL,fill=spec$value)+theme_panel()+theme(axis.text.x=element_text(angle=90,hjust=1)))
  if (type == "gating") {
    colour <- if("ASMA_gate"%in%names(d))"ASMA_gate" else "cell_type"
    plots <- lapply(spec$values,function(v)point_plot(d,v,"Histone3_asinh",colour,v,size=.12))
    if("UMAP1"%in%names(d))plots <- c(plots,list(point_plot(d,"UMAP1","UMAP2",colour,size=.12)))
    return(wrap_plots(plots,ncol=length(plots)))
  }
  if (type == "density_qc") {
    if(spec$qc_kind=="scatter") {
      plots <- list()
      for(x in c("mean_delaunay_neighbors_80um","mean_distance_5nn_um"))for(y in c("mean_glycolysis_score","mean_tca_oxphos_score")) {
        stat <- suppressWarnings(cor.test(d[[x]],d[[y]],method="spearman",exact=TRUE))
        plots[[length(plots)+1L]] <- ggplot(d,aes(.data[[x]],.data[[y]],colour=Progression_skin,size=n_met_hi_fib))+geom_point()+
          scale_colour_manual(values=colours_for(d$Progression_skin))+labs(x=display_labels(x),y=display_labels(y),subtitle=sprintf("rho = %.2f; p = %.3g",stat$estimate,stat$p.value))+theme_panel()
      }
    } else plots <- lapply(c("n_analyzed_cells","analyzed_cells_per_mm2"),function(v)box_panel(d,"Progression_skin",v,title=display_labels(v),test="kruskal"))
    return(wrap_plots(plots,ncol=2))
  }
  if (type == "model_frequency") return(ggplot(d,aes(scenario_label,frequency,fill=predicted_ec_state_like))+geom_col(position="dodge")+scale_y_continuous(labels=scales::percent)+labs(x=NULL,y="Predicted state frequency",fill=NULL)+theme_panel())
  if (type == "model_spatial") return(ggplot(d,aes(x_centroid,y_centroid,colour=transition))+geom_point(size=1.3)+facet_wrap(~donor,scales="free")+scale_y_reverse()+labs(x="X (um)",y="Y (um)",colour="State transition")+theme_panel())
  if (type == "seahorse") {
    endpoint <- d[d$plot=="Respiration endpoint",]; time <- d[d$plot!="Respiration endpoint",]
    require_packages_panel(c("lme4","lmerTest","pbkrtest"))
    models <- lapply(split(endpoint,endpoint$measure),function(dd) {
      dd$group <- factor(dd$group,levels=c("Healthy","Progressive SSc"));dd$EC_line<-factor(dd$EC_line)
      fit <- lmerTest::lmer(value~group+EC_line+(1|donor),data=dd)
      result <- anova(fit,ddf="Kenward-Roger")
      data.frame(measure=dd$measure[1],p=result["group","Pr(>F)"])
    })
    stats <- do.call(rbind,models)
    endpoint$group<-factor(endpoint$group,levels=c("Healthy","Progressive SSc"))
    means<-endpoint%>%group_by(group,measure)%>%summarise(mean=mean(value),sem=sd(value)/sqrt(n()),.groups="drop")
    p1 <- ggplot(means,aes(group,mean,fill=group))+geom_col(width=.65,alpha=.65)+geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),width=.2)+
      geom_point(data=endpoint,aes(group,value),inherit.aes=FALSE,position=position_jitter(width=.12,height=0,seed=1),size=1.3)+
      facet_wrap(~measure,scales="free_y")+scale_fill_manual(values=colours_for(endpoint$group))+theme_panel()+labs(x=NULL,y="OCR")+
      geom_text(data=stats,aes(x=1.5,y=Inf,label=sprintf("p = %.3g",p)),inherit.aes=FALSE,vjust=1.5,size=3)
    a <- aggregate(time$value,list(group=time$group,time=time$time_minutes),function(x)c(mean=mean(x),se=sd(x)/sqrt(length(x))))
    a$mean <- a$x[,"mean"]; a$se <- a$x[,"se"]
    p2 <- ggplot(a,aes(time,mean,colour=group))+geom_line()+geom_point()+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=.4)+scale_colour_manual(values=colours_for(a$group))+labs(x="Time (min)",y="OCR")+theme_panel()
    return(p2/p1)
  }
  stop("Unknown panel renderer: ",type)
}
freeze_layout_units <- function(grob) {
  # Resolve text-dependent dimensions while the measuring device is active.
  # Flexible panel dimensions retain their relative units.
  resolve <- function(units, direction) {
    if(!length(units))return(units)
    converted<-lapply(seq_along(units),function(i){
      u<-units[i]
      if(any(grepl("null|npc",as.character(u))))u else if(direction=="width")grid::convertWidth(u,"inches") else grid::convertHeight(u,"inches")
    })
    do.call(grid::unit.c,converted)
  }
  if(inherits(grob,"gtable")) {
    grob$grobs<-lapply(grob$grobs,freeze_layout_units)
    grob$widths<-resolve(grob$widths,"width")
    grob$heights<-resolve(grob$heights,"height")
  } else if(inherits(grob,"gTree") && length(grob$children)) {
    grob$children<-do.call(grid::gList,lapply(grob$children,freeze_layout_units))
  }
  grob
}
render_panel <- function(figure, panel, data_dir, output_dir, formats=c("pdf","png")) {
  key <- paste(figure,panel,sep="/"); spec <- panel_registry[[key]]
  if (is.null(spec)) stop("Unknown figure panel: ",key)
  d <- read_panel(spec,data_dir)
  set.seed(1)
  plot <- make_original_panel(d,spec,key,data_dir)
  if(is.null(plot))plot <- make_panel(d,spec)
  if(inherits(plot,"patchwork"))plot<-plot & theme(text=element_text(family="Arial")) else if(inherits(plot,"ggplot"))plot<-plot+theme(text=element_text(family="Arial"))
  directory <- file.path(output_dir,paste0("Figure_",figure));dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  prefix <- file.path(directory,paste0("Panel_",panel))
  # Build one layout on the export device, then reuse it for every format.
  # Separate ggplot builds can inherit different text metrics from prior devices.
  layout_device <- tempfile(fileext=".png")
  grDevices::png(layout_device,width=spec$width,height=spec$height,units="in",res=200,type="cairo",pointsize=12,family="Arial")
  plot_grob <- tryCatch({
    grid::grid.newpage()
    set.seed(1)
    layout<-if(inherits(plot,c("grob","gtable")))plot else if(inherits(plot,"patchwork"))patchwork::patchworkGrob(plot) else ggplot2::ggplotGrob(plot)
    freeze_layout_units(layout)
  },finally={grDevices::dev.off();unlink(layout_device)})
  for (format in formats) {
    filename <- paste0(prefix,".",format)
    if (format=="pdf") grDevices::cairo_pdf(filename,width=spec$width,height=spec$height,pointsize=12,family="Arial") else if(format=="png") grDevices::png(filename,width=spec$width,height=spec$height,units="in",res=200,type="cairo",pointsize=12,family="Arial") else stop("Supported formats: pdf, png")
    tryCatch({grid::grid.newpage();set.seed(1);grid::grid.draw(plot_grob)},finally=grDevices::dev.off())
    if (file.info(filename)$size < 500) stop("Empty plot export: ",filename)
  }
  invisible(list(figure=figure,panel=panel,rows=nrow(d),table=spec$table,type=spec$type))
}
panel_main <- function(figure,panel,args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)<2)stop("Usage: Rscript figure_regeneration/Figure_<id>/Panel_<id>.R DATA_DIR OUTPUT_DIR")
  result<-render_panel(figure,panel,args[1],args[2])
  writeLines(capture.output(sessionInfo()),file.path(args[2],paste0("Figure_",figure),paste0("Panel_",panel,"_sessionInfo.txt")))
  invisible(result)
}
