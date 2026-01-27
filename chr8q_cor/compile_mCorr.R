compile_mCorr <- function(DEGs, p = 0.05, FDR.features = 0.05, 
                         n.dot.features = 10, 
                         estCol = "Spearman.est", 
                         pCol = "Spearman.p",
                         qCol = "Spearman.q") {
  ## create heatmap data frames
  # extract feature names for each omics type
  types <- names(DEGs)
  for (i in 1:length(types)) {
    DEGs[[types[i]]]$feature_name <- colnames(DEGs[[types[i]]])[[1]]
    colnames(DEGs[[types[i]]])[[1]] <- "feature"
  }
  
  # collapse DEG results across omics types
  DEG.df <- as.data.frame(data.table::rbindlist(DEGs, use.names = TRUE, idcol = "type"))
  DEG.df$minusLogP <- -log(DEG.df[,pCol], base = 10)
  DEG.df$minusLogFDR <- -log(DEG.df[,qCol], base = 10)
  DEG.df$sig <- FALSE
  DEG.df[!is.na(DEG.df[,pCol]) & !is.na(DEG.df[,qCol]) & 
           DEG.df[,pCol] < p & DEG.df[,qCol] < FDR.features, ]$sig <- TRUE
  
  # summarize results for each feature
  DEG.df$rankVal <- DEG.df[,estCol]
  DEG.df$P.Value <- DEG.df[,pCol]
  mean.DEG.df <- plyr::ddply(na.omit(DEG.df), .(feature), summarize,
                             mean_rankVal = mean(rankVal),
                             sd_rankVal = sd(rankVal),
                             Fisher_p = metap::sumlog(na.omit(P.Value))$p,
                             types = paste0(type, collapse = ", "),
                             N_types = length(unique(type)),
                             N_sig = length(sig[sig]),
                             sig_types = paste0(type[sig], collapse = ", "))
  if (length(unique(mean.DEG.df$Fisher_p)) > 1) {
    mean.DEG.df$adj_Fisher_p <- 
      qvalue::qvalue(mean.DEG.df$Fisher_p, pi0=1)$qvalues
  } else {
    mean.DEG.df$adj_Fisher_p <- NA
  }
  
  # order results by rankVal for dot plot
  mean.DEG.df <- dplyr::arrange(mean.DEG.df, desc(mean_rankVal))
  
  # reduce plot data down to top results
  sig.DEG.df <- mean.DEG.df[mean.DEG.df$N_sig > 0, ]
  
  if (nrow(sig.DEG.df) >= n.dot.features) {
    top.DEG.df <- 
      sig.DEG.df %>% dplyr::slice_max(abs(mean_rankVal), n = n.dot.features)
  } else {
    top.DEG.df <- 
      mean.DEG.df %>% dplyr::slice_max(abs(mean_rankVal), n = n.dot.features)
  }
  dot.df <- na.omit(DEG.df[DEG.df$feature %in% top.DEG.df$feature, ])
  dot.df$qVal <- dot.df[,qCol]
  if (any(dot.df$qVal == 0)) {
    dot.df[dot.df$qVal == 0,]$qVal <- 0.0001
  }
  
  ## create venn diagram
  # compile significant results for each type in list
  venn.list <- list()
  for (i in 1:length(types)) {
    venn.list[[types[i]]] <- DEG.df[DEG.df$type == types[i] &
                                      DEG.df$sig, ]$feature
  }
  
  # generate venn diagram
  venn.plot <- ggvenn::ggvenn(venn.list) # only displays first 4 types
  
  ## create dot plot
  # set theme
  bg.theme <- ggplot2::theme(
    legend.background = element_rect(), legend.position = "top",
    legend.text = element_text(size = 14),
    legend.key = element_blank(),
    legend.title = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16, colour = "black"),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16, colour = "black"),
    plot.title = element_text(
      lineheight = .8, face = "bold", size = 36
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black")
  )
  
  # generate dot plot
  dot.plot <- ggplot2::ggplot(
    na.omit(dot.df),
    ggplot2::aes(
      x = type, y = feature, color = rankVal,
      size = -log10(qVal)
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = mean.DEG.df[
      mean.DEG.df$feature %in% top.DEG.df$feature, ]$feature) +
    viridis::scale_color_viridis() +
    bg.theme +
    ggplot2::labs(
      x = "Input Data",
      y = "Feature Set",
      color = estCol, size = "-log(FDR)"
    ) + geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21)
  
  ## run correlations
  # extract Log2FC, -logFDR values across omics types
  Log2FC.df <- reshape2::dcast(na.omit(DEG.df), feature ~ type,
                               value.var = estCol, fill = NA
  )
  minusLogFDR.df <- reshape2::dcast(na.omit(DEG.df), feature ~ type,
                                    value.var = "minusLogFDR", fill = NA
  )
  
  # convert from data frame to numeric matrix
  rownames(Log2FC.df) <- Log2FC.df$feature
  Log2FC.mat <- as.matrix(Log2FC.df[, 2:ncol(Log2FC.df)])
  
  # create correlation matrix
  corr.mat <- stats::cor(Log2FC.mat)
  
  # plot correlation matrix
  corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
  
  ## compile outputs
  outputs <- list(
    results = DEG.df,
    mean.results = mean.DEG.df,
    est.df = Log2FC.df,
    minusLogFDR.df = minusLogFDR.df,
    venn.diagram = venn.plot,
    dot.plot = dot.plot,
    corr = corr.mat,
    corr.matrix = corr.mat.plot
  )
  return(outputs)
}
