assign_subgroup = function(group1, group2){
  group1 <<- group1
  group2 <<- group2
  key <<- paste0(paste(group1, collapse = "_"), "&", paste(group2, collapse = "_"))
  subgroupsce <<- inputData[, colData(inputData)$hippo_cluster%in% c(group1,group2)]
  subgroupsce_lognorm <<- inputData_lognorm[, colData(inputData)$hippo_cluster %in% c(group1,group2)]
  subgroupsce_integrated <<- inputData_integrated[, colData(inputData)$hippo_cluster %in% c(group1,group2)]
  cellgroup <<- ifelse(colData(subgroupsce)$hippo_cluster%in% group1, 1, 2)
  genemean1 = rowMeans(subgroupsce@assays@data$counts[,cellgroup==1])
  genemean2 = rowMeans(subgroupsce@assays@data$counts[,cellgroup==2])
  subgroupsce@metadata$log2mean <<- log2(genemean1*genemean2)/2
  subgroupsce@metadata$log2meandiff <<- log2(abs(genemean1-genemean2))
}

FC_mean_plot = function(xvar, xlimits, xLabel){
  p = list()
  for (i in 1:length(dflist)){
    if(is.null(dflist[[i]][[key]]$genes)){
      x = xvar[match(rownames(dflist[[i]][[key]]), names(xvar))]
    }else{
      x = xvar[match(dflist[[i]][[key]]$genes,names(xvar))]
    }
    if(is.null(dflist[[i]][[key]]$log2FC)){
      y = dflist[[i]][[key]]$avg_log2FC
    }else{
      y = dflist[[i]][[key]]$log2FC
    }
    p[[i]] = ggplot(na.omit(data.frame(x = x, y = y, hits = dflist[[i]][[key]]$hits)), aes(x = x, y = y, colour = hits)) +
      geom_point(alpha = 0.5, size = 0.5) +
      geom_hline(yintercept=log2(1.5),linetype=2) +
      geom_hline(yintercept=-log2(1.5),linetype=2) +
      ggtitle(titlelist[i]) + theme_minimal() + xlim(xlimits) + ylim(c(-4,4)) + 
      scale_color_manual(values = c("gray", "blue")) +
      xlab(xLabel) + ylab("log2 Fold Change") + theme(legend.position = "bottom")
  }
  return(ggarrange(plotlist = p,common.legend = TRUE, legend = "right", ncol = 4, nrow = 2))
}

mean_meandiff_plot = function(){
  p = list()
  for (i in 1:length(dflist)){
    if(is.null(dflist[[i]][[key]]$genes)){
      x = subgroupsce@metadata$log2mean[match(rownames(dflist[[i]][[key]]), rownames(subgroupsce))]
      y = subgroupsce@metadata$log2meandiff[match(rownames(dflist[[i]][[key]]), rownames(subgroupsce))]
    }else{
      x = subgroupsce@metadata$log2mean[match(dflist[[i]][[key]]$genes,rownames(subgroupsce))]
      y = subgroupsce@metadata$log2meandiff[match(dflist[[i]][[key]]$genes,rownames(subgroupsce))]
    }
    p[[i]] = ggplot(na.omit(data.frame(x = x, y = y, hits = dflist[[i]][[key]]$hits)), aes(x = x, y = y, colour = hits)) +
      geom_point(alpha = 0.5, size = 0.5) +
      ggtitle(titlelist[i]) + theme_minimal() + xlim(c(-6,6)) + ylim(c(-15,6))+ 
      scale_color_manual(values = c("gray", "blue")) + 
      xlab("log2 mean") + ylab("log2 mean difference") + theme(legend.position = "bottom")
  }
  return(ggarrange(plotlist = p,common.legend = TRUE, legend = "right", ncol = 3, nrow = 4))
}


create_heatmap = function(){
  matscglmm = subgroupsce@assays@data$counts[gene_ind, ]
  #row annotation
  annotation_REvariation = data.frame(REvariation = pois_glmm_df[[key]]$REvariation[match(gene_ind,pois_glmm_df[[key]]$genes)])
  rownames(annotation_REvariation) = rownames(matscglmm)
  #col annotation
  donorlist = unique(subgroupsce@colData$donor)
  annotation_df = data.frame(donor = subgroupsce@colData$donor, group = as.factor(colData(subgroupsce)$hippo_cluster))
  rownames(annotation_df) = colnames(matscglmm)
  annotation_df = annotation_df[with(annotation_df, order(donor, group)), ] 
  
  # cell level
  pheatmap(matscglmm[,rownames(annotation_df)], annotation_col = annotation_df, annotation_row = annotation_REvariation, cluster_rows=F, cluster_cols=F, show_colnames = F, color=colorRampPalette(c("navy", "white", "red"))(10), main = fig.name, breaks = seq(0,5, length.out = 10))
  # donor level
  raw_count = melt(t(matscglmm), value.name = "count")
  colnames(raw_count) = c("cell", "gene", "count")
  raw_count$group = colData(subgroupsce)$hippo_cluster
  raw_count$donor = subgroupsce@colData$donor
  pseudobulk_mean = aggregate(count ~ donor + group + gene , data = raw_count, FUN = mean)
  pseudobulk_mean$cell_donor = paste0(pseudobulk_mean$group, "_", pseudobulk_mean$donor)
  pseudobulk_mean = reshape(pseudobulk_mean[,c("cell_donor","gene", "count")], idvar = "gene", timevar = c("cell_donor"), direction = "wide")
  rownames(pseudobulk_mean) = pseudobulk_mean[,1]
  pseudobulk_mean[,1] = NULL
  
  annotation_df = data.frame(group = c(rep(as.character(group1), 5), rep(as.character(group2), 5)), donor = rep(unique(raw_count$donor),2))
  rownames(annotation_df) = colnames(pseudobulk_mean)
  annotation_df = annotation_df[with(annotation_df, order(donor, group)), ] 
  
  pheatmap(pseudobulk_mean[, rownames(annotation_df)], annotation_col = annotation_df, annotation_row = annotation_REvariation, cluster_rows=F, cluster_cols=F, show_colnames = F, color=colorRampPalette(c("navy", "white", "red"))(10), main = fig.name, breaks = seq(0,5, length.out = 10))
}
