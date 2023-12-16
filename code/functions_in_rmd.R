assign_subgroup = function(group1, group2){
  group1 <<- group1
  group2 <<- group2
  key <<- paste0(paste(group1, collapse = "_"), "&", paste(group2, collapse = "_"))
  subgroupsce <<- inputData[, inputData$hippo_cluster %in% c(group1,group2)]
#  subgroupsce_lognorm <<- inputData_lognorm[, colData(inputData)$hippo_cluster %in% c(group1,group2)]
  subgroupsce_integrated <<- inputData_integrated[, inputData$hippo_cluster %in% c(group1,group2)]
  cellgroup <<- ifelse(subgroupsce$hippo_cluster%in% group1, 1, 2)
  genemean1 = rowMeans(counts(subgroupsce)[,cellgroup==1])
  genemean2 = rowMeans(counts(subgroupsce)[,cellgroup==2])
  subgroupsce@metadata$log2mean <<- log2(genemean1*genemean2)/2
  subgroupsce@metadata$log2meandiff <<- log2(abs(genemean1-genemean2))
}

assign_celltype = function(celltype){
  group1 <<- "ctrl"
  group2 <<- "stim"
  key <<- celltype
  subgroupsce <<- sce[, sce$cluster_id%in% celltype]
  cellgroup <<- ifelse(subgroupsce$group_id%in% group1, 1, 2)
  genemean1 = rowMeans(counts(subgroupsce)[,cellgroup==1])
  genemean2 = rowMeans(counts(subgroupsce)[,cellgroup==2])
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
    p[[i]] = ggplot(na.omit(data.frame(x = x, y = y, DEGs = dflist[[i]][[key]]$hits)), aes(x = x, y = y, colour = DEGs)) +
      geom_point(alpha = 0.5, size = 0.5) +
      geom_hline(yintercept=log2(1.5),linetype=2) +
      geom_hline(yintercept=-log2(1.5),linetype=2) +
      ggtitle(titlelist[i]) + theme_minimal() + xlim(xlimits) + ylim(c(-4,4)) + 
      scale_color_manual(values = c("gray", "blue")) +
      xlab(xLabel) + ylab("Log2 Fold Change") + theme(legend.position = "bottom")
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
    p[[i]] = ggplot(na.omit(data.frame(x = x, y = y, DEGs = dflist[[i]][[key]]$hits)), aes(x = x, y = y, colour = DEGs)) +
      geom_point(alpha = 0.5, size = 0.5) +
      ggtitle(titlelist[i]) + theme_minimal() + xlim(c(-6,6)) + ylim(c(-15,6))+ 
      scale_color_manual(values = c("gray", "blue")) + 
      xlab("Log2 mean") + ylab("Log2 mean difference") + theme(legend.position = "bottom")
  }
  return(ggarrange(plotlist = p,common.legend = TRUE, legend = "right", ncol = 3, nrow = 4))
}

variation_analysis = function(){
    
  subset_ind = inputData$hippo_cluster%in%c(group1,group2)
  variation_raw = data.frame()
  variation_vst = data.frame()
  variation_cpm = data.frame()
  variation_integrated = data.frame()
  donor = subgroupsce$donor
  
  for(gene in intersect(rownames(subgroupsce_integrated), rownames(vstcounts))){
    
    y = log2(assays(subgroupsce)$counts[gene,] + 1)
    y_vst = vstcounts[gene,subset_ind]
    y_cpm = log2(assays(subgroupsce)$cpm[gene,] + 1)
    y_int = pmax(assays(subgroupsce_integrated)$counts[gene,],0)
    
    s_raw = summary(lm(y ~ donor + cellgroup))
    s_vst = summary(lm(y_vst ~ donor + cellgroup))
    s_cpm = summary(lm(y_cpm ~ donor + cellgroup))
    s_int = summary(lm(y_int ~ donor + cellgroup))
    
    variation_raw = rbind(variation_raw, 
                          data.frame(donor = var(model.matrix(s_raw)[,2:5]%*% coef(s_raw)[2:5]), 
                                     celltype = var(model.matrix(s_raw)[,6] * coef(s_raw)[6]), 
                                     res = var(s_raw$residuals)))
    
    variation_vst = rbind(variation_vst, 
                              data.frame(donor = var(model.matrix(s_vst)[,2:5]%*% coef(s_vst)[2:5]), 
                                         celltype = var(model.matrix(s_vst)[,6] * coef(s_vst)[6]), 
                                         res = var(s_vst$residuals)))
    
    variation_cpm = rbind(variation_cpm, 
                          data.frame(donor = var(model.matrix(s_cpm)[,2:5]%*% coef(s_cpm)[2:5]), 
                                     celltype = var(model.matrix(s_cpm)[,6] * coef(s_cpm)[6]), 
                                     res = var(s_cpm$residuals)))
    
    variation_integrated = rbind(variation_integrated, 
                                 data.frame(donor = var(model.matrix(s_int)[,2:5]%*% coef(s_int)[2:5]), 
                                            celltype = var(model.matrix(s_int)[,6] * coef(s_int)[6]), 
                                            res = var(s_int$residuals)))
  }
  variation_raw = variation_raw/rowSums(variation_raw)
  variation_vst = variation_vst/rowSums(variation_vst)
  variation_cpm = variation_cpm/rowSums(variation_cpm)
  variation_integrated = variation_integrated/rowSums(variation_integrated)
  variation_raw$data = "Raw UMI"
  variation_vst$data = "VST"
  variation_cpm$data = "CPM"
  variation_integrated$data = "Integrated"
  
  top500 = order(variation_raw$res)[1:500]
  v_list = list(variation_raw[top500,], 
                variation_vst[top500,],
                variation_cpm[top500,],
                variation_integrated[top500,])

  for (i in 1:length(v_list)){
    v_list[[i]]$quantile_donor = as.factor(as.numeric(cut_number(v_list[[i]]$donor,4))*100/4)
    v_list[[i]]$quantile_celltype = as.factor(as.numeric(cut_number(v_list[[i]]$celltype,4))*100/4)
    v_list[[i]]$quantile_res = as.factor(as.numeric(cut_number(v_list[[i]]$res,4))*100/4)
  }
  combined_df <<- bind_rows(v_list)
  melt_combined <<- reshape2::melt(combined_df, id.vars = c("data", "quantile_res"), measure.vars = c("donor","celltype"), value.name = "variation", variable.name = "source")
  wide_combined <<- reshape(data=cbind(index = c(rep(1:500,4), rep(501:1000,4)), melt_combined[,-2]), idvar = c("index", "source"),
                          v.names = "variation",
                          timevar = "data",
                          direction="wide")
  colnames(wide_combined)[3:6] <<- c("Raw UMI", "VST", "CPM", "Integrated")
}

value_comparison = function(value = c("pval", "log2FC")){
  p = list()
  for (i in 2:length(dflist)){
    if (nrow(dflist[[i]][[key]])==0){next}
    if (value == "pval"){
      xlimit = c(0,1)
      ylimit = c(0,1)
      if (is.null(dflist[[i]][[key]]$pval)){
        x = dflist[[i]][[key]]$p_val
      }else
        x = dflist[[i]][[key]]$pval
    }
    if (value == "log2FC"){
      xlimit = c(-3,3)
      ylimit = c(-3,3)
      if (is.null(dflist[[i]][[key]]$log2FC)){
        x = dflist[[i]][[key]]$avg_log2FC
      }else
        x = dflist[[i]][[key]]$log2FC
    }
    y = dflist[[1]][[key]][match(dflist[[i]][[key]]$genes, dflist[[1]][[key]]$genes),value]
    p[[i-1]] = ggplot(na.omit(data.frame(x = x, y = y)), aes(x = x, y = y)) +
      geom_point(alpha = 0.5, size = 0.5) + theme_minimal() + xlim(xlimit) + ylim(ylimit)+ 
      xlab(titlelist[i]) + ylab(titlelist[1])
  }
  figure = ggarrange(plotlist = p,nrow = 2, ncol = 4)
  
  # Annotate the figure by adding a common labels
  return(annotate_figure(figure,
                  top = paste(value, "comparison",key)))
}
