assign_subgroup = function(group1, group2){
  group1 <<- group1
  group2 <<- group2
  cur_key <<- paste0(paste(group1, collapse = "_"), "&", paste(group2, collapse = "_"))
  subgroupsce <<- sce_fallopian_tubes[, sce_fallopian_tubes$hippo_cluster %in% c(group1,group2)]
  subgroup_Seurat <<- Seurat_fallopian_tubes[, sce_fallopian_tubes$hippo_cluster %in% c(group1,group2)]
  cellgroup <<- ifelse(subgroupsce$hippo_cluster%in% group1, 1, 2)
  genemean1 = rowMeans(counts(subgroupsce)[,cellgroup==1])
  genemean2 = rowMeans(counts(subgroupsce)[,cellgroup==2])
  subgroupsce@metadata$log2mean <<- log2(genemean1*genemean2)/2
  subgroupsce@metadata$log2meandiff <<- log2(abs(genemean1-genemean2))
}

assign_celltype = function(celltype){
  group1 <<- "ctrl"
  group2 <<- "stim"
  cur_key <<- celltype
  subgroupsce <<- sce[, sce$cluster_id%in% celltype]
  subgroup_Seurat <<- Seurat_sce[, Seurat_sce$cluster_id %in% celltype]
  cellgroup <<- ifelse(subgroupsce$group_id%in% group1, 1, 2)
  genemean1 = rowMeans(counts(subgroupsce)[,cellgroup==1])
  genemean2 = rowMeans(counts(subgroupsce)[,cellgroup==2])
  subgroupsce@metadata$log2mean <<- log2(genemean1*genemean2)/2
  subgroupsce@metadata$log2meandiff <<- log2(abs(genemean1-genemean2))
}

FC_mean_plot = function(xvar, xlimits, xLabel){
  p = list()
  for (i in 1:length(dflist)){
    if(is.null(dflist[[i]][[cur_key]]$genes)){
      x = xvar[match(rownames(dflist[[i]][[cur_key]]), names(xvar))]
    }else{
      x = xvar[match(dflist[[i]][[cur_key]]$genes,names(xvar))]
    }
    if(is.null(dflist[[i]][[cur_key]]$log2FC)){
      y = dflist[[i]][[cur_key]]$avg_log2FC
    }else{
      y = dflist[[i]][[cur_key]]$log2FC
    }
    p[[i]] = ggplot(na.omit(data.frame(x = x, y = y, DEGs = dflist[[i]][[cur_key]]$hits)), aes(x = x, y = y, colour = DEGs)) +
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
    if(is.null(dflist[[i]][[cur_key]]$genes)){
      x = subgroupsce@metadata$log2mean[match(rownames(dflist[[i]][[cur_key]]), rownames(subgroupsce))]
      y = subgroupsce@metadata$log2meandiff[match(rownames(dflist[[i]][[cur_key]]), rownames(subgroupsce))]
    }else{
      x = subgroupsce@metadata$log2mean[match(dflist[[i]][[cur_key]]$genes,rownames(subgroupsce))]
      y = subgroupsce@metadata$log2meandiff[match(dflist[[i]][[cur_key]]$genes,rownames(subgroupsce))]
    }
    p[[i]] = ggplot(na.omit(data.frame(x = x, y = y, DEGs = dflist[[i]][[cur_key]]$hits)), aes(x = x, y = y, colour = DEGs)) +
      geom_point(alpha = 0.5, size = 0.5) +
      ggtitle(titlelist[i]) + theme_minimal() + xlim(c(-6,6)) + ylim(c(-15,6))+ 
      scale_color_manual(values = c("gray", "blue")) + 
      xlab("Log2 mean") + ylab("Log2 mean difference") + theme(legend.position = "bottom")
  }
  return(ggarrange(plotlist = p,common.legend = TRUE, legend = "right", ncol = 3, nrow = 4))
}

variation_analysis = function(){
    
  subset_ind = sce_fallopian_tubes$hippo_cluster%in%c(group1,group2)
  variation_raw = data.frame()
  variation_vst = data.frame()
  variation_cpm = data.frame()
  variation_integrated = data.frame()
  donor = subgroupsce$donor

  for(gene in intersect(rownames(subgroup_Seurat@assays$integrated$data), 
                        rownames(vstcounts))){
    y = log2(assays(subgroupsce)$counts[gene,] + 1)
    y_vst = vstcounts[gene,subset_ind]
    y_cpm = log2(assays(subgroupsce)$cpm[gene,] + 1)
    y_int = as.numeric(subgroup_Seurat@assays$integrated$data[gene,])
    meta_data = as.data.frame(colData(subgroupsce))
    
    a_raw = anova(lm(y ~ donor + cellgroup))
    a_vst = anova(lm(y_vst ~ donor + cellgroup))
    a_cpm = anova(lm(y_cpm ~ donor + cellgroup))
    a_int = anova(lm(y_int ~ donor + cellgroup))
    
    variation_raw = rbind(variation_raw, 
                          data.frame(gene = gene, 
                                     donor = a_raw$`Sum Sq`[1]/sum(a_raw$`Sum Sq`),
                                     celltype = a_raw$`Sum Sq`[2]/sum(a_raw$`Sum Sq`),
                                     res = a_raw$`Sum Sq`[3]/sum(a_raw$`Sum Sq`)))
    variation_vst = rbind(variation_vst, 
                          data.frame(gene = gene, 
                                     donor = a_vst$`Sum Sq`[1]/sum(a_vst$`Sum Sq`),
                                     celltype = a_vst$`Sum Sq`[2]/sum(a_vst$`Sum Sq`),
                                     res = a_vst$`Sum Sq`[3]/sum(a_vst$`Sum Sq`)))
    variation_cpm = rbind(variation_cpm, 
                          data.frame(gene = gene, 
                                     donor = a_cpm$`Sum Sq`[1]/sum(a_cpm$`Sum Sq`),
                                     celltype = a_cpm$`Sum Sq`[2]/sum(a_cpm$`Sum Sq`),
                                     res = a_cpm$`Sum Sq`[3]/sum(a_cpm$`Sum Sq`)))
    variation_integrated = rbind(variation_integrated, 
                                 data.frame(gene = gene, 
                                            donor = a_int$`Sum Sq`[1]/sum(a_int$`Sum Sq`),
                                            celltype = a_int$`Sum Sq`[2]/sum(a_int$`Sum Sq`),
                                            res = a_int$`Sum Sq`[3]/sum(a_int$`Sum Sq`)))
  }

  top500 = order(variation_raw$res)[1:500]
  v_list = list(variation_raw[top500,], 
                variation_vst[top500,],
                variation_cpm[top500,],
                variation_integrated[top500,])
  data_name = c("Raw UMI", "VST", "CPM", "Integrated")
  for (i in 1:length(v_list)){
    v_list[[i]]$data = data_name[i]
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
    if (nrow(dflist[[i]][[cur_key]])==0){next}
    if (value == "pval"){
      xlimit = c(0,1)
      ylimit = c(0,1)
      if (is.null(dflist[[i]][[cur_key]]$pval)){
        x = dflist[[i]][[cur_key]]$p_val
      }else
        x = dflist[[i]][[cur_key]]$pval
    }
    if (value == "log2FC"){
      xlimit = c(-3,3)
      ylimit = c(-3,3)
      if (is.null(dflist[[i]][[cur_key]]$log2FC)){
        x = dflist[[i]][[cur_key]]$avg_log2FC
      }else
        x = dflist[[i]][[cur_key]]$log2FC
    }
    y = dflist[[1]][[cur_key]][match(dflist[[i]][[cur_key]]$genes, dflist[[1]][[cur_key]]$genes),value]
    p[[i-1]] = ggplot(na.omit(data.frame(x = x, y = y)), aes(x = x, y = y)) +
      geom_point(alpha = 0.5, size = 0.5) + theme_minimal() + xlim(xlimit) + ylim(ylimit)+ 
      xlab(titlelist[i]) + ylab(titlelist[1])
  }
  figure = ggarrange(plotlist = p,nrow = 2, ncol = 4)
  
  # Annotate the figure by adding a common labels
  return(annotate_figure(figure,
                  top = paste(value, "comparison",cur_key)))
}
