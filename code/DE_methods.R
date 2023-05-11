hippo_mean_DE = function(sce, cellgroup1, cellgroup2){
  tmpcount1 = sce@assays@data$counts[,cellgroup1]
  tmpcount2 = sce@assays@data$counts[,cellgroup2]
  if(length(cellgroup1) == 1){
    tmpmean1 = mean(tmpcount1)
  }else{
    tmpmean1 = Matrix::rowMeans(tmpcount1)
  }
  if(length(cellgroup2) == 1){
    tmpmean2 = mean(tmpcount2)
  }else{
    tmpmean2 = Matrix::rowMeans(tmpcount2)
  }
  df = data.frame(genes = rownames(sce@assays@data$counts))
  df$meandiff = abs(tmpmean1 - tmpmean2)
  pooled_mean = (length(cellgroup1)*tmpmean1 + length(cellgroup2)*tmpmean2)/(length(cellgroup1)+length(cellgroup2))
  df$sd = sqrt(pooled_mean*(1/length(cellgroup1)+1/length(cellgroup2)))
  df$z = df$meandiff/df$sd
  df = df[order(df$z, decreasing = TRUE), ]
  df$pval = pnorm(-df$z)
  df$BHpval = p.adjust(df$pval, method = "BH")
  return(df)
}

hippo_prop_DE = function(sce, cellgroup1, cellgroup2){
  tmpcount1 = sce@assays@data$counts[,cellgroup1]
  tmpcount2 = sce@assays@data$counts[,cellgroup2]
  if(length(cellgroup1) == 1){
    tmpprop1 = mean(tmpcount1==0)
  }else{
    tmpprop1 = Matrix::rowMeans(tmpcount1==0)
  }
  if(length(cellgroup2) == 1){
    tmpprop2 = mean(tmpcount2==0)
  }else{
    tmpprop2 = Matrix::rowMeans(tmpcount2==0)
  }
  df = data.frame(genes = rownames(sce@assays@data$counts))
  pooled_prop = (length(cellgroup1)*tmpprop1 + length(cellgroup2)*tmpprop2)/(length(cellgroup1)+length(cellgroup2))
  df$propdiff = abs(tmpprop1 - tmpprop2)
  df$sd = sqrt(pooled_prop*(1-pooled_prop)*(1/length(cellgroup1)+1/length(cellgroup2)))
  df$z = df$propdiff/df$sd
  df = df[order(df$z, decreasing = TRUE), ]
  df$pval = pnorm(-df$z)
  df$BHpval = p.adjust(df$pval, method = "BH")
  return(df)
}

hippo_poisson_glm_DE = function(sce, cellgroups, repgroups = NA, freq_expressed = 0.05){
  countdf = data.frame(cellgroups = cellgroups, repgroups = repgroups)
  pval = rep(NA,nrow(sce@assays@data$counts))
  mu = rep(NA,nrow(sce@assays@data$counts))
  beta_cellgroups = rep(NA,nrow(sce@assays@data$counts))
  status = rep("done", nrow(sce@assays@data$counts))
  Rsquared = rep(NA,nrow(sce@assays@data$counts))
  for(i in 1:nrow(sce@assays@data$counts)){
    countdf$count = round(pmax(sce@assays@data$counts[i,],0))
    if (mean(countdf$count!=0, na.rm = TRUE) <= freq_expressed) {
      if(mean(countdf$count!=0, na.rm = TRUE) == 0){
        status[i] = "zero mean"
        next
      }else{
        status[i] = "lowly expressed"
        next
      }
    }
    if (is.na(repgroups[1])){
      gm = stats::glm(count~cellgroups, family = poisson, data = countdf)
    }else{
      gm = stats::glm(count~cellgroups + repgroups, family = poisson, data = countdf)
    }
    if(!gm$converged){
      status[i] = "not converge"
      next}
    gm = summary(gm)
    pval[i] = gm$coefficients["cellgroups",4]
    mu[i] = gm$coefficients["(Intercept)", 1]
    beta_cellgroups[i] = gm$coefficients["cellgroups",1]
    Rsquared[i] = with(gm, 1 - deviance/null.deviance)
  } 
  df = data.frame(genes = rownames(sce@assays@data$counts))
  df$mu = mu
  df$beta_cellgroups = beta_cellgroups
  df$log2FC = log2(exp(beta_cellgroups))
  df$pval = pval
  df$status = status
  df$BH = p.adjust(df$pval, method = "BH")
  df$Rsquared = Rsquared
  return(df)
}


hippo_poisson_glmm_DE = function(sce, cellgroups, repgroups, freq_expressed = 0.05){
  countdf = data.frame(cellgroups = as.factor(cellgroups), repgroups = as.factor(repgroups))
  pval = rep(NA,nrow(sce@assays@data$counts))
  mu = rep(NA,nrow(sce@assays@data$counts))
  beta_cellgroups = rep(NA,nrow(sce@assays@data$counts))
  status = rep("done", nrow(sce@assays@data$counts))
  res_square = rep(NA,nrow(sce@assays@data$counts))
  REvariation = rep(NA,nrow(sce@assays@data$counts))
  FEvariation = rep(NA,nrow(sce@assays@data$counts))
  RESvariation = rep(NA,nrow(sce@assays@data$counts))
  
  for(i in 1:nrow(sce@assays@data$counts)){
    countdf$count = round(pmax(sce@assays@data$counts[i,],0))
    if (mean(countdf$count!=0, na.rm = TRUE) <= freq_expressed) {
      if(mean(countdf$count!=0, na.rm = TRUE) == 0){
        status[i] = "zero mean"
        next
      }else{
        status[i] = "lowly expressed"
        next
      }
    }
    gm = tryCatch(summary(MASS::glmmPQL(count~cellgroups, random = ~1|repgroups, family = poisson, data = countdf, verbose = FALSE)),
                error = function(e){NULL})
    if (is.null(gm)){
      status[i] = "not converge"
      next
    }
    gm_null = tryCatch(summary(MASS::glmmPQL(count~1, random = ~1|repgroups, family = poisson, data = countdf, verbose = FALSE)),
                  error = function(e){NULL})
    pval[i] = gm$tTable[2 ,"p-value"]
    res_square[i] = gm$sigma^2
    mu[i] = gm$coefficients$fixed[1]
    beta_cellgroups[i] = gm$coefficients$fixed[2]
    rsquared = tryCatch(r.squaredGLMM(gm, gm_null, pj2014 = T), error = function(e){NULL})
    if (!is.null(rsquared)){
      REvariation[i] = rsquared[1,2] - rsquared[1,1]
      FEvariation[i] = rsquared[1,1]
      RESvariation[i] = 1-rsquared[1,2]
    }
  } 
  df = data.frame(genes = rownames(sce@assays@data$counts))
  df$mu = mu
  df$beta_cellgroups = beta_cellgroups
  df$log2FC = log2(exp(beta_cellgroups))
  df$sigma_square = res_square
  df$pval = pval
  df$status = status
  df$BH = p.adjust(df$pval, method = "BH")
  df$REvariation = REvariation
  df$FEvariation = FEvariation
  df$RESvariation = RESvariation
  return(df)
}

hippo_binomial_glmm_DE = function(sce, cellgroups, repgroups, freq_expressed = 0.05){
  countdf = data.frame(cellgroups = as.factor(cellgroups), repgroups = as.factor(repgroups))
  pval = rep(NA,nrow(sce@assays@data$counts))
  beta_cellgroups = rep(NA,nrow(sce@assays@data$counts))
  mu = rep(NA,nrow(sce@assays@data$counts))
  status = rep("done", nrow(sce@assays@data$counts))
  res_square = rep(NA,nrow(sce@assays@data$counts))
  for(i in 1:nrow(sce@assays@data$counts)){
    countdf$count = 1*(sce@assays@data$counts[i,]>0)
    if (mean(countdf$count) <= freq_expressed) {
      if(mean(countdf$count) == 0){
        status[i] = "zero mean"
        next
      }else{
        status[i] = "lowly expressed"
        next
      }
    }
    genemean_bygroup = aggregate(count ~ cellgroups, data = countdf, FUN = mean)
    if (genemean_bygroup$count[1] == genemean_bygroup$count[2]){
      status[i] = "no difference between groups"
      next
    }
    gm = tryCatch(summary(MASS::glmmPQL(count~cellgroups, random = ~1|repgroups, family = binomial, data = countdf, verbose = FALSE, niter = 50)),
                  error = function(e){NULL})
    if (is.null(gm)){
      status[i] = "not converge"
      next
    }
    pval[i] = gm$tTable[2 ,"p-value"]
    beta_cellgroups[i] = gm$coefficients$fixed[2]
    mu[i] = gm$coefficients$fixed[1]
    res_square[i] = gm$sigma^2
  } 
  df = data.frame(genes = rownames(sce@assays@data$counts))
  df$beta_cellgroups = beta_cellgroups
  df$log2FC = log2(exp(beta_cellgroups))
  df$mu = mu
  df$sigma_square = res_square
  df$pval = pval
  df$status = status
  df$BH = p.adjust(df$pval, method = "BH")
  return(df)
}

pseudobulk_deseq2 = function(sce, cellgroups, repgroups){
  raw_count = melt(t(round(pmax(sce@assays@data$counts,0))), value.name = "count")
  colnames(raw_count) = c("cell", "gene", "count")
  raw_count$cellgroups = rep(as.factor(cellgroups), nrow(raw_count)/length(cellgroups))
  raw_count$repgroups = rep(as.factor(repgroups), nrow(raw_count)/length(repgroups))
  pseudobulk_count = aggregate(count ~ repgroups + cellgroups + gene , data = raw_count, FUN = sum)
  coldata = unique(pseudobulk_count[c("repgroups","cellgroups")])
  pseudobulk_count$cell_rep = paste0(pseudobulk_count$cellgroups, "_", pseudobulk_count$repgroups)
  pseudobulk_count = reshape(pseudobulk_count[,c("cell_rep","gene", "count")], idvar = "gene", timevar = c("cell_rep"), direction = "wide")
  rownames(pseudobulk_count) = pseudobulk_count[,1]
  pseudobulk_count[,1] = NULL
  dds = DESeqDataSetFromMatrix(countData = pseudobulk_count,
                               colData = coldata,
                               design= ~ repgroups + cellgroups)
  dds = DESeq(dds)
  df = results(dds)[c("pvalue", "padj", "baseMean", "log2FoldChange")]
  colnames(df) = c("pval", "BH", "basemean", "log2FC")
  return(df)
}

raw2TPM = function(raw_count){
  x = raw_count/nrow(raw_count)
  return(t(t(x)*1e6/colSums(x)))
}

MAST_DE = function(sce, cellgroups, repgroups, freq_expressed = 0.05){
  sca = FromMatrix(log2(edgeR::cpm(round(pmax(sce@assays@data$counts,0)))+1), 
                      data.frame(cellgroups = as.factor(cellgroups), repgroups = as.factor(repgroups)), 
                      data.frame(gene = rownames(sce@assays@data$counts)))
  expressed_genes = freq(sca) > freq_expressed
  sca = sca[expressed_genes,]
  colData(sca)$cdr = colSums(assay(sca)>0)/dim(sca)[1]
  zlmCellgroups = zlm( ~ cellgroups + cdr + repgroups, sca)
  lrt = lrTest(zlmCellgroups, "cellgroups")
  df = data.frame(genes = rownames(sca))
  df$coef_cellgroups = zlmCellgroups@coefC[,2]
  df$log2FC = log2(exp(df$coef_cellgroups))
  df = merge(df, data.frame(genes = names(lrt[,3,3]), pval = lrt[,3,3]), by = "genes", all.x = TRUE)
  df$BH = p.adjust(df$pval, method = "BH")
  return(df)
}

identifyhits = function(BH, log2FC, log2mean, log2meandiff = -Inf, 
                        pvalucutoff = 0.05, log2FCcutoff = log2(1.5), 
                        log2meancutoff = -2.25, log2meandiffcutoff = -1, newcriteria = F){
  if(newcriteria){
    hits = BH<pvalcutoff & abs(log2FC)>log2FCcutoff & (log2mean > log2meancutoff | log2meandiff > log2meandiffcutoff)
  }else{
    hits = BH<pvalcutoff & abs(log2FC)>log2FCcutoff
  }
  
  hits = ifelse(is.na(BH), NA, hits)
  return(hits)
}

