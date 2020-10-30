

### objective function for progression model
# parameters: h, t_diab, t_npdr, eta, K
getRMSE_filter = function(p, x, y, df.rank, fitControl, file.out){
  # check
  print(p)
  t_diab = p[1]
  t_npdr = p[2]
  eta = p[3]
  K = p[4]
  n = p[5]
  # subset data
  if (t_diab < t_npdr){
    top.feats = df.rank[prog == 'progAB'][order(rank)][1:n, ensemblID]
  } else {
    top.feats = df.rank[prog == 'progBA'][order(rank)][1:n, ensemblID]
  }
  x.top = x[, top.feats]
  # assign y values for regression
  y.num = meta %>%
    mutate(tmp = disease_group_detailed) %>%
    mutate(tmp = gsub('Control', 1, tmp)) %>%
    mutate(tmp = gsub('Diabetic', t_diab, tmp)) %>%
    mutate(tmp = gsub('NPDR [+] DME', 4, tmp)) %>%
    mutate(tmp = gsub('NPDR', t_npdr, tmp)) %>%
    .$tmp %>%
    as.numeric
  # carry out cross-validation
  spls.fit = caret::train(
    x.top,
    y.num,
    method = "spls",
    preProcess = c("center", "scale"),
    trControl = fitControl,
    tuneGrid = data.frame(
      eta = eta, K = K, kappa = 0.5))
  # return mean RMSE
  meanRMSE = mean(spls.fit$resample$RMSE)
  write(paste0(c(p, meanRMSE), collapse ='\t'), file=file.out, append=T)
  print(meanRMSE)
  return(meanRMSE)
}



#' 
#' Compare List of Genesets
#'
#' This function uses simple hypergeomtric tests to determine the
#' overlap between sets in a list of sets.
#'  
#' @param geneSets  list of gene sets
#' @param bkgSet set of background genes
#' 
#' @return object$setInfo, object$setTable, object$matrix$...
#' 
compareSets = function(geneSets, bkgSet, setInfoIn = NULL, minSize=10){
  # remove all genes not in background
  geneSets = lapply(geneSets, function(geneSet)
    intersect(geneSet, bkgSet))
  # remove sets with size < minSize
  geneSets = geneSets[!(lapply(geneSets, length) < minSize)]
  # statistics on individual gene sets
  setInfoOut = data.table(
    setName=names(geneSets),
    setSize=sapply(geneSets, function(x) length(x))
  )
  # merge with setInfoIn
  if (!is.null(setInfoIn)){
    setInfoOut = merge(setInfoIn, setInfoOut, by='setName')
  }
  # prepare matrices
  mat = matrix(nrow=length(geneSets), ncol=length(geneSets))
  colnames(mat) = names(geneSets)
  rownames(mat) = names(geneSets)
  mat.p = mat.n = mat.exp = mat.fold = mat
  # comparison of genesets TODO: parallelize
  for (i in 1:length(geneSets)){
    for (j in 1:length(geneSets)){
      nameA = names(geneSets)[i]
      nameB = names(geneSets)[j]
      setA = geneSets[[i]]
      setB = geneSets[[j]]
      setAB = intersect(setA, setB)
      p = phyper(
        length(setAB), 
        length(setA), 
        length(bkgSet) - length(setA), 
        length(setB), lower.tail = F)
      # write matrices
      mat.p[i, j] = p
      mat.n[i, j] = length(setAB)
      mat.exp[i, j] = length(setA) / length(bkgSet) * length(setB)
      mat.fold = mat.n / mat.exp
    }
  }
  # combine matrices to data table
  df.p = getUpper(mat.p, value.name = 'pVal')
  df.n = getUpper(mat.n, value.name = 'observed')
  df.exp = getUpper(mat.exp, value.name = 'expected')
  df.fold = getUpper(mat.fold, value.name = 'fold')
  setTable = Reduce(function(x, y) merge(x, y, by=c('Var1', 'Var2')), list(df.exp, df.n, df.fold, df.p))
  setTable$setSize1 = setInfoOut[match(setTable$Var1, setInfoOut$setName), setSize]
  setTable$setSize2 = setInfoOut[match(setTable$Var2, setInfoOut$setName), setSize]
  # multiple testing
  setTable$fdr = p.adjust(setTable$pVal, method = 'BH')
  # return
  out = list()
  out$setInfo = setInfoOut
  out$setTable = setTable
  out$matrix$pVal = mat.p
  out$matrix$observed = mat.n
  out$matrix$expected = mat.exp
  out$matrix$fold = mat.fold
  out
}

# upper triangle of matrix to dataframe
getUpper = function(m, value.name = 'value'){
  rowCol = expand.grid(rownames(m), colnames(m))
  df = rowCol[as.vector(upper.tri(m,diag=F)),]
  df[, value.name] = m[upper.tri(m,diag=F)]
  data.table(df)
}

# function for CPM plots
plotCPM = function(genes){
  # gather data
  df.CPMplot = lapply(genes, function(gene)
    iowa.metadata %>%
      tibble::rownames_to_column('sampleID') %>%
      dplyr::select(sampleID, sample_site, disease_group_dme) %>%
      mutate(ensemblID = gene) %>%
      mutate(geneID = df.map[OriginalName == gene, geneID]) %>%
      mutate(expression = t(iowa.full[gene, sampleID])) %>%
      magrittr::set_colnames(
        c('sampleID', 'Tissue', 'Disease',
          'ensemblID', 'geneID', 'logCPM'))) %>%
    bind_rows %>%
    mutate(Disease = factor(Disease, 
                            levels = c('Control', 'Diabetic', 'NPDR', 'NPDR/PDR + DME'))) %>%
    data.table  
  # plot
  ggplot(df.CPMplot, aes(x=Disease, y=logCPM, fill=Tissue)) +
    geom_boxplot(outlier.size = 1, alpha=0.7) +
    geom_jitter(position=position_jitterdodge(), size=0.65) +
    scale_fill_manual(values = c('#ed7d31', '#5b9bd5')) +
    labs(y='CPM [log2]') +
    ggpubr::grids() +
    theme(
      axis.text.x = element_text(angle=90, hjust=1, vjust=0.4),
      axis.title.x = element_blank(),
      text=element_text(size=6),
      plot.background = element_rect(fill = "transparent", color = NA), 
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major = element_line(),
      panel.grid.minor = element_line(),
      rect = element_rect(fill = "transparent"),
      legend.key.size = unit(0.75,"line"),
      legend.position = 'bottom') +
    facet_wrap(~geneID, ncol=5)}





# function for confounder correction
compensate.confounders <- function(dataset = NULL, confounders = NULL, parallel = FALSE, n.cluster = 20, pval.thresh = 0.05, verbosity = TRUE, offset2zero = TRUE){
  function.env <- new.env()
  library(parallel)
  if(is.null(dataset) | is.null(confounders)){
    stop("Dataset or counfounder table not provided!")
  }
  cat('running koljas version of cC\n')
  
  # check for variables (genes) which are 0 in all samples
  # NOTE: does this even make sense if not operating on counts data?
  # NOTE: possibly move this step outside of function
  cat("Looking for variables that are 0 in all samples...\n", fill = TRUE)
  all.zero.vars <- apply(X = dataset, 
                         MARGIN = 2, 
                         FUN = function(gene){
                           sum(gene) == 0 
                         })
  
  if(sum(all.zero.vars) > 0){
    
    cat(paste(sum(all.zero.vars), "variables have been removed:\n"))
    if(verbosity){
      cat(paste(colnames(dataset)[which(all.zero.vars)],"\n",sep = ""))
    }
    dataset <- dataset[, !all.zero.vars]
    
  }else{
    
    cat("No variables were removed", fill = TRUE)
  }
  
  
  # ...
  if(parallel) {
    cat("Parallel run set up...\n", fill = TRUE)
    
    cl <- makeCluster(n.cluster, envir = function.env)
    clusterExport(cl, list("confounders", "dataset", "pval.thresh"), envir = function.env)
    
    data.corr <- parApply(cl, 
                          X = dataset,
                          MARGIN = 2,
                          FUN = function(var){
                            # creates a dataframe with gene expression and confounder info
                            var.dat<- data.frame(var, confounders)
                            # creates formular for confunder (as required by the linear model)
                            confounder.formula <- paste("var ~", paste(colnames(var.dat)[-1], collapse = " + "))
                            # fits linear model using confounder formula and gene specific data
                            model <- lm(confounder.formula, var.dat)
                            # analyse model fit using anova
                            an <- anova(model)
                            # collect p values from model fit
                            pvals <- as.matrix(subset(an, subset = rownames(an) != "Residuals", select = "Pr(>F)"))
                            # if any p-value is below significance threshold ...
                            if(sum(pvals <= pval.thresh)){
                              # create confounder formula for significant regressors
                              confounder.formula.sig <- paste("var ~", paste(rownames(pvals)[which(pvals <= pval.thresh)], collapse = " + "))
                              # fit model once again using only significant regressors
                              model.sig <- lm(confounder.formula.sig, var.dat)
                              # correcting for confounding effect? Does this also work for multiple variables?
                              data.out <- residuals(model.sig) + coef(model.sig)["(Intercept)"]
                            }else{
                              data.out <- var.dat$var
                            }
                            data.out
                          })
    stopCluster(cl)
    
    
  }else{
    cat("Non-parallel run set up...")
    data.corr <- apply(X = as.data.frame(dataset),
                       MARGIN = 2,
                       FUN = function(var){
                         
                         var.dat<- data.frame(var, confounders)
                         
                         confounder.formula <- paste("var ~", paste(colnames(var.dat)[-1], collapse = " + "))
                         model <- lm(confounder.formula, var.dat)
                         an <- anova(model)
                         pvals <- as.matrix(subset(an, subset = rownames(an) != "Residuals", select = "Pr(>F)"))
                         if(sum(pvals <= pval.thresh) > 0){
                           confounder.formula.sig <- paste("var ~", paste(rownames(pvals)[which(pvals <= pval.thresh)], collapse = " + "))
                           model.sig <- lm(confounder.formula.sig, var.dat)
                           data.out <- residuals(model.sig) + coef(model.sig)["(Intercept)"]
                         }else{
                           data.out <- var.dat$var
                         }
                       })
    if(class(data.corr) == "list"){
      
      data.corr <- do.call("cbind", data.corr)
    }
    
  }
  
  cat(paste0('Number of corrected rows: ',
      sum(colSums(dataset) != colSums(data.corr)), '\n'))
  
  # rownames(data.corr) <- rownames(dataset)
  
  if(offset2zero & min(data.corr) < 0){
    
    data.corr <- data.corr - min(data.corr)
    
  }
  
  cat("Correction done!")
  
  return(data.corr)
}