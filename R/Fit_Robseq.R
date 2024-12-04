#### Main Function ####
#' Robseq Function
#'
#' @param features A dataframe with gene expression counts in the row and samples in the column.
#' @param metadata A dataframe with information on the subjects such as disease status, gender and etc.
#' @param norm.method The normalization method to be used. The user can choose from 5 different methods such as TMM, RLE, CPM, Upper quartile and Qauntile. ‘TMM’ by default.
#' @param expVar The name of the variable on which the differential expression will be evaluated. If the user provides no name then the metadata should have a column named as exposure which should have information on things such as disease status, treatment conditions or etc. ‘Exposure’ by default.
#' @param coVars The names of the covariates/confounders that needs to adjusted for in the differential expression analysis. 'NULL' by default
#' @param filter If true, only genes with sufficiently large counts are retained. 'FALSE' by default
#' @param parallel If true, the analysis will be performed on multiple cores with faster runtimes.'FALSE' by default
#' @param ncores The number of cores on which the analysis will be serially performed. The user needs to specify this only when parallel = TRUE.'1' by default
#' @param verbose If true, it will print the progress report. 'FALSE' by default
#'
#' @return output
#'
#'
#' @examples
#'\dontrun{Robseq(features, metadata)}
#'
#'@export
robust.dge <- function(features,
                       metadata,
                       norm.method = 'RLE',
                       expVar = 'Exposure',
                       coVars = NULL,
                       filter = FALSE,
                       parallel = FALSE,
                       ncores = 1,
                       verbose = FALSE){

  #### Creating Regression Pre-requisites ####

  start.time <- Sys.time()
  if(is.null(coVars)){
    regData <- metadata[, c(expVar), drop = FALSE]
  }else{
    regData <- metadata[, c(expVar, coVars)]
  }
  regData[sapply(regData, is.character)] <- lapply(regData[sapply(regData, is.character)], as.factor)
  formula <- as.formula(paste("expr ~ ", paste(colnames(regData), collapse = "+")))

  #### Normalizing Expression counts ####

  if(norm.method == 'TMM'){
    norm.y = suppressMessages(tmm_norm(features, metadata))
    norm.y = log2(norm.y + 0.5)
  }else if(norm.method == 'RLE'){
    norm.y = suppressMessages(rle_norm(features, metadata))
    norm.y = log2(norm.y + 0.5)
  }else if(norm.method == 'CPM'){
    norm.y = data.frame(cpm(features, log = TRUE, prior.count = 1))
  }else if(norm.method == 'Quantile'){
    norm.y = suppressMessages(quan_norm(features, metadata))
    norm.y = log2(norm.y + 0.5)
  }else if(norm.method == 'UQuantile'){
    norm.y = suppressMessages(uqrt_norm(features, metadata))
    norm.y = log2(norm.y + 0.5)
  }

  #### Apply model ####

  if(parallel){
    if(verbose){
      cl <- parallel::makeCluster(ncores)
      doSNOW::registerDoSNOW(cl)
      packages <- c('MASS', 'dfadjust')
      exports <- c('perGene.mod')
      pb <- txtProgressBar(max = nrow(norm.y), style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      res <- foreach(j = 1:nrow(norm.y), .combine = rbind,
                     .packages = packages, .options.snow = opts,
                     .export = exports) %dopar% {
                       expr <- as.numeric(norm.y[j, ])
                       tmpfit <- perGene.mod(expr = expr,
                                             formula = formula,
                                             regData = regData,
                                             expVar = expVar)
                       return(tmpfit)
                     }
      close(pb)
      stopCluster(cl)
    }else{
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      packages <- c('MASS', 'dfadjust')
      exports <- c('perGene.mod')
      res <- foreach(j = 1:nrow(norm.y), .combine = rbind, .packages = packages, .export = exports) %dopar% {
                       expr <- as.numeric(norm.y[j, ])
                       tmpfit <- perGene.mod(expr = expr,
                                             formula = formula,
                                             regData = regData,
                                             expVar = expVar)
                       return(tmpfit)
                     }
      stopCluster(cl)
    }
  }else{
    if(verbose){
      res <- pbapply(norm.y, 1, function(x) perGene.mod(expr = x,
                                                        formula = formula,
                                                        regData = regData,
                                                        expVar = expVar))
    }else{
      res <- apply(norm.y, 1, function(x) perGene.mod(expr = x,
                                                      formula = formula,
                                                      regData = regData,
                                                      expVar = expVar))
    }
    res <- do.call('rbind', res)
  }

  #### Summarizing Output ####

  Genes <- rownames(norm.y)
  output <- data.frame(Genes = Genes, res)
  output$adjPval <- p.adjust(output$Pval, method = 'BH')
  output <- output[order(output$adjPval), ]
  stop.time <- Sys.time()
  Time.min = round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  finres = list(Method = 'Robseq',
                res = output,
                features = features,
                Time.min = Time.min)
  return(finres)
}
