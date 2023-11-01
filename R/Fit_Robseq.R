#### Main Function ####
#' Robseq Function
#'
#' @param features Gene Expressions
#' @param metadata Native Tissue or Tumor
#' @param norm.method Normalizing method
#' @param expVar 'Exposure' by default
#' @param coVars 'NULL' by default
#' @param parallel 'FALSE' by default
#' @param ncores '1' by default
#'
#' @return output
#'
#'
#' @examples
#'\dontrun{Robseq(features, metadata)}
#'
#'@export
Robseq <- function(features,
                   metadata,
                   norm.method = 'tmm',
                   expVar = 'Exposure',
                   coVars = NULL,
                   parallel = FALSE,
                   ncores = 1){

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

  if(norm.method == 'tmm'){
    norm.y = suppressMessages(tmm_norm(features, metadata))
    norm.y = log2(norm.y + 0.5)
  }else if(norm.method == 'rle'){
    norm.y = suppressMessages(rle_norm(features, metadata))
    norm.y = log2(norm.y + 0.5)
  }else{
    norm.y = data.frame(cpm(features, log = TRUE, prior.count = 1))
  }

  #### Apply model ####

  if(parallel){
    cl <- parallel::makeCluster(ncores)
    doSNOW::registerDoSNOW(cl)
    packages <- c('MASS', 'dfadjust')
    exports <- c('perGene.mod', 'expVar', 'coVars')
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
    res <- pbapply(norm.y, 1, function(x) perGene.mod(expr = x,
                                                      formula = formula,
                                                      regData = regData,
                                                      expVar = expVar))
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
