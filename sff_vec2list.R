#' Transform vector of results of sff_Mediation to list
#' 
#' Provides quick plotting and extraction of estimated functions from the bootstrap output of sff_Mediation.
#' @param result a numeric vector with results of sff_Mediation
#' @param result_names a character vector with names of functions. The names must have the words: 'afunction', 'bsurface', 'gfunction', 'd1function', 'd2function', 'absurface', 'abline'
#' @param stats a string 'all' indicates that all functions will be plotted and/or returned. \code{stats} should be one of: 'afunction', 'bsurface', 'gfunction', 'd1function', 'd2function', 'absurface', 'abline'
#' @param plot a logical scalar to indicate if results should be plotted
#' @param ask a logical scalar to indicate if R should ask before switching to the next plot
#' @param returns a logical scalar to indicate whether output be produced. If \code{stat} is 'all', then a list with all functions and surfaces is returned. Else, a vector (for functions) or matrix (for surfaces) is returned.
#' @param surface a list with two entries (s, t). s are the rows (time parameter for mediator) and t are columns (time parameter for outcome) to extract from surface. Default is to extract the whole surface

sff_vec2list <- function(result, len=NULL,result_names=names(result), stats="all", plot=FALSE, ask=TRUE, returns=TRUE, surface=NULL){
  results = list()
  
  results$afunction  = result[grep('^afunction_', result_names)]
  if(is.null(len)){
    len = length(results$afunction)
  }
  
  if(stats=="all"){
  results$absurface  = t(matrix(result[grep("^absurface_",result_names)], len, len))
  results$abfunction = result[grep('^abfunction_', result_names)]
  results$bsurface   = t(matrix(result[grep("^bsurface_",result_names)], len, len))
  results$gfunction  = result[grep('^gfunction_', result_names)]
  results$d1function = result[grep('^d1function_', result_names)]
  results$d2function = result[grep('^d2function_', result_names)]
  
  # Plot results
  if(plot==TRUE){
    mipar = par()$mfrow
    par(mfrow=c(2,1),ask = ask)
    
    tfine = seq(0,1, length.out=len)
    plot(tfine, results$d1function, type="l", main="'Delta1' function")
    plot(tfine, results$afunction, type="l", main="'Alpha' function")
    
    plot(tfine, results$d2function, type="l", main="'Delta2' function")
    plot(tfine, results$gfunction, type="l", main="'Gamma' function")
    
    par(mfrow=c(1,1))
    image(t(results$bsurface), col  = gray((0:32)/32), main="'Beta' function")
    image(t(results$absurface), col  = gray((0:32)/32), main="'Alpha*Beta' surface")
    plot(tfine, results$abfunction, type="l", main="'Alpha*Beta' integral")
    par(mfrow=mipar, ask=FALSE)
  }
  }else if(stats %in% c('afunction', 'bsurface', 'gfunction', 'd1function', 'd2function', 
                       'absurface', 'abfunction')){
    if(length(grep('surface', stats)) == 0){
      results = result[grep(paste0('^',stats), result_names)]
      
      if(plot==TRUE){
        tfine = seq(0,1, length.out=len)
        plot(tfine, results, type="l", main=stats)
      }
    }else{
      if(is.null(surface)){
      results = t(matrix(result[grep(paste0('^',stats),result_names)], len, len))
      
      if(plot==TRUE){
        tfine = seq(0,1, length.out=len)
        image(t(results), col  = gray((0:32)/32), main=stats)
      }}else{
        if(surface$s == "all") surface$s=1:len
        if(surface$t == "all") surface$t=1:len
        
        results = t(matrix(result[grep(paste0('^',stats),result_names)], len, len))
        results = results[surface$s, surface$t]
        
        if(plot==TRUE & length(dim(results)) == 0){
          tfine = seq(0,1, length.out=len)
          plot(tfine, results, type="l", main=stats)
        }else{
          if(plot==TRUE & length(dim(results)) == 2){
            tfine = seq(0,1, length.out=len)
            image(t(results), col  = gray((0:32)/32), main=stats)  
        }}
      }
    }
    
  }else{
    stop("stats is not one of: 'afunction', 'bsurface', 'gfunction', 'd1function', 'd2function', 
         'absurface', 'abline'")
  }
  
  if(returns==TRUE)return(results)
}
  