#' Pad numbers
#'  @param x numeric sequence
#'  

num_pad <- function(o) {
  max_size = nchar(as.character(max(o)))
  string = rep("", length(o))
  for(i in 1:length(o)){
    if(nchar(as.character(o[i])) < max_size){
      dif = max_size - nchar(as.character(o[i]))
      string[i] = paste0(paste(rep("0",dif), collapse=""), o[i])
    }else{
      string[i] = o[i]
    }
  }
  return(string)
}