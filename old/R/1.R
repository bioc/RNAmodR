.get_weights <- function(file){
  args <- RNAmodRargs(file)
  weights <- as.numeric(strsplit(getParam(args, "RiboMeth", "weights"),";")[[1]])
  names(weights) <- as.numeric(strsplit(getParam(args, "RiboMeth", "weights_rel_pos"),";")[[1]])
  return(weights)
}


x$around <- unlist(lapply(seq_len(nrow(x)), function(i){
  
}))

x <- data.frame(position = ks01$position,
                counts = ks01$counts,
                ratio = unlist(lapply(seq_len(nrow(ks01)),function(i){ifelse((i-1) == 0,NA,ks01[i,]$counts/ks01[i-1,]$counts)})),
                ratio2 = unlist(lapply(seq_len(nrow(ks01)),function(i){ifelse((i-1) == 0,NA,ks01[i-1,]$counts/ks01[i,]$counts)})),
                scoreA = .aggregate_score_A(ks01,ks01$position,.get_weights("inst/extdata/RiboMeth.param")),
                scoreB = .aggregate_score_B(ks01,ks01$position,.get_weights("inst/extdata/RiboMeth.param")))