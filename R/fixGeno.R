#' Fix potential genotype error
#'
#' fix potential genotype error
#'
#' @param x a vector
#' @param fix.size define the neighbor size to fix the error
#'
#' @return a vector
#' @export
#'

fixGeno <- function(x,fix.size=10){

    wind.geno.rle <- rle(x)
    error.id <- which(wind.geno.rle$lengths < fix.size)
    for(i in error.id){
        left.id <- sum(wind.geno.rle$lengths[1:i]) - wind.geno.rle$lengths[i]
        right.id <- sum(wind.geno.rle$lengths[1:i])
        if(i==1){
            x[(left.id+1):right.id] <- x[right.id+1]
        }else if(x[left.id]==x[right.id+1]){
            x[(left.id+1):right.id] <- x[left.id]}
        }
  return(x)

}
