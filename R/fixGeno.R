#' Fix potential genotype error
#'
#' fix potential genotype error
#'
#' @param t binmap object
#' @param fix.size define the neighbor size to fix the error
#'
#' @return a vector
#' @export
#'

fixGeno <- function(t,fix.size=10){
    #print(t)
    #print(class(t))
    x <- t$code
    #is.na(x)
    wind.geno.rle <- rle(x)
    error.id <- which(wind.geno.rle$lengths < fix.size)
    for(i in error.id){
        left.id <- sum(wind.geno.rle$lengths[1:i]) - wind.geno.rle$lengths[i]
        right.id <- sum(wind.geno.rle$lengths[1:i])
        if(i==1&right.id==length(x)){
            next
        }else if(i==1&right.id<length(x)&is.na(x[right.id+1])){
            next
        }else if(i==1&right.id<length(x)&is.na(x[right.id+1])==FALSE){
            x[(left.id+1):right.id] <- x[right.id+1]
        }else if(right.id==length(x)&is.na(x[left.id])){
            next
        }else if(right.id==length(x)&is.na(left.id)==FALSE){
            x[(left.id+1):right.id] <- x[left.id]
        }else if(is.na(x[left.id])&is.na(x[right.id+1])){
            next
        }else if(is.na(x[left.id])&is.na(x[right.id+1])==FALSE){
            x[(left.id+1):right.id]<- x[right.id]
        }else if(is.na(x[left.id])==FALSE&is.na(x[right.id+1])){
            x[(left.id+1):right.id]<- x[left.id]
        }else if(x[left.id]==x[right.id+1]){
            x[(left.id+1):right.id] <- x[left.id]}
    }
    t$fix_code <- x
    return(t)

}
