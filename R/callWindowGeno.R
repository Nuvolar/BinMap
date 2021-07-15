#' Call genotype
#'
#'get genotype of a chrom of single sample
#'
#' @param x data.frame
#' @param window.type "number" ,"kb" two types of call window default is number
#' @param window.size  window size dafalut is 15
#' @param low default is 0.2
#' @param high default is 0.8
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @return data.frame
#' @export
#'
#'
callWindowGeno<-function(x,
                         window.type="number",window.size=15,
                         low=0.2,high=0.8){
    if (window.type=="number"){
        #print(window.type=="number")
        wind_sum<- dplyr::mutate(x,group=(as.numeric(rownames(x))-1) %/% window.size) %>%
        dplyr::group_by(.data$group) %>%
        dplyr::summarise(start=min(.data$POS),end=max(.data$POS),
                         code_sum=if(sum(is.na(.data$code))>=floor(window.size*2/3)){
                             NA
                         }else{
                             round(sum(.data$code,na.rm = TRUE)*window.size/(window.size-sum(is.na(.data$code))))
                         }) %>%
        dplyr::ungroup()
        wind_geno<-wind_sum
        code<-lapply(wind_sum$code_sum,FUN = function(x){
            if (is.na(x)==TRUE){
                NA
            }else if(x < (window.size*2*low)){
                0
            }else if(x > (window.size*2*high)){
                2
            }else{
                1
            }
        })
        wind_geno$code<-unlist(code)

    }else if(window.type=="kb"){

        wind_sum<-dplyr::mutate(x,group=as.numeric(x[,2]) %/% (window.size*1000)) %>%
        dplyr::group_by(.data$group) %>%
        dplyr::summarise(start=min(.data$POS),end=max(.data$POS),code_sum=ifelse(sum(is.na(.data$code))>=floor(length(.data$code)*2/3),NA,
                                                                                 round(sum(.data$code,na.rm = T)/(length(.data$code)-sum(is.na(.data$code))),digits = 3)
                                                                                 )) %>%
        dplyr::ungroup()
        wind_geno<-wind_sum
        code<-lapply(wind_sum$code_sum,FUN = function(x){
            if (is.na(x)==TRUE){
                NA
            }else if(x < (2*low)){
                0
            }else if(x > (2*high)){
                2
            }else{
                1
            }
        })
        wind_geno$code<-unlist(code)
    }
    return(wind_geno)
}
