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
        wind_sum<- dplyr::mutate(x,group=(as.numeric(rownames(x))-1) %/% window.size) %>%
        dplyr::group_by(.data$group) %>%
        dplyr::summarise(start=min(.data$POS),end=max(.data$POS),code_sum=ifelse(sum(is.na(.data$code))>=floor(window.size*2/3),NA,sum(.data$code)*window.size/(window.size-sum(is.na(.data$code))))) %>%
        dplyr::ungroup()
        wind_geno<-dplyr::mutate(wind_sum,code=ifelse(is.na(.data$code_sum),NA,
                                                      ifelse(.data$code_sum < (window.size*2*low), 0,
                                                             ifelse(.data$code_sum > (window.size*2*high), 2, 1 )
                                                             )
                                                      )
                                 )

        return(wind_geno)
    }else if(window.type=="kb"){

        wind_sum<-dplyr::mutate(x,group=as.numeric(x[,2]) %/% (window.size*1000)) %>%
        dplyr::group_by(.data$group) %>%
        dplyr::summarise(start=min(.data$POS),end=max(.data$POS),code_sum=ifelse(sum(is.na(.data$code))>=floor(length(.data$code)*2/3),NA,sum(.data$code)/(length(.data$code)-sum(is.na(.data$code))))) %>%
        dplyr::ungroup()
        wind_geno<-dplyr::mutate(wind_sum,code=ifelse(is.na(.data$code_sum),NA,
                                                      ifelse(.data$code_sum < (2*low), 0,
                                                             ifelse(.data$code_sum > (2*high), 2, 1 )
                                                             )
                                                      )
                                 )

        return(wind_geno)
    }

}
