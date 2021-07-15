#' Choose output file type
#'
#' @param x binmap object
#' @param outputfile output file name
#' @param filetype filetype

#' @importFrom utils write.table

#' @return result file
#' @export
#'

outfiletype<-function(x,outputfile,filetype){
  if (filetype=="qtl"){
    sample_name<-colnames(x)
    id<-x$SNP
    CHR<-x$CHR
    temp_map<-rbind(t(id),t(CHR))
    rownames(temp_map)<-c("id","")
    utils::write.table(temp_map,file = paste(outputfile,"qtl","csv",sep = "."),append = FALSE,sep = ",",quote = FALSE,row.names = TRUE, col.names =FALSE)

    for(i in sample_name[5:ncol(x)]){

      a<-lapply(x[[i]],FUN = function(y){
        if (is.na(y)==TRUE){
          "NA"
        }else if(y==0){
          "A"
        }else if(y==1){
          "H"
        }else if(y==2){
          "B"
        }
      })
      temp_new_map<-as.data.frame(rbind(t(id),unlist(a)))
      rownames(temp_new_map)<-c("id",i)
      write.table(temp_new_map[2,],file = paste(outputfile,"qtl","csv",sep = "."),sep = ",",append = TRUE,row.names = TRUE,col.names = FALSE,quote = FALSE)
    }

  }else{
      cat("This file type not support")
  }

}
