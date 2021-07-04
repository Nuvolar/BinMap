#' Data Input
#'
#' Read a file in table format from vcf/ped/bed
#'
#' @param inputfile input file name
#' @param datatype inputfile type. support vcf/ped/bed
#' @param father support a character vector or a number that father index of sample
#' @param mother support a character vector or a number that father index of sample
#' @param outputfile output file name
#' @param screening a logical value. use plink to screening
#' @param maf if screening == T ,it useful
#' @param geno if screening == T ,it useful
#' @param mind if screening == T ,it useful
#' @param hwe if screening == T ,it useful
#'
#' @importFrom data.table fread
#' @importFrom utils write.table
#' @return data.frame
#' @export
#'
#'

readgeno<-function(inputfile, datatype,outputfile,
                   screening=FALSE,
                   maf=NA,
                   geno=NA,
                   mind=NA,
                   hwe=NA,
                   father = "1/1",
                   mother = "0/0"){

    p_sys<-system("plink", intern = TRUE)

    if (datatype!="bed"&datatype!="ped"&datatype!="vcf"){
        cat(paste(datatype,"not support"))
        return(NULL)
        }


    if ((datatype=="bed"|datatype=="ped")&attributes(p_sys)$status != 10){
        cat("  please make sure plink soft in your enviroment path    ")
        return(NULL)
    }

    datatype_cmd<-""
    if (datatype=="bed"){
        datatype_cmd<-paste("plink --bfile", inputfile,sep = " ")
    }else if(datatype=="ped"){
        datatype_cmd<-paste("plink --file", inputfile,sep = " ")
    }else if(datatype=="vcf"){
        datatype_cmd<-paste("plink --vcf", paste0(inputfile,".vcf"),"--vcf-half-call missing ",sep = " ")
    }

    fit_cmd<-""
    if (screening == TRUE){
        if (attributes(p_sys)$status != 10){
            cat("  please make sure plink soft in your enviroment path    ")
            return(NULL)
        }
        maf_cmd  <- ifelse ( class(maf)  == "numeric" & maf  <=1.0 & maf  >=0.0, paste("--maf" , maf ,sep=" "),"")
        geno_cmd <- ifelse ( class(geno) == "numeric" & geno <=1.0 & geno >=0.0, paste("--geno", geno,sep=" "),"")
        mind_cmd <- ifelse ( class(mind) == "numeric" & mind <=1.0 & mind >=0.0, paste("--mind", mind,sep=" "),"")
        hwe_cmd  <- ifelse ( class(hwe)  == "numeric" & hwe  <=1.0 & hwe  >=0.0, paste("--hwe" , hwe ,sep=" "),"")
        fit_cmd  <- paste(maf_cmd,geno_cmd,mind_cmd,hwe_cmd)
    }

    ####
    if (screening==TRUE&(length(father)>1|length(mother)>1)){

        total_cmd<-paste(datatype_cmd," --recode vcf-iid --out",paste0(outputfile,".temp"),sep = " ")
        system(total_cmd,intern = TRUE)

        geno<-data.table::fread(paste(outputfile,"temp","vcf",sep = "."),header=T,sep="\t")
        if (length(father)!=1&length(father)!=nrow(geno)) return("father input is wrong")
        if (length(mother)!=1&length(mother)!=nrow(geno)) return("mother input is wrong")

        if (length(father)==1&class(father)=="numeric") father<-geno[,(father+10)]
        if (length(mother)==1&class(mother)=="numeric") mother<-geno[,(mother+10)]

        if (length(father)==nrow(geno)) father<-as.character(lapply(father, FUN = function(x)substr(x,1,3)))
        if (length(mother)==nrow(geno)) mother<-as.character(lapply(mother, FUN = function(x)substr(x,1,3)))

        geno<-cbind(geno,father,mother)
        colnames(geno)<-c(colnames(geno)[1:(ncol(geno)-2)],"TempFather","TempMother")
        write.table(geno,file=paste(outputfile,"vcf",sep="."),sep="\t",append = FALSE,row.names = FALSE,col.names = TRUE,quote = FALSE)
        father<-(ncol(geno)-10);mother<-(ncol(geno)-9)

        total_cmd<-paste("plink --vcf", paste0(outputfile,".vcf"),fit_cmd,"--recode vcf-iid --out",outputfile,sep = " ")
        system(total_cmd,intern = TRUE)
    }else if (screening==TRUE&length(father)==1&length(mother)==1){
        total_cmd<-paste(datatype_cmd,fit_cmd," --recode vcf-iid --out",outputfile,sep = " ")
        system(total_cmd,intern = TRUE)
    }

    # read vcf file
    ifelse(screening==FALSE,geno<-data.table::fread(paste(inputfile,"vcf",sep = "."),header=T,sep="\t"),geno<-data.table::fread(paste(outputfile,"vcf",sep = ".")))
    sample_name<-colnames(geno[,10:length(geno)])
    temp_geno<-lapply(geno[,10:length(geno)], FUN = function(x)substr(x,1,3))
    temp_geno<-as.data.frame(temp_geno)

    if (length(father)==1&class(father)=="numeric") father<-temp_geno[,father]
    if (length(mother)==1&class(mother)=="numeric") mother<-temp_geno[,mother]

    geno<-cbind(geno[,c(1,2,4,5)],temp_geno,father,mother)

    colnames(geno)<-c("CHR","POS","P1","P2",sample_name,"father","mother")
    rm(temp_geno)

    code_num<-ncol(geno)
    write.table(geno[0,1:(code_num-2)],file=paste(outputfile,"temp","txt",sep="."),sep="\t",append = FALSE,row.names = FALSE,col.names = TRUE,quote = FALSE)

    apply(geno,MARGIN = 1,FUN=function(x){
        f_value<-x[code_num-1]
        m_value<-x[code_num]
        if(f_value=="0/0"){
            temp_p<-x[3]
            x[3]<-x[4];x[4]<-temp_p

        if (f_value!=m_value&(f_value=="1/1"|f_value=="0/0")&(m_value=="1/1"|m_value=="0/0")){
            x<-as.data.frame(t(x))
            x[,5:code_num]<-lapply(x[,5:code_num], FUN = function(a){
                if (a==f_value) {
                    a<-2
                }else if(a==m_value){
                    a<-0
                }else if(a=="0/1"){
                    a<-1
                }else if(a=="./."|a=="./0"|a=="./1"){
                    a<-NA
                }else{
                    print("there are something wrong")
                }
            })
        }

            write.table(x[,1:(code_num-2)],file=paste(outputfile,"temp","txt",sep="."),sep="\t",append = TRUE,row.names = FALSE,col.names = FALSE,quote = FALSE)
        }
    })
    geno<-data.table::fread(paste(outputfile,"temp","txt",sep="."),header=T,sep="\t")
    return(geno)
}




