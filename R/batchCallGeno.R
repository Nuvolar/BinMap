#' Genotyping all sample by chromsome
#'
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
#' @param window.type "number" ,"kb" two types of call window default is number
#' @param window.size  window size dafalut is 15
#' @param low default is 0.2
#' @param high default is 0.8
#' @param fix a logical value ,default is T
#' @param fix.size define the neighbor size to fix the error
#'
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' @importFrom dplyr filter
#' @importFrom data.table fread
#' @importFrom utils write.table
#' @return files of relust
#' @export
#'
#'
batchCallGeno<-function(inputfile, datatype,outputfile,
                        screening=FALSE,
                        maf=NA,geno=NA,mind=NA,hwe=NA,
                        father = "1/1",mother = "0/0",
                        window.type="number",window.size=15,
                        low=0.2,high=0.8,
                        fix=TRUE,fix.size=10
                        ){

    snp_data<-readgeno(inputfile,datatype,outputfile,screening,maf,geno,mind,hwe,father,mother)

    snp_df<-tidyr::gather(snp_data,key = "line",value = "code",-c("CHR","POS","P1","P2"))

    chrom<-unique(snp_data$CHR)

    sample_name<-colnames(snp_data)[5:length(snp_data)]

    rm(snp_data)

    for (i in sample_name){
        for (chr in chrom){
            snp_ind_chr<-dplyr::filter(snp_df,.data$CHR==chr&.data$line==i)
            wind_geno<-callWindowGeno(snp_ind_chr,window.type,window.size,low,high)

        if (fix==TRUE)wind_geno$fix_code<-fixGeno(wind_geno$code,fix.size)

        wind_geno$CHR<-chr
        wind_geno$ind<-i

        ifelse(i==sample_name[1]&chr==chrom[1],write.table(wind_geno,file = paste(outputfile,"wind_geno","txt",sep = "."),sep='\t', row.names = FALSE, col.names =TRUE, quote =FALSE,append = FALSE),
               write.table(wind_geno,file = paste(outputfile,"wind_geno","txt",sep = "."),sep='\t', row.names = FALSE, col.names =FALSE, quote =FALSE,append = TRUE))
        }

    }
    rm(snp_df)
    wind_geno<-data.table::fread(paste(outputfile,"wind_geno","txt",sep = "."),header=T,sep="\t")
    ifelse(fix==TRUE,bin_geno<-tidyr::spread(wind_geno[,c(1,2,3,6,7,8)],key="ind",value="fix_code"),bin_geno<-tidyr::spread(wind_geno[,c(1,2,3,5,6,7)],key="ind",value="code"))
    rm(wind_geno)

    temp_bin_geno<-bin_geno[0,]
    for(chr in chrom){
        snp_chr<-dplyr::filter(bin_geno,.data$CHR==chr)
        for (i in 1:(nrow(snp_chr)-1)){
            if (i==i) temp_bin_geno<-rbind(temp_bin_geno,snp_chr[1,])
            temp_value<-snp_chr[i,6:ncol(snp_chr)]==snp_chr[i+1,6:ncol(snp_chr)]
            if (sum(temp_value)==length(sample_name)){
                value_end<-snp_chr[i+1,3]
                value_snp<-paste(chr,temp_bin_geno[nrow(temp_bin_geno),2],value_end,sep = "_")

                temp_bin_geno[nrow(temp_bin_geno),3]<-value_end
                temp_bin_geno[nrow(temp_bin_geno),4]<-value_snp

            }else{
                temp_bin_geno<-rbind(temp_bin_geno,snp_chr[i+1,])
            }

        }
    }
    write.table(temp_bin_geno,file = paste(outputfile,"bin_geno","txt",sep = "."),sep='\t', row.names = FALSE, col.names =TRUE, quote =FALSE,append = FALSE)

}
