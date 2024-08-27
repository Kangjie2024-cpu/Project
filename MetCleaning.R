
library(stringi)
library(e1071)
library(openxlsx)
library(reshape2)
library(dplyr)
suppressMessages(library(MetAnalyzer))
library(MetCleaning)
library(stringr)
options(stringsAsFactors = FALSE)

arg = commandArgs(T)

# --- set input ---- #
path = arg[1]
setwd(path)

print("xxxxxxxxxxxxxxxxx metcleaning!")

data_func <- function(ION = ION){
          ion <- tolower(ION)
          #setwd(paste0(path,"/",ION,"/report"))
          #setwd(paste0(path,"/",ION,"/"))
          setwd(path)
          if (ION == "POS"){
            polar <- "positive"
            file_idresults = "pos-del-iso.csv"
          }else{
            polar <- "negative"
            file_idresults = "neg-del-iso.csv"
          }
          #file_idresults <- paste0(path,"/",ION,"/idresults.csv")
         print("a")
          idresults <- read.csv(file_idresults,sep=",",header=TRUE,encoding='UTF-8',check.names=F,stringsAsFactors = FALSE)

          #colnames(idresults) <- sub("mzmed","mz",colnames(idresults ))
          colnames(idresults) = lapply(colnames(idresults),function(x){ifelse(x=="mzmed","mz",x)}) %>% unlist()
          #colnames(idresults) <- sub("rtmed","rt",colnames(idresults ))
          colnames(idresults) = lapply(colnames(idresults),function(x){ifelse(x=="rtmed","rt",x)}) %>% unlist()
          colnames(idresults) <- str_replace_all(colnames(idresults),".HILIC.*","")  #HILIC常用情况
          colnames(idresults) <- str_replace_all(colnames(idresults),paste0("_",ION,".*"),"")  #_NEG常用情况
          colnames(idresults) <- str_replace_all(colnames(idresults),paste0("_",ion,".*"),"")  #_neg常用情况
          colnames(idresults)[1:(ncol(idresults)-15)] <- str_replace_all(colnames(idresults)[1:(ncol(idresults)-15)],"\\.","-")  #HILIC常用情况
         print("b")
        if(file.exists(paste0(path,"/newname.xlsx"))){
          newname <- read.xlsx(paste0(path,"/newname.xlsx"), rowNames = F, check.names = F)
          newname <- newname[,1:2]
          Name <- melt(newname, id = colnames(newname)[1])
          inputname <- dplyr::select(Name, colnames(Name)[1], value) %>% filter(value != "")
          names(inputname) <- c("name", "group")
          inputname <- inputname[order(inputname$group),] #01.07  #保证同一组在表单排序在一起
          }else{
            samplename <- colnames(idresults)[grep("-\\d+$", colnames(idresults))]
            inputname <- data.frame(name = samplename, group = str_remove(samplename, "-\\d+$"))
          }
          if ("QC" %in% inputname$group){
            inputname <- rbind(filter(inputname, group != "QC"), filter(inputname, group == "QC"))
          }
          taball <- unique(inputname$group)
            print("c")
        #MetCleaning
        shape_func<- function(data,ION){
            tabnoQC <- taball[taball != "QC"]
            if (nrow(data) <= 20000) {
              threshold <- 0.3
            }else if (nrow(data) > 20000 & nrow(data) <= 30000){
              threshold <- 0.5
            }else if (nrow(data) > 30000 & nrow(data) <= 50000){
              threshold <- 0.66
            }else{
              threshold <- 0.8
            }
            k <- 1
            tmp <- data
            sample_num <-length(inputname[inputname$group != "QC" ,]$name)*threshold
            for (x in 1:nrow(data)){
              m <- 0
              sample_2 <- data[x,colnames(data) %in% tabnoQC]
              if ( data$hits.forward_zhumetlib[x] == "" & data$hits.forward_metlinlib[x] == "" & data$Ref_Compound[x] == "" & sum(sample_2) <= sample_num ){
                m <- m+1
              } 
              if (m==0){
                tmp[k,] <- data[x,]
                k <- k+1
              }
            }
            data <- tmp[1:(k-1),]
            return(data)
          }
        print("d")


        amide_ms1  <- shape_func(idresults,ION = ION)
        # if( nrow(amide_ms1) > 15000) {
        #   amide_ms1  <- shape_func(idresults,ION = ION,threshold = (2/3))
        # }


        amide_ms1  <- amide_ms1[,colnames(amide_ms1) != "X"]
        amide_path <- paste0(path,"/amide ",ion," identification")
        dir.create(amide_path)
        setwd(amide_path)



        write.table(amide_ms1,paste0("amide ",ion," ms1.csv"),sep = ",",na = "",row.names = F)
        file.copy(paste0(path,"/sample.information.csv"),amide_path)
        MetCleaning(data=paste0("amide ",ion," ms1.csv"), polarity = polar, mz.tolerance=30,rt.tolerance=60, qc.outlier.filter=FALSE, subject.outlier.filter=FALSE,met.plot = F)

        #data_filter for neg-dele-iso.csv, pos-dele-iso.csv
        filter_func <- function(data,ION){
                data <- data[,colnames(data) != "polarity"]
                #colnames(data) <- sub("mz","mzmed",colnames(data))
                #colnames(data) <- sub("rt","rtmed",colnames(data))
                colnames(data) = lapply(colnames(data) ,function(x){ifelse(x=="mz","mzmed",x)}) %>% unlist()
                colnames(data) = lapply(colnames(data) ,function(x){ifelse(x=="rt","rtmed",x)}) %>% unlist()
                data$name <- str_replace_all(data$name,paste0("_",ION),"")
                data$isotopes <- str_replace_all(data$isotopes,'\\[\\d+\\]',"")   #isotopes:[M]—&""
                data <- subset(data,data$isotopes == "[M]+"|data$isotopes == ""|data$isotopes == "[M]-")
                
                # 空值处理
                for (q in 1:length(inputname$name)) {
                    data[,inputname$name[q]] <- ifelse(data[,inputname$name[q]]==0,NA,data[,inputname$name[q]])
                }

                write.table(data,paste0(path,"/",ION,"_before_impute.csv"),sep = ",",na = "",row.names = F)
                # 样本种如果有一个空值，就去掉这一行
                if(FALSE){
                    k <- 1
                    tmp <- data
                    for (x in 1:nrow(data)){
                        m <- 0
                        sample_1 <- data[x,colnames(data) %in% inputname$name]
                        if (any(is.na(sample_1))){
                            m <- m+1
                        }
                        if (m==0){
                            tmp[k,] <- data[x,]
                            k <- k+1
                        }
                    }
                    data <- tmp[1:(k-1),]
                }
                # 用knn方法对数据进行NA值填充
                
                if(FALSE){
                    library(impute)
                    #knndf = data.frame()
                    knnlist = list()
                    #data=preProcess(data,method="knnImpute",k=5)
                    samplegroup = unique(inputname$group[which(inputname$group!="QC")])
                    samplename = c()
                    for(eachgroup in samplegroup){
                        eachgroupsample = data[,colnames(data) %in% inputname[inputname$group==eachgroup,]$name]
                        knngroupample = impute.knn(as.matrix(eachgroupsample))$data
                        #knndf=cbind(knndf,as.data.frame(knngroupample))
                        knnlist[[eachgroup]] = knngroupample
                        samplename = c(samplename,colnames(knngroupample))
                    }
                    knndf = data.frame(knnlist)
                    colnames(knndf) = samplename
                    for(eachcol in colnames(knndf)){
                        # names = strsplit(eachcol,"\\.")[[1]]
                        # groupname = names[1]
                        # samname = paste(names[2:length(names)],collapse="-")
                        # print(samname)
                        data[,eachcol] = knndf[,eachcol]
                        }

                }

                if(TRUE){
                  library(impute)
                  # 数据进行KNN填充        
                  samplegroup = unique(inputname$group)
                  eachgroupsample = data[,colnames(data) %in% inputname[inputname$group %in% samplegroup,]$name]
                  knngroupample = t(impute.knn(as.matrix(t(eachgroupsample)),rowmax = 0.9,colmax = 0.99)$data)
                  for(eachcol in colnames(knngroupample)){
                      data[,eachcol] = knngroupample[,eachcol]
                      }
                }



                QC <-data[,colnames(data) %in% inputname[inputname$group== "QC",]$name]
                QC[QC =="N/A"] <- NA
                QC[QC =="NA"] <- NA
                # 计算RSD(RSD=STDEV/Average,stdev=sd,average=mean)
                RSD <- function(x){
                  m <- mean(x,na.rm=TRUE)
                  s <- sd(x,na.rm=TRUE)
                  RSD <- s/m
                  return(RSD=RSD)
                }
                QC_RSD <- apply(QC,1,RSD)
                data$RSD<- QC_RSD
                loc <- grep("isotopes",colnames(data))
                data <- data.frame(data[,1:(loc-1)],data[,colnames(data) %in% inputname$name],RSD =data[,ncol(data)],data[,loc:(loc+14)],check.names = F)
                data <- data %>% dplyr::arrange(RSD,desc(isotopes))
                write.table(data,paste0(path,"/",ion,"-data_before_FilterByRsd.csv"),sep = ",",na = "",row.names = F)
                FilterByRsd = function(data){
                    # nnum = nrow(data)
                    # less05 = sum(data$RSD>0.5)
                    # lassnum = nnum - less05
                    # if(lassnum<=5000){
                    #     data = data
                    # }else{
                    #     data = data[data$RSD<0.5,]
                    # }
                    data = data[data$RSD<0.5,]
                    return(data)
                }
                # 根据nrow进行Rsd的过滤
                data = FilterByRsd(data)


                QC_stat1<- data[(data$isotopes == "[M]+"|data$isotopes == "[M]-") & data$RSD < 0.5,colnames(data) %in% inputname[inputname$group== "QC",]$name]
                QC_stat1_sum <- apply(QC_stat1,2,sum,na.rm =T)
                QC_stat<- data[data$isotopes == "[M]+"|data$isotopes == "[M]-",colnames(data) %in% inputname[inputname$group== "QC",]$name]
                QC_stat_sum <- apply(QC_stat,2,sum,na.rm =T)
                percent30 <- max(QC_stat1_sum/QC_stat_sum)
                feature_numb <- nrow(data)
                meta_numb <- nrow(data[data$hits.forward_zhumetlib != "",])
                record <- c()
                record <- c("ION",ION)
                record <- rbind(record,c("percent30", percent30))
                record <- rbind(record,c("feature_numb",feature_numb))
                record <- rbind(record,c("meta_numb",meta_numb))
                report <- paste0(path,"/POS_report_info.txt")

                return(data)
            }
        file_data_after_pre <- paste0(amide_path,"/data_after_pre.csv")
        data_after_pre <- read.csv(file_data_after_pre,sep=",",header=TRUE,encoding='UTF-8',check.names=F,stringsAsFactors = FALSE)
        colnames(data_after_pre)[grep("^QC\\d+$",colnames(data_after_pre))] <- str_replace_all(colnames(data_after_pre)[grep("^QC\\d+$",colnames(data_after_pre))] ,"^QC","QC-")
        dele_iso <- filter_func(data_after_pre,ION = ION)
        write.table(dele_iso,paste0(path,"/",ion,"-dele-iso.csv"),sep = ",",na = "",row.names = F)
}

data_func(ION ="POS")
data_func(ION ="NEG")
setwd(path)
unlink("amide pos identification",recursive = T)
unlink("amide neg identification",recursive = T)
