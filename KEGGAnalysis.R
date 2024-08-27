suppressMessages(library(RCurl))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
library(httr)

args <- commandArgs(T)
path <- as.character(args[1])
inputfile <- as.character(args[2])
species <- as.character(args[3])
# inputfile <- "all_list.txt"
# species <- "mtu"
#
setwd(path)

data <- read.table(file = inputfile, sep = "\t", header = F, stringsAsFactors = F)
names(data) <- "uniprotid"

#################################################
#将uniprot id 通过Uniprot api转化为kegg id
#################################################

uniprot2keggid <- function(id){
  e_mail <- 'kgu@sibs.ac.cn' #Enter your email address here to help hosts troubleshoot issues
  us_er <- paste0('R ', e_mail)
  acc1 <- 'ACC' #Enter the starting accession id (e.g. 'ACC' for Uniprot accession)
  acc2 <- 'KEGG_ID' #Enter the target accession id (e.g. 'HPA_ID' for Human Protein Atlas accession)
  fmt <- 'tab' #Enter format (e.g. 'tab' for TSV)
  qry <- "A0A0B4J2F2" #Enter query term
  qry <- paste(id, collapse = ",")
  user_agent(us_er)
  r <- POST('https://www.uniprot.org/id-mapping/', body = list(from= acc1, to = acc2, format = fmt, query = qry), encode = "form")
  
  res <- content(r, type = "text")
  res <- read_tsv(res)
  names(res) <- c("uniprotid", "keggid")
  res$keggid <- as.character(res$keggid)
  
  return(res)
}


#################################################
#将uniprot id 通过Uniprot api转化为entrez id
#################################################
uniprot2entrez <- function(id){
  e_mail <- 'kgu@sibs.ac.cn' #Enter your email address here to help hosts troubleshoot issues
  us_er <- paste0('R ', e_mail)
  acc1 <- 'ACC' #Enter the starting accession id (e.g. 'ACC' for Uniprot accession)
  acc2 <- 'P_ENTREZGENEID' #Enter the target accession id (e.g. 'HPA_ID' for Human Protein Atlas accession)
  fmt <- 'tab' #Enter format (e.g. 'tab' for TSV)
  qry <- "A0A0B4J2F2" #Enter query term
  qry <- paste(id, collapse = ",")
  user_agent(us_er)
  r <- POST('https://www.uniprot.org/uploadlists/', body = list(from= acc1, to = acc2, format = fmt, query = qry), encode = "form")
  
  res <- content(r, type = "text")
  res <- read_tsv(res)
  names(res) <- c("uniprotid", "geneid")
  res$geneid <- as.character(res$geneid)
  
  return(res)
}

################################################
#通过kegg api 将uniprot id 转化为kegg id
################################################
uniprot2kegg <- function(id, species){
  url <- paste0("http://rest.kegg.jp/conv/", species, "/uniprot")
  url.exists(url)
  d <- debugGatherer()
  pro2keggid <- read.table(url)
  #pro2keggid <- read_tsv(pro2keggid, col_names = F)
  names(pro2keggid) <- c("uniprotid", "keggid")

  pro2keggid$uniprotid <- unlist(lapply(as.character(pro2keggid$uniprotid), function(x){
    strsplit(x, ":")[[1]][2]
  }))

  res <- pro2keggid[!duplicated(pro2keggid$uniprotid),]
  res <- dplyr::filter(res, uniprotid %in% id)
  
  return(res)

}

################################################
#通过kegg api 将entrez id 转化为kegg id
################################################
entrez2kegg <- function(pro2gene, species){
  url <- paste0("http://rest.kegg.jp/conv/", species, "/ncbi-geneid")
  url.exists(url)
   gene2keggid<- read.table(url)
   gene2keggid$V1<-gsub("ncbi-geneid:","",gene2keggid$V1)
  #d <- debugGatherer()
  #gene2keggid <- getURL(url=url, debugfunction = d$update, verbose = TRUE)
  #gene2keggid <- read_tsv(gene2keggid, col_names = F)
  names(gene2keggid) <- c("geneid", "keggid")
  
  #gene2keggid$geneid <- unlist(lapply(gene2keggid$geneid, function(x){
  #  strsplit(x, ":")[[1]][2]
  #}))
  
  gene2keggid <- gene2keggid[!duplicated(gene2keggid$geneid),]
  
  result <- left_join(pro2gene, gene2keggid, by = "geneid")
  result <- dplyr::filter(result, !is.na(keggid))
  
  return(result)
}

##################################################################################

#res <- uniprot2keggid(id = data$uniprotid) %>%
#  filter(grepl(pattern = paste0(species, ":"), keggid)) %>%
#  distinct(uniprotid, .keep_all = TRUE)
#
#nomapping_id <- data$uniprotid[!data$uniprotid %in% res$uniprotid]

#if (length(nomapping_id) > 0){
#  res2 <- uniprot2entrez(nomapping_id)
#  if (nrow(res2)== 0){
    pro2kegg <- uniprot2kegg(id = data$uniprotid, species = species)
#    df_result <- bind_rows(res, pro2kegg)
#  }else{
#    pro2kegg <- entrez2kegg(pro2gene = res2, species = species)
#    df_result <- bind_rows(res, pro2kegg)
#    df_result <- df_result[,1:2]
#  }
#}else{
  df_result <- pro2kegg
#}

df_result2 <- left_join(data, df_result, by = "uniprotid")

#######################################################################

if (sum(is.na(df_result$keggid)) > length(df_result$keggid)/10*7.5){
  stop(paste0("\n", paste0(rep("#", 100), collapse = ""), "\n#\n#\n#The proportion of mapping to keggid is below 25%, maybe you should try to use KAAS !!!\n#\n#\n", paste0(rep("#", 100), collapse = ""), "\n\n"))
}

write.table(df_result2, file = "all.ko", sep = "\t", col.names = F, row.names = F, quote = F, na = "")
