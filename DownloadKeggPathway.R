suppressMessages(library(dplyr))
suppressMessages(library(RCurl))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(plyr))

args <- commandArgs(T)

species <- as.character(args[1])
path <- as.character(args[2])
setwd(path)
pwd <- getwd()

ko2name_url <- paste0("http://rest.kegg.jp/list/", species)
url.exists(ko2name_url)
ko2name <- read.table(ko2name_url, check.names = F, sep = "\t", comment.char = "%",quote = "")
#ko2name <- read_tsv(ko2name, col_names = F)
returnname <- function(x){
  if(grepl(";",x[4])){
    a <- unlist(strsplit(x[4],";"))[1]
    return(unlist(strsplit(a,","))[1])
  }else{
    return("")
  }
}
returndef <- function(x){
  if(grepl(";",x[4])){
    return(unlist(strsplit(x[4],";"))[2])
  }else{
    return(x[4])
  }
}
n <- apply(ko2name,1,returnname)
d <- apply(ko2name,1,returndef)
d <- str_sub (d,2)
ko2name <- data.frame(ko2name[,1],n,d,stringsAsFactors = F, check.names = F)
names(ko2name) <- c("keggid","Name","Definition")


pathway2name_url <- paste0("http://rest.kegg.jp/list/pathway/", species)
url.exists(pathway2name_url)

pathway2name <-  read.table(pathway2name_url, check.names = F, sep = "\t", comment.char = "%",quote = "")
#pathway2name <- read_tsv(pathway2name, col_names = F)
names(pathway2name) <- c("pathwayid", "pathwayname")

pathway2name$pathwayid <- str_replace(pathway2name$pathwayid, "path:", "")
if(species!="ko"){
  pathway2name$pathwayname <- sub("(.*)-.*", "\\1", pathway2name$pathwayname)
}

keggid2pathway_url <- paste0("http://rest.kegg.jp/link/pathway/", species)
url.exists(keggid2pathway_url)
keggid2pathway <- read.table(keggid2pathway_url, check.names = F, sep = "\t", comment.char = "%",quote = "")
#keggid2pathway <- read_tsv(keggid2pathway, col_names = F)
names(keggid2pathway) <- c("keggid", "pathwayid")

keggid2pathway$pathwayid <- str_replace(keggid2pathway$pathwayid, "path:", "")
if(species == "ko"){
  keggid2pathway <- keggid2pathway[grepl("ko",keggid2pathway$pathwayid,fixed = T),]
}

resultmed <- merge(keggid2pathway,pathway2name,,by="pathwayid",all.x = TRUE)
resultpathway <- merge(ko2name,resultmed,by="keggid",all.y = TRUE)
names(resultpathway)[1] <- "Kegg_id"
if(species == "ko"){
  resultpathway$Kegg_id <- sub("ko:","",resultpathway$Kegg_id,fixed=T)
}
resultpathway <- filter(resultpathway, !(grepl("(^[a-zA-Z]+01100$)",pathwayid)|grepl("(^[a-zA-Z]+01110$)",pathwayid)|grepl("(^[a-zA-Z]+01120$)",pathwayid)|grepl("(^[a-zA-Z]+01130$)",pathwayid)|grepl("(^[a-zA-Z]+01200$)",pathwayid)|grepl("(^[a-zA-Z]+01210$)",pathwayid)|grepl("(^[a-zA-Z]+01212$)",pathwayid)|grepl("(^[a-zA-Z]+01230$)",pathwayid)|grepl("(^[a-zA-Z]+01220$)",pathwayid)))
write.table(resultpathway,paste0(species,".pathway"),sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)

