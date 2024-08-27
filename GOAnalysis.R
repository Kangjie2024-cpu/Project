suppressMessages(library(dplyr))
suppressMessages(library(GO.db))
library(openxlsx)

######################FUNCTION######################

GO2TERM <- function(ontology){
  term <- toTable(GOTERM)
  x <- term[term$Ontology == ontology, ]
  return(x)
}



GOID2ALL <- function(annotfile, ontology){
  goterm <- GO2TERM(ontology)
  
  go2proid <- filter(annotfile, go_id %in% goterm$go_id)[c(2,1)]
  #上面这步会去掉注释到secondary的那部分。。。注意！
  
  secondary_goid <- filter(goterm[,-1], !is.na(Secondary))
  aa <- filter(annotfile, go_id %in% secondary_goid$Secondary)
  
  if (ontology == "BP") {
    tb1 <- toTable(GOBPANCESTOR)
  } else if (ontology == "MF") {
    tb1 <- toTable(GOMFANCESTOR)
  } else {
    tb1 <- toTable(GOCCANCESTOR)
  }
  colnames(tb1) <- c("go_id", "ancestor")
  
  go2ancestor <- filter(tb1, go_id %in% unique(go2proid$go_id) & ancestor != "all")
  
  tmp_go2proid <- inner_join(go2ancestor, go2proid, by = "go_id")
  names(tmp_go2proid)[2] <- "go_id"
  
  all_go2proid <- rbind(tmp_go2proid[,2:3], go2proid) %>% unique()
  
  return(all_go2proid)
}

#获取topGO实例
get_GOdata <- function(annotfile, ontology, prolist){
  goid2pro_df <- GOID2ALL(annotfile, ontology)
  pro2go <- split(goid2pro_df[,1], goid2pro_df[,2])
  pronames <- names(pro2go)
  pro <- factor(as.integer(pronames %in% prolist))
  names(pro) <- pronames
  GOdata <- new("topGOdata", allGenes = pro, ontology = ontology, 
                annot = annFUN.gene2GO, gene2GO = pro2go)
  return(GOdata)
}


get_golevel <- function(ontology, allGO){
  switch(
    ontology,
    BP = {
      level1 <- "GO:0008150"
      Children <- GOBPCHILDREN
    },
    MF = {
      level1 <- "GO:0003674"
      Children <- GOMFCHILDREN
    },
    CC = {
      level1 <- "GO:0005575"
      Children <- GOCCCHILDREN
    }
  )
  
  alreadylevel <- c()
  alreadylevel[level1] <- 1
  lev <- 1
  
  while(length(alreadylevel) < length(allGO)){
    index <- names(alreadylevel)[which(alreadylevel == lev)]
    index <- index[index %in% allGO]
    
    child_go <- mget(index, Children, ifnotfound = NA)
    child_go <- unlist(lapply(child_go, function(x){
      x[!x %in% names(alreadylevel)]
    }))
    child_go <- unique(child_go[child_go %in% allGO])
    
    lev <- lev + 1
    alreadylevel[child_go] <- lev
  }
  
  data <- data.frame(level=as.numeric(alreadylevel), go_id=names(alreadylevel), stringsAsFactors = F)
  return(data)
}



get_level2 <- function(go2pro, level, ontology, len){
  switch(
    ontology,
    BP = {
      type <- "biological_process"
    },
    MF = {
      type <- "molecular_function"
    },
    CC = {
      type <- "cellular_component"
    }
  )
  
  goterm <- unique(GO2TERM(ontology)[,2:3])
  
  lev2 <- inner_join(level, goterm, by = "go_id") %>% filter(level == 2)
  lev2$Type <- type
  
  lev2_len <- filter(go2pro, go_id %in% lev2$go_id) %>% group_by(go_id) %>% summarise(Seqs = n())
  lev2_len <- mutate(lev2_len, per = Seqs/len*100)
  
  lev2_data <- inner_join(lev2, lev2_len, by = "go_id")
  lev2_data[c(1,2)] <- lev2_data[c(2,1)]
  names(lev2_data)[c(1,2)] <- names(lev2_data)[c(2,1)]
  
  return(lev2_data)
}



get_gofile <- function(go2pro, level, ontology){
  switch(
    ontology,
    BP = {
      type <- "Biological Process"
    },
    MF = {
      type <- "Molecular Function"
    },
    CC = {
      type <- "Cellular Component"
    }
  )
  
  go2pro_list <- split(go2pro[,2], go2pro[,1])
  
  go2pro_len <- sapply(go2pro_list, length)
  go2pro_sequences <- sapply(go2pro_list, function(x){
    base::paste(x, collapse = ",")
  })
  
  goterm <- unique(GO2TERM(ontology)[,2:3])
  
  go_df <- inner_join(level, goterm, by = "go_id")
  go_df$Type <- type
  go_df$Seqs_Num <- go2pro_len[go_df$go_id]
  go_df$Sequences <- go2pro_sequences[go_df$go_id]
  
  return(go_df)
}



enricher <- function(BGgo, DEPgo, ontology, BGnum, DEPnum){
  
  switch(
    ontology,
    BP = {
      level1 <- "GO:0008150"
      category <- "P"
    },
    MF = {
      level1 <- "GO:0003674"
      category <- "F"
    },
    CC = {
      level1 <- "GO:0005575"
      category <- "C"
    }
  )
  
  DEPgo <- DEPgo[,c(-1,-4)] %>% filter(!go_id %in% level1)
  
  BGgo <- filter(BGgo, go_id %in% DEPgo$go_id)
  BGgo_list <- split(BGgo[,2], BGgo[,1])
  BG_go2pro_len <- sapply(BGgo_list, length)
  BG_go2pro_seq <- sapply(BGgo_list, function(x){
    base::paste(x, collapse = ",")
  })
  
  BG_df <- data.frame(go_id=names(BGgo_list), Ref=BG_go2pro_len, RefSeqs=BG_go2pro_seq, stringsAsFactors = F, row.names = NULL)
  
  enrich_data <- inner_join(DEPgo, BG_df, by = "go_id")
  enrich_data <- mutate(enrich_data, Category=category, TestAll=DEPnum, RefAll=BGnum, Test_per=Seqs_Num/DEPnum, Ref_per=Ref/BGnum)
  names(enrich_data)[3] <- "Test"
  enrich_data$Over_Under <- ifelse(enrich_data$Test_per >= enrich_data$Ref_per, "over", "under")
  enrich_data <- enrich_data[c(1,2,7,3,5,8,9,10,11,12,4,6)]
  
  p.value <- phyper(enrich_data$Test-1,enrich_data$Ref, enrich_data$RefAll-enrich_data$Ref,enrich_data$TestAll,lower.tail = FALSE)
  #  FDR <- p.adjust(p.value,method="fdr",n=length(p.value))
  #  richFactor <- enrich_data$Test/enrich_data$Ref
  
  #  enrich_data <- base::cbind(enrich_data, p.value, FDR, richFactor)
  enrich_data <- base::cbind(enrich_data, p.value)
  
  return(enrich_data)
}

###################DATA INPUT###################

args <- commandArgs(T)
annotfile <- as.character(args[1])
DEPfile <- as.character(args[2])

annot <- read.table(file = annotfile, sep = "\t", stringsAsFactors = F)
annot <- unique(annot)
DEP <- read.table(file = DEPfile, sep = "\t",stringsAsFactors = F)

names(annot) <- c("pro_id", "go_id")
DEP_annot <- filter(annot, pro_id %in% DEP$V1)


###################CREATE file for annotationfile.pl###################

ontology <- base::rbind(GO2TERM("BP"), GO2TERM("MF"), GO2TERM("CC"))[,-1] %>% filter(go_id %in% annot$go_id) 

proid_go_group <- left_join(annot, unique(ontology[,1:3]), by = "go_id")

outputfile <-  paste0(strsplit(annotfile, "_")[[1]][1], ".annot")

write.table(proid_go_group, file = outputfile, sep = "\t",row.names = F, col.names = F, quote = F)


###################DEP GO classification annotation###################

bp2proid <- GOID2ALL(annotfile = DEP_annot, ontology = "BP")
mf2proid <- GOID2ALL(annotfile = DEP_annot, ontology = "MF")
cc2proid <- GOID2ALL(annotfile = DEP_annot, ontology = "CC")

bplevel <- get_golevel(allGO = unique(bp2proid$go_id), ontology = "BP")
mflevel <- get_golevel(allGO = unique(mf2proid$go_id), ontology = "MF")
cclevel <- get_golevel(allGO = unique(cc2proid$go_id), ontology = "CC")

bplevel2 <- get_level2(go2pro = bp2proid, level = bplevel, ontology = "BP", len = length(DEP$V1))
mflevel2 <- get_level2(go2pro = mf2proid, level = mflevel, ontology = "MF", len = length(DEP$V1))
cclevel2 <- get_level2(go2pro = cc2proid, level = cclevel, ontology = "CC", len = length(DEP$V1))

golevel2_df <- bind_rows(bplevel2, mflevel2, cclevel2)

write.table(golevel2_df, file = "GOLevel2.txt", sep = "\t",row.names = F, quote = F)


###################BP/MF/CC ANALYSIS###################

DEP_bp <- get_gofile(go2pro = bp2proid, level = bplevel, ontology = "BP")
DEP_mf <- get_gofile(go2pro = mf2proid, level = mflevel, ontology = "MF")
DEP_cc <- get_gofile(go2pro = cc2proid, level = cclevel, ontology = "CC")

bp_classification <- DEP_bp
names(bp_classification)[1:2] <- c("Level", "GO_ID")
mf_classification <- DEP_mf
names(mf_classification)[1:2] <- c("Level", "GO_ID")
cc_classification <- DEP_cc
names(cc_classification)[1:2] <- c("Level", "GO_ID")

write.table(bp_classification, file = "BP.txt", sep = "\t", row.names = F, quote = F)
write.table(mf_classification, file = "MF.txt", sep = "\t", row.names = F, quote = F)
write.table(cc_classification, file = "CC.txt", sep = "\t", row.names = F, quote = F)

#wb <- createWorkbook()
#modifyBaseFont(wb, fontSize = 11, fontName = "Arial")
#
#addWorksheet(wb, sheetName = "BP")
#addWorksheet(wb, sheetName = "MF")
#addWorksheet(wb, sheetName = "CC")
#
#header_style <- createStyle(textDecoration = "bold", halign = "left")
#
#addStyle(wb, sheet = "BP", rows = 1, cols = 1:ncol(bp_classification), style = header_style)
#addStyle(wb, sheet = "MF", rows = 1, cols = 1:ncol(mf_classification), style = header_style)
#addStyle(wb, sheet = "CC", rows = 1, cols = 1:ncol(cc_classification), style = header_style)
#
#writeData(wb, sheet = "BP", x = bp_classification)
#writeData(wb, sheet = "MF", x = mf_classification)
#writeData(wb, sheet = "CC", x = cc_classification)
#
#setColWidths(wb, sheet = "BP", cols = 2:4, widths = c(12,30,18))
#setColWidths(wb, sheet = "MF", cols = 2:4, widths = c(12,30,18))
#setColWidths(wb, sheet = "CC", cols = 2:4, widths = c(12,30,18))

#saveWorkbook(wb, "GO.xlsx", overwrite = T)



###################Enrich ANALYSIS###################

BGfile <- as.character(args[3])

if (is.na(BGfile)){
  #  saveWorkbook(wb, "GO.xlsx", overwrite = T)
  cat("No enrich analysis, exit!\n")
  q()
}

BG <- read.table(file = BGfile, sep = "\t", stringsAsFactors = F)

BG_annot <- filter(annot, pro_id %in% BG$V1)

BG_bp2proid <- GOID2ALL(annotfile = BG_annot, ontology = "BP")
BG_mf2proid <- GOID2ALL(annotfile = BG_annot, ontology = "MF")
BG_cc2proid <- GOID2ALL(annotfile = BG_annot, ontology = "CC")

enrich_bp <- enricher(BGgo = BG_bp2proid, DEPgo = DEP_bp, ontology = "BP", BGnum = as.numeric(length(BG$V1)), DEPnum = as.numeric(length(DEP$V1)))
enrich_mf <- enricher(BGgo = BG_mf2proid, DEPgo = DEP_mf, ontology = "MF", BGnum = as.numeric(length(BG$V1)), DEPnum = as.numeric(length(DEP$V1)))
enrich_cc <- enricher(BGgo = BG_cc2proid, DEPgo = DEP_cc, ontology = "CC", BGnum = as.numeric(length(BG$V1)), DEPnum = as.numeric(length(DEP$V1)))

enrich_go <- base::rbind(enrich_bp, enrich_mf, enrich_cc)
names(enrich_go)[c(1,11)] <- c("GO_id", "TestSeqs")

enrich_go <- arrange(enrich_go, p.value) %>% arrange(Over_Under)

enrich_go$FDR <- p.adjust(enrich_go$p.value,method="fdr",n=length(enrich_go$p.value))
enrich_go$richFactor <- enrich_go$Test/enrich_go$Ref

write.table(enrich_go, file = "Enrichment.txt", sep = "\t", row.names = F, quote = F)

#addWorksheet(wb, sheetName = "Enrichment")
#addStyle(wb, sheet = "Enrichment", rows = 1, cols = 1:ncol(enrich_go), style = header_style)
#
#cell_style <- createStyle(fgFill = "yellow")
#np <- filter(enrich_go, p.value <= 0.05 & Over_Under == "over") %>% nrow()
#addStyle(wb, sheet = "Enrichment", rows = 2:np, cols = 2, style = cell_style)
#addStyle(wb, sheet = "Enrichment", rows = 2:np, cols = 10, style = cell_style)
#addStyle(wb, sheet = "Enrichment", rows = 2:np, cols = 13, style = cell_style)
#
#writeData(wb, sheet = "Enrichment", x = enrich_go)
#setColWidths(wb, sheet = "Enrichment", cols = c(1,2), widths = c(12,30))
#
#saveWorkbook(wb, "GO.xlsx", overwrite = T)

