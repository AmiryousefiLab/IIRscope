# TODO: check different trnfixer/trnFixer versions and see which one performs better.
#       and join them.


#' Get gene names
#' 
#' @param gb gb file to get the names from.
#' @param type either gene, tRNA, mRNA or rRNA.
#' 
#' @return vector with the gene names.
gene.name<-function(gb, type){
  if(type=="gene"){
    t <- gb[grep("  gene  ", gb)+1]
    if(length(t) > 0){
      for (i in 1:length(t)){
        crude<-strsplit(t[i], " ")[[1]]
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  gene  ", gb)+2][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  gene  ", gb)+3][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  gene  ", gb)+4][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        crude.name<-crude[which(crude!=rep("", length(crude)))]
        t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
      }
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("gene ", gb)+1][na]), ""), ""))
      }
    }
  }
  else if(type=="tRNA"){
    t <- gb[grep(" tRNA  ", gb)+1]
    if(length(t) > 0){
      for (i in 1:length(t)){
        crude<-strsplit(t[i], " ")[[1]]
        if(length(grep("=", crude))==0){
          t[i] <- gb[grep("tRNA  ", gb)+2][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  tRNA  ", gb)+3][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  tRNA  ", gb)+4][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        crude.name<-crude[which(crude!=rep("", length(crude)))]
        t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
      }
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("tRNA ", gb)+1]), ""), ""), "\n")
      }
    }
  }
  else if(type=="rRNA"){
    t <- gb[grep("rRNA  ", gb)+1]
    if(length(t) > 0){
      for (i in 1:length(t)){
        crude<-strsplit(t[i], " ")[[1]]
        if(length(grep("=", crude))==0){
          t[i] <- gb[grep("rRNA  ", gb)+2][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  rRNA  ", gb)+3][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  rRNA  ", gb)+4][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        crude.name<-crude[which(crude!=rep("", length(crude)))]
        t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
      }
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("gene ", gb)+1][na]), ""), ""))
      }
    }
  }
  else if(type=="mRNA"){
    t <- gb[grep("mRNA  ", gb)+1]
    if(length(t) > 0){
      for (i in 1:length(t)){
        crude<-strsplit(t[i], " ")[[1]]
        if(length(grep("=", crude))==0){
          t[i] <- gb[grep("mRNA  ", gb)+2][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  mRNA  ", gb)+3][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  mRNA  ", gb)+4][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        crude.name<-crude[which(crude!=rep("", length(crude)))]
        t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
      }
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("gene ", gb)+1][na]), ""), ""))
      }
    }
  }
  else {
    stop("The type should be defined as either gene or mRNA or rRNA or tRNA")
  }
  t<-t[!is.na(t)]
  return(t)
}

#' Get gene coordinates
#' 
#' @param gb gb file to get the genes from.
#' 
#' @return he crude table extracted from gb file for genes.
#' 
gb.gene.cor<- function(gb){
  gene<- gb[grep("  gene  ", gb)]
  trna<- gb[grep("  tRNA ", gb)]
  rrna<- gb[grep("  rRNA ", gb)]
  t<-c(gene, trna, rrna)
  m<- matrix(0, length(t), 2)
  for (i in 1:length(t)){
    if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]==","){#check the last character if its ',' and add the next one from gb file to that
      #previous code--> if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]=="," && strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])-1]==")")
      s<-t[i]
      t[i]<-paste(t[i], gsub(" ", "", gb[which(t[i]==gb)+1]), "")
    }
    t[i]<- gsub(", ", ",", t[i])
    t[i]<- gsub(") ", ")", t[i])
    if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]==","){#check again if the last characted if its ',' and add the next one from gb file to that
      #previous code--> if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]=="," && strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])-1]==")")
      t[i]<-paste(t[i], gsub(" ", "", gb[which(s==gb)+2]), "")
    }
    t[i]<- gsub(", ", ",", t[i])
    t[i]<- gsub(") ", ")", t[i])
    m[i,2]<-strsplit(t[i], " ")[[1]][length(strsplit(t[i], " ")[[1]])]
  }
  names<-  c(gene.name(gb, "gene"), gene.name(gb, "tRNA"), gene.name(gb, "rRNA") )
  if(length(t)!=length(names)){
    stop("Error: Check the gene names of the gb file")
  }
  m[,1]<-names
  return(m)
}

#' Get genes coordinates
#' 
#' @param gb gb file to get the genes from.
#' 
#' @return hpolished table of genes and their coordinates, with the gene names,
#' start, end of first part and the second parts (first to fifth columns).
#' The orders are reflected with the swapping of values.
#' 
gene.cordinates<- function(gb){
  gb.gene.cor.out<-gb.gene.cor(gb)
  l<- length(gb.gene.cor.out[,1])
  m<- matrix(0, l, 5)
  m[,1]<- gb.gene.cor.out[,1]
  for (i in 1:l){
    if(length(strsplit(gb.gene.cor.out[i,2], ",")[[1]])==2){
      m[i,2]<-sub("join\\(", "", strsplit(gb.gene.cor.out[i,2], ",")[[1]][1])
      m[i,2]<-sub("order\\(", "", m[i,2])
      m[i,4]<-sub(")", "", strsplit(gb.gene.cor.out[i,2], ",")[[1]][2])
      cord<-m[i,2]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]]
      }
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
      cord<-m[i,4]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(4,5)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        m[i,c(4,5)]<- strsplit(cord, "\\..")[[1]]
      }
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
    }
    else {
      m[i,2]<-gb.gene.cor.out[i,2]
      cord<-m[i,2]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        # print(cord) TODO fix problem
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]]
      }
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
    }
  }
  m<- gsub(")", "", m)
  m<- gsub("\\(", "", m)
  m<- gsub("<", "", m)
  m<- gsub(">", "", m)
  return(m)
}


#' Generate gene table from parsed gb file
#'
#' @param gb A list containing parsed GB file information. It is generated by
#' function \code{\link[genbankr]{parseGenBank}}
#' @param genome a DNAstring object. It contains the genome sequence.
#'
#' @return a data frame. It contains information of genes.
#' @importFrom magrittr %>%
#' @import dplyr
#' 
geneTableParsed <- function(gb, genome){
  feature <- vector(mode = "list", length = length(gb$FEATURES))
  names(feature) <- names(gb$FEATURES)
  for (i in names(feature)){
    feature[[i]] <- as.data.frame(gb$FEATURES[[i]])
  }
  type <- lapply(feature, "[[", "type")
  cols <- c("start", "end", "strand",
            "type", "gene", "pseudo", "product")
  info <- NULL
  for(i in 1:length(feature)){
    tmp <- feature[[i]]
    miscol <- cols[!cols %in% colnames(tmp)]
    df <- data.frame(matrix(rep(NA, length(miscol) * nrow(tmp)),
                            nrow = nrow(tmp), ncol = length(miscol)),
                     stringsAsFactors = FALSE)
    colnames(df) <- miscol
    tmp <- cbind.data.frame(tmp, df)
    tmp <- tmp[, which(colnames(tmp) %in% cols)]
    info <- rbind.data.frame(info, tmp)
  }
  
  
  # gene
  info$gene[is.na(info$gene)] <- info$product[is.na(info$gene)]
  info$pseudo[is.na(info$pseudo)] <- FALSE
  
  info$gene[grepl(".*([0-9\\.]+)S.*", info$gene)] <-
    rrnFixer(info$gene[grepl(".*([0-9\\.]+)S.*", info$gene)])
  info$gene[grepl("^trn.*", info$gene, ignore.case=TRUE)] <-
    trnFixer(info$gene[grepl("^trn.*", info$gene, ignore.case=TRUE)])
  
  
  gene_table <- info %>%
    dplyr::filter(type %in% c("gene", "tRNA", "rRNA")) %>%
    dplyr::select(start, end, strand, gene, pseudo) %>%
    stats::na.omit() %>%
    unique() %>%
    dplyr::mutate(chr = rep("chr1", n()))
  # for (i in 1:nrow(gene_table)){
  #   if (gene_table$strand[i] == "-" ){
  #     tmp <- gene_table$start[i]
  #     gene_table$start[i] <- gene_table$end[i]
  #     gene_table$end[i] <- tmp
  #   }
  # }
  # gene_table <- select(gene_table, chr, start, end, gene)
  
  # remove duplicated tRNA and rRNA
  gene_table <- gene_table[order(gene_table[, "start"], -gene_table[, "end"]), ]
  gene_table <- gene_table[!duplicated(gene_table[, c("start", "strand", "gene")]),]
  gene_table <- gene_table[!duplicated(gene_table[, c("end", "strand", "gene")]),]
  
  # codon usage
  cds <- info[which(info$type == "CDS"),]
  cds_cu <- codonUsage(cds, genome)
  gene_table <- dplyr::left_join(gene_table, cds_cu, by = c("gene", "strand",
                                                            "start"))
  
  # gc content per gene
  gene_table <- gc_count_gene(genome, gene_table)
  return(gene_table)
}

#' Generate gene table from GenBankRecord
#'
#' @param gb Formalclass GenBankRecord. It is generated by function
#' \code{\link[genbankr]{readGenBank}}
#' @param genome a DNAstring object. It contains the genome sequence.
#'
#' @return a data frame. It contains information of genes.
#' @importFrom magrittr %>%
#' @import dplyr
#' 
geneTableRead <- function(gb, genome){
  genes <- as.data.frame(genbankr::genes(gb))
  
  genes$gene[is.na(genes$gene)] <- genes$gene_id[is.na(genes$gene)]
  
  if (!"pseudo" %in% colnames(genes)){
    genes$pseudo <- rep(FALSE, nrow(genes))
  }
  features <- as.data.frame(genbankr::otherFeatures(gb))
  
  features <- features %>%
    dplyr::mutate(pseudo = rep(FALSE, n())) %>%
    dplyr::filter(type %in% c("rRNA", "tRNA")) %>%
    dplyr::filter(!gene %in% genes$gene)
  
  features$gene[is.na(features$gene)] <- features$product[is.na(features$gene)]
  
  if (nrow(features) != 0){
    gene_table <- genes %>%
      select(start, end, gene, strand, pseudo) %>%
      rbind.data.frame(select(features, start, end, gene, strand, pseudo)) %>%
      mutate(chr = rep("chr1", n())) %>%
      stats::na.omit() %>%
      select(chr, start, end, gene, strand, pseudo) %>%
      unique()
  } else {
    gene_table <- genes %>%
      select(start, end, gene, strand, pseudo) %>%
      mutate(chr = rep("chr1", n())) %>%
      stats::na.omit() %>%
      select(chr, start, end, gene, strand, pseudo) %>%
      unique()
  }
  
  gene_table$gene[grepl(".*([0-9\\.]+)S.*", gene_table$gene)] <-
    rrnFixer(gene_table$gene[grepl(".*([0-9\\.]+)S.*", gene_table$gene)])
  gene_table$gene[grepl("^trn.*", gene_table$gene, ignore.case=TRUE)] <-
    trnFixer(gene_table$gene[grepl("^trn.*", gene_table$gene, ignore.case=TRUE)])
  
  # codon usage
  cds <- as.data.frame(genbankr::cds(gb))
  cds_cu <- codonUsage(cds, genome)
  if (is.null(cds_cu)){
    gene_table$cu_bias <- rep(NA, nrow(gene_table))
  } else {
    gene_table <- dplyr::left_join(gene_table, cds_cu,
                                   by = c("gene", "strand", "start"))
  }
  
  # gc content per gene
  gene_table <- gc_count_gene(genome, gene_table)
  return(gene_table)
}

rrnFixer <- function(rRNA){
  rRNA <- sub("[a-zA-Z]*([0-9\\.]*)S.*", "\\1", rRNA)
  rRNA <- paste("rrn", rRNA, sep = "")
}

trnfixer <- function(Genelist){
  for (i in 1:length(Genelist[,1])){
    
    if (Genelist[i,1] %in% c("trnK-UUU", "tRNA-Lys")){
      Genelist[i,1]<- "trnK"
    }
    
    if (Genelist[i,1] %in% c("trnN-UUA", "tRNA-Asn")){
      Genelist[i,1]<- "trnN"
    }
    
    if (Genelist[i,1] %in% c("trnT-UGA", "tRNA-Thr")){
      Genelist[i,1]<- "trnT"
    }
    
    if (Genelist[i,1] %in% c("trnR-UCU", "tRNA-Arg")){
      Genelist[i,1]<- "trnR"
    }
    
    if (Genelist[i,1] %in% c("trnM-UAC", "tRNA-Met")){
      Genelist[i,1]<- "trnM"
    }
    
    if (Genelist[i,1] %in% c("trnI-UAA", "tRNA-Ile")){
      Genelist[i,1]<- "trnI"
    }
    
    if (Genelist[i,1] %in% c("trnQ-GUU", "tRNA-Gln")){
      Genelist[i,1]<- "trnQ"
    }
    
    if (Genelist[i,1] %in% c("trnH-GUA", "tRNA-His", "trnH-GUG")){
      Genelist[i,1]<- "trnH"
    }
    
    if (Genelist[i,1] %in% c("trnP-GGG", "tRNA-Pro")){
      Genelist[i,1]<- "trnP"
    }
    
    if (Genelist[i,1] %in% c("trnE-CUC", "tRNA-Glu")){
      Genelist[i,1]<- "trnE"
    }
    
    if (Genelist[i,1] %in% c("trnD-CUA", "tRNA-Asp")){
      Genelist[i,1]<- "trnD"
    }
    
    if (Genelist[i,1] %in% c("trnA-CGC", "tRNA-Ala")){
      Genelist[i,1]<- "trnA"
    }
    
    if (Genelist[i,1] %in% c("trnG-CCC", "tRNA-Gly")){
      Genelist[i,1]<- "trnG"
    }
    
    if (Genelist[i,1] %in% c("trnV-CAC", "tRNA-Val")){
      Genelist[i,1]<- "trnV"
    }
    
    if (Genelist[i,1] %in% c("trnY-AUA", "tRNA-Tyr")){
      Genelist[i,1]<- "trnY"
    }
    
    if (Genelist[i,1] %in% c("trnS-AGA", "tRNA-Ser")){
      Genelist[i,1]<- "trnS"
    }
    
    if (Genelist[i,1] %in% c("trnW-ACC", "tRNA-Trp")){
      Genelist[i,1]<- "trnW"
    }
    
    if (Genelist[i,1] %in% c("trnC-ACA", "tRNA-Cys")){
      Genelist[i,1]<- "trnC"
    }
    
    if (Genelist[i,1] %in% c("trnL-AAU", "tRNA-Leu")){
      Genelist[i,1]<- "trnL"
    }
    
    if (Genelist[i,1] %in% c("trnF-AAA", "tRNA-Phe")){
      Genelist[i,1]<- "trnF"
    }
  }
  return(Genelist)
}

trnFixer <- function(tRNA) {
  #tRNA <- gene_table$gene[grepl("^trn.*", gene_table$gene, ignore.case=TRUE)]
  tRNA <- sub("-", "", tRNA)
  tRNA <- sub("^tRNA", "trn", tRNA)
  aa_table <- rbind(c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu",
                      "Gln", "Gly", "His", "He", "Leu", "Lys",
                      "Met", "Phe", "Pro", "Ser", "Thr", "Trp",
                      "Tyr", "Val"),
                    c("A", "R", "N", "D", "C", "E", "Q", "G", "H",
                      "I", "L", "K", "M", "F", "P", "S", "T", "W",
                      "Y", "V"))
  
  for (i in 1:ncol(aa_table)){
    tRNA <- sub(aa_table[1, i], aa_table[2, i], tRNA)
  }
  tRNA <- sub("(trnf*[A-Z]).*", "\\1", tRNA)
  #gene_table$gene[grepl("^trn.*", gene_table$gene, ignore.case=TRUE)] <- tRNA
  #return(gene_table)
  return(tRNA)
}

trnCut <- function(GL){
  name<- GL[,1]
  for (i in 1:length(name)){
    if (length(strsplit(as.character(name[i]), split="")[[1]])>5){
      name[i]<-paste(strsplit(as.character(name[i]), split="")[[1]][1:4], collapse = "")
    }
  }
  GL[,1]<- as.character(name)
  return(GL)
}

codonUsage <- function(cds, genome){
  # Forward strand
  cds_cu_f <- NULL
  cds_f <- cds[which(cds$strand == "+"),]
  tmp <- cds_f
  if (nrow(tmp) != 0){
    cds_seq_f <- Biostrings::DNAStringSet(genome, start = cds_f$start[1],
                                          end = cds_f$end[1])
    if (nrow(tmp) > 1){
      i <- 2
      repeat{
        if (tmp$gene[i] == tmp$gene[i - 1]) {
          cds_seq_f[[i - 1]] <- c(cds_seq_f[[i - 1]],
                                  Biostrings::subseq(genome,
                                                     start = cds_f$start[i],
                                                     end = cds_f$end[i]))
          t <- tmp[i, ]
          tmp <- tmp[-i, ]
        } else {
          cds_seq_f <- append(cds_seq_f, Biostrings::DNAStringSet(genome,
                                                                  start = cds_f$start[i],
                                                                  end = cds_f$end[i]))
          t <- tmp[i, ]
          i <- i + 1
        }
        if (identical(t, cds_f[nrow(cds_f),])) {
          break()
        } else {
          t <- NULL
        }
      }
    }
    
    names(cds_seq_f) <- tmp$gene
    cds_cu_f <- coRdon::codonTable(cds_seq_f)
    cds_cu_f <- as.vector(coRdon::MILC(cds_cu_f))
    cds_cu_f <- data.frame(cu_bias = cds_cu_f, gene = tmp$gene,
                           start = tmp$start,
                           strand = rep("+", nrow(tmp)),
                           stringsAsFactors = FALSE)
  }
  
  
  # Reverse strand
  cds_cu_r <- NULL
  cds_r <- cds[which(cds$strand == "-"),]
  tmp <- cds_r
  if (nrow(tmp) != 0){
    cds_seq_r <- Biostrings::DNAStringSet(genome, start = cds_r$start[1],
                                          end = cds_r$end[1])
    if (nrow(tmp) > 1){
      i <- 2
      repeat{
        if (tmp$gene[i] == tmp$gene[i - 1]) {
          cds_seq_r[[i - 1]] <- c(cds_seq_r[[i - 1]],
                                  Biostrings::subseq(genome,
                                                     start = cds_r$start[i],
                                                     end = cds_r$end[i]))
          t <- tmp[i, ]
          tmp <- tmp[-i, ]
        } else {
          cds_seq_r <- append(cds_seq_r, Biostrings::DNAStringSet(genome,
                                                                  start = cds_r$start[i],
                                                                  end = cds_r$end[i]))
          t <- tmp[i, ]
          i <- i + 1
        }
        if (identical(t, cds_r[nrow(cds_r),])) {
          break()
        } else {
          t <- NULL
        }
      }
    }
    names(cds_seq_r) <- tmp$gene
    cds_seq_r <- Biostrings::reverseComplement(cds_seq_r)
    cds_cu_r <- coRdon::codonTable(cds_seq_r)
    cds_cu_r <- as.vector(coRdon::MILC(cds_cu_r))
    cds_cu_r <- data.frame(cu_bias = cds_cu_r, gene = tmp$gene,
                           start = tmp$start,
                           strand = rep("-", nrow(tmp)),
                           stringsAsFactors = FALSE)
  }
  
  if (!is.null(cds_cu_f) & !is.null(cds_cu_r)){
    cds_cu <- rbind.data.frame(cds_cu_f, cds_cu_r)
    colnames(cds_cu) <- c("cu_bias", "gene", "start", "strand")
  } else if(!is.null(cds_cu_f)){
    cds_cu <- cds_cu_f
    colnames(cds_cu) <- c("cu_bias", "gene", "start", "strand")
  } else if(!is.null(cds_cu_r)){
    cds_cu <- cds_cu_r
    colnames(cds_cu) <- c("cu_bias", "gene", "start", "strand")
  } else {
    cds_cu <- NULL
  }
  
  return(cds_cu)
}

#' Reverse SSC section in the plot.
#'
#' @param g gene coordinates for SSC region.
#' @param SSC_start number: SSC start
#' @param SSC_end number: SSC end
#'
#' @return g edited to be reversed.
#' 
SSCrev<- function(g, SSC_start, SSC_end){
  g<- as.data.frame(g)
  g[, 2]<- as.numeric(g[,2])
  g[, 3]<- as.numeric(g[,3])
  for (i in 1:length(g[,1])){
    new_start <- SSC_end - as.numeric(g[i, 'end']) + SSC_start
    new_end <- SSC_end - as.numeric(g[i, 'start']) + SSC_start
    
    if(as.numeric(g[i,'start']) >= SSC_start && as.numeric(g[i,'end']) <= SSC_end){
      g[i, 2]<- new_start
      g[i, 3]<- new_end
      if(g[i,'strand'] == '+'){
        g[i, 'strand'] <- '-'
      } else {
        g[i,"strand"] <- '+'
      }
    }
    
    #If the gene is starting from SSC then it should also be fixed
    
    if(as.numeric(g[i,'start']) < SSC_start && as.numeric(g[i,'end']) > SSC_start){
      # We save the part outside SSC
      g <- rbind(g, c(g[i,"gene"], g[i,'start'], SSC_start, 0, 0, g[i,'strand']))

      g[i, 2]<- new_start
      g[i, 3]<- SSC_end
      if(g[i,'strand'] == '+'){
        g[i, 'strand'] <- '-'
      } else {
        g[i,"strand"] <- '+'
      }
    } else if(as.numeric(g[i,'start']) < SSC_end && as.numeric(g[i,'end']) > SSC_end){
      # We save the part outside SSC
      g <- rbind(g, c(g[i,"gene"], SSC_end, g[i,'end'], 0, 0, g[i,'strand']))

      g[i, 2]<- SSC_start
      g[i, 3]<- new_end
      if(g[i,'strand'] == '+'){
        g[i, 'strand'] <- '-'
      } else {
        g[i,"strand"] <- '+'
      }
    }
  }
  g <- as.matrix(g)
  return(g)
}
