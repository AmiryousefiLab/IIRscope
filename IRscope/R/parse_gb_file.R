# TODO: check different functions versions and see which one performs better.

#' Getting data from gb file
#' 
#' Extracting the needed information from a gb file.
#' 
#' @param gbfile A character. It is the GI or GenBank accession for the
#' interested specie's chloroplast genome files.
#' @param local.file A logical value. If it is \code{TRUE}, the value passed to
#' \code{gbfile} should be a local GB file.
#' @return list with:
#' \itemize{
#'   \item gb: read gb file.
#'   \item genome: A DNAString object. It is a class defined in \code{Biostrings}
#'    package. It contains the whole genome sequence of a specie from 
#'    the given gb file.
#'   \item l: length of the genome.
#'   \item gene_table: table with the genes and positions.
#' }
get_data_from_gb <- function(gbfile, local.file = FALSE){
  my_env <- new.env(parent = emptyenv())
  my_env$gb <- NULL
  my_env$genome <- NULL
  my_env$l <- NULL
  my_env$gene_table <- NULL # TODO: use this instead of the other way to genes
  if (local.file){
    tryCatch({
      my_env$gb <- genbankr::readGenBank(gbfile)
      my_env$genome <- genbankr::getSeq(my_env$gb)[[1]]
      my_env$l<- Biostrings::nchar(my_env$genome)
      my_env$gene_table <- geneTableRead(my_env$gb, my_env$genome)
    }, error = function(e){
      tryCatch({
        # TODO: parseGenBank doesn't detect ORIGIN without trailing spaces, 
        #       this is a temporary solution
        tx = readLines(gbfile)
        tx = gsub(pattern = "ORIGIN", replace = "ORIGIN ", x = tx)
        writeLines(tx, con=gbfile)
        
        my_env$gb <- genbankr::parseGenBank(gbfile)
        my_env$genome <- my_env$gb$ORIGIN[[1]] 
        my_env$l <- Biostrings::nchar(my_env$genome)
        my_env$gene_table <- geneTableParsed(my_env$gb, my_env$genome)
      }, error = function(e){
        my_env$gb <- read.gb(gbfile)
        my_env$genome <- get_genome(my_env$gb)
        my_env$l <- Biostrings::nchar(my_env$genome)
        my_env$gene_table <- gene.cordinates(my_env$gb)
      })
    })
    
  } else {
    tryCatch({
      my_env$gb <- genbankr::readGenBank(text = fetch.gb(gbfile))
      my_env$genome <- genbankr::getSeq(my_env$gb)
      my_env$l <- Biostrings::nchar(my_env$genome)
      my_env$gene_table <- geneTableRead(my_env$gb, my_env$genome)
    }, error = function(e){
      my_env$gb <- genbankr::parseGenBank(text = fetch.gb(gbfile))
      my_env$genome <- my_env$gb$ORIGIN[[1]]
      my_env$l <- Biostrings::nchar(my_env$genome)
      my_env$gene_table <- geneTableParsed(my_env$gb, my_env$genome)
    })
  }
  return(list(gb = my_env$gb, genome = my_env$genome, l = my_env$l, 
              gene_table = my_env$gene_table))
}

# TODO: join next three functions (get_genome and FasExtract). Similar behaviour 
#       but different return. See which one performs better.

#' # Questions:
#' # Is it necessary we check column number of ORIGIN field in function "get_genome"?
#' 
#' Genome extracter
#'
#' Extracting the fasta format chloroplast genome from a GeneBank File
#' @param gb vector of characters. Genebank file returned by
#' \code{\link{fetch.gb}} or \code{\link{read.gb}} function.
#' @return The genome of the sequence that is deposited at the end of the GeneBank file in DNAString object format
get_genome <- function(gb){
  fasta <- gb[(grep("ORIGIN", gb)+1):length(gb)]
  while(fasta[length(fasta)]=="") {
    fasta <- fasta[1:length(fasta)-1]
  }
  while(fasta[length(fasta)]=="//") {
    fasta <- fasta[1:length(fasta)-1]
  }
  fas <- ""
  for (i in 1:length(fasta)){
    sort.let<- sort(unique(c(grep("c", strsplit(fasta[i], " ")[[1]]),
                             grep("a", strsplit(fasta[i], " ")[[1]]), 
                             grep("t", strsplit(fasta[i], " ")[[1]]), 
                             grep("g", strsplit(fasta[i], " ")[[1]]))))
    try(if(length(sort.let)==!6) stop("Check the gb file; the columns of ORIGIN should be 6"))
    fasta[i]<- paste0(strsplit(fasta[i], " ")[[1]][sort.let[1]],
                      strsplit(fasta[i], " ")[[1]][sort.let[2]],
                      strsplit(fasta[i], " ")[[1]][sort.let[3]],
                      strsplit(fasta[i], " ")[[1]][sort.let[4]],
                      strsplit(fasta[i], " ")[[1]][sort.let[5]],
                      strsplit(fasta[i], " ")[[1]][sort.let[6]])
    fas<-paste0(fas, fasta[i])
  }
  fas<-gsub("NA", "", fas)
  genome <- strsplit(fas, "")[[1]]
  genome <- paste0(genome, collapse = '')
  DNAString(genome)
}

#' Genome extracter
#'
#' Extracting the fasta format chloroplast genome from the GeneBank File
#' @param gb vector of characters. Genebank file returned by
#' \code{\link{fetch.gb}} or \code{\link{read.gb}} function.
#' @return The genome of the sequence that is deposited at the end of the 
#' GeneBank file in fasta format
FasExtract<- function(gb){
  fasta<-gb[(grep("ORIGIN", gb)+1):length(gb)]
  while(fasta[length(fasta)]=="") {
    fasta<- fasta[1:length(fasta)-1]
  }
  while(fasta[length(fasta)]=="//") {
    fasta<- fasta[1:length(fasta)-1]
  }
  fas<-""
  for (i in 1:length(fasta)){
    sort.let<- sort(unique(c(grep("c", strsplit(fasta[i], " ")[[1]]),grep("a", strsplit(fasta[i], " ")[[1]]), grep("t", strsplit(fasta[i], " ")[[1]]), grep("g", strsplit(fasta[i], " ")[[1]]))))
    try(if(length(sort.let)==!6) stop("Check the gb file; the columns of ORIGIN should be 6"))
    fasta[i]<- paste(strsplit(fasta[i], " ")[[1]][sort.let[1]],strsplit(fasta[i], " ")[[1]][sort.let[2]],strsplit(fasta[i], " ")[[1]][sort.let[3]],strsplit(fasta[i], " ")[[1]][sort.let[4]],strsplit(fasta[i], " ")[[1]][sort.let[5]],strsplit(fasta[i], " ")[[1]][sort.let[6]], sep="")
    fas<-paste(fas, fasta[i], sep="")
  }
  #fasta[length(fasta)]<- gsub("NA", "", fasta[length(fasta)])
  fas<-gsub("NA", "", fas)
  strsplit(fas, "")[[1]]
}


# Alternate version
#' get_genome <- function(gb){
#' 
#'   fasta<-gb[(grep("ORIGIN", gb)+1):length(gb)]
#' 
#'   # Remove empty lines
#'   fasta <- fasta[-which(fasta == "" | fasta == "//")]
#'   # while(fasta[length(fasta)]=="") {
#'   #   fasta<- fasta[1:length(fasta)-1]
#'   # }
#'   #
#'   # while(fasta[length(fasta)]=="//") {
#'   #   fasta<- fasta[1:length(fasta)-1]
#'   # }
#' 
#'   # Remove position indexes
#'   fasta <- substring(fasta, 11)
#'   # Extract all letters from fasta
#'   fasta <- gsub("[^a-zA-Z\\-]", "", fasta)
#'   fasta <- Reduce(c, strsplit(fasta, ""))
#'   fasta <- paste(fasta, collapse = "")
#'   fasta <- Biostrings::DNAString(fasta)
#'   fasta <- rdnFixer(fasta)
#'   return(fasta)
#' }
 

#' Randomly fill nucleotides
#'
#' This fuction change liters "u", "r", "y", "s", "w", "k", "m", "b", "d", "h",
#' "v", "n" or "-" in genome sequence (randomly) into corresponding nucleotide
#' "a", "t", "c" or "g".
#'
#' @param genome A DNAString object. It is a class defined in \code{Biostrings}
#' package. It contains the whole genome sequence of a specie.
#' @return A large character vector. It countains the nucleotide sequence with
#' all letters except "a", "t", "c", "g" changed.
rdnFixer<- function(gb){
  seq<- FasExtract(gb)
  seq[which(seq=="u")]<-sample(c("t"), length(which(seq=="u")), TRUE)
  seq[which(seq=="r")]<-sample(c("a", "g"), length(which(seq=="r")), TRUE)
  seq[which(seq=="y")]<-sample(c("c", "t"), length(which(seq=="y")), TRUE)
  seq[which(seq=="s")]<-sample(c("c", "g"), length(which(seq=="s")), TRUE)
  seq[which(seq=="w")]<-sample(c("a", "t"), length(which(seq=="w")), TRUE)
  seq[which(seq=="k")]<-sample(c("g", "t"), length(which(seq=="k")), TRUE)
  seq[which(seq=="m")]<-sample(c("c", "a"), length(which(seq=="m")), TRUE)
  seq[which(seq=="b")]<-sample(c("c", "g", "t"), length(which(seq=="b")), TRUE)
  seq[which(seq=="d")]<-sample(c("a", "g", "t"), length(which(seq=="d")), TRUE)
  seq[which(seq=="h")]<-sample(c("c", "a", "t"), length(which(seq=="h")), TRUE)
  seq[which(seq=="v")]<-sample(c("c", "a", "g"), length(which(seq=="v")), TRUE)
  seq[which(seq=="n")]<-sample(c("c", "g", "t", "a"), length(which(seq=="n")), TRUE)
  seq[which(seq=="-")]<-sample(c("c", "g", "t", "a"), length(which(seq=="-")), TRUE)
  return(seq)
}

# Alternate version
#' rdnFixer<- function(genome){
#'   seq <- Biostrings::toString(genome)
#'   seq <- unlist(strsplit(seq, ""))
#'   seq[which(seq=="U")]<-sample(c("T"), length(which(seq=="U")), TRUE)
#'   seq[which(seq=="R")]<-sample(c("A", "G"), length(which(seq=="R")), TRUE)
#'   seq[which(seq=="Y")]<-sample(c("C", "T"), length(which(seq=="Y")), TRUE)
#'   seq[which(seq=="S")]<-sample(c("C", "G"), length(which(seq=="S")), TRUE)
#'   seq[which(seq=="W")]<-sample(c("A", "T"), length(which(seq=="W")), TRUE)
#'   seq[which(seq=="K")]<-sample(c("G", "T"), length(which(seq=="K")), TRUE)
#'   seq[which(seq=="M")]<-sample(c("C", "A"), length(which(seq=="M")), TRUE)
#'   seq[which(seq=="B")]<-sample(c("C", "G", "T"), length(which(seq=="B")), TRUE)
#'   seq[which(seq=="D")]<-sample(c("A", "G", "T"), length(which(seq=="D")), TRUE)
#'   seq[which(seq=="H")]<-sample(c("C", "A", "T"), length(which(seq=="H")), TRUE)
#'   seq[which(seq=="V")]<-sample(c("C", "A", "G"), length(which(seq=="V")), TRUE)
#'   seq[which(seq=="N")]<-sample(c("C", "G", "T", "A"), length(which(seq=="N")), TRUE)
#'   seq[which(seq=="-")]<-sample(c("C", "G", "T", "A"), length(which(seq=="-")), TRUE)
#'   seq <- paste(seq, collapse = "")
#'   seq <- Biostrings::DNAString(seq)
#'   return(seq)
#' }

#' Extract specie's name from GB file
#'
#' @param gb a character. It contains the the read GenBank file.
#' @return a character contains the specie's name.
#' @export
sp.name<- function(gb){
  sp <- paste(strsplit(gb[2], " ")[[1]][3], strsplit(gb[2], " ")[[1]][4], sep=" ")
  if(length(sp) < 2){
    sp <- gb[(grep('ORGANISM', gb))]
    sp <- sub('.*ORGANISM  ', '', sp)
  }
  if(length(sp) < 2){
    sp <- gb[(grep('DEFINITION', gb))]
    sp <- sub('.*DEFINITION  ', '', sp)
    sp <- sub(' (plastid|chloroplast).*', '', sp)
  }
  return(sp)
}

# Alternate version
#' sp.name <- function(definition){
#'   # if (text){
#'     # sp <- gsub("(DEFINITION\\ \\ )", "", definition[2], perl = TRUE)
#'     # sp <- sub("(\\w+\\s+\\w+).*", "\\1", sp, perl = TRUE)
#'   # } else {
#'     sp <- sub("(\\w+\\s+\\w+).*", "\\1", definition, perl = TRUE)
#'   # }
#'   return(sp)
#' }
#' 

#' Transform the dogma file gene list to one suited for the functions.
#'
#' @param dogma a character. It contains the the read dogma file.
#' @return a list with the genes.
GnlBuilder<- function(dogma){
  l<- length(dogma[,1])
  list<- cbind(as.character(rep(0, l)),as.character(rep(0, l)),as.character(rep(0, l)),as.character(rep(0, l)),as.character(rep(0, l)))
  for (i in 1:l){
    if (as.character(dogma[i,4])=="-"){
      p<-dogma[i, 2]
      dogma[i,2]<- dogma[i,1]
      dogma[i,1]<- p
    }
    list[i,1]<-as.character(dogma[i,3])
    list[i,2]<-as.character(dogma[i,1])
    list[i,3]<-as.character(dogma[i,2])
  }
  return(list)
}

#' Transform the dogma file gene list to one suited for the functions with the
#' gene namesfixed.
#'
#' @param dogma a character. It contains the the read dogma file.
#' @return a list with the genes with the gene names fixed.
trnDogma<- function(dogma){
  name<- as.character(dogma[,3])
  for (i in 1:length(name)){
    if (length(strsplit(as.character(name[i]), split="")[[1]])>5){
      name[i]<-paste(strsplit(as.character(name[i]), split="")[[1]][1:4], collapse = "")
    }
  }
  dogma[,3]<- as.character(name)
  return(dogma)
}
