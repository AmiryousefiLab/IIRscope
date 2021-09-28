#' # Questions:
#' # Is it necessary we check column number of ORIGIN field in function "FasExtract"?
#' 
#' 
#' #' Genome extracter
#' #'
#' #' Extracting the fasta format chloroplast genome from the GeneBank File
#' #' @param gb a large character vector returned by function
#' #' \code{\link{fetch.gb}} or \code{\link{read.gb}}
#' #'
#' #' @return The genome of the sequence that is deposited at the end of the
#' #' GeneBank file in fasta format
#' #' @export
#' FasExtract <- function(gb){
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
#' 
#' #' Randomly fill nucleotides
#' #'
#' #' This fuction change liters "u", "r", "y", "s", "w", "k", "m", "b", "d", "h",
#' #' "v", "n" or "-" in genome sequence (randomly) into corresponding nucleotide
#' #' "a", "t", "c" or "g".
#' #'
#' #' @param genome A DNAString object. It is a class defined in \code{Biostrings}
#' #' package. It contains the whole genome sequence of a specie.
#' #'
#' #' @return A large character vector. It countains the nucleotide sequence with
#' #' all letters except "a", "t", "c", "g" changed.
#' #' @export
#' #'
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


FasExtract<- function(gb){
  #' Genome extracter
  #'
  #' Extracting the fasta format chloroplast genome from the GeneBank File
  #' @param file Name of the GeneBank file of a sequence on the working reposotory
  #' @return The genome of the sequence that is deposited at the end of the GeneBank file in fasta format
  #' @export
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


sp.name<- function(gb){
  paste(strsplit(gb[2], " ")[[1]][3], strsplit(gb[2], " ")[[1]][4], sep=" ")
}

# TODO change to this one better
#' #' Extract specie's name from GB file
#' #'
#' #' @param definition a charactor. It contain the definition field from the
#' #' GenBank file.
#' #'
#' #' @return a character contains the specie's name.
#' #' @export
#' #'
#' sp.name <- function(definition){
#'   # if (text){
#'     # sp <- gsub("(DEFINITION\\ \\ )", "", definition[2], perl = TRUE)
#'     # sp <- sub("(\\w+\\s+\\w+).*", "\\1", sp, perl = TRUE)
#'   # } else {
#'     sp <- sub("(\\w+\\s+\\w+).*", "\\1", definition, perl = TRUE)
#'   # }
#'   return(sp)
#' }

