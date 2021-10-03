# TODO: change functions to accept either one IRinfo format or the other 
#       (gbfiles/dogma) and not be repeated.
# TODO: pass onto functions the global variables or save them somewhere else.

#### This section contains the functions related to the gbfiles input ####

#' Generates all needed information to plot the IR region given gbfiles
#'
#' @param gbfiles named list with the gbfiles. The parameters of the list are:
#' \itemize{
#'   \item gbfile: string with the accession number os path with the gbfile.
#'   \item local.file: boolean TRUE if it is a pth with a local gbfile or 
#'         FALSE if accession number.
#' }
#' @param Sfiles list of booleans. Indicates if the file needs a SSC reversion
#' @param progress progress bar from shiny to show how many species have been 
#'        analyzed so far.
#' @return None. The variables set here are used as global variables (l, IRList, 
#'         IndelList, FasList, nuclw, spnames, GeneList)
#' @export
IRs<- function(gbfiles, Sfiles, progress){
  l<<- length(gbfiles)  # number of species
  FasList<<-  list()    # list with the fasta genomes
  IRList<<-   list()    # list with the IR infos
  IndelList <<- list()  # list with the indel info
  GeneList<<- list()    # list with the genes 
  spnames<<-  list()    # list with the species name
  nuclw<<- numeric(l)   # vector to scale the plot later
  
  for (i in 1:l){
    # extracts data from gb file to get the ir region and the indel
    gb.info <- gbfiles[[i]]
    data <- get_data_from_gb(gb.info$gbfile, gb.info$local.file)
    
    res <- get_ir(data$genome)
    IRList[[i]] <<- res$ir_table
    if(is.null(res$indel_table)){
      IndelList[[i]] <<- NA
    } else {
      IndelList[[i]] <<- res$indel_table
    }
    
    # gets data in the needed format to get the genes. 
    #   TODO: change gene.cordinates, sp.name... to receive the same kind of data as above.
    if(gb.info$local.file) {
      gb.data <- read.gb(gb.info$gbfile)
    } else {
      gb.data <- fetch.gb(gb.info$gbfile)
    }
    GeneList[[i]] <<-gene.cordinates(gb.data)
    
    # if ssc reversion was ticked, it reverse that part in the gene list
    rev <- Sfiles[[i]]
    if (rev){
      GeneList[[i]]<<- SSCrev(GeneList[[i]], (IRList[[i]][1]+IRList[[i]][3]), IRList[[i]][2])
    }
    spnames[[i]]<<- sp.name(gb.data)
    
    FasList[[i]]<<- rdnFixer(gb.data)
    nuclw[i]<<- 100/IRList[[i]][4]
    
    # getting genelist in dogma format
    GeneList[[i]]<<- trnfixer(GeneList[[i]])
    GeneList[[i]]<<- trnCut(GeneList[[i]])
    
    # when one specie is analysed, progress bar grows
    progress$set(value = i, 
                 detail = paste0('This may take a while... ', i, '/', l))
  }
}

#### Transforming data from gbfile section to dogma and viceversa ####

#' Transform IRInfo format to the one used in gbfiles functions
#'
#' @param IRDinp vector with ira ini, ira end, irb ini, irb end, len genome.
#' @return vector with ira ini, ira end, len ir, len genome.
toIRinfo <- function(IRDinp){
  if(IRDinp[2] - IRDinp[1] != IRDinp[4] - IRDinp[3]){
    stop("IR lens are not the same")
  }
  return(c(IRDinp[1], IRDinp[3], IRDinp[2] - IRDinp[1], IRDinp[4]))
}

#' Transform IRInfo format to the one used in dogma functions
#'
#' @param IRinfo vector with ira ini, ira end, len ir, len genome.
#' @return vector with ira ini, ira end, irb ini, irb end, len genome
toIRDinp <- function(IRinfo){
  return(c(IRinfo[1], IRinfo[1]+IRinfo[3], IRinfo[2], IRinfo[2]+IRinfo[3], IRinfo[4]))
}


#### This section contains the functions related to the dogma input ####


#' Generates all needed information to plot the IR region given fasta and gene 
#' info files
#'
#' @param dgfiles list with the dogma files (files with gene info) paths.
#' @param fastafiles list with the fasta files paths.
#' @param irfiles list with ir information provided.
#' @param nfiles list of species names.
#' @param progress progress bar from shiny to show how many species have been 
#'        analyzed so far.
#' @return None. The variables set here are used as global variables (IRListDinp, 
#'         IndelList, FasList, GeneList, nuclw, spnames)
#' @export
IRsD<- function(dgfiles, fastafiles, irfiles, nfiles, progress){
  l<<- length(dgfiles)  # number of species
  FasList<<-  list()    # list with the fasta genomes
  IRListDinp<<- list()    # list with the IR infos
  IndelList <<- list()  # list with the indel info
  GeneList<<- list()    # list with the genes 
  spnames<<-  list()    # list with the species name
  nuclw<<- numeric(l)   # vector to scale the plot later
  
  for (i in 1:l){
    dg <<- dgfiles[[i]]
    FasList[[i]] <<- fastafiles[[i]]
    GeneList[[i]] <<- dgfiles[[i]]
    spnames[[i]] <<- nfiles[[i]]
    
    # if no ir info is provided
    if (irfiles[[i]][1]==999){
      # gets genome and transforms the IRinfo into the one used in the functions 
      # for dogma files.
      #   TODO: change functions to accept either one IRinfo format or the other
      #         and not be repeated.
      genome <- readDNAStringSet(FasList[[i]])[[1]]
      res <- get_ir(genome)
      sss <- res$ir_table
      IRListDinp[[i]] <<- c(sss[1], sss[1]+sss[3], sss[2], sss[2]+sss[3], sss[4])
      if(is.null(res$indel_table)){
        IndelList[[i]] <<- NA
      } else {
        IndelList[[i]] <<- res$indel_table
      }
    } else { # the ir info is provided
      IRListDinp[[i]] <<- irfiles[[i]]
      IndelList[[i]] <<- NA
    }
    nuclw[i] <<- 100/length(FasList[[1]])
    GeneList[[i]] <<- GeneList[[i]]
    
    progress$set(value = i, 
                 detail = paste0('This may take a while... ', i, '/', l))
  }
}