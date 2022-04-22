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
    genome_length <- Biostrings::nchar(data$genome)
    
    IRList[[i]] <<- get_ir_positions(res$ir_table, genome_length)
    if(is.null(res$indel_table)){
      IndelList[[i]] <<- NA
    } else {
      IndelList[[i]] <<- res$indel_table
    }
    
    # GeneList[[i]] <<- gene.cordinates(gb.data)
    
    # if ssc reversion was ticked, it reverses that part in the gene list
    #   TODO: join reverse regions functions to work as one.
    rev <- Sfiles[[i]]
    gene_table <- data$gene_table
    gene_table$start <- (gene_table$start - 1) %% genome_length # Needed to inclue start inside the gene.
    
    if(is.data.frame(data$gene_table)){
      if (rev){
        data_rev <<- convert_region(ir_table = res$ir_table, genome_length, 
                                         region = "SSC", genome = data$genome, 
                                         gene_table = gene_table, 
                                         indel_table = res$indel_table)
        gene_table <- data_rev$gene_table
        
        if(!is.null(res$indel_table)){
          IndelList[[i]] <<- data_rev$indel_table
        }
      }
      GeneList[[i]] <<- toGeneList(gene_table)
    } else {
      # If it's not a dataframe it has to use another function
      GeneList[[i]] <<- unique(toGeneListWithStrand(gene_table))
      if (rev){
        GeneList[[i]]<<- SSCrev(GeneList[[i]], (IRList[[i]][1]+IRList[[i]][3]), IRList[[i]][2])
      }
    }

    # gets data in the needed format to get more information 
    #   TODO: change sp.name, rdnFixer... to receive the same kind of data as above or skip them.
    if(gb.info$local.file) {
      gb.data <- read.gb(gb.info$gbfile)
    } else {
      gb.data <- fetch.gb(gb.info$gbfile)
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

# TODO: use gene_table everywhere instead of the old matrix. Adapt code to not 
#       use these two 0's columns anymore in the matrix.

#' Transform gene_table data.frame to a matrix with less info used in the old functions
#'
#' @param gene_table data.frame with the gene information. The columns used in this
#'  function are 1 (start), 2 (end), 4 (gene name) and 6 (strand).
#' @return matrix with the following columns: gene_name, start, end, 0's, 0's and strand.
#' (these 0's columns are there for convenience with the old format).
toGeneList <- function(gene_table){
  # We transform the dataframe into a matrix with only the data we are interested:
  #   gene_name, start and end.
  m <- as.matrix(subset(gene_table, select=c(4,1,2)))
  # We make sure data types are correct
  m[,1] <- as.character(m[,1])
  m[,2] <- as.integer(m[,2])
  m[,3] <- as.integer(m[,3])
  #m[,4] <- mapvalues(m[,4], from = c('+', '-'), to = c(1,0))
  
  # We add one empty columns for convenience with the previous format ()
  m <- cbind(m, replicate(2, numeric(NROW(m))))
  # and the information about the strand
  m <- cbind(m, gene_table[[3]])
  
  return(m)
}

# TODO: use gene_table everywhere instead of the old matrix. Adapt code to not 
#       use these two 0's columns anymore in the matrix.

#' Transform gene_table data.frame to a matrix with less info used in the old functions
#'
#' @param gm gene matrix with the gene information.
#' @return matrix with the following columns: gene_name, start, end, 0's, 0's and strand.
#' (these 0's columns are there for convenience with the old format).
toGeneListWithStrand <- function(gm){
  colnames(gm) <- c("gene", "start", "end", "start1", "end1")
  df <- as.data.frame(gm)
  
  # We separate those columns where the 4th and 5th columns have values
  seconds <- which(df$start1 != 0 | df$end1 !=0)
  seconds_len <- length(seconds)
  if(seconds_len > 0){
    for(i in 1:seconds_len){
      df[nrow(df) + 1,] <- c(df[seconds[i],"gene"], df[seconds[i],"start1"], df[seconds[i],"end1"], 0, 0)
      df[seconds[i],'start1'] <- 0
      df[seconds[i],'end1'] <- 0
    }  
  }
  
  # We order the values and add the strand column
  df_len <- nrow(df)
  df$strand <- rep(NA, df_len)
  for(i in 1:df_len){
    if(as.numeric(df[i,'start']) < as.numeric(df[i,'end'])){ # if they are in order the strand is -
      df[i,'strand'] <- '-'
    } else { # else we reverse the order, and the strand is +
      aux <- df[i,'start']
      df[i,'start'] <- df[i, 'end']
      df[i, 'end'] <- aux
      df[i,'strand'] <- '+'
    }
  }
  return(as.matrix(df))
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
      genome_length <- Biostrings::nchar(genome)
      
      sss <<- get_ir_positions(res$ir_table, genome_length)

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