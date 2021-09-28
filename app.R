library(shiny)
library(shinyjs)
library(shinyvalidate)
library(shinyWidgets)
library(shinyalert)
library(purrr)

library(BiocManager)
options(repos = BiocManager::repositories())

## TODO cambiar a libreria
source('IRscope.R') #source('/home/dmcarmen/Desktop/IRscope_Chloroplot/IRscope/IRscope.R')
pkgload::load_all(path = "Chloroplot") #pkgload::load_all(path = "/home/dmcarmen/Desktop/IRscope_Chloroplot/IRscope/Chloroplot")

## AUXILIAR FUNCTIONS TODO move
rowInputGB <- function(x) {
  fluidRow(column(width = 2, 
                  textInput(inputId = paste0("acc",x), label = "Accession No.", 
                            placeholder = paste0(x," Accession ..."))),
           column(width = 3,
                  fileInput(inputId = paste0("data",x), label = "GeneBank File",
                            placeholder = "or GB file")),
           # TODO
           # column(width = 2, actionButton(paste0("removefile",x), "Remove", icon = icon("trash"))),
           column(width = 5, br(), verbatimTextOutput(paste0("sp",x), placeholder = TRUE)), 
           column(width = 1, strong("SSC"), br(), 
                  checkboxInput(paste0("S",x), label = icon("transfer", lib = "glyphicon"), value = FALSE)), 
           column(width = 1, br())
  )
}

rowInputManual <- function(x) {
  fluidRow(br(), 
           column(width = 3, 
                  fileInput(inputId = paste0("dogma",x), label = "Annotations", 
                            placeholder = paste0(x,"th Annotation ..."))),
           column(width = 3, 
                  fileInput(inputId = paste0("fas",x), label = "Genome Fasta", 
                            placeholder = "And Fasta file")), 
           column(width = 6, 
                  textInput(inputId = paste0("IR",x),label = "IR info (optional)", 
                            placeholder = "IRb end, IRb start, IRa end, IRa start")), 
           column(width = 2, br())
  )
}

# # TEST
# # Make a palette of 40 colors
# colors <- rainbow(40, alpha = NULL)
# # Mirror the rainbow, so we cycle back and forth smoothly
# colors <- c(colors, rev(colors[c(-1, -40)]))
# # TEST
# TODO 
# urlToGetAllOpenBugs = "https://api.github.com/repos/jquery/jquery/issues?state=open&labels=bug";

ui <- fluidPage(
  br(),
  useShinyalert(),
  HTML('<!DOCTYPE html>
       <html>
       <head>
       <style type="text/css">
       h1 {color:#267A43;}
       h4 {color:black;}
       h3 {color:navy;}
       p {color:navy;}
       </style>
       </head>
       <body>
       
       <div style="opacity:0.1;position:absolute;right:16px;left:16px;height:178px;border-radius:300px;background-color:#40B3DF"></div>
       <div style="padding:20px;border-radius:300px;border:20px solid #1E8449";">
       <div style="opacity:0.3;position:absolute;background-color:#8E44AD"></div>
       
       <h1 style="font-family:Copperplate Gothic Bold;letter-spacing:8px;text-align:center;">IRSCOPE</h1>
       <h4 style = "font-family:Charlesworth;text-align:center;color:navy;"><em> Tool for visualizing the junction sites of the chloroplast genome </em> </h4>

       </div>
       </body>
       </html>
       '),
  
  # tags$h1(style = "font-family:Copperplate Gothic Bold", "IRSCOPE"),
  #tags$h4(style = "font-family:Charlesworth", tags$em("Tools for analyzing and visualizing the chloroplast genomes.")),
  #tags$hr(),
  tags$br(),
  tags$head(
    tags$head(tags$style("#StatusGB{color: white;
                         background: lightslategrey;
                         white-space: pre-wrap;}")),
    # # TEST TODO
    # tags$script('$(document).ready(function () {
    #   $.getJSON(urlToGetAllOpenBugs, function (allIssues) {
    #     $("div").append("found " + allIssues.length + " issues</br>");
    #     $.each(allIssues, function (i, issue) {
    #       $("div")
    #       .append("<b>" + issue.number + " - " + issue.title + "</b></br>")
    #       .append("created at: " + issue.created_at + "</br>")
    #       .append(issue.body + "</br></br></br>");
    #     });
    #   });
    # });'),
    # # TEST
    ),
  
  #navlistPanel(widths = c(1, 11), well = FALSE,
  # tabPanel("Home",column(width = 6, h4(tags$em("Welcome to IRscope...")), br(), p(tags$strong("Chloroscope"),  "is a suit for analyzing and visualazing the chloroplast genomes. The tools can be found in the 'Tools' tab. Currently only the IRScope is a functional program of the website and the more information related to that can be found in the related pages."))),
  #tabPanel("Tools",
  fluidRow(column(width = 1, ""), column(width = 11, offset = 1,
                                         tabsetPanel(
                                           tabPanel(strong("Home"), br(),
                                                    column( width=6, h4(tags$em("Welcome to IRscope...")), br(), tags$p(strong("IRscope"),"is a tool for  visualizing the genes on the boundaries of the junction sites of the chloroplast genome. The input files can be uploaded using the 'Files upload' tab. After the submission, the program will start searching the", tags$em("Inverted Repeat"), "regions. Note that the species selected will lead to a consensus radius for each junction site for the genes to be plotted. This indicates that if the input files are not representing closely related species, the program may lead to unsatisfactory result. Finally note that the input files can be either manual annotations, GeneBank accession numbers, or related GB files of the chloroplast genomes.")
                                                    )),
                                           tabPanel("Files upload", br(), fluidRow(column(width = 11, br(), p("You need to provide information related to at least two species for the analysis to proceed. After you have successfully provided your input files, press the 'Submit' button. This will activate its following box. Be patient till you are communicated through that box on how to proceed for downloading your output. Also note that the computations may take up to five minutes. Please be patient. For running the analysis, one can use either of the following cases."), strong("1. GB Files"), p("In this method, you can provide the GB files or the accession number of the related species to obtain the IR plot. This function is active under the 'GB File' tab"),  strong("2. Manual Files"),p("This method is useful for the times when  simple tabular format of the annotations, the genome sequences of the targeted species, and possibly their corresponding inverted repeat regions boundary coordinates are available. This function is active under 'Manual Files' tab."),  p("Further instructions for each case is given specifically under relevant tab."))),
                                                    br(), br(),br(),
                                                    tabsetPanel(
                                                      tabPanel("GB File", br(), br(), 
                                                               p("In this section you can upload either your GeneBank files or provide their corresponding accession numbers if applicable. You can denote if you want to depict the SSC region in the reverse order by ticking the 'SSC', box for a corresponding input file.",
                                                                 br(), br(), "For example, you can test the program with these accession numbers:"), tags$em("NC_007898, NC_008096"), br(), br(),
                                                               
                                                               numericInput("nGB", "Number of GB files (min 1, max 20):", 2, min = 1, max = 20),
                                                               uiOutput("dynUIGB"), # columns displayed depending on the numericInput above
                                                               column(width=10,verbatimTextOutput("stat"), actionButton("Gbsub", "Submit"), 
                                                                      br(), br(), p("After submission the section below will turn innactive, meaning that analysis is in process. Be patient until you are prompted to download your plot."),
                                                                      verbatimTextOutput("StatusGB", placeholder = TRUE), br(), verbatimTextOutput("StatusGB2", placeholder = FALSE), 
                                                                      radioButtons("file_type", "Please choose a format to download the plot.",
                                                                                   choices = c("pdf", "png", "jpeg","bmp", "tiff"), inline=TRUE),
                                                                      downloadButton(outputId = "downloadData", label = "Download"),
                                                                      # plotOutput("Plot"),
                                                                      ),
                                                               ),
                                                      tabPanel("Manual Files", br(), 
                                                               p("In case your GB files are of low quality or otherwise not available and still you would like to obtain the plot, you can provide your annotations as a simple tabular format (primary", 
                                                                 tags$a("DOGMA", href="https://dogma.ccbb.utexas.edu"), "output) beside their genomes in fasta format and proceed to obtain the plot. In case you also have the values for the junction sites coordinates, you can provide them on the 'IR info' section as JLB, JSB, JSA, and JLA respectively. This decreases the processing time dramatically and would be useful in cases where for example the IR regions are not identical.", 
                                                                 br(), br(), "A minimal example of a annotation input can be like", 
                                                                 downloadLink('downloadDogma', 'this'), "file."),
                                                               
                                                               numericInput("nManual", "Number of GB files (min 1, max 20):", 2, min = 1, max = 20),
                                                               uiOutput("dynUIManual"), # columns displayed depending on the numericInput above
                                                               column(width=10, actionButton("Go", "Submit"), br(), br(), p("After submission the section below will turn innactive, meaning that analysis is in process. Be patient until you are prompted to download your plot."), verbatimTextOutput("Status", placeholder = TRUE), br(), verbatimTextOutput("Status2", placeholder = FALSE), downloadButton("Down","Download file")))
                                                    )),
                                           tabPanel("FAQ", br(), column(width = 1, ""), br(), column(width = 6, h4(strong("Q:"), tags$em(strong("Why the gene names are not displayed properly?"))), p(strong("A:"),"Check your input file(s), IRscope is based on the gene names specified in them."), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("Why some of the genes are extended to the large portion of the tracks?"))), p(strong("A:"),"This is due to the bad annotation implemented in the input file(s). Often this is related to the genes with introns or trans-spliced genes like rps12. Check and fix your input files."), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("I have a well annotated genome but not getting the plot?"))), p(strong("A:"),"IRscope can handle the redundant base pairs but if the sequence is of a very poor quality the program may fail in detecting the inverted regions and hence not produce any output."), br(),
                                                                                                     # h4(strong("Q:"), tags$em(strong("I see genes overplotting on each other, how can I fix this?"))), p(strong("A:"),"This is mostly because at least one of your species differs significantly from the rest so that it creates too large radius for finding the genes in the vicinity of the junction site and consequently causes the populated genes in these respective areas. Try reducing the species sampling to a smaller group of more closely related species."), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("There is a warning about the spacing around the junction not being in scale. Why?"))), p(strong("A:"),"If at least one of your species differs significantly from the rest so that it creates too large radius for finding the genes in the vicinity of the junction site, the radius in that junction will adapt to each of the species best, instead of using the same radius. Try reducing the species sampling to a smaller group of more closely related species."), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("Can I modify and run the codes on my own computer?"))), p(strong("A:"),"Yes, please download the file in the 'About' section, run the whole script in your R session. This will set up a pseudo web interface on your computer which resemles and works exactly like the online version. You can further tune the functions as you wish." ), br(),
                                                                                                     # h4(strong("Q:"), tags$em(strong("My analysis runs slow. Can I enhance its speed?"))), p(strong("A:"),"Yes, after downloading the file in the 'About' section, alter the 'parallel' arguments of the function 'IRinfo' into 'TRUE'. Then run the whole code on your R console. This will launch an instance with default 4 CPUs. If your machine is equipped with more cores, you can increase the number of CPUs with the 'nCPU' in 'p.d' subfunction inside the 'IRinfo' function."), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("What should be the format of the annotation files in the 'Manual Files' section?"))), p(strong("A:"),"This should be a plain text format (.txt) of the four tab separated columns, as start of the gene, end of the gene, its name, and the direction of it (+ or -). Note that the file needs to be without header." ), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("My species names are not plotted correctly, why?"))), p(strong("A:"),"The genome names in the 'GB Files' are read from the 'Organism' line of the file while the first space separated text of the genomes first line in the manual part is a determinant of species name in this section. In the 'Manual Files' section, the names are read from the first line after the '>' of the fasta file uploaded. Please edit them if they do not follow your expectations and rerun." ), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("My analysis took excessively long time and even after that the plot is of an overall low quality, why is that?"))), p(strong("A:"),"The IRscope is primarily optimized for the angiosperms and it often turns reliable results for other seed plants as well. However, further departure from the Embryophyta will extend this program to its limits. In such cases, you may want to consider use of the manual section with providing the the IR coordinates." ), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("How can I cite the program"))), p(strong("A:"), "The paper describing this web app in detail with the title", tags$em("IRscope: An online program to visualize the junction sites of chloroplast genomes"), "is", tags$a(href="https://doi.org/10.1093/bioinformatics/bty220", "published"), "in", tags$em("Bioinformatics"), "journal. Please, cite the paper when you use the program."), br()
                                                                                                     
                                           )),
                                           tabPanel("About", br(), column(width = 6, p("IRscope and its all dependencies are coded in R. The detailed instruction on how to use each subfuction is provided as well. Also the accession numbers that used to test and train the program are listed at the end of the file.", "Click", downloadLink('downloadDat', 'here'), "to download the file after which, you can run all the codes (apart from accession numbers) into your R Console to obtain this app on your PC. Here are two example outputs of the program:" , downloadLink('downloadPIC1', 'example1'), "and", downloadLink('downloadPIC2', 'example2'), "."))),
                                           tabPanel("Contact", br(), column(width = 6, p("The paper describing this web app in detail with the title", tags$em("IRscope: An online program to visualize the junction sites of chloroplast genomes"), "is", tags$a(href="https://doi.org/10.1093/bioinformatics/bty220", "published"),"in", tags$em("Bioinformatics"), "journal. Please cite that papers in all your correspondence. In case of practical issues please consult the FAQ sections first or email us at", tags$em("ali(dot)amiryousefi@helsinki(dot)fi."))))
                                         )
  )
  ))

server <- function(input, output, session) {
  # # TEST TODO #########################
  # pos <- 0L
  # 
  # # Returns a hex color string, e.g. "#FF0073"
  # nextColor <- function() {
  #   # Choose the next color, wrapping around to the start if necessary
  #   pos <<- (pos %% length(colors)) + 1L
  #   colors[[pos]]
  # }
  # 
  # observe({
  #   # Send the next color to the browser
  #   session$sendCustomMessage("background-color", nextColor())
  #   
  #   # Update the color every 100 milliseconds
  #   invalidateLater(100)
  # })
  # # TEST #

  ###Input validator section########################################
  iv <- InputValidator$new()
  iv$add_rule("nGB", sv_between(1, 20))
  iv$add_rule("nManual", sv_between(1, 20))
  iv$enable()
  
  ###Dynamic UI section#############################################
  # TODO guardar proceso si eso (nombres de gb y asi, se renderean de nuevo cada vez que sale el numero)
  output$dynUIGB <- renderUI({
    if (is.numeric(input$nGB) & input$nGB >= 1 & input$nGB <= 20) {
      row_idx <- seq_len(input$nGB)
      row_idx %>% map(~ rowInputGB(.x))
    }
  })
  
  output$dynUIManual <- renderUI({
    if (is.numeric(input$nManual) & input$nManual >= 1 & input$nManual <= 20) {
      row_idx <- seq_len(input$nManual)
      row_idx %>% map(~ rowInputManual(.x))
    }
  })
  
  ###Remove file section TODO commented bc it is not working ################
  # values <- reactiveValues(
  #   upload_state = NULL
  # )
  # 
  # listenRemoveFile <- reactive({
  #   removeButtonsInput <- list()
  #   for (i in seq_len(input$nGB)){
  #     if(is.null(input[[paste0('removefile', i)]])){
  #       removeButtonsInput <- append(removeButtonsInput, 0)
  #     } else {
  #       removeButtonsInput <- append(removeButtonsInput, input[[paste0('removefile', i)]])
  #     }
  #   }
  #   removeButtonsInput
  # })
  # 
  # observeEvent(listenRemoveFile(), {
  #   for (i in seq_len(input$nGB)){
  #     
  #     if(!is.null(input[[paste0('removefile', i)]])) {
  #       if(input[[paste0('removefile', i)]] == 1){
  #         shinyjs::reset(paste0("data", i))
  #         input[[paste0("data", i)]] <- NULL
  #       }
  #     }
  #   }
  # })
  # 
  # observeEvent(input$removefile, {
  #   shinyjs::reset("gb_file")
  #   values$upload_state <- 'reset'
  # })
  
  
  
  ###GB File section#############################################
  ###Status update for the GB file submission
  output$StatusGB <- renderText({
    if(input$Gbsub){
      paste("Thank you for your Manual submission at", paste0(date(), "!"),
            "\nYou may now download your plot via the tab below.", sep= " ")
    }
  })
  
  
  ####Calculating the IR values with IRs in a pseudo text function in the GB file section
  output$StatusGB2 <- renderText({
    if(input$Gbsub){
      tryCatch({
        rm(list = c('FasList', 'IRList', 'GeneList', 'spnames', 'l', 'nuclw'))
      })
      
      gbfiles <- isolate(gb())
      l <- length(gbfiles)
      
      progress <- Progress$new(session, min=0, max=l)
      on.exit(progress$close())
      progress$set(message = 'Calculation in progress',
                   detail = paste0('This may take a while... 0/', l))
      
      dist<- isolate(IRs(gbfiles, SRev(), progress))
      return(dist)
    }
  })
  
  # TODO
  # output$Plot <- renderPlot({
  #   if(input$Gbsub){
  #     showNotification("Generating plot...", duration = NULL, id = "message")
  #     IRs2()
  #     showNotification("Plot generated", duration = 2, id = "message")
  #   }
  # })
  
  
  ###Downloading the result with the Download botton in the GB files section
  output$downloadData <- downloadHandler(
    filename = function(){
      return(paste('IR', input$file_type, sep="."))
    },
    content = function(file){
      showNotification("Downloading plot...", duration = NULL, id = "message")
      
      # dev.new(width, height)
      
      if (input$file_type  == "pdf") {
        grDevices::pdf(file, width=8.3, height=(l+2)*8.3/12)
      } else {
        do.call(input$file_type, args=list(filename=file,
                                           width = 8.3, height = (l+2)*8.3/12,
                                           units = "in", res = 300))
        
      }
      IRs2(file=file)
      # if(IRs2(file=file)){
      #   shinyalert("Warning!", "The spacing around the junction is not in scale.
      #              Refer to the FAQ to see why.", type = "warning")
      # }
      while (!is.null(dev.list())) dev.off()
      removeNotification(id = "message")
    }
  )
  
  ###Reactive list of the objects needed for the calculation of this section i.e. GB files
  gb <- reactive({
    gbFiles <- list()
    
    for (i in seq_len(input$nGB)) {
      acci <- input[[paste0("acc", i)]]
      if(acci != "") {
        gbFiles[[i]] <- list(gbfile = acci,
                             local.file = FALSE)
      }
    }
    
    for (i in seq_len(input$nGB)) {
      datai <- input[[paste0("data", i)]]
      inFile <- datai
      if (is.data.frame(inFile)){
        gbFiles[[i]] <- list(gbfile = inFile$datapath,
                             local.file = TRUE)
      }
    }
    
    gbFiles[which(gbFiles!="NULL")]
  })
  
  SRev <- reactive({
    SFiles<- list()
    
    for (i in seq_len(input$nGB)) {
      if (input[[paste0("acc", i)]] != "" || is.data.frame(input[[paste0("data", i)]])) {
        SFiles[[i]]<- input[[paste0("S", i)]]
      }
    }
    SFiles[which(SFiles!="NULL")]
  })
  
  
  
  ###Dogma File section#############################################
  
  ###Status file for the submission of the files
  output$Status<- renderText({
    if(input$Go){
      paste("Thank you for your Manual submission at", paste0(date(), "!"),
            "You may now download your plot via the tab below.", sep= " ")
    }
  })
  
  
  
  ###Calculating the IR with the IRsD in a pseudo text function in the Dogma section
  output$Status2<- renderText({
    if(input$Go){
      dist <- isolate(IRsD(df(), ff(),  ii(), nf(), file=file))
      dist
    }
  })
  
  ###Downloading the result with the Download botton in the DOGMA files section
  output$Down <- downloadHandler(
    filename = "IR.pdf",
    content = function(file) {
      #mboard(datagb(), file = file)
      isolate(IRsD2(df(), file=file))
      #jpeg(file,  width=8.3, height=8.9, units="in", res=100)
      #plot(rnorm(100))
      #dev.off()
    }
  )
  
  
  ###Reactive list of objects for the calculation of this section
  ###df for dogma files
  df <- reactive({
    dFiles<- list()
    
    for (i in seq_len(input$nManual)) {
      inFilei <- input[[paste0("dogma", i)]]
      if (is.data.frame(inFilei)) {
        dFiles[[i]]<-GnlBuilder(trnDogma(read.table(inFile$datapath)))
      }
    }
    
    dFiles[which(dFiles!="NULL")]
  })
  
  ###ff for fasta files
  ff <- reactive({
    fFiles<- list()
    
    for (i in seq_len(input$nManual)) {
      inFilei <- input[[paste0("fas", i)]]
      if (is.data.frame(inFilei)) {
        fFiles[[i]]<-rdnFixerD(as.character(read.fasta(inFile$datapath)[[1]]))
      }
    }
    
    fFiles[which(fFiles!="NULL")]
  })
  
  ###ii for IR info
  ii <- reactive({
    iFiles<- list()
    
    for (i in seq_len(input$nManual)) {
      inFilei <- input[[paste0("fas", i)]]
      
      if (input[[paste0("IR", i)]] == "" && is.data.frame(inFilei)) {
        iFiles[[i]]<- 999
      }
      else if(input[[paste0("IR", i)]] != ""){
        iFiles[[i]]<- c(as.numeric(unlist(strsplit(input[[paste0("IR", i)]], split = ","))) , length(as.character(read.fasta(inFilei$datapath)[[1]])))
      }
    }
    
    iFiles[which(iFiles!="NULL")]
  })
  
  ###nf for the names of the species
  nf<- reactive({
    nFiles<- list()
    
    for (i in seq_len(input$nManual)) {
      inFilei <- input[[paste0("fas", i)]]
      if (is.data.frame(inFilei)) {
        nFiles[[i]]<-gsub("_", " ", names(read.fasta(inFilei$datapath)))
      }
    }
    
    nFiles[which(nFiles!="NULL")]
  })
  
  
  
  #output$markdown <- renderUI({
  #  HTML(markdown::markdownToHTML(knit('Rmarkdown.Rmd', quiet = TRUE)))
  #})
  
  
  # doog <- readLines("dogma.txt")
  # 
  # output$downloadDogma <- downloadHandler(###Downloading the Dogma file
  #   filename = function() {
  #     paste("annotation-", Sys.Date(), ".txt", sep="")
  #   },
  #   content = function(file) {
  #     writeLines(doog, file)
  #   }
  # )
  # 
  # pic1 <- readJPEG("IR1.jpg")
  # pic2 <- readJPEG("IRprevioustest.jpg")
  # 
  # output$downloadPIC1 <- downloadHandler(###Downloading the first PIC
  #   filename = function() {
  #     paste("IR1example-", Sys.Date(), ".jpg", sep="")
  #   },
  #   content = function(file) {
  #     writeJPEG(pic1, file)
  #   }
  # )
  
  ###Downloading the first PIC
  output$downloadPIC2 <- downloadHandler(
    filename = function() {
      paste("IR2example-", Sys.Date(), ".jpg", sep="")
    },
    content = function(file) {
      writeJPEG(pic2, file)
    }
  )
  
  
  ###Downloading the R codes
  ### TODO IRscope.R es parte no todo
  data <- readLines("IRscope.R")
  
  output$downloadDat <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".R", sep="")
    },
    content = function(file) {
      writeLines(data, file)
    }
  )
  
  
  ### Updating the section areas of the browsed files in GB section
  
  ## TODO ahora 1:20 pero mejor seq_len(input$nGB) pero pide reactive
  lapply(1:20, function(i) {
    output[[paste0("sp",i)]] <- renderText({
      if (input[[paste0("acc",i)]] != "") {
        print(head(fetch.gb(input[[paste0("acc",i)]]))[2])
      }
      else {
        frame <- input[[paste0("data",i)]]
        if (is.data.frame(frame)) {
          print(readLines(frame$datapath)[2])
        }
      }
    })
  })
  
}

shinyApp(ui, server)