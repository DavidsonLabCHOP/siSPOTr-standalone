library(shiny)
library(DT)
library(shinythemes)
library(Biostrings)
library(data.table)
library(stringr)
library(shinyjs)
library(writexl)
library(shinyalert)

# UI
ui <- fluidPage(
  shinyalert::useShinyalert(),  # Initialize shinyalert
  theme = shinytheme("flatly"),  # Sleek modern theme
  titlePanel(
    div(
      tags$img(src = "logo.png", style = "max-width: 75%; height: auto;", alt = "Logo"),
      h1("siSPOTr: CRISPR Guide RNA Analysis", style = "text-align: center; margin-top: 10px;")
    )
  ),
  navbarPage(
    "I want to:",
    id = "nav",
    
    # Tab 1: Introduction
    tabPanel(
      title = tagList(icon("info-circle"), "Learn about siSPOTr"),
      mainPanel(
        div(
          # Paragraph Above the Image
          p("RNA interference (RNAi) serves as a powerful and widely used gene silencing tool for basic biological research and is being developed as a therapeutic avenue to suppress disease-causing genes. However, the specificity and safety of RNAi strategies remains under scrutiny because small inhibitory RNAs (siRNAs) induce off-target silencing. Currently, the tools available for designing siRNAs are biased toward efficacy as opposed to specificity. Prior work from our laboratory and others’ supports the potential to design highly specific siRNAs by limiting the promiscuity of their seed sequences (positions 2–8 of the small RNA), the primary determinant of off-targeting. Here, a bioinformatic approach to predict off-targeting potentials was established using publically available siRNA data from more than 50 microarray experiments. With this, we developed a specificity-focused siRNA design algorithm and accompanying online tool which, upon validation, identifies candidate sequences with minimal off-targeting potentials and potent silencing capacities. This tool offers researchers unique functionality and output compared with currently available siRNA design programs. Furthermore, this approach can greatly improve genome-wide RNAi libraries and, most notably, provides the only broadly applicable means to limit off-targeting from RNAi expression vectors."),
          
          # Image
          tags$img(
            src = "si-intro2.png",
            style = "max-width: 75%; height: auto;",
            alt = "Guide RNA Process Schematic"
          ),
          
          # Paragraph Below the Image
          p("Diagram of on- and off-target silencing by siRNAs. (A) Cartoon depicting a siRNA duplex designed to exhibit proper strand-biasing [i.e. strong G-C (blue) and weak A/G-U (red) binding at the respective 5′ and 3′ ends of the sense strand] and contain a low off-targeting potential seed (green highlight). Upon loading into RISC, the antisense strand may direct on-target silencing (intended) and off-target silencing (unintended). (B) Schematic highlighting the relationship between the frequencies of seed complement binding sites in the 3′-UTRome and the off-targeting potential for siRNAs."),
          
          # Citations Section
          div(
            h4("Resources"),
            tags$ul(
              tags$li(
                "Boudreau, Ryan L et al. “siSPOTR: a tool for designing highly specific and potent siRNAs for human and mouse.",
                tags$a(href = "https://pmc.ncbi.nlm.nih.gov/articles/PMC3592398/", "Read more", target = "_blank")
              ),
              tags$li(
                "To see more exciting work at the Davidson Lab ",
                tags$a(href = "https://www.thedavidsonlab.com/", "click here!", target = "_blank")
              )
            )
          )
        )
      )
    ),
    
    
    # Tab 2: Guide RNA Finder
    tabPanel(
      title = tagList(icon("dna"), "Find Best Guide RNAs"),
      sidebarLayout(
        sidebarPanel(
          h4("Input Options"),
          selectInput(
            "species", "Select Organism",
            choices = c("human", "mouse"),
            selected = "human"
          ),
          numericInput(
            "gsize", "Guide RNA Size (nt)",
            value = 22,
            min = 19, max = 22,
            step = 1
          ),
          fileInput(
            "fastaFile1", "Upload FASTA File",
            accept = c(".fasta", ".fa")
          ),
          actionButton(
            "runButton1", "Run Analysis", 
            class = "btn-primary"
          )
        ),
        mainPanel(
          h4("Results"),
          dataTableOutput("resultTable1"),
          downloadButton("downloadTable1", "Download Results as Excel"),
          actionButton("infoTable1", "Explain Column Headers", icon = icon("info-circle"))
        )
      )
    ),
    
    # Tab 3: Off-Target Assessment
    tabPanel(
      title = tagList(icon("crosshairs"), "Assess Off-Target Risk"),
      sidebarLayout(
        sidebarPanel(
          h4("Input Options"),
          selectInput(
            "species2", "Select Organism",
            choices = c("human", "mouse"),
            selected = "human"
          ),
          checkboxInput(
            "revcomp", "These are target sites (anti-sense)",
            value = TRUE
          ),
          fileInput(
            "fastaFile2", "Upload FASTA File",
            accept = c(".fasta", ".fa")
          ),
          actionButton(
            "runButton2", "Run Off-Target Analysis", 
            class = "btn-warning"
          )
        ),
        mainPanel(
          h4("Results"),
          dataTableOutput("resultTable2"),
          downloadButton("downloadTable2", "Download Results as Excel"),
          actionButton("infoTable2", "Explain Column Headers", icon = icon("info-circle")),
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  
  
  #Function 1 takes gene sequence and spits out guide RNAs
  processFASTA1 <- function(fastaFile1) {
    target.mRNA <- readDNAStringSet(fastaFile1$datapath)
    names(target.mRNA) <- make.names(names(target.mRNA), unique = TRUE)
    POTS.file <- "POTS.txt"
    POTS.dt <- fread(POTS.file)
    guide.size <- input$gsize #22 # Probably limit to 19-22. Or warn if different.
    target.species <- input$species # For seed off-target POTS score. Only Human and Mouse right now.
    # allow.seed.gu.wobble <- TRUE 
    gc.min <- 0.2
    gc.max <- 0.7
    n.gc.pass <- 2 # Hard-code this for interface. 
    is.pol3 <- TRUE # Will need to do the pol3 check once guides are put into scaffolds (to-do)
    
    # Functions ----
    # Splits target mRNA into guide-sized windows
    get.all.target.sites.f <- function(seqs, size=22){
      seq.list <- lapply(seq_along(seqs),
                         function(i){
                           t.seq <- seqs[[i]]
                           seq.len <- nchar(t.seq)
                           nm <- names(seqs)[[i]]
                           v <- trim(Views(t.seq, start = seq(seq.len), width = size))
                           v <- v[width(v)==size]
                           names(v) <- paste0(nm, "_p_", start(v), "_", end(v))
                           return(RNAStringSet(v))
                         })
      do.call(c, seq.list)
    }
    # Takes the target sites and generates appropriate guide and passenger,injecting G:U wobbles for biasing, where appropriate.
    # Also adds, POTS, gc% and flag for passing filters.
    characterize.and.bias.guides.f <- function(target.seq, POTS=POTS.dt, n.pass.gc=n.gc.pass, species=target.species){
      # Help bias guides first
      guides.rna <- reverseComplement(target.seq)
      guide.names <- names(target.seq)
      guides.rna.right <- as.character(subseq(guides.rna, start=2, end=width(guides.rna)))
      guides.rna.1 <- subseq(guides.rna, start=1, width=1) %>% ifelse(.=="C", "u", .)
      guide.rna.out <- paste0(guides.rna.1, guides.rna.right)
      # Now passenger across from the 5' end of guide
      passenger.rna <- reverseComplement(subseq(guides.rna, 1, width=guide.size-2))
      pass.width <- width(passenger.rna)
      passenger.rna.guideP1P2 <- as.character(subseq(passenger.rna, 
                                                     start=pass.width-1, 
                                                     end=pass.width)) %>% gsub("C", "u", .)
      passenger.rna.out <- paste0(subseq(passenger.rna, 1, pass.width-2),
                                  passenger.rna.guideP1P2,
                                  "uu")
      dt <- data.table(seq.id=guide.names,
                       target.seq=as.character(target.seq),
                       guide.seq=guide.rna.out,
                       pass.seq=passenger.rna.out)
      # Calc GC%, and other features from user-specified 
      pass.5p.search <- ifelse(n.pass.gc<2, 2, n.pass.gc) # Will search the first 2 bases of the passenger strand for GC, unless user specifies more than 2
      dt[, `:=`(Seed=subseq(guide.seq, start=2, end=8),
                GC=stringr::str_count(guide.seq, "[GC]")/nchar(guide.seq),
                passenger.bias=stringr::str_count(subseq(pass.seq, 1, width=pass.5p.search), "[GC]"))]
      dt.merged <- merge.data.table(dt, POTS, by="Seed")
      dt.merged[, c("target.name", "target.pos.tmp"):=tstrsplit(seq.id, split="_p_")]
      dt.merged[, c("start.target", "end.target"):=tstrsplit(target.pos.tmp, split="_", type.convert = TRUE)]
      
      output.order <- c("seq.id", "guide.seq", "pass.seq", "target.name", "start.target", "end.target", "GC", "passenger.bias", "POTS.HSA", "POTS.MMU", "SPS")
      # Clean and output
      if(tolower(species)=="mouse" | tolower(species)=="mmu"){
        setorder(dt.merged, POTS.MMU)
      } else{
        setorder(dt.merged, POTS.HSA)
      }
      return(dt.merged[, .SD, .SDcols=output.order])
    }
    
    data1 <- get.all.target.sites.f(target.mRNA) %>% 
      characterize.and.bias.guides.f(target.seq=.)
    
    # Pass target.sites all to shiny. Will set filters here, but may want to be reactive ----
    
    #target.sites.all[between(GC, gc.min, gc.max) & passenger.bias>=n.gc.pass]
    return(data1)
    
  }
  # Function 2 takes guide RNAs and assesses off-target potential
  processFASTA2 <- function(fastaFile2) {
    POTS.file <- "POTS.txt"
    POTS.dt <- fread(POTS.file)
    # User-supplied:
    is.guide <- TRUE # Have user select if the inputs are Guides (Antisense) or target sites
    # Guide input:
    # THis is a user-interface thing. I'd recommend allowing file upload, or
    # direct pasting of 22-mers. Pasting could be one per line, or fasta
    # For now. Hard-coded and example set of guides
    guide.fl <- "data/test_guides.fa" 
    guide.in <- readBStringSet(fastaFile2$datapath)
    if(is.guide==FALSE){
      guides <- RNAStringSet(reverseComplement(guide.in))
    } else{
      guides <- RNAStringSet(guide.in)
    }
    if(length(names(guides))==0){
      names(guides) <- paste0("Guide_", seq_along(guides))
    }
    seeds <- subseq(guides, start = 2, end = 8)
    guides.dt <- data.table(seq.id=names(guides),
                            seq.in=as.character(guide.in),
                            guide.seq=as.character(guides),
                            Seed=as.character(seeds))
    data2 <- merge.data.table(guides.dt, POTS.dt, on="Seed", all.x=TRUE)
    #guide.dt.annot
    return(data2)
  }
  
  # Store the processed data as reactive values
  data1 <- reactiveValues()
  data2 <- reactiveValues()
  
  # Event handler for Analysis 1 Run button
  observeEvent(input$runButton1, {
    if (!is.null(input$fastaFile1)) {
      data1$table <- processFASTA1(input$fastaFile1)
    }
  })
  
  # Event handler for Analysis 2 Run button
  observeEvent(input$runButton2, {
    if (!is.null(input$fastaFile2)) {
      data2$table <- processFASTA2(input$fastaFile2)
    }
  })
  
  # Render the Analysis 1 data table in the UI
  output$resultTable1 <- renderDataTable({
    if (!is.null(data1$table)) {
      data1$table
    }
  })
  
  # Render the Analysis 2 data table in the UI
  output$resultTable2 <- renderDataTable({
    if (!is.null(data2$table)) {
      data2$table
    }
  })
  # Add download handler for Analysis 1
  output$downloadTable1 <- downloadHandler(
    filename = function() {
      paste("Guide_RNAs_", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      if (!is.null(data1$table)) {
        writexl::write_xlsx(data1$table, path = file)
      }
    }
  )
  
  # Add download handler for Analysis 2
  output$downloadTable2 <- downloadHandler(
    filename = function() {
      paste("Off_Target_Risk_", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      if (!is.null(data2$table)) {
        writexl::write_xlsx(data2$table, path = file)
      }
    }
  )
  
  # Column explanations for Tab 2
  observeEvent(input$infoTable1, {
    shinyalert(
      title = "Column Headers Explanation (T2):",
      text = paste(
        "<b>seq.id:</b> Sequence Identifier: A unique identifier for each sequence.<br>",
        "<b>guide.seq:</b> Guide Sequence: The RNA sequence used for targeting.<br>",
        "<b>pass.seq:</b> Passenger Sequence: The complementary sequence paired with the guide.<br>",
        "<b>GC:</b> GC Content: The percentage of guanine and cytosine in the sequence.<br>",
        "<b>POTS.HSA:</b> POTS Score (Human): Predicted off-target score for humans.<br>",
        "<b>POTS.MMU:</b> POTS Score (Mouse): Predicted off-target score for mice.",
        sep = ""
      ),
      html = TRUE,
      type = "info"
    )
  })
  
  # Column explanations for Tab 3
  observeEvent(input$infoTable2, {
    shinyalert(
      title = "Column Headers Explanation (T3):",
      text = paste(
        "<b>seq.id:</b> Sequence Identifier: A unique identifier for each sequence.<br>",
        "<b>guide.seq:</b> Guide Sequence: The RNA sequence used for targeting.<br>",
        "<b>pass.seq:</b> Passenger Sequence: The complementary sequence paired with the guide.<br>",
        "<b>GC:</b> GC Content: The percentage of guanine and cytosine in the sequence.<br>",
        "<b>POTS.HSA:</b> POTS Score (Human): Predicted off-target score for humans.<br>",
        "<b>POTS.MMU:</b> POTS Score (Mouse): Predicted off-target score for mice.",
        sep = ""
      ),
      html = TRUE,
      type = "info"
    )
  })
  
}

# Run the Shiny app
shinyApp(ui, server)
