

###############################################################################################################################
#App
###############################################################################################################################

library(shiny)
library(shinyWidgets)
library(shinycssloaders)



ui <- fluidPage(
  
  useSweetAlert(),
  
  navbarPage("EXIS", id = "navbar",
             ################################################################################################################################
             #get GEO
             ################################################################################################################################
             
             tabPanel("GEO accession", value = "panel1",
                      
                      
                      
                      mainPanel(
                        
                        textInput(inputId = "getGEO", 
                                  label = "Enter GEO accession", 
                                  value = "GSE36980"), 
                        
                        
                        actionBttn(inputId = "downloadGEO", 
                                   label = "Next",
                                   style = "jelly",
                                   color = "primary"),
                        
                        actionBttn(inputId = "infopanel1", 
                                   label = "Help",
                                   style = "jelly",
                                   color = "success")
                        
                        
                        
                      )
                      
                      
                      
             ),
             
             
             ################################################################################################################################
             #get meta data
             ################################################################################################################################
             
             tabPanel("Meta data", value = "panel2",
                      
                      
                      
                      sidebarPanel(
                        
                        
                        awesomeCheckbox(inputId = "dependentselect", 
                                        label = "Dependent samples", 
                                        value = FALSE),
                        
                        uiOutput("groups"),
                        
                        uiOutput("pairs"),
                        
                        actionBttn(inputId = "meta.ok", 
                                   label = "Next",
                                   style = "jelly",
                                   color = "primary")
                        
                      ),
                      
                      mainPanel(
                        downloadButton("downloadmeta", "Download table"),
                        
                        tableOutput(outputId = "grouping") %>% withSpinner(color="#0dc5c1")
                      )
                      
                      
                      
             ),
             
             ################################################################################################################################
             #get isoform expression
             ################################################################################################################################
             
             tabPanel("Isoform expression", value = "panel3",
                      
                      sidebarPanel(
                        
                        selectInput(inputId = "brainarray",
                                    label = "Probe set definition",
                                    choices = c("Exon-specific probe sets" = "exon", "Transcript-specific probe sets" = "transcript"),
                                    selected = "transcript"),
                        
                        uiOutput("robustout"),
                        
                        uiOutput("chipout"),
                        
                        uiOutput("organismout"),
                        
                        uiOutput("annotationout"),
                        
                        uiOutput("versionout"),
                        
                        awesomeCheckbox(inputId = "outlier",
                                        label = "Keep all samples",
                                        value = TRUE,
                                        status = "danger"),
                        
                        uiOutput("outliersout"),
                        
                        awesomeCheckbox(
                          inputId = "backgroundcor",
                          label = "Background correction", 
                          value = TRUE),
                        
                        awesomeCheckbox(
                          inputId = "quantilenorm",
                          label = "Quantile normalization", 
                          value = TRUE),
                        
                        actionBttn(inputId = "ann.ok",
                                   label = "Calculate",
                                   style = "fill",
                                   color = "danger"),
                        
                        
                        uiOutput("proceedann")
                        
                      ),
                      
                      mainPanel(
                        
                        downloadButton("downloadexpr", "Download expression data"),
                        
                        plotOutput("boxplot") %>% withSpinner(color="#0dc5c1"),
                        
                        plotOutput("normboxplot") %>% withSpinner(color="#0dc5c1"),
                        
                        plotOutput("hist") %>% withSpinner(color="#0dc5c1"),
                        
                        plotOutput("normhist") %>% withSpinner(color="#0dc5c1")
                        
                      )
                      
             ),
             
             
             
             ################################################################################################################################
             #PCA
             ################################################################################################################################
             
             tabPanel("PCA", value = "panel4",
                      
                      sidebarPanel(
                        selectInput(inputId = "hpca", 
                                    label = "Horizontal axis",
                                    choices = c("PC1","PC2","PC3", "PC4", "PC5", "PC6", "PC7", "PC8"),
                                    selected = "PC1"),
                        
                        selectInput(inputId = "vpca", 
                                    label = "Vertical axis",
                                    choices = c("PC1","PC2","PC3", "PC4", "PC5", "PC6", "PC7", "PC8"),
                                    selected = "PC2"),
                        
                        
                        actionBttn(inputId = "pca.ok", 
                                   label = "Next",
                                   style = "jelly",
                                   color = "primary")
                      ),
                      
                      mainPanel(
                        plotlyOutput("pca") %>% withSpinner(color="#0dc5c1")
                      )
                      
                      
             ),
             
             ################################################################################################################################
             #differential isoform expression
             ################################################################################################################################
             
             tabPanel("Differential expression analysis", value = "panel5",
                      
                      sidebarPanel(
                        
                        awesomeCheckbox(inputId = "allexonsin",
                                        label = "Genome-wide analysis",
                                        value = TRUE),
                        
                        uiOutput("geneorlistout"),
                        
                        uiOutput("geneselectionout"), 
                        
                        uiOutput("exonselectionout"),
                        
                        uiOutput("exonsonlyout"),
                        
                        actionBttn(inputId = "diffexpr.ok", 
                                   label = "Calculate", 
                                   style = "fill", 
                                   color = "danger"),
                        
                        uiOutput("proceedtovol")
                        
                      ),
                      
                      mainPanel(
                        uiOutput("comparisonout"),
                        
                        uiOutput("exonboxplotout"),
                        
                        uiOutput("mappingplotout"),
                        
                        dataTableOutput("finaltable") %>% withSpinner(color="#0dc5c1"),
                        
                        downloadButton("downloadfinal", "Download table"),
                        
                        uiOutput("uiexprboxplot"),
                        
                        plotlyOutput("probemapping") %>% withSpinner(color="#0dc5c1")
                      )
                      
             ),
             
             
             
             ################################################################################################################################
             #volcano plot
             ################################################################################################################################
             
             tabPanel("Volcano Plot", value = "panel6",
                      
                      sidebarPanel(
                        prettyRadioButtons(
                          inputId = "raworfdr",
                          label = NULL, 
                          choices = c("Raw P-value", "FDR")
                        ),
                        
                        numericInput(
                          inputId = "pthreshold",
                          label = "P threshold",
                          value = 0.05
                        ),
                        
                        numericInput(
                          inputId = "logfcthreshold",
                          label = "logFC threshold",
                          value = 1
                        )
                        
                      ),
                      
                      mainPanel(
                        uiOutput("comparisonout1"),
                        
                        plotlyOutput("volcano") %>% withSpinner(color="#0dc5c1"),
                        
                        dataTableOutput("voltable") %>% withSpinner(color="#0dc5c1")
                      )
             )
             
  )
)








server <- function(input, output, session){
  
  ################################################################################################################################
  #get GEO
  ################################################################################################################################
  
  hideTab("navbar", target = "panel2")
  hideTab("navbar", target = "panel3")
  hideTab("navbar", target = "panel4")
  hideTab("navbar", target = "panel5")
  hideTab("navbar", target = "panel5")
  hideTab("navbar", target = "panel6")
  
  observeEvent(input$downloadGEO, {
    showTab("navbar", target = "panel2")
  })
  
  
  #download data
  
  gset <- eventReactive(input$downloadGEO, {
    withProgress(message = "Downloading data.....", value = 0, {
      getGEO(input$getGEO, GSEMatrix =TRUE, getGPL = FALSE)
    })
  })
  
  
  observeEvent(input$infopanel1, {
    sendSweetAlert(
      session = session,
      title = "Information",
      text = "Enter the GEO accession number. Find your number at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi",
      type = "info"
    )
  })
  
  
  #go to next tab
  
  observeEvent(input$downloadGEO,{
    updateNavbarPage(session, "navbar",
                     selected = "panel2")
  })
  
  
  
  
  
  
  
  ################################################################################################################################
  #get meta
  ################################################################################################################################
  
  
  #select pairing columns
  
  output$groups <- renderUI({
    
    checkboxGroupInput(inputId = "groupselect", 
                       label = "Select grouping variable(s)", 
                       choices = colnames(get_grouping(gset())),
                       selected = auto_group(gset()))
    
    
  })
  
  
  
  output$pairs <- renderUI({
    
    if (length(input$dependentselect) > 0){
      if (input$dependentselect == TRUE) {
        checkboxGroupInput(inputId = "pairselect", 
                           label = "Select pairing variable(s)", 
                           choices = colnames(get_grouping(gset())))
      } 
    }
    
  })
  
  
  
  #get meta data
  meta <- reactive({
    
    
    if (length(input$dependentselect) > 0){
      if (input$dependentselect == FALSE) {
        groups <- get_grouping(gset())
        meta <- get_meta(gset = gset(), grouping_column = groups[,input$groupselect])
        return(meta)
        
      }
    }
    
    
    
    
    if (length(input$dependentselect) > 0){
      if (input$dependentselect == TRUE & (length(input$groupselect) > 0)) {
        groups <- get_grouping(gset())
        meta <- get_meta(gset = gset(), grouping_column = groups[,input$groupselect], pairing_column = groups[,input$pairselect])
        return(meta)
        
      }
    }
    
  })
  
  
  #make table of meta data
  output$grouping <- renderTable({
    
    meta()
    
  })
  
  
  #download table of meta data
  output$downloadmeta <- downloadHandler(
    filename = "meta data",
    content = function(file){
      write.table(meta(), file)
    }
  )
  
  
  observeEvent(input$meta.ok, {
    updateNavbarPage(session, "navbar",
                     selected = "panel3")
    
    
  })
  
  
  observeEvent(input$meta.ok, {
    showTab("navbar", target = "panel3")
  })
  
  
  ################################################################################################################################
  #get isoform expression
  ################################################################################################################################
  
  #change input
  
  output$chipout <- renderUI({
    
    if (length(input$brainarray) > 0){
      if (input$brainarray == "exon" | input$brainarray == "transcript") {
        selectInput(inputId = "chipin",
                    label = "Chiptype",
                    choices = c("hugene10st", "huex10st"),
                    selected = get_chiptype(gset()))
      }
    }
    
  })
  
  
  output$organismout <- renderUI({
    
    if (length(input$brainarray) > 0){
      if (input$brainarray == "exon" | input$brainarray == "transcript") {
        selectInput(inputId = "organismin",
                    label = "Organism",
                    choices = "hs",
                    selected = get_organism(gset()))
      }
    }
    
  })
  
  
  output$versionout <- renderUI({
    
    if (length(input$brainarray) > 0){
      
      if (input$brainarray == "exon" | input$brainarray == "transcript"){
        selectInput(inputId = "versionin",
                    label = "Brainarray version",
                    choices = 25,
                    selected = 25)
      }
    }
  })
  
  output$robustout <- renderUI({
    
    if (length(input$brainarray) > 0){
      
      if (input$brainarray == "transcript"){
        awesomeCheckbox(inputId = "robustin",
                        label = "Robust probe sets only",
                        value = TRUE,
                        status = "success")
      }
    }
  })
  
  
  output$outliersout <- renderUI({
    
    if (length(input$outlier) > 0){
      
      if(input$outlier == FALSE){
        samples <- meta()[,1]
        pickerInput(inputId = "outliersin", 
                    label = "Select samples to be removed", 
                    choices = meta()[,1],
                    multiple = TRUE)
        
      }
    }
  })
  
  
  
  
  #Make reactive variables from input
  
  
  brainarray <- eventReactive(input$ann.ok, {
    input$brainarray
  })
  
  robust <- eventReactive(input$ann.ok, {
    if (length(input$robustin) > 0){
      if(input$robustin == TRUE){
        robust <- "robust"
      }
      if(input$robustin == FALSE){
        robust <- "all"
      }
      return(robust)
    }
  })
  
  
  
  chiptype <- eventReactive(input$ann.ok, {
    
    if(length(input$chipin) > 0){
      chiptype = input$chipin
    }
    
    if(length(input$chipin) < 1) {
      chiptype = NULL
    }
    
    return(chiptype)
  })
  
  
  
  organism <- eventReactive(input$ann.ok, {
    
    if(length(input$organismin) > 0){
      organism = input$organismin
    }
    
    if(length(input$chipin) < 1) {
      organism = NULL
    }
    
    return(organism)
  })
  
  
  
  
  version <- eventReactive(input$ann.ok, {
    if(length(input$versionin) > 0){
      version = input$versionin
    }
    
    if(length(input$versionin) < 1){
      version = 25
    }
    
    return(version)
  })
  
  
  
  outliers <- eventReactive(input$ann.ok, {
    if (input$outlier == TRUE){
      outliers = NULL
    }
    
    if (input$outlier == FALSE){
      if (length(input$outliersin) < 1){
        outliers = NULL
      }
    }
    
    if (input$outlier == FALSE){
      if (length(input$outliersin) > 0){
        outliers = input$outliersin
      }
    }
    
    return(outliers)
  })
  
  
  
  #Read CEL files
  data1 <- eventReactive(input$ann.ok, {
    
    withProgress(message = "Reading CEL files.....", value = 0, {
      
      if (brainarray() == "exon") {
        
        data1 <- readcels(gset = gset(), 
                          chiptype = chiptype(), 
                          organism = organism(), 
                          version = version(),
                          annotation = "ense",
                          robust = TRUE,
                          outliers = outliers())
      }
      
      if (brainarray() == "transcript") {
        
        if (robust() == "all") {
          data1 <- readcels(gset = gset(), 
                            chiptype = chiptype(), 
                            organism = organism(), 
                            version = version(),
                            annotation = "enst",
                            robust = FALSE,
                            outliers = outliers())
        }
        
        if (robust() == "robust") {
          data1 <- readcels(gset = gset(), 
                            chiptype = chiptype(), 
                            organism = organism(), 
                            version = version(),
                            annotation = "enst",
                            robust = TRUE,
                            outliers = outliers())
        }
        
      }
      return(data1)
    })
  })
  
  
  #get Exon-annotation file
  annotated <- eventReactive(input$ann.ok, {
    
    if (brainarray() == "transcript") {
      annotated <- NULL
    }
    
    
    if (brainarray() == "exon") {
      
      if (file.exists(paste(chiptype(), organism(), version(), "OfficialNewProbeSets2.txt", sep = "_"))) {
        annotated <- read.table(file = paste(chiptype(), organism(), version(), "OfficialNewProbeSets2.txt", sep = "_"), header = TRUE)
        
      }
      
      if (!file.exists(paste(chiptype(), organism(), version(), "OfficialNewProbeSets2.txt", sep = "_"))) {
        print("annotation file not found")
        
      }
      
    }
    
    
    return(annotated)
  })
  
  
  
  
  #RMA background correction, normalization, and summarization
  data.expr <- eventReactive(input$ann.ok, {
    
    withProgress(message = "Normalizing the data.....", value = 0.2, {
      
      data.rma <- affy::rma(data1(), normalize = input$quantilenorm, background = input$backgroundcor)
      data.expr <- exprs(data.rma)
      data.expr <- data.expr[rownames(data.expr) != "nonsense",]
      
    })
    
    return(data.expr)
  })
  
  
  
  
  #Get sample names
  samples <- eventReactive(input$ann.ok, {
    samples <- fuzzyjoin::fuzzy_inner_join(as.data.frame(colnames(data.expr())), meta(), by = c("colnames(data.expr())" = "GEO ID"), match_fun = str_detect)
    samples <- samples[,2]
    return(samples)
  })
  
  
  output$boxplot <- renderPlot(NULL)
  output$normboxplot <- renderPlot(NULL)
  output$hist <- renderPlot(NULL)
  output$normhist <- renderPlot(NULL)
  
  
  
  observeEvent(input$ann.ok, {
    
    #raw boxplot
    output$boxplot <- renderPlot({
      if(length(data1())>0){
        withProgress(message = "Making quality plots.....", value = 0.6, {
          par(mar=c(10,2,1,1))
          boxplot(data1(),which='pm', col = "red", names = samples(), las = 2, main = "Boxplot of raw data", ylab = "log intensity")
        })
      }
    })
    
    
    #normalized boxplot
    output$normboxplot <- renderPlot({
      if(length(data.expr())>0){
        
          par(mar=c(10,2,1,1))
          boxplot(data.expr(), col = "blue", names = samples(), las = 2, main = "Boxplot of normalized data", ylab = "log intensity")
        
      }
    })
    
    #raw histogram
    output$hist <- renderPlot({
      if(length(data())>0){
        
          par(mar=c(10,2,1,1))
          hist(data1(),lwd=2,which='pm',ylab='Density',xlab='Log2 intensities',main='Density plot of raw data')
        
      }
      
    })
    
    #normalized histogram
    output$normhist <- renderPlot({
      if(length(data.expr())>0){
        
          par(mar=c(10,2,1,1))
          plotDensity(data.expr(),lwd=2,ylab='Density',xlab='Log2 intensities',main='Density plot of normalized data')
        
      }
      
    })
    
    #get proceed button
    output$proceedann <- renderUI({
      actionBttn(inputId = "ann.proceed",
                 label = "Next",
                 style = "jelly",
                 color = "primary")
    })
    
    #Download data expr
    output$downloadexpr <- downloadHandler(
      filename = "data expr",
      content = function(file){
        write.table(data.expr(), file, sep = "\t")
      }
    )
    
  })
  
  
  
  observeEvent(input$ann.proceed, {
    updateNavbarPage(session, "navbar",
                     selected = "panel4")
    
  })
  
  observeEvent(input$ann.proceed, {
    showTab("navbar", target = "panel4")
  })
  
  
  ################################################################################################################################
  #PCA
  ################################################################################################################################
  
  data.PC <- eventReactive(input$ann.proceed, {
    prcomp(t(data.expr()),scale.=TRUE)
  })
  
  pc.x <- reactive({
    
    switch(input$hpca,
           "PC1" = 1,
           "PC2" = 2,
           "PC3" = 3,
           "PC4" = 4,
           "PC5" = 5,
           "PC6" = 6,
           "PC7" = 7,
           "PC8" = 8)
    
  })
  
  pc.y <- reactive({
    switch(input$vpca,
           "PC1" = 1,
           "PC2" = 2,
           "PC3" = 3,
           "PC4" = 4, 
           "PC5" = 5,
           "PC6" = 6,
           "PC7" = 7,
           "PC8" = 8)
    
  })
  
  
  #make reactive for input
  
  output$pca <- renderPlotly({
    pca.plot(data.PC(), meta(), pc.x(), pc.y())
  })
  
  
  observeEvent(input$pca.ok, {
    updateNavbarPage(session, "navbar",
                     selected = "panel5")
    
  })
  
  observeEvent(input$pca.ok, {
    showTab("navbar", target = "panel5")
  })
  
  ################################################################################################################################
  #differential isoform expression
  ################################################################################################################################
  
  output$comparisonout <- renderUI({
    pickerInput(
      inputId = "comparisonin",
      label = "Select comparisons of interest",
      choices = get_contrasts(meta()),
      options = list(
        style = "btn-primary")
    )
  })
  
  output$exonsonlyout <- renderUI({
    if (brainarray() == "exon") {
      awesomeCheckbox(inputId = "filter.exons",
                      label = "Unique exons only",
                      value = TRUE,
                      status = "danger")
    }
    
  })
  
  
  output$geneorlistout <- renderUI({
    
    if (length(input$allexonsin) > 0){
      
      if (input$allexonsin == FALSE) {
        
        prettyRadioButtons(
          inputId = "geneorlistin",
          label = "Select transcripts based on gene(s) or select transcripts manually?",
          choices = c("Gene(s) of interest" = "gene", 
                      "Manual selection" = "manual"))
      }
    }
  })
  
  output$exonboxplotout <- renderUI({
    
    if (length(input$allexonsin) > 0){
      
      if (input$allexonsin == FALSE) {
        
        awesomeCheckbox(
          inputId = "exonboxplotin",
          label = "Boxplots",
          status = "danger",
          value = FALSE)
      }
    }
  })
  
  
  
  
  output$geneselectionout <- renderUI({
    
    if (length(input$allexonsin) > 0){
      
      if (input$allexonsin == FALSE) {
        
        if (length(input$geneorlistin) > 0){
          
          if (input$geneorlistin == "gene") {
            
            textAreaInput(
              inputId = "geneselectin",
              label = "Enter gene(s) of interest")
          }
        }
      }
    }
  })
  
  
  
  
  output$exonselectionout <- renderUI({
    
    if (length(input$allexonsin) > 0){
      
      if (input$allexonsin == FALSE) {
        
        if (length(input$geneorlistin) > 0){
          
          if (input$geneorlistin == "manual") {
            
            textAreaInput(
              inputId = "exonselectionin",
              label = "Enter transcripts of interest")
          }
        }
      }
    }
  })
  
  
  
  genes <- eventReactive(input$diffexpr.ok, {
    
    if (length(input$allexonsin) > 0){
      if (input$allexonsin == TRUE) {
        genes <- NULL
      }
    }
    
    if (length(input$geneorlistin) > 0 ){
      if (input$geneorlistin == "gene") {
        genes <- as.vector(str_split(input$geneselectin, "\n", simplify = TRUE))
      }
      
      if(input$geneorlistin == "manual") {
        genes <- NULL
      }
    }
    
    return(genes)
    
  })
  
  
  transcripts <- eventReactive(input$diffexpr.ok, {
    
    if (length(input$allexonsin) > 0){
      if (input$allexonsin == TRUE) {
        exons <- NULL
      }
    }
    
    if(length(input$geneorlistin) > 0){
      if(input$geneorlistin == "manual"){
        exons <- as.vector(str_split(input$exonselectionin, "\n", simplify = TRUE))
      }
      if(input$geneorlistin == "gene"){
        exons <- NULL
      }
    }
    
    return(exons)
    
  })
  
  
  #get statistics
  
  select.top.table <- eventReactive(input$diffexpr.ok,{
    
    data.expr <- data.expr()
    top.table <- diff_expr(data.expr, meta(), comparisons = get_contrasts(meta()))
    
    
    if (brainarray() == "transcript"){
      select.top.table <- transcript_selection(top.table, gene = genes(), transcripts = transcripts())
      
    }
    
    if (brainarray() == "exon"){
      select.top.table <- exon_selection1(top.table, annotated = annotated(), gene = genes(), transcripts = transcripts(), unique_exons = input$filter.exons, genome_wide = input$allexonsin)
      
    }
    
    return(select.top.table)
  })
  
  
  
  #select comparison
  table.choice <- reactive({
    input$comparisonin
  })
  
  
  #print table
  output$finaltable <- renderDataTable({
    select.top.table()[[table.choice()]]
  })
  
  #Make boxplot
  
  exonboxplot <- reactive({
    input$exonboxplotin
  })
  
  output$exprboxplot <- renderPlot(NULL)
  
  output$exprboxplot <- renderPlot({
    if (brainarray() == "exon"){
      if (length(exonboxplot()) > 0){
        if (exonboxplot() == TRUE){
          makeBoxplots(contrast = table.choice(), 
                       meta = meta(), 
                       data.expr = data.expr(), 
                       annotated = annotated(), 
                       gene = genes(), 
                       transcripts = transcripts(), 
                       unique_exons = input$filter.exons, 
                       genome_wide = input$allexonsin)
        }
      }
    }
    
    
  })
  
  output$exprboxplot <- renderPlot({
    if (brainarray() == "transcript"){
      if (length(exonboxplot()) > 0){
        if (exonboxplot() == TRUE){
          makeBoxplots1(
            contrast = table.choice(), 
            meta = meta(), 
            data.expr = data.expr(), 
            select.top.table = select.top.table())
          
        }
      }
    }
    
  })
  
  
  output$uiexprboxplot = renderUI({
    if (length(exonboxplot()) > 0) {
      if (exonboxplot() == TRUE){
        plotOutput("exprboxplot") %>% withSpinner(color="#0dc5c1")
      }
    }
  })
  
  #Probe mapping plot
  output$mappingplotout <- renderUI({
    
    if (length(input$allexonsin) > 0){
      
      if (input$allexonsin == FALSE) {
        
        if (length(input$geneorlistin) > 0){
          if (input$geneorlistin == "gene"){
            awesomeCheckbox(
              inputId = "mappingplotin",
              label = "Probe mapping",
              status = "danger",
              value = FALSE)
          }
          
        }
        
      }
    }
  })
  
  mappingplot <- reactive({
    input$mappingplotin
  })
  
  output$probemapping <- renderPlotly(NULL)
  
  output$probemapping <- renderPlotly({
    
    if (brainarray() == "transcript"){
      if (length(mappingplot()) > 0){
        if (mappingplot() == TRUE){
          ggplotly(probemapping.enst(chiptype(), organism(), version(), genes()))
          
        }
      }
    }
          
  })
  
  
  #Download table
  output$downloadfinal <- downloadHandler(
    filename = "diff_expr",
    content = function(file){
      write.table(select.top.table()[[table.choice()]], file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
  
  #Go to next tab
  observeEvent(input$diffexpr.ok, {
    output$proceedtovol <- renderUI({
      actionBttn(inputId = "diffexpr.proceed", 
                 label = "Next",
                 style = "jelly",
                 color = "primary")
    })
  })
  
  observeEvent(input$diffexpr.proceed, {
    updateNavbarPage(session, "navbar",
                     selected = "panel6")
    
  })
  
  observeEvent(input$diffexpr.proceed, {
    showTab("navbar", target = "panel6")
  })
  
  ################################################################################################################################
  #volcano plot
  ################################################################################################################################
  
  
  output$comparisonout1 <- renderUI({
    pickerInput(
      inputId = "comparisonin1",
      label = "Select comparisons of interest",
      choices = get_contrasts(meta()),
      options = list(
        style = "btn-primary")
    )
  })
  
  
  table.choice1 <- reactive({
    input$comparisonin1
  })
  
  
  p.choice <- reactive({
    input$raworfdr
  })
  
  
  p.threshold <- reactive({
    input$pthreshold
  })
  
  logFC.threshold <- reactive({
    input$logfcthreshold
  })
  
  
  
  
  
  
  
  output$volcano <- renderPlotly({
    
    if (length(input$comparisonin1) > 0){
      
      plotdata <- select.top.table()[[table.choice1()]]
      
      
      if (p.choice() == "Raw P-value"){
        plotdata$color[(plotdata$logFC < logFC.threshold() & plotdata$logFC > (-1 * logFC.threshold())) | plotdata$P.value > p.threshold()] = "darkgrey"
        plotdata$color[plotdata$logFC < (-1 * logFC.threshold()) & plotdata$P.value <= p.threshold()] = "blue"
        plotdata$color[plotdata$logFC >= logFC.threshold() & plotdata$P.value <= p.threshold()] = "red"
        
        Ensembl.ID = plotdata[,1]
        
        volcano <- ggplot(plotdata, aes(x = logFC, y = -log10(P.value), key = Ensembl.ID)) +
          geom_point(color = plotdata$color) +
          geom_hline(yintercept = -log10(p.threshold()), color = "grey", linetype = "dotted", size = 0.5) +
          geom_vline(xintercept = c(-1 *logFC.threshold(), logFC.threshold()), color = "grey", linetype = "dotted", size = 0.5) +
          ggtitle(table.choice1()) +
          theme_light()
      }
      
      
      if (p.choice() == "FDR"){
        plotdata$color[(plotdata$logFC < logFC.threshold() & plotdata$logFC > (-1 * logFC.threshold())) | plotdata$FDR > p.threshold()] = "darkgrey"
        plotdata$color[plotdata$logFC < (-1 * logFC.threshold()) & plotdata$FDR <= p.threshold()] = "blue"
        plotdata$color[plotdata$logFC >= logFC.threshold() & plotdata$FDR <= p.threshold()] = "red"
        
        Ensembl.ID = plotdata[,1]
        
        volcano <- ggplot(plotdata, aes(x = logFC, y = -log10(FDR), key = Ensembl.ID)) +
          geom_point(color = plotdata$color) +
          geom_hline(yintercept = -log10(p.threshold()), color = "grey", linetype = "dotted", size = 0.5) +
          geom_vline(xintercept = c(-1 *logFC.threshold(), logFC.threshold()), color = "grey", linetype = "dotted", size = 0.5) +
          ggtitle(table.choice1()) +
          theme_light()
      }
      
      
      ggplotly(volcano)
      
    }
  })
  
  
  
  output$voltable <- renderDataTable({
    
    if (length(input$comparison1) > 0){
      voldata <- select.top.table()[[table.choice1()]]
      
      if (p.choice() == "Raw P-value"){
        voltable <- voldata %>%
          filter(P.value <= p.threshold()) %>%
          filter(abs(logFC) >= logFC.threshold())
      }
      
      if (p.choice() == "FDR"){
        voltable <- voldata %>%
          filter(FDR <= p.threshold()) %>%
          filter(abs(logFC) >= logFC.threshold())
      }
      
      return(voltable)
      
    }
    
  })
  
  
  
}


shinyApp(ui = ui, server = server)





