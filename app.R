library(shiny)
library(shinythemes)
library(tidyverse)
library(viridis)


ui <- navbarPage(title = 'Methods',
                 fluid = TRUE,
                 tabPanel('qRT-PCR',
                          sidebarPanel(tags$h3("Parameters"),
                                       fileInput(inputId = "upload",
                                                 label = "Upload your data (results.csv or results.txt)",
                                                 accept = c(".csv")),
                                       textInput(inputId = "endogenousControl",
                                                 label = "Endogenous Control Gene",
                                                 placeholder = "GAPDH"),
                                       textInput(inputId = "referenceSample",
                                                 label = "Reference Condition",
                                                 placeholder = "Reference Condition"),
                                       textInput(inputId = "title",
                                                 label = "Plot title",
                                                 placeholder = "Plot title"),
                                       textInput(inputId = "subtitle",
                                                 label = "Plot subtitle",
                                                 placeholder = "Plot subtitle"),
                                       textInput(inputId = "xlabel",
                                                 label = "X-axis label",
                                                 placeholder = "X-axis label"),
                                       textInput(inputId = "ylabel",
                                                 label = "Y-axis label",
                                                 placeholder = "Y-axis label"),
                                       actionButton(inputId = "submit",
                                                    label = "Submit"),
                                       
                                       
                          ),
                          mainPanel(
                            includeHTML("qPCR.html"),
                            #br(),
                            #br(),
                            #br(),
                            hr(),
                            #br(),
                            #br(),
                            #br(),
                                downloadButton(outputId = "downloadData", label = "Download Data"),
                                #br(),
                                #br(),
                                #br(),
                                #br(),
                                #br(),
                                tableOutput("finalData"),
                            hr(),
                                downloadButton(outputId = "downloadPlot", label = "Download Plot"),
                                radioButtons(inputId = "var3", label = NULL, choices = list("png", "pdf", "svg", "tiff", "jpeg"), inline=TRUE),
                                plotOutput("barPlot")
                            
                          )
                 )
      )

server <- function(input, output, session) {
  
  data <- reactive({
    req(input$upload)
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
           csv = vroom::vroom(input$upload$datapath, delim = ","),
           validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })
  
  final <- reactive({a <- data() %>% group_by(Target, Sample, Biological.Set.Name) %>% summarise(mean.Cq = mean(Cq))
                    b <- a %>% filter(Target == input$endogenousControl) %>% rename(ref.Cq = mean.Cq)
                    c <- a %>% filter(Target != input$endogenousControl)
                    d <- left_join(c, b, by = "Sample")
                    e <- d %>% mutate(del.Cq = mean.Cq - ref.Cq)
                    f <- e %>% select(Target.x, Sample, Biological.Set.Name.x, del.Cq)
                    g <- f %>% filter(Biological.Set.Name.x == input$referenceSample) %>% group_by(Target.x) %>% summarise(mean.del.Cq = mean(del.Cq))
                    h <- f %>% filter(Biological.Set.Name.x != input$referenceSample)
                    i <- left_join(f, g, by = "Target.x") %>% mutate(del.del.Cq = del.Cq - mean.del.Cq)
                    final <- i %>% select(Target.x, Sample, Biological.Set.Name.x, del.del.Cq) %>% group_by(Target.x, Biological.Set.Name.x) %>% summarise(mean.dd.Cq = mean(del.del.Cq), SD = sd(del.del.Cq)) %>% mutate(Fold_Change = 2^-(mean.dd.Cq))
                    finalPlot <- final %>% rename(Target = Target.x, Condition = Biological.Set.Name.x)
  })
  
  output$finalData <- renderTable({
    input$submit
    isolate({
      final()
    })
  })
  
  plot <- reactive({
    ggplot(final(), aes(x = Condition, y = Fold_Change, fill = Target)) +
      geom_bar(stat = "identity", width = 0.5, position = position_dodge(0.6), color = "#434343") +
      geom_errorbar(aes(ymin = Fold_Change - SD, ymax = Fold_Change + SD), width=.2, position = position_dodge(0.6))+
      theme_bw()+
      scale_fill_viridis(discrete = TRUE)+
      labs(title = input$title,
           subtitle = input$subtitle,
           x = input$xlabel,
           y = input$ylabel,
           fill = "Target gene",
           caption = NULL)+
      theme(plot.title=element_text(size=25),
            axis.text=element_text(size=12, color = "#434343"),
            axis.title=element_text(size=18),
            legend.text=element_text(size=10),
            legend.title=element_text(size=15),
            legend.position = c(0.2, 0.8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    })
  
  output$barPlot <- renderPlot({
    input$submit
    isolate({
    #here goes the ggplot command and remove plot()
    plot()
    })
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(final(), file)
    }
  )
  
  output$downloadPlot <- downloadHandler(
    # Specify the file name
    filename = function() {
      paste("plot", input$var3, sep = ".")
    },
    # Open the device and create the plot and close the device
    content = function(file) {
      ggsave(file, plot = plot(), device = input$var3, dpi = 600)
    }
  )
  
}

shinyApp(ui, server)
