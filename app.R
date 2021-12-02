library(shiny)
library(shinythemes)
library(tidyverse)


ui <- navbarPage(title = 'Methods',
                 tabPanel('qRT-PCR',
                          sidebarPanel(tags$h3("Parameters"),
                                       fileInput(inputId = "upload",
                                                 label = "Upload your data (results.csv or results.txt)",
                                                 accept = c(".csv")),
                                       textInput(inputId = "endogenousControl",
                                                 label = "Endogenous Control gene",
                                                 placeholder = "Endogenous Control"),
                                       textInput(inputId = "referenceSample",
                                                 label = "Reference Sample",
                                                 placeholder = "Reference Sample"),
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
                            downloadButton(outputId = "downloadData", label = "Download"),
                           # 'Head',
                           # tableOutput("head"),
                           # 'Plot',
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
  #output$head <- renderTable({
  #  head(data())
  #})
  
  final <- reactive({a <- data() %>% group_by(Target, Sample, Biological.Set.Name) %>% summarise(mean.Cq = mean(Cq))
                    b <- a %>% filter(Target == input$endogenousControl) %>% rename(ref.Cq = mean.Cq)
                    c <- a %>% filter(Target != input$endogenousControl)
                    d <- left_join(c, b, by = "Sample")
                    e <- d %>% mutate(del.Cq = mean.Cq - ref.Cq)
                    f <- e %>% select(Target.x, Sample, Biological.Set.Name.x, del.Cq)
                    g <- f %>% filter(Biological.Set.Name.x == input$referenceSample) %>% group_by(Target.x) %>% summarise(mean.del.Cq = mean(del.Cq))
                    h <- f %>% filter(Biological.Set.Name.x != input$referenceSample)
                    i <- left_join(f, g, by = "Target.x") %>% mutate(del.del.Cq = del.Cq - mean.del.Cq)
                    final <- i %>% select(Target.x, Sample, Biological.Set.Name.x, del.del.Cq) %>% group_by(Target.x, Biological.Set.Name.x) %>% summarise(mean.dd.Cq = mean(del.del.Cq), sd.Cq = sd(del.del.Cq)) %>% mutate(fc = 2^-(mean.dd.Cq))
                    finalPlot <- final %>% rename(Target = Target.x, Biological.Set.Name = Biological.Set.Name.x)
  })
  output$barPlot <- renderPlot({
    input$submit
    isolate({
    ggplot(final(), aes(x = Biological.Set.Name, y = fc, fill = Target)) +
      geom_bar(stat = "identity", width = 0.5, position = position_dodge(0.6), color = "#434343") +
      geom_errorbar(aes(ymin = fc - sd.Cq, ymax = fc - sd.Cq), width=.2, position = position_dodge(0.6))+
      theme_bw()+
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
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(final(), file)
    }
  )
  
}

shinyApp(ui, server)
