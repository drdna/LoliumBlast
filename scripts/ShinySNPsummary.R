library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- 1. Data Parsing Function ---
parse_snp_data <- function(filepath) {
  lines <- readLines(filepath)
  data_start <- which(lines == "DATA")
  if(length(data_start) == 0) stop("File must contain a 'DATA' header.")
  
  raw_data <- read.table(text = lines[(data_start + 1):length(lines)], 
                         sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  
  processed <- raw_data %>%
    rename(Contig = V1, Position = V2, RefAlt = V3, Count = V4) %>%
    unite("StrainList", V5:last_col(), sep = " ", na.rm = TRUE) %>%
    mutate(StrainList = trimws(StrainList)) %>%
    separate_rows(StrainList, sep = " ") %>%
    filter(StrainList != "")
  
  return(processed)
}

# --- 2. UI Definition ---
ui <- fluidPage(
  titlePanel("Genomic SNP Track Viewer"),
  hr(),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("file", "1. Upload SNP Data"),
      uiOutput("contig_selector"),
      hr(),
      h4("SNP Counts (Current View)"),
      tableOutput("snp_counts_table")
    ),
    
    mainPanel(
      width = 9,
      plotOutput("snp_plot", height = "600px"),
      br(),
      wellPanel(
        uiOutput("pos_slider")
      )
    )
  )
)

# --- 3. Server Logic ---
server <- function(input, output, session) {
  
  dataset <- reactive({
    req(input$file)
    parse_snp_data(input$file$datapath)
  })
  
  output$contig_selector <- renderUI({
    df <- dataset()
    selectInput("selected_contig", "2. Select Contig:", choices = sort(unique(df$Contig)))
  })
  
  output$pos_slider <- renderUI({
    req(input$selected_contig)
    df <- dataset() %>% filter(Contig == input$selected_contig)
    
    # Setting step to 1 for maximum precision
    sliderInput("pos_range", "3. Adjust Position Range:",
                min = min(df$Position), 
                max = max(df$Position),
                value = c(min(df$Position), max(df$Position)),
                step = 1, 
                width = "100%",
                sep = "")
  })
  
  # Reactive filtering for both plot and table
  filtered_data <- reactive({
    req(dataset(), input$selected_contig, input$pos_range)
    dataset() %>%
      filter(Contig == input$selected_contig,
             Position >= input$pos_range[1],
             Position <= input$pos_range[2])
  })
  
  # Sidebar Table: Count SNPs per strain
  output$snp_counts_table <- renderTable({
    req(filtered_data())
    filtered_data() %>%
      group_by(Strain = StrainList) %>%
      summarise(SNPs = n()) %>%
      arrange(desc(SNPs))
  }, striped = TRUE, spacing = "s", width = "100%")
  
  # Main Plot
  output$snp_plot <- renderPlot({
    df_plot <- filtered_data()
    req(nrow(df_plot) > 0)
    
    ggplot(df_plot, aes(x = Position, y = StrainList)) +
      geom_segment(aes(x = input$pos_range[1], xend = input$pos_range[2], 
                       y = StrainList, yend = StrainList), 
                   color = "grey92") +
      geom_point(color = "firebrick", size = 3.5, alpha = 0.8) +
      theme_minimal() +
      labs(x = "Genomic Position (bp)", y = NULL) +
      theme(
        axis.text.y = element_text(face = "bold", size = 12, family = "mono"),
        panel.grid.major.y = element_blank()
      )
  })
}

shinyApp(ui, server)
