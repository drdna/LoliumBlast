library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- 1. Data Parsing Function ---
parse_snp_data <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  # Locate the DATA block
  data_start <- which(trimws(lines) == "DATA")
  if(length(data_start) == 0) return(NULL)
  
  # Read tab-separated columns 1-4, then everything else as a single string
  # We use fill=TRUE and a high number of col.names to catch all space-separated strains
  raw_data <- read.table(text = lines[(data_start + 1):length(lines)], 
                         sep = "\t", stringsAsFactors = FALSE, fill = TRUE,
                         col.names = c("Contig", "Position", "RefAlt", "Count", "StrainsRaw"))
  
  # Process the variable-length strain list
  processed <- raw_data %>%
    # Convert 'StrainsRaw' (and any accidental extra columns) into a single space-separated string
    unite("AllStrains", StrainsRaw:last_col(), sep = " ", na.rm = TRUE) %>%
    mutate(AllStrains = trimws(AllStrains)) %>%
    # Expand: Create one row per strain for each SNP position
    separate_rows(AllStrains, sep = "\\s+") %>%
    filter(AllStrains != "") %>%
    rename(Strain = AllStrains) %>%
    mutate(Position = as.numeric(Position))
  
  return(processed)
}

# --- 2. UI Definition ---
ui <- fluidPage(
  titlePanel("Genomic SNP Explorer"),
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
      tabsetPanel(
        tabPanel("Track View", 
                 plotOutput("snp_plot", height = "600px")),
        tabPanel("SNP Density", 
                 plotOutput("density_plot", height = "400px"),
                 br(),
                 numericInput("binwidth", "Density Window Size (bp):", 
                              value = 10000, min = 1, step = 100))
      ),
      br(),
      wellPanel(
        uiOutput("pos_slider")
      )
    )
  )
)

# --- 3. Server Logic ---
server <- function(input, output, session) {
  
  # Reactive dataset
  dataset <- reactive({
    req(input$file)
    parse_snp_data(input$file$datapath)
  })
  
  # Dynamic UI: Contig Selector
  output$contig_selector <- renderUI({
    req(dataset())
    selectInput("selected_contig", "2. Select Contig:", 
                choices = sort(unique(dataset()$Contig)))
  })
  
  # Dynamic UI: Bottom Position Slider
  output$pos_slider <- renderUI({
    req(input$selected_contig)
    df_sub <- dataset() %>% filter(Contig == input$selected_contig)
    
    sliderInput("pos_range", "3. Adjust Position Range (Fine Increment):",
                min = min(df_sub$Position), 
                max = max(df_sub$Position),
                value = c(min(df_sub$Position), max(df_sub$Position)),
                step = 1, # Finer increment as requested
                width = "100%",
                sep = "")
  })
  
  # Filtered data based on UI inputs
  filtered_data <- reactive({
    req(dataset(), input$selected_contig, input$pos_range)
    dataset() %>%
      filter(Contig == input$selected_contig,
             Position >= input$pos_range[1],
             Position <= input$pos_range[2])
  })
  
  # Sidebar Table: Count SNPs per strain in the current window
  output$snp_counts_table <- renderTable({
    req(nrow(filtered_data()) > 0)
    filtered_data() %>%
      group_by(Strain) %>%
      summarise(SNPs = n()) %>%
      arrange(desc(SNPs))
  }, striped = TRUE, spacing = "s", width = "100%")
  
  # Tab 1: Track View (One strain per line)
  output$snp_plot <- renderPlot({
    df_plot <- filtered_data()
    req(nrow(df_plot) > 0)
    
    ggplot(df_plot, aes(x = Position, y = Strain)) +
      # Grey horizontal guides
      geom_segment(aes(x = input$pos_range[1], xend = input$pos_range[2], 
                       y = Strain, yend = Strain), 
                   color = "grey95") +
      # The SNP points
      geom_point(color = "firebrick", size = 3, alpha = 0.8) +
      theme_minimal() +
      labs(x = "Genomic Position (bp)", y = NULL, 
           title = paste("SNP Track:", input$selected_contig)) +
      theme(
        axis.text.y = element_text(face = "bold", family = "mono", size = 11),
        panel.grid.major.y = element_blank()
      )
  })
  
  # Tab 2: SNP Density (Ignores multiple counts per site)
  output$density_plot <- renderPlot({
    df_unique <- filtered_data() %>% distinct(Position)
    req(nrow(df_unique) > 0)
    
    ggplot(df_unique, aes(x = Position)) +
      geom_histogram(binwidth = input$binwidth, fill = "midnightblue", color = "white") +
      theme_minimal() +
      labs(
        title = "Unique SNP Density (Collapsed across strains)",
        x = "Genomic Position (bp)",
        y = "Unique SNPs per Window"
      ) +
      scale_x_continuous(limits = input$pos_range)
  })
}

shinyApp(ui, server)