# app.R
library(shiny)
library(glue)
library(DT)
library(reticulate)
library(shinydashboard)
library(promises)
library(future)
library(shinyWidgets)

plan(multisession)

# ---- Python Imports ----
os        <- import("os")
glob      <- import("glob")
re        <- import("re")
subprocess<- import("subprocess")
argparse  <- import("argparse")
multi     <- import("multiprocessing")
shutil    <- import("shutil")

# Ensure Python path
reticulate::py_run_string(
  "import sys; import os; sys.path.append(f'{os.getcwd()}')"
)

Sys.setenv(PYTHONUNBUFFERED = "1")


options(shiny.maxRequestSize = 100*1024^2)


# ---- UI ----
ui <- dashboardPage(
  dashboardHeader(
    title      = HTML("<em>Pyricularia oryzae</em><br>MiniMonsterPlex Tool"),
    titleWidth = 350
  ),
  dashboardSidebar(
    width = 350,
    sidebarMenu(
      menuItem("Inputs",            tabName = "inputs",  icon = icon("upload")),
      menuItem("Alignment Summary", tabName = "summary", icon = icon("table"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .form-group { margin-bottom: 15px; }
      .btn-primary { background-color: #007bff; border-color: #007bff; }
      .btn-primary:hover { background-color: #0056b3; border-color: #004085; }
    "))),
    tabItems(
      # ---- Inputs tab ----
      tabItem(tabName = "inputs",
              fluidRow(
                # Left: placeholder for uploads
                box(
                  width = 6, title = "Upload Data", status = "primary", solidHeader = TRUE,
                  uiOutput("file_upload_ui")
                ),
                # Right: project selector + metadata filename
                box(
                  width = 6, title = "Settings", status = "primary", solidHeader = TRUE,
                  selectizeInput(
                    "project_id", "Project ID",
                    choices = NULL, selected = "",
                    options = list(create = TRUE, placeholder = "Enter Project ID")
                  ),
                  textInput("metadata.file", "Metadata filename", value = "metadata.csv")
                )
              ),
              fluidRow(
                column(12, align = "center",
                       actionButton(
                         "myButton", "Run MiniMonsterPlex",
                         class = "btn-primary",
                         style = "width:250px; margin-top:20px;"
                       )
                )
              )
      ),
      
      # ---- Alignment Summary tab ----
      tabItem(tabName = "summary",
              fluidRow(
                box(
                  width = 12, title = "Alignment Summary Table", status = "info", solidHeader = TRUE,
                  DTOutput("alignment_table")
                )
              )
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  
  
  setwd(getwd())
  base_dir     <- getwd()
  projects_dir <- file.path(base_dir, "Projects")
  dir.create(projects_dir, showWarnings = FALSE)
  
  # Create necessary directories as soon as project ID is selected
  observeEvent(input$project_id, {
    req(nzchar(input$project_id))  # Ensure project_id is not empty
    
    project_dir     <- file.path(projects_dir, input$project_id)
    fastq_dir       <- file.path(project_dir, "newFastq")
    output_dir      <- file.path(project_dir, "output")
    metadata_dir    <- file.path(project_dir, "metadata")
    
    dir.create(fastq_dir,    recursive = TRUE, showWarnings = FALSE)
    dir.create(output_dir,   recursive = TRUE, showWarnings = FALSE)
    dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
    
    message("Project directories created for: ", input$project_id)
    message("  - newFastq: ", fastq_dir)
    message("  - output: ", output_dir)
    message("  - metadata: ", metadata_dir)
  })
  
  
  # Populate project dropdown
  observe({
    existing <- list.dirs(projects_dir, full.names = FALSE, recursive = FALSE)
    updateSelectizeInput(session, "project_id",
                         choices = existing,
                         selected = "",
                         server = TRUE)
  })
  
  # Dynamically render fileInputs only when project_id is non-empty
  output$file_upload_ui <- renderUI({
    req(nzchar(input$project_id))
    tagList(
      fileInput("fastq_files",     "Upload FASTQ files:", multiple = TRUE),
      fileInput("metadata_upload", "Upload metadata file")
    )
  })
  
  # Copy FASTQ files
  observe({
    req(nzchar(input$project_id), input$fastq_files)
    fastq_dir <- file.path(projects_dir, input$project_id, "newFastq")
    dir.create(fastq_dir, recursive = TRUE, showWarnings = FALSE)
    for (i in seq_along(input$fastq_files$datapath)) {
      file.copy(
        input$fastq_files$datapath[i],
        file.path(fastq_dir, input$fastq_files$name[i]),
        overwrite = TRUE
      )
    }
    message("FASTQ files copied to: ", fastq_dir)
  })
  
  # Copy metadata file
  observe({
    req(nzchar(input$project_id), input$metadata_upload)
    meta_dir <- file.path(projects_dir, input$project_id, "metadata")
    dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(
      input$metadata_upload$datapath,
      file.path(meta_dir, input$metadata_upload$name),
      overwrite = TRUE
    )
    message("Metadata copied to: ", file.path(meta_dir, input$metadata_upload$name))
  })
  
  # Run analysis
  alignment_trigger <- reactiveVal(0)
  
  # Reactive expression for tracking project ID
  project_id_reactive <- reactive({
    input$project_id
  })
  
  observeEvent(input$myButton, {
    req(nzchar(input$project_id), input$metadata.file)
    
    project_id <- input$project_id
    project_dir   <- file.path(projects_dir, project_id)
    fastq_dir     <- file.path(project_dir, "testFastq")
    metadata_dir  <- file.path(project_dir, "metadata")
    metadata_file <- file.path(metadata_dir, input$metadata.file)
    
    dir.create(fastq_dir,    recursive = TRUE, showWarnings = FALSE)
    dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
    
    default_meta <- file.path(base_dir, input$metadata.file)
    if (file.exists(default_meta)) {
      file.copy(default_meta, metadata_file, overwrite = TRUE)
      message("Copied default metadata to: ", metadata_file)
    }
    
    print(paste("trying to start main", project_id, input$metadata.file))
    reticulate::source_python("MiniMonsterPlex_shiny.py")
    main(project_id, input$metadata.file)
    message("Analysis complete for: ", project_id)
    
    # Trigger table refresh
    alignment_trigger(alignment_trigger() + 1)
  })
  
  # Render the alignment summary table reactively
  output$alignment_table <- renderDT({
    req(nzchar(input$project_id))
    
    project_id <- input$project_id
    summary_path <- file.path(projects_dir, project_id, "output/alignment_summary.csv")
    
    if (file.exists(summary_path)) {
      df <- read.csv(summary_path)
      colnames(df) <- c("sample ID", "# reads", "# aligned", "percent")
      # Remove duplicates based on "sample ID", keeping the latest occurrence
      df <- df[!duplicated(df[["sample ID"]], fromLast = TRUE), ]
      df <- df[order(df[["sample ID"]]), ]
      datatable(
        df, extensions = 'Buttons', rownames = FALSE,
        options = list(dom='Bfrtip', buttons=c('copy','csv','excel','pdf','print'))
      )
    } else {
      datatable(
        data.frame(Message = "No alignment summary found for this project."),
        rownames = FALSE
      )
    }
  })
  
  
  
}

shinyApp(ui, server)
