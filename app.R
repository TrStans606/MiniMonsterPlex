library(shiny)
library(glue)
options(shiny.maxRequestSize = 100*1024^2)
runMiniMonsterPlex <- function(output.folder,metadata.file,input.folder,isolate.list,
                               isolate.file,host.list,host.file){
  # If the isolate list is empty, set it to NULL; otherwise, split it into a list
  if (isolate.list == "") {
    isolate.list = NULL
  } else {
    isolate.list = unlist(strsplit(isolate.list, " "))
  }
  
  # If the isolate file is empty, set it to NULL
  if (isolate.file == "") {
    isolate.file = NULL
  }
  
  # If the host list is empty, set it to NULL; otherwise, split it into a list
  if (host.list == "") {
    host.list = NULL
  } else {
    host.list = unlist(strsplit(host.list, " "))
  }
  
  # If the host file is empty, set it to NULL
  if (host.file == "") {
    host.file = NULL
  }
  
  # Construct the command options based on the provided isolate and host lists/files
  if (is.null(isolate.list)) {
    isolate.list.command = " "
  } else {
    isolate.list.command = glue(' -i {isolate.list}')
  }
  if (is.null(isolate.file)) {
    isolate.file.command = " "
  } else {
    isolate.file.command = glue(' -il {isolate.file}')
  }
  if (is.null(host.list)) {
    host.list.command = " "
  } else {
    host.list.command = glue(' -hf {host.list}')
  }
  if (is.null(host.file)) {
    host.file.command = " "
  } else {
    host.file.command = glue(' -hfl {host.file}')
  }
  
  # Construct the final command to run the MiniMonsterPlex.py script with the provided options
  command = glue('python MiniMonsterPlex.py -o {output.folder} -m {metadata.file} -f {input.folder} {isolate.list.command} {isolate.file.command} {host.list.command} {host.file.command}')
  
  # Execute the constructed command
  system(command)
  
}


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("MiniMonsterPlex"),
  
  fileInput("fastq_files", "Upload the data you want anaylazed:",multiple = TRUE),
  fileInput("metadata_upload", "Upload a metadata file"),
  fileInput("isolate_upload", "Upload a isolate file"),
  fileInput("hosts_upload", "Upload a hosts file"),
  actionButton("copyBtn", "Upload Files"),
  
  textInput("output.folder", "Enter output folder here", value = "output",placeholder = "Must not already exist"),
  textInput("metadata.file", "Enter metadata_file here", value = "metadata.csv"),
  textInput("input.folder", "Enter input folder here", value = "fastq"),
  textInput("isolate.list", "Enter a space seperated list of isolates you want included here. All will be included by default.", value = ""),
  textInput("isolate.file", "Enter a file containing a new line list of isolates you want included here. All will be included by default.", value = ""),
  textInput("host.list", "Enter a space seperated list of hosts you want included here. All will be included by default.", value = ""),
  textInput("host.file", "Enter a file containing a new line list of hosts you want included here. All will be included by default.", value = ""),
  actionButton("myButton", "Run MiniMonsterPlex")
)

server <- # Define server logic required to draw a histogram
  function(input, output, session) {
    observeEvent(input$copyBtn, {
      req(input$fastq_files)
      fastq_upload <- input$fastq_files
      for(i in 1:length(fastq_upload$datapath)){
        file.copy(fastq_upload$datapath[i], file.path(getwd(), 'fastq',fastq_upload$name[i]))
      }
      req(input$metadata_upload)
      metadata_file <- input$metadata_upload
      file.copy(metadata_file$datapath, file.path(getwd(), metadata_file$name))
      
      req(input$isolate_upload)
      isolates_file <- input$isolate_upload
      file.copy(isolates_file$datapath, file.path(getwd(), isolates_file$name))
      
      req(input$hosts_upload)
      hosts_file <- input$hosts_upload
      file.copy(hosts_file$datapath, file.path(getwd(), hosts_file$name))
    })
    observeEvent(input$myButton, {
      runMiniMonsterPlex(input$output.folder,input$metadata.file,input$input.folder,input$isolate.list,input$isolate.file,input$host.list,input$host.file)
    })

    
  }

shinyApp(ui=ui, server = server)