# CODE THAT WILL GIVE GRAPH OF GO OUTPUT:
library(shiny) # Loads Shiny web app library
library(ggplot2) # Load ggplot2 for data visualization



parse_go_output <- function(go_output) {
  # Map the single letters to the protein structure type
  structure_map <- list('H' = 'Helix', 'E' = 'Sheet', 'T' = 'Turn', 'C' = 'Coil')
  # Initialize and empty data frame
  parsed_structure <- data.frame(Type = character(), Start = integer(), End = integer())
  
  # Splits the go output into individual characters
  structure <- unlist(strsplit(go_output, ""))
  # Start position of the first structure segment
  start <- 1
  # Adds the first character (current character being evaluated)
  last_char <- structure[1]  # Initialize with the first character
  
  # Loop through the structure array
  for (i in seq_along(structure)) {
    if (structure[i] != last_char || i == length(structure)) {
      if (start < i) {
        # New row added to parsed_structure, if a new structure begins
        # Type of structure (helix, turn, etc.) looked up in the structure map using "last_char"
        parsed_structure <- rbind(parsed_structure, 
                                  data.frame(Type = structure_map[[last_char]], Start = start, End = i - 1))
      }
      start <- i # Push the counter ahead. Update the start position for the new segment.
      last_char <- structure[i] # Update the last character to the current one.
    }
  }
  
  # Ensure the last character sequence is added
  if (start <= length(structure)) {
    parsed_structure <- rbind(parsed_structure, 
                              data.frame(Type = structure_map[[last_char]], Start = start, End = length(structure)))
  }
  
  return(parsed_structure) # Return the parsed structure
}

# Define a function to get the correct command based on the OS
get_command <- function() {
  os_type <- Sys.info()["sysname"]
  command <- ""
  
  if (os_type == "Windows") {
    command <- "CSB_PSP.exe"  # Assuming the executable has .exe extension on Windows
  } else {
    command <- "./CSB_PSP"  # Assuming the executable does not need an extension on Unix/Linux/macOS
  }
  
  return(command)
}

# Define UI
ui <- fluidPage(
  titlePanel("Protein Secondary Structure Prediction"),
  sidebarLayout(
    sidebarPanel(
      textInput("input_string", "Enter Protein Sequence:"),
      actionButton("process_button", "Process")
    ),
    mainPanel(
      verbatimTextOutput("output_text"),
      plotOutput("structure_plot")  # Output for the plot
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  observeEvent(input$process_button, {
    command <- get_command()
    go_output <- system2(command, args = input$input_string, stdout = TRUE)
    
    # Call the Go program and capture its output
    # go_output <- system2("./CSB_PSP", args = input$input_string, stdout = TRUE)
    print(go_output)  # Debug: print the output to R console
    # go_output <- system2("./CSB_PSP", args = c("-sequence", input$input_string), stdout = TRUE)
    
    # Parse Go output into a structured format
    parsed_structure <- parse_go_output(go_output)
    
    # Render text output
    output$output_text <- renderText({
      paste("Processed output:", go_output)
    })
    

    
    # Render plot output
    output$structure_plot <- renderPlot({
      if (nrow(parsed_structure) > 0) {
        ggplot(parsed_structure, aes(x = Start - 0.5, xend = End + 0.5, y = 0, yend = 0, color = Type)) +
          geom_segment(size = 10) +  # Increase size for better visibility
          scale_color_manual(values = c("Helix" = "red", "Sheet" = "green", "Turn" = "blue", "Coil" = "black")) +
          theme_minimal() +
          labs(title = "Secondary Structure Prediction", x = "Position", y = "", color = "") +
          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
      }
    })   
    
    print(parsed_structure)  # Added this line in the server function to see the DataFrame before plotting
    
  })
}

# Run the application 
shinyApp(ui, server)


