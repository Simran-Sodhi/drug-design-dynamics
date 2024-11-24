# Check if shiny is installed; if not, install it
if (!require("shiny")) {
  install.packages("shiny")
}

# Check if ggplot2 is installed; if not, install it
if (!require("ggplot2")) {
  install.packages("ggplot2")
}

library(shiny)
library(bslib)
library(ggplot2)

# First, set your working directory to source file location.
# you should stand at currentDir:
## |-- currentDir
## |--------- app.R
## |--------- output
## |-------------------- simulations.csv
## |--------- monteCarlo
## |--------- machineLearning


# format of output file (can also see /simulations.csv/ in github):
## BindingEnergy
## -22632847.078032
## 6574893379.883727
## -15532783.512819
## -50041433.277568


# Define UI
ui <- page_sidebar(
  title = "Simulation of Protein-ligand Interaction",
  
  sidebar = sidebar(
    fileInput("Protein", "Upload Protein", accept = c(".pdb", ".mol2")),
    fileInput("Ligands", "Upload Ligands", accept = c(".pdb", ".mol2"), 
              multiple = TRUE),

    radioButtons("method", "Simulation Method:",
                 c("Monte Carlo" = "mc",
                   "Machine Learning" = "ml")),
    
    actionButton("runGoCode", "Run Simulation")
  ),
  
  navset_card_underline(
    title = "Visualizations",
    nav_panel("Plot", plotOutput("plot")),
    nav_panel("Table", tableOutput("table"))
  )
)

server <- function(input, output) {
  observeEvent(input$runGoCode, {
    req(input$Protein, input$Ligands)  # Ensure both files are uploaded
    
    proteinFile <- input$Protein$datapath
    ligandFiles <- input$Ligands$datapath
    ligandFileNames <- input$Ligands$name
    
    simulationMethod <- input$method
    
    if (simulationMethod == "mc") {
      # Compile the Go program
      compile_result <- 
        system("cd monteCarlo && go build -o monteCarlo main.go", intern = TRUE)
      print(compile_result)
      
      # Run the compiled Go program
      run_result <- 
        system(paste("cd monteCarlo && ./monteCarlo", proteinFile, 
                     paste(ligandFiles, collapse = " ")), intern = TRUE)
      print(run_result)
    } else if (simulationMethod == "ml") {
      
      compile_result <- 
        system("cd machineLearning && go build -o machineLearning main.go", 
               intern = TRUE)
      print(compile_result)
      
      # Run the compiled Go program
      run_result <- 
        system(paste("cd machineLearning && ./machineLearning", proteinFile, 
                     paste(ligandFiles, collapse = " ")), intern = TRUE)
      print(run_result)
    }
    
    
    # Read results from the output file: output/simulations.csv
    outputFile <- "output/simulations.csv"
    if (file.exists(outputFile)) {
      results <- read.csv(outputFile, header = TRUE)
      print(results)
      
      # Add ligandfile names as a column
      results$LigandFile <- ligandFileNames
      
      # Render results in the table
      output$table <- renderTable({
        results
      })
      
      # Render scatter plot of the results
      output$plot <- renderPlot({
        ggplot(results, aes(x = LigandFile, y = BindingEnergy)) +
          geom_point() +
          labs(
            title = "Binding Energy of Multiple Ligands",
            x = "Ligand",
            y = "Binding Energy"
          ) +
          theme_minimal()
      })
    } else {
      warning("Output file not found.")
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)