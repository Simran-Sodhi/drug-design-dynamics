# Remember to specify the path to the Python executable 
# to ensure compatibility with the required environment at line 241

# Check if shiny is installed; if not, install it
if (!require("shiny")) {
  install.packages("shiny")
}

# Check if ggplot2 is installed; if not, install it
if (!require("ggplot2")) {
  install.packages("ggplot2")
}

# Check if bslib is installed; if not, install it
if (!require("bslib")) {
  install.packages("bslib")
}

# Check if dplyr is installed; if not, install it
if (!require("dplyr")) {
  install.packages("dplyr")
}

# Check if plotly is installed; if not, install it
if (!require("plotly")) {
  install.packages("plotly")
}

# Check if shinyFiles is installed; if not, install it
if (!require("shinyFiles")) {
  install.packages("shinyFiles")
}

# Check if bio3d is installed; if not, install it
if (!require("bio3d")) {
  install.packages("bio3d")
}

# Check if r3dmol is installed; if not, install it
if (!require("r3dmol")) {
  install.packages("r3dmol")
}


# Remember to check the packages
library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(plotly)
library(shinyFiles)
library(bio3d)
library(r3dmol)


# Define UI
ui <- page_sidebar(
  title = "Simulation of Protein-ligand Interaction",
  
  sidebar = sidebar(
    width = 350, 
    
    radioButtons("method", "Simulation Method:",
                 c("Metropolis" = "mc",
                   "Machine Learning" = "ml")),
    
    # Placeholder for dynamic file inputs
    uiOutput("dynamic_file_inputs"),
    
    actionButton("runGoCode", "Run Simulation")
  ),
  
  # Conditional rendering for metropolis
  conditionalPanel(
    condition = "input.method == 'mc'",
    navset_card_underline(
      nav_panel("Plot", plotlyOutput("plot", height = "500px")),
      nav_panel("Structure",
                tags$video(
                  src = "output/results.mp4", 
                  type = "video/mp4",
                  controls = "controls",     # Add playback controls
                  width = "100%",            
                  height = "600px"           
                ))
    )
  ),
  
  # Conditional rendering for Machine Learning
  conditionalPanel(
    condition = "input.method == 'ml'",
    navset_card_underline(
      nav_panel("Evaluation", tableOutput("evaluation_results")), 
      nav_panel("Feature Importance", uiOutput("feature_importance")) 
    )
  )
)

server <- function(input, output, session) {
  outputDir <- file.path(getwd(), "output")  # Fixed output directory
  
  # Ensure the output directory exists
  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
  }
  
  # Expose the 'output' directory to the web
  addResourcePath("output", outputDir)
  
  # Enable directory selection
  shinyDirChoose(input, "ligandDir",
                 roots = c(MCdata = file.path(getwd(), "MCdata")),
                 session = session)
  
  # Dynamic file input rendering based on the selected method
  output$dynamic_file_inputs <- renderUI({
    if (input$method == "mc") {
      tagList(
        fileInput("Protein", "Upload Protein (.pdb/.mol2):",
                  accept = c(".pdb", ".mol2")),
        shinyDirButton("ligandDir", "Select Ligand Directory",
                       "Choose directory")
        )
    } else if (input$method == "ml") {
      tagList(
        fileInput("PLI", "Upload PLI data (.csv):",
                  accept = c(".csv"), multiple = TRUE),
        selectInput( 
          "features", 
          "Select features:", 
          list("electrostatic",
               "polar_solvation",
               "non_polar_solvation",
               "vdW"), 
          multiple = TRUE 
        ),
        selectInput( 
          "models", 
          "Select models:", 
          list("RandomForest", "DecisionTree",
               "XGBoost", "LightGBM", "SVM"), 
          multiple = TRUE 
        )
      )
    }
  })
  
  observeEvent(input$runGoCode, {
    simulationMethod <- input$method
    
    if (simulationMethod == "mc") {
      
      req(input$Protein, input$ligandDir)  # Ensure both files are uploaded
      
      # Get protein file
      proteinFile <- input$Protein$datapath

      # Get all ligand files in the directory
      ligandDir <- parseDirPath(c(MCdata = file.path(getwd(), "MCdata")),
                                input$ligandDir)
      ligandFiles <- list.files(ligandDir, pattern = "\\.mol2$",
                                full.names = TRUE)
      ligandFileNames <- list.files(ligandDir, pattern = "\\.mol2$",
                                    full.names = FALSE)
      
      # Ensure ligand files exist
      if (length(ligandFiles) == 0) {
        showNotification("No .mol2 files found in the selected directory.", type = "error")
        return()
      }
    
      # Compile Go program
      compile_result <- 
        system("cd metropolis && go build -o metropolis main.go", intern = FALSE)
      # intern = FALSE if we want to print go output to R console
      
      # print(compile_result)
      
      # Run Go program
      run_result <- 
        system(paste("cd metropolis && ./metropolis", proteinFile, 
                     paste(ligandFiles, collapse = " ")), intern = FALSE)
      print(run_result)
      
      # Read outputs
      outputFile <- file.path(outputDir, "simulations.csv")
      if (file.exists(outputFile)) {
        results <- read.csv(outputFile, header = TRUE)
        
        print(results)
        
        # Add ligandFileNames as a column
        results$LigandFile <- ligandFileNames
        
        
        # Render scatter plot of the results
        output$plot <- renderPlotly({
          # Add a column to identify the point with the lowest BindingEnergy
          results <- results %>%
            mutate(
              LigandName = gsub("\\_ligand.mol2$", "", LigandFile),
              is_lowest = ifelse(BindingEnergy == min(BindingEnergy), 
                                 "Minimum Energy", "Other"),
              hover_label = 
                paste("Ligand:", LigandName, "<br>BindingEnergy:",
                      round(BindingEnergy, 4)))
          
          # Filter the data to retain only the "Lowest" point
          results_filtered <- results %>%
            filter(is_lowest == "Minimum Energy")
          
          p <- ggplot(results, aes(x = LigandFile, y = BindingEnergy,
                                   color = is_lowest, text = hover_label, group = 1)) +
            geom_point(data = results_filtered, 
                       aes(color = is_lowest, text = hover_label), 
                       size = 1.3, show.legend = FALSE) +
            geom_line(aes(group = 1), color = "black", size = 0.5) +
            labs(
              title = "Binding Energy of Multiple Ligands",
              x = "Ligand",
              y = "Binding Energy"
            ) +
            scale_color_manual(values = c("Minimum Energy" = "red")) +  
            # Set color for points
            theme_bw() +
            theme(axis.text.x = element_blank(),
                  axis.text.y = element_text(angle = 45),
                  plot.title = element_text(hjust = 0.5)) 
          
          ggplotly(p, tooltip = "text")
        })
      } else {
        warning("Output file not found.")
      }
      
    } else if (simulationMethod == "ml") {
      req(input$PLI, input$features, input$models)
      
      # Run Python script
      python_command <- paste(
        "python",  # Specify the path to executable Python
        file.path(getwd(), "machineLearning", "ml.py"),
        "--data", paste(input$PLI$datapath, collapse=" "),
        "--features", paste(input$features, collapse = " "),
        "--models", paste(input$models, collapse=" "),
        "--output_Dir", shQuote(outputDir)
      )
      
      # call the ML program
      system(python_command, intern = TRUE) 
      # if intern = FALSE, show python output at R console
      
      # Display evaluation results
      output$evaluation_results <- renderTable({
        result_evaluation <- file.path(outputDir, "evaluation_results.csv")
        if (file.exists(result_evaluation)) {
          read.csv(result_evaluation)
        } else {
          cat("Result file does not exist:", result_evaluation, "\n")
          NULL
        }
      })
      
      # Display feature importance results
      # Dynamically load feature importance for each model
      output$feature_importance <- renderUI({
       model_tables <- lapply(input$models, function(model) {
          result_fi <- file.path(outputDir, paste0("feature_importance_", model, ".csv"))
          if (file.exists(result_fi)) {
            tableOutput(paste0("feature_table_", model))
          } else {
            cat("Result file does not exist:", result_fi, "\n")
            NULL
          }
        })
        
        do.call(tagList, model_tables)
      })
      
      # Render tables for each model's feature importance
      lapply(input$models, function(model) {
        output[[paste0("feature_table_", model)]] <- renderTable({
          result_fi <- file.path(outputDir, 
                                 paste0("feature_importance_", model, ".csv"))
          if (file.exists(result_fi)) {
            read.csv(result_fi)
          } else {
            NULL
          }
        })
      })
    }
})
}

# Run the application
shinyApp(ui = ui, server = server)