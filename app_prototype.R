# Load required libraries
library(shiny)
library(ggthemes)
library(foreach)
library(dplyr)
library(matrixStats)
library(readr)
library(tidyr)

# Load training data and betas
training_data <- readRDS("data/training_data.Rds")
training_betas <- readRDS("data/training_betas.Rds")

tau_seq <- round(seq(0.31, 0.52, 0.01), digits = 2)
n_tau <- length(tau_seq)

# Helper functions
standardize <- function(x) {
  return((x - mean(x)) / sd(x))
}

decision_test <- function(beta_train, X_train, X_test){
  X_train <- as.matrix(X_train)
  X_test <- as.matrix(X_test)
  N_train <- nrow(X_train)
  N_test <- nrow(X_test)
  beta0 <- beta_train[1]
  beta_vals <- beta_train[-1]
  decision_out <- vector(length = nrow(X_test))
  for(i in 1:N_test){
    x_prod <- c(X_test[i,] %*% t(X_train))
    decision_out[i] <- beta0 + sum(beta_vals*x_prod)
  }
  return(decision = decision_out)
}

# Define UI
ui <- fluidPage(
  titlePanel("Parkinson's disease treatment for new patients naive to dopaminergic therapy."),
  p("Upload patient data to determine the optimal first line treatment for Parkinson's Disease between LDOPA or DRA. Please follow these steps:"),
  tags$ol(
    tags$li("Upload your CSV file using the file input below."),
    tags$li("Data should have one row per patient with baseline covariates in columns ordered by patients' sex, BMI (kg), age (years), time from diagnosis (years), MDS-UPDRS Part 2, and MDS-UPDRS Part 3."),
    tags$li("Click the 'Process Data' button to analyze your data."),
    tags$li("View the results in the table that will appear."),
    tags$li("Use the 'Download Results' button to save the results as a CSV file.")
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload your CSV file"),
      actionButton("process", "Process Data"),
      downloadButton("downloadData", "Download Results")
    ),
    mainPanel(
      h3("Probability that LDOPA is the optimal treatment. Larger values favor LDOPA, smaller values favor DRA."),
      tableOutput("results")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  user_data <- reactive({
    req(input$file)
    read_csv(input$file$datapath)
  })
  
  get_decision <- function(user_data) {
    decisions <- foreach(tau_ind = 1:n_tau, .combine = rbind) %:%
      foreach(imp_ind = 1:10, .combine = rbind) %:%
      foreach(fold_ind = 1:100, .combine = rbind) %do% {
        train_dat <- training_data[[imp_ind]][[fold_ind]]
        train_means <- colMeans(train_dat[,-1])
        train_sds <- colSds(as.matrix(train_dat[,-1]))
        user_data_std <- t(apply(user_data[,-1], 1, function(x) (x - train_means) / train_sds))
        X_train <- as.matrix(mutate(as.data.frame(train_dat), across(-1, standardize)))
        X_test <- cbind(user_data[,1], user_data_std)
        decision <- sign(decision_test(training_betas[[imp_ind]][[tau_ind]][[fold_ind]], X_train, X_test))
        data.frame(id = 1:nrow(user_data),
                   tau = tau_seq[tau_ind],
                   imp_ind = imp_ind,
                   fold_ind = fold_ind,
                   decision = decision)
      }
    
    mean_decisions <- decisions %>%
      group_by(id, tau) %>%
      summarize(mean_decision = mean(decision == -1), .groups = 'drop')
    
    pivot_wider(mean_decisions, id_cols = id, names_from = tau, values_from = mean_decision)
  }
  
  results <- eventReactive(input$process, {
    withProgress(message = 'Processing data, please wait...', {
      get_decision(user_data())
    })
  })
  
  output$results <- renderTable({
    results()
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("levodopa_probability_results_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(results(), file, row.names = FALSE)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)