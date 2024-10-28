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
ref_tau_seq <- round(seq(0.31, 0.52, 0.01), digits = 2)

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
  tags$head(
    tags$style(HTML("
      .highlight {
        background-color: yellow !important;
      }
    "))
  ),
  titlePanel("Parkinson's disease treatment for new patients naive to dopaminergic therapy."),
  p("Welcome to the First-Line therapy for Parkinson's Disease (FLIP) tool. This tool is designed to provide data-driven, patient-specific decisions for optimal first-line therapy for patients with early Parkinson's Disease. For further details on the methodology used for the decision rule, please refer to LINK TO PAPER."),
  tabsetPanel(
    tabPanel("Baseline Covariates and Treatment Effects",
             h3("Baseline Covariates"),
             p("The following table summarizes the baseline covariates in the training data. These data were obtained through the CALM-PD trial."),
             htmlOutput("baseline_html"),
             h3("Treatment effects"),
             p("This table summarizes the estimates of the treatment effects found in the training data. The benefit outcome in the analysis was the reduction in MDS-UPDRS part 3 scores and the risk outcome was the proportion of patients in treatment groups that experienced a dopaminergic event."),
             htmlOutput("tmle_html"),
    ),
    tabPanel("Patient Data and Results",
             h3("Input Patient Data and View Results"),
             p("Choose how you want to input patient data:"),
             radioButtons("input_type", "Input Type:",
                          choices = c("Upload CSV" = "csv", "Manual Input" = "manual")),
             conditionalPanel(
               condition = "input.input_type == 'csv'",
               fileInput("file", "Upload your CSV file"),
               p("Data should have one row per patient with baseline covariates in columns ordered by patients' sex, BMI (kg), age (years), time from diagnosis (years), MDS-UPDRS Part 2, and MDS-UPDRS Part 3.")
             ),
             conditionalPanel(
               condition = "input.input_type == 'manual'",
               numericInput("sex", "Sex (1 for male, 0 for female)", value = NULL),
               numericInput("bmi", "BMI (kg)", value = NULL),
               numericInput("age", "Age (years)", value = NULL),
               numericInput("time_from_diagnosis", "Time from diagnosis (years)", value = NULL),
               numericInput("updrs_part2", "MDS-UPDRS Part 2", value = NULL),
               numericInput("updrs_part3", "MDS-UPDRS Part 3", value = NULL)
             ),
             fluidRow(
               column(4, numericInput("tau_min", "Minimum tau", value = 0.31, min = 0.31, max = 0.52, step = 0.01)),
               column(4, numericInput("tau_max", "Maximum tau", value = 0.52, min = 0.31, max = 0.52, step = 0.01)),
               column(4, numericInput("tau_increment", "Tau increment", value = 0.01, min = 0.01, max = 0.21, step = 0.01))
             ),
             actionButton("process", "Process Data"),
             downloadButton("downloadData", "Download Results"),
             hr(),
             h3("Results"),
             p("Probability that LDOPA is the optimal treatment. Larger values favor LDOPA, smaller values favor DRA."),
             tableOutput("results")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  output$baseline_html <- renderUI({
    HTML(readLines("data/calm_baseline_covariates.html"))
  })
  
  output$tmle_html <- renderUI({
    HTML(readLines("data/calm_tmle_estimates.html"))
  })
  
  user_data <- reactive({
    if (input$input_type == "csv") {
      req(input$file)
      read_csv(input$file$datapath)
    } else {
      req(input$sex, input$bmi, input$age, input$time_from_diagnosis, input$updrs_part2, input$updrs_part3)
      data.frame(
        sex = input$sex,
        bmi = input$bmi,
        age = input$age,
        time_from_diagnosis = input$time_from_diagnosis,
        updrs_part2 = input$updrs_part2,
        updrs_part3 = input$updrs_part3
      )
    }
  })
  
  tau_seq <- reactive({
    req(input$tau_min, input$tau_max, input$tau_increment)
    round(seq(input$tau_min, input$tau_max, by = input$tau_increment), digits = 2)
  })
  
  get_decision <- function(user_data, tau_sequence) {
    decisions <- foreach(tau = tau_sequence, .combine = rbind) %:%
      foreach(imp_ind = 1:10, .combine = rbind) %:%
      foreach(fold_ind = 1:100, .combine = rbind) %do% {
        tau_ind <- which.min(abs(round(seq(0.31, 0.52, 0.01), digits = 2) - tau))
        train_dat <- training_data[[imp_ind]][[fold_ind]]
        train_means <- colMeans(train_dat[,-1])
        train_sds <- colSds(as.matrix(train_dat[,-1]))
        user_data_std <- t(apply(user_data[,-1], 1, function(x) (x - train_means) / train_sds))
        X_train <- as.matrix(mutate(as.data.frame(train_dat), across(-1, standardize)))
        X_test <- cbind(user_data[,1], user_data_std)
        decision <- sign(decision_test(training_betas[[imp_ind]][[tau_ind]][[fold_ind]], X_train, X_test))
        data.frame(id = 1:nrow(user_data),
                   tau = tau,
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
      get_decision(user_data(), tau_seq())
    })
  })
  
  output$results <- renderTable({
    req(results())
    res <- results()
    
    # Function to apply highlighting
    highlight_first_ge_0.5 <- function(row) {
      indices <- which(row >= 0.5)
      if (length(indices) > 0) {
        first_index <- min(indices)
        row[first_index] <- sprintf('<span class="highlight">%.3f</span>', row[first_index])
      }
      return(row)
    }
    
    # Apply highlighting to each row
    res[-1] <- t(apply(res[-1], 1, highlight_first_ge_0.5))
    
    res
  }, sanitize.text.function = function(x) x, escape = FALSE)
  
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