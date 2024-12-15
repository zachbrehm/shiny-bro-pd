# Load required libraries
library(shiny)
library(ggplot2)
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
      .highlight-blue {
        background-color: #77acd3 !important;
      }
      .highlight-orange {
        background-color: #f69035 !important;
      }
      .title-panel {
        background-color: #00205B;
        color: white;
        padding: 15px;
        margin-bottom: 20px;
        position: relative;
      }
      .logo {
        position: absolute;
        bottom: 10px;
        right: 10px;
        max-height: 50px;
        max-width: 100px;
      }
    "))
  ),
  div(class = "title-panel",
      titlePanel("Parkinson's disease treatment for new patients naive to dopaminergic therapy.")
  ),
  HTML("<p><b>Welcome to the First-Line therapy for Parkinson's Disease (FLIP) tool</b>.<br>This tool is designed to provide data-driven, patient-specific decisions for optimal first-line therapy for patients with early Parkinson's Disease. This tool should be considered a prototype only and not be used in practice at this time.</p>"),
  p("Please refer to the tables below to see relevant summary statistics from the training data. To input patient data and view results, please navigate to the second tab."),
  tabsetPanel(
    tabPanel("Baseline Covariates and Treatment Effects",
             h3("Baseline Covariates"),
             p("The following table summarizes the baseline covariates in the training data. These data were obtained through the CALM-PD trial."),
             htmlOutput("baseline_html"),
             h3("Treatment effects"),
             p("This table summarizes the estimates of the treatment effects found in the training data. The benefit outcome in the analysis was the reduction in MDS-UPDRS part 3 scores and the risk outcome was the proportion of patients in treatment groups that experienced a dopaminergic event."),
             p("Patients in both treatment groups saw an improvement in motor function on average. However patients taking levodopa typically experienced a greater improvement than those taking DRA."),
             p("This additional benefit came at the cost of a higher likelihood of experiencing a dopaminergic event."),
             htmlOutput("tmle_html"),
    ),
    tabPanel("Patient Data and Results",
             h3("Input Patient Data and View Results"),
             p("Choose how you want to input patient data, either by uploading a csv file for multiple patients or manually inputting data for a single patient."),
             p("If uploading a csv, data should have one row per patient with baseline covariates in columns ordered by patients' sex, BMI (kg), age (years), time from diagnosis (years), MDS-UPDRS Part 2, and MDS-UPDRS Part 3."),
             radioButtons("input_type", "Input Type:",
                          choices = c("Manual Input" = "manual", "Upload CSV" = "csv")),
             conditionalPanel(
               condition = "input.input_type == 'manual'",
               h3("Patient data"),
               fluidRow(
                 column(2, 
                        numericInput("sex", "Sex (1 for male, 0 for female)", value = 1),
                        numericInput("bmi", "BMI (kg)", value = 28)
                 ),
                 column(2, 
                        numericInput("age", "Age (years)", value = 60),
                        numericInput("time_from_diagnosis", "Time from diagnosis (years)", value = 1.5)
                 ),
                 column(2, 
                        numericInput("updrs_part2", "MDS-UPDRS Part 2", value = 10),
                        numericInput("updrs_part3", "MDS-UPDRS Part 3", value = 30)
                 )
               )
             ),
             conditionalPanel(
               condition = "input.input_type == 'csv'",
               fileInput("file", "Upload your CSV file")
             ),
             hr(),
             h3("Risk tolerance range"),
             fluidRow(
               column(2, numericInput("tau_min", "Minimum risk tolerance", value = 0.35, min = 0.31, max = 0.52, step = 0.01)),
               column(2, numericInput("tau_max", "Maximum risk tolerance", value = 0.50, min = 0.31, max = 0.52, step = 0.01)),
               column(2, numericInput("tau_increment", "Risk increment", value = 0.05, min = 0.01, max = 0.21, step = 0.01))
             ),
             actionButton("process", "Process Data"),
             downloadButton("downloadData", "Download Results"),
             hr(),
             h3("Results"),
             #p("Results table displays the proportion of decision outcomes for each risk tolerance that assign levodopa to the patient."),
             #p("Estimates that are below 0.5 are shaded blue to indicate that DRA is the preferred treatment, while values above 0.5 are shaded orange to indicate levodopa is the preferred treatment"),
             #p("Estimates whose standard errors overlap 0.5 are uncertain and not shaded."),
             tableOutput("results"),
             hr(),
             h3("Patient-specific Graph"),
             p("To view a plot of the estimates and their error bars for a specific patient, please select their ID below."),
             selectInput("patient_id", "Select Patient ID", choices = NULL),
             plotOutput("patient_graph")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
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
        tau_ind <- which.min(abs(ref_tau_seq - tau))
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
    
    mean_decisions <- decisions |> 
      group_by(id, imp_ind, tau) |> 
      summarize(imp_mean = mean(decision == -1),
                imp_se = sqrt(imp_mean * (1 - imp_mean) / n()),
                .groups = 'drop') |> 
      group_by(id, tau) |>
      summarize(
        p_mean = mean(imp_mean),
        p_lower = if(all(imp_mean == 1)) 1 else max(0, quantile(imp_mean, 0.025)),
        p_upper = if(all(imp_mean == 0)) 0 else min(1, quantile(imp_mean, 0.975)),
        .groups = 'drop'
      )
    
    results_long <- mean_decisions
    results_long
  }
  
  results_long <- eventReactive(input$process, {
    withProgress(message = 'Processing data, please wait...', {
      get_decision(user_data(), tau_seq())
    })
  })
  
  output$results <- renderTable({
    req(results_long())
    res_long <- results_long()
    
    # Ensure all necessary columns exist
    required_cols <- c("id", "tau", "p_mean", "p_lower", "p_upper")
    missing_cols <- setdiff(required_cols, names(res_long))
    if (length(missing_cols) > 0) {
      stop("Missing columns in results: ", paste(missing_cols, collapse = ", "))
    }
    
    # Select and rename columns
    res <- res_long %>%
      select(ID = id, 
             `Risk Tolerance` = tau, 
             `Probability of LDOPA` = p_mean, 
             `Lower CI Bound` = p_lower, 
             `Upper CI Bound` = p_upper)
    
    # Round numeric columns to 3 decimal places
    res %>% mutate(across(c(`Probability of LDOPA`, `Lower CI Bound`, `Upper CI Bound`), ~round(., 3)))
  }, digits = 3)
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("levodopa_probability_results_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(results(), file, row.names = FALSE)
    }
  )
  
  # Update patient ID choices based on results
  observe({
    req(results_long())
    patient_ids <- unique(results_long()$id)
    updateSelectInput(session, "patient_id", choices = patient_ids)
  })
  
  # Placeholder for patient-specific graph
  output$patient_graph <- renderPlot({
    req(input$patient_id)
    dat <- filter(results_long(), id == input$patient_id) |> 
      mutate(
        pick = case_when(
          p_mean < 0.5 & p_upper < 0.5 ~ "DRA",
          p_mean > 0.5 & p_lower > 0.5 ~ "LDPA",
          TRUE ~ "Undecided"
        ),
        pick = factor(pick, levels = c("DRA", "LDPA", "Undecided"))
      )
    
    ggplot(dat) + 
      geom_point(aes(x = tau, y = p_mean, color = pick), size = 3) + 
      geom_errorbar(aes(x = tau, ymin = pmax(p_lower, 0), ymax = pmin(p_upper, 1), color = pick), width = 0.01) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) + 
      scale_color_manual(values = c("DRA" = "#4f81af", "LDPA" = "#d45b21", "Undecided" = "#000000")) + 
      theme_minimal() + 
      labs(
        title = paste0("Summary of treatment assignments for patient ", input$patient_id),
        subtitle = "Values below 0.5 favor DRA and values above 0.5 favor Levodopa",
        x = "Risk tolerance", 
        y = "Probability of assigning Levodopa",
        color = "Decision"
      ) +
      theme(legend.position = "bottom")
  })
}

# Run the app
shinyApp(ui = ui, server = server)