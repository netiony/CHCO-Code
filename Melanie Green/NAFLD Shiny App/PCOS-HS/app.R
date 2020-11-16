#-------------------------------------------------------------------------------
#
# This Shiny app takes inputs based on a LRM developed by Laura Pyle and returns 
# the likelihood of developing NAFLD. Requires BMI percentile (or height, weight
# and age to calculate it), a waist measurement, ALT, and SHBG. 
#
# v 2.0 
# Tim Vigers and Laura Pyle 
# 9/5/18
#
#-------------------------------------------------------------------------------

# Load the libraries.
library(shiny)
library(AGD)

ui <- fluidPage(
  helpText("ROC analysis was used to identify the cutoff of the probability that maximized the Youden Index. A threshold of 43% provides 82% sensitivity and 69% specificity in a training cohort and 91% sensitivity and 70% specificity in a validation cohort. This model has so far only been validated in adolescent girls with PCOS."),
  # Put inputs in a grid at the bottom of the page.   
  fluidRow(
    column(6,
           # BMI percentile input.      
           numericInput(inputId = "bmi.perc",
                        label = "BMI Percentile",
                        value = 0,min = 0,max = 100),
           # Height input.
           numericInput(inputId = "height",
                        label = "Height",
                        value = 0,min = NA,max = NA),
           # Option for height measurement units.
           radioButtons(inputId = "height_units",
                        label = NULL,
                        c("cm" = 1, "in" = 2),
                        inline = TRUE),
           # Weight input.
           numericInput(inputId = "weight",
                        label = "Weight",
                        value = 0,min = NA,max = NA),
           # Weight measurement units.
           radioButtons(inputId = "weight_units",
                        label = NULL,
                        c("kg"  = 1,"lbs" = 2),
                        inline = TRUE),
           # Age input.
           numericInput(inputId = "age",
                        label = "Age (years)",
                        value = 0,min = NA,max = NA)),
    # New column.
    column(width = 6,
           # Waist measurement input.      
           numericInput(inputId = "Waist_cm",
                        label = "Waist",
                        value = 0,min = NA,max = NA),
           # Waist measurement units.      
           radioButtons(inputId = "Waist_units",
                        label = NULL,
                        c("cm" = 1,"in" = 2),
                        inline = TRUE),
           # ALT input.      
           numericInput(inputId = "ALT",
                        label = "ALT",
                        value = 0,min = NA,max = NA),
           # SHBG input.      
           numericInput(inputId = "SHBG",
                        label = "SHBG",
                        value = 0,min = NA,max = NA),
           # SHBG unit selection.
           radioButtons(inputId = "SHBG_units",
                        label = "SHBG Units",
                        c("nmol/L" = 1, "ug/dL" = 2),
                        inline = TRUE),
           # Run button.
           actionButton(inputId = "run",
                        label = "Run")
    )),
  # Output panel on righthand side. 
  mainPanel(
    titlePanel("NAFLD Probability (%)"),
    textOutput("highALT"),
    textOutput("highprobability"),
    textOutput("lowprobability"),
    tags$head(tags$style("#highALT{color: red;
                         font-size: 50px;
                         font-style: italic;
                         }"
                         )
    ),
    tags$head(tags$style("#highprobability{color: red;
                         font-size:200px;
                         font-style: italic;
                         }"
                         )
    ),
    tags$head(tags$style("#lowprobability{color: green;
                         font-size: 200px;
                         font-style: italic;
                         }"
                         )
  )
    )
    )
server <- function(input, output, session) {
  # Calculate probability based on reduced model.  
  prob <- eventReactive(input$run, {
    # Define coefficients.
    intercept <- 25.18765308
    bmicoeff <- -0.34114772
    waistcoeff <- 0.06148684
    ALTcoeff <- 0.09373780
    SHBGcoeff <- -0.07953961
    # Convert units.    
    if (input$height_units == 2) {
      height <- input$height * 2.54
    } else {
      height <- input$height
    }
    if (input$weight_units == 2) {
      weight <- input$weight * 0.45359237
    } else {
      weight <- input$weight
    }
    if (input$Waist_units == 2) {
      Waist_cm <- input$Waist_cm * 2.54
    } else {
      Waist_cm <- input$Waist_cm
    }
    if (input$SHBG_units == 2) {
      SHBG <- input$SHBG / 34.6741
    } else {
      SHBG <- input$SHBG
    }
    # Calculate BMI percentile from height, weight, and age.
    if (input$height != 0 && input$weight != 0 && input$age != 0) {
      bmi <- weight / ((height / 100)^2)
      bmiperc <- 100*round(pnorm(y2z(y=bmi,x=input$age,sex="F",ref=cdc.bmi)),4)
    } else {
      bmiperc <- input$bmi.perc
    }
    (exp(intercept+
           (bmicoeff * bmiperc)+
           (waistcoeff * input$Waist_cm)+
           (ALTcoeff * input$ALT)+
           (SHBGcoeff * SHBG)))/
      (1 + exp(intercept+
                 (bmicoeff * bmiperc)+
                 (waistcoeff * input$Waist_cm)+
                 (ALTcoeff * input$ALT)+
                 (SHBGcoeff * SHBG)))
  })
  # Reactive high ALT warning.
  highALT <- eventReactive(input$run,{
    if (!is.null(prob()) && input$ALT != "" && input$ALT >= 44) {
      print("ALT at or above 44, further screening recommended.")
    } else {
      print("")
    }
  })
  # Reactive green output for low probability (below the 0.43 cutoff).  
  lowprob <- eventReactive(input$run,{
    if (!is.null(prob()) && as.numeric(prob()) < 0.43 && input$ALT < 44) {
      round((prob() * 100),digits = 2)
    }
  })
  # Reactive red output for high probability (above the 0.43 cutoff).
  highprob <- eventReactive(input$run,{
    if (!is.null(prob()) && as.numeric(prob()) >= 0.43 && input$ALT < 44) {
      round((prob() * 100),digits = 2)
    }
  })
  # Print low probability.
  output$lowprobability <- renderText({
    lowprob()
  })
  # Print high probability.  
  output$highprobability <- renderText({
    highprob()
  })
  # Print high ALT warning.
  output$highALT <- renderText({
    highALT()
  })
}
shinyApp(ui, server)