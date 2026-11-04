### SHINY APP: LM-BASED REFERENCE INTERVALS ###
# RIs for any combination of age, sex and menopausal status, from the LM constants written by 05_reference_intervals
# (lm_constants.csv); online at https://nicogo.shinyapps.io/ricalculation_app/

library(shiny)
library(dplyr)

constants <- read.csv("lm_constants.csv")


# UI
ui <- fluidPage(
  titlePanel("Metabolite Reference Interval Calculator"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("age", "Age:", value = 50, min = 40, max = 76),               # age range of the reference population
      selectInput("sex", "Sex:",
                  choices = c("Male", "Female")),
      uiOutput("menopause_input")
    ),
    mainPanel(
      h3("Amino Acids"),
      tableOutput("amino_acids"),
      h3("Acylcarnitines"),
      tableOutput("acylcarnitines")
    )
  )
)

# Server
server <- function(input, output, session) {

  calculateReferenceIntervals <- function(age, menopause, sex) {

    constants %>%
      mutate(prediction = Intercept +
               age * coeff_age +
               (sex == "Female") * coeff_sex +
               (menopause == "Post-Menopause") * coeff_menopausal_status +
               (sex == "Female") * age * coeff_age_sex +
               (menopause == "Post-Menopause") * age * coeff_age_menopausal_status,
             `lower border` = pmax(round(prediction + lower_border, digits = 2), 0),
             `upper border` = round(prediction + upper_border, digits = 2),
             unit = ifelse(grepl("carnitine", metabolite, ignore.case = TRUE), "[nM]", "[μM]")) %>%
      select(metabolite, unit, `lower border`, `upper border`)
  }

  # Render the reference intervals based on the input values
  output$amino_acids <- renderTable({
    req(input$menopause)
    calculateReferenceIntervals(input$age, input$menopause, input$sex) %>%
      filter(unit == "[μM]")
  })
  output$acylcarnitines <- renderTable({
    req(input$menopause)
    calculateReferenceIntervals(input$age, input$menopause, input$sex) %>%
      filter(unit == "[nM]")
  })

  # Generate menopause status input dynamically
  output$menopause_input <- renderUI({
    if (input$sex == "Male") {
      selectInput("menopause", "Menopausal Status:", choices = c("Male"), selected = "Male")
    } else {
      selectInput("menopause", "Menopausal Status:", choices = c("Pre-Menopause", "Post-Menopause"))
    }
  })
}

# Run the app
shinyApp(ui = ui, server = server)
