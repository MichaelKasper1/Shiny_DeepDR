
Pred_Gen <- function(mut_data, exp_data, model){
  sample_final <- list(as.matrix(mut_data), as.matrix(exp_data))
  predictions <- predict(model, sample_final)
  return(data.frame(predictions))
}

ui <- fluidPage(
  "Hello, world!"
)
server <- function(input, output, session) {
}
shinyApp(ui, server)
