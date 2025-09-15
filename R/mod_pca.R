#' module_pca UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_module_pca_ui <- function(id) {
  ns <- NS(id)
  tagList(
 
  )
}
    
#' module_pca Server Functions
#'
#' @noRd 
mod_module_pca_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_module_pca_ui("module_pca_1")
    
## To be copied in the server
# mod_module_pca_server("module_pca_1")
