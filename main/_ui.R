library(shiny)
library(bslib)
library(gargoyle)

source("modules/module_dashboard.R")
source("modules/module_overview.R")
source("modules/module_sidebar_overview.R")

ui <- page_fillable(
  title = "MOVIDApp",
  theme = bs_theme(version = 5, brand = "_brand.yml"),
  page_navbar(
    title = "MOVIDA",
    navbar_options = navbar_options(
      underline = TRUE
    ),
    nav_panel(
      title = "Overview",
      uiOutput("selected_feature"), # Placeholder for the sidebar UI
      layout_sidebar(
        sidebar = sidebar(uiOutput("sidebar_ui"), width = 450), # dynamically rendered sidebar
        mod_overview_ui("overview")
      )
    ),
    nav_panel(
      title = "Dashboard",
      mod_dashboard_ui("dashboard")
    ),
    nav_panel(
      title = "Plot Factory",
      navset_pill(
        nav_panel(
          title = "PCA",
          uiOutput("1")
        ),
        nav_panel(
          title = "Enrichment map",
          uiOutput("2")
        ),
        nav_panel(
          title = "Expression or abundance",
          uiOutput("2")
        ),
        nav_panel(
          title = "Heatmap",
          uiOutput("2")
        ),
        nav_panel(
          title = "From selected data",
          uiOutput("3")
        )
      )
    ),
    nav_panel(
      title = "Biological processes",
      uiOutput("bio_maps")
    ),
    nav_panel(
      title = "Multiomics",
      uiOutput("multiomics_dashboard")
    ),
  )
)
