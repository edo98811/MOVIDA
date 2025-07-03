

ui <- page_fillable(
  title = "MOVIDApp",
  theme = bs_theme(version = 5, brand = "inst/extdata/_brand.yml"),
  page_navbar(
    title = "MOVIDA",
    navbar_options = navbar_options(
      underline = TRUE
    ),
    nav_panel(
      title = "Overview",
      uiOutput("selected_feature"), # Placeholder for the sidebar UI
      layout_sidebar(
        sidebar = sidebar(uiOutput("sidebar_ui"), width = 460, position = "right"), # dynamically rendered sidebar
        mod_overview_ui("overview")
      )
    ),
    nav_panel(
      title = "Dashboard",
      mod_dashboard_ui("dashboard")
    ),
    nav_panel(
      title = "Plot Factory",
      mod_plotting_ui("plotting")
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
