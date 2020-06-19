##----------------------------------------------------------------------------##
## Tab: Color management.
##----------------------------------------------------------------------------##

tab_color_management <- tabItem(
  tabName = "color_management",
  box(
    title = tagList(
      boxTitle("Samples"),
      cerebroInfoButton("colors_for_samples_info")
    ),
    status = "primary",
    solidHeader = TRUE,
    width = 6,
    collapsible = TRUE,
    uiOutput("colors_for_samples_UI")
  ),
  box(
    title = tagList(
      boxTitle("Clusters"),
      cerebroInfoButton("colors_for_clusters_info")
    ),
    status = "primary",
    solidHeader = TRUE,
    width = 6,
    collapsible = TRUE,
    uiOutput("colors_for_clusters_UI")
  )
)
