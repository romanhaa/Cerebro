##----------------------------------------------------------------------------##
## Tab: Trajectory.
##----------------------------------------------------------------------------##

# elements needed
# - user interface
# - plot 1: trajectory
# - plot 2: meta data over pseudotime

tab_trajectory <- tabItem(
  tabName = "trajectory",
  uiOutput("trajectory_UI")
)

