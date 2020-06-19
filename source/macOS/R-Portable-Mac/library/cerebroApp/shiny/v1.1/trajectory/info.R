##----------------------------------------------------------------------------##
## Panel: Trajectory.
##----------------------------------------------------------------------------##

trajectory_projection_info <- list(
  title = "Trajectory",
  text = p("This plot shows cells projected into trajectory space, colored by the specified meta info, e.g. sample or cluster. The path of the trajectory is shown as a black line. Specific to this analysis, every cell has a 'pseudotime' and a transcriptional 'state' which corresponds to its position along the trajectory path.")
)

trajectory_density_info <- list(
  title = "Distribution along pseudotime",
  text = p("This plot shows the distribution of the variable selected above to color cells by along pseudotime. If this is a categorical variable, e.g. 'sample' or 'cluster', you will see a density plot. In contrast, if you have selected a continuous variable, e.g. nUMI or nGene, cells will be colored by the state they belong to.")
)

trajectory_number_of_cells_by_state_info <- list(
  title = "Number of cells by state",
  text = p("This table shows how many cells are assigned to each state. If the variable selected to color the cells in the projection is categorical, e.g. 'sample' or 'cluster', the number of cells in each state is also split into the subgroups.")
)

states_by_sample_info <- list(
  title = "Composition of states by sample",
  text = p("Percentage bar plot representation of the composition of states by sample. Allows to see which samples contribute most strongly to each state. Samples can be removed from the plot by clicking on them in the legend.")
)

states_by_cluster_info <- list(
  title = "Composition of states by cluster",
  text = p("Percentage bar plot representation of the composition of states by cluster. Allows to see which clusters contribute most strongly to each state. Clusters can be removed from the plot by clicking on them in the legend.")
)

states_by_cell_cycle_seurat_info <- list(
  title = "Composition of states by cell cycle (Seurat)",
  text = p("Cell cycle distribution by sample using the method embedded in the Seurat framework. For each cell, it calculates scores for both G2M and S phase based on lists of genes (see 'Analysis info' tab on the left) and assigns the cell cycle phase on the basis of these scores.")
)

states_nUMI_info <- list(
  title = "Number of transcripts by state",
  text = p("Violin plot of the number of transcripts (UMIs) found in each state.")
)

states_nGene_info <- list(
  title = "Number of expressed genes by state",
  text = p("Violin plot of the number of expressed genes found in each state.")
)
