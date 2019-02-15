source("http://bioconductor.org/biocLite.R")

workflowInstall <- function(pkg, ...)
{
    repos <- c(biocinstallRepos(), 
               sprintf("%s//bioconductor.org/packages/%s/workflows", 
                       BiocInstaller:::.protocol(), 
                       BiocInstaller::biocVersion()))
    
    install.packages(pkg, repos=repos, ...)
}
