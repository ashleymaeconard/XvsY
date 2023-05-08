# install all of the CRAN Packages
install.packages('devtools', repos = "http://cran.us.r-project.org")
install.packages('cowplot', repos = "http://cran.us.r-project.org")
remotes::install_version('ggplot2', version='3.3.6')
# install.packages('ggplot2', repos = "http://cran.us.r-project.org")
install.packages('stringr', repos = "http://cran.us.r-project.org")
install.packages('scales', repos = "http://cran.us.r-project.org")
install.packages('ggfortify', repos = "http://cran.us.r-project.org")
install.packages('dplyr', repos = "http://cran.us.r-project.org")
install.packages('VennDiagram', repos = "http://cran.us.r-project.org")

options(BioC_mirror = "http://bioconductor.org", repos = "http://cran.r-project.org")

# First Install all Bioconductor R Packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.Dm.eg.db")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("annotate")
