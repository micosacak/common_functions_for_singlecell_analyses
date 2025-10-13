
{
# Set options to avoid interactive prompts
options(repos = "https://cloud.r-project.org/") # Set a default CRAN mirror
options(BiocManager.check_repositories = FALSE) # Avoid BiocManager prompts
options(install.packages.check.source = "no") # Skip source checks
options(warn = -1) # Suppress warnings

# Install CRAN packages
install.packages(c(
  "BiocManager", "devtools", "doMC", "Seurat", "remotes",
  "ggraph", "Signac", "R.utils", "ggseqlogo", "pheatmap",
  "viridis", "ggforce", "scCustomize", "NMF", "anndata"
), dependencies = TRUE, ask = FALSE)

# Install Bioconductor packages
BiocManager::install(c(
  "GenomeInfoDb", "GenomicRanges", "IRanges", "Rsamtools",
  "BiocGenerics", "GenomicFeatures", "AnnotationDbi",
  "SingleCellExperiment", "AnnotationHub", "TFBSTools",
  "JASPAR2024", "JASPAR2020", "biomaRt",
  "BSgenome.Drerio.UCSC.danRer11", "ComplexHeatmap",
  "motifmatchr", "chromVAR", "Nebulosa", "cicero",
	"scDblFinder","pcaMethods"
), update = FALSE, ask = FALSE)

# Install GitHub packages
devtools::install_github("caleblareau/BuenColors", upgrade = "never")
devtools::install_github("buenrostrolab/FigR", upgrade = "never")
devtools::install_github("jokergoo/circlize", upgrade = "never")
devtools::install_github("jokergoo/ComplexHeatmap", upgrade = "never")
devtools::install_github("jinworks/CellChat", upgrade = "never")
devtools::install_github('kharchenkolab/pagoda2', upgrade = "never")
devtools::install_github("velocyto-team/velocyto.R", upgrade = "never")

# Install Seurat-related packages from GitHub
remotes::install_github("satijalab/seurat-wrappers", ref = "seurat5", quiet = TRUE, upgrade = "never")
remotes::install_github("satijalab/seurat-data", ref = "seurat5", quiet = TRUE, upgrade = "never")
remotes::install_github("Moonerss/scrubletR", upgrade = "never")
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", force = TRUE, upgrade = "never")

# Reset warnings to default
options(warn = 0)
}
