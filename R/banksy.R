#!/usr/bin/env Rscript
# =============================================================================
# BANKSY Spatial Transcriptomics Clustering Analysis
# =============================================================================
#
# DESCRIPTION:
#   This script performs BANKSY (Building Aggregates with a Neighborhood Kernel
#   and Spatial Yardstick) analysis on spatial transcriptomics data. BANKSY
#   unifies cell typing and tissue domain segmentation by incorporating spatial
#   neighborhood information into the clustering process.
#
# CLINICAL CONTEXT:
#   Spatial transcriptomics enables gene expression profiling while preserving
#   tissue architecture. This analysis identifies spatially coherent cell
#   populations and tissue domains, which is valuable for:
#   - Tumor microenvironment characterization
#   - Tissue architecture analysis
#   - Cell type identification with spatial context
#
# AUTHOR:
#   Vladimir Kovacevic
#
# VERSION: 1.0.0
#
# REFERENCES:
#   Singhal et al. (2024) Nature Genetics. BANKSY unifies cell typing and
#   tissue domain segmentation for scalable spatial omics data analysis.
#   https://prabhakarlab.github.io/Banksy/
#
# USAGE:
#   Rscript banksy.R --input <path_to_h5ad> --output <output_dir> [options]
#
# EXAMPLES:
#   # Basic analysis with default parameters
#   Rscript banksy.R --input data/sample.h5ad --output results/
#
#   # Cell-typing focused analysis (lambda=0.2)
#   Rscript banksy.R --input data/sample.h5ad --output results/ --lambda 0.2
#
#   # Domain segmentation (lambda=0.8)
#   Rscript banksy.R --input data/sample.h5ad --output results/ --lambda 0.8
#
#   # Multiple lambda values for comparison
#   Rscript banksy.R --input data/sample.h5ad --output results/ --lambda 0.2,0.8
#
# =============================================================================

# =============================================================================
# 1. COMMAND-LINE ARGUMENT PARSING
# =============================================================================

#' Parse command-line arguments using base R
#'
#' This function provides robust argument parsing without external dependencies,
#' ensuring portability across different R environments (local, Docker, HPC).
#'
#' @return Named list of parsed arguments with defaults applied
parse_arguments <- function() {
    # Get raw command line arguments
    raw_args <- commandArgs(trailingOnly = TRUE)

    # Define default values for all parameters
    defaults <- list(
        input        = NULL,
        output       = NULL,
        sample_name  = NULL,
        lambda       = "0.2,0.8",
        k_geom       = "15,30",
        resolution   = 0.8,
        n_pcs        = 20,
        qc_lower     = 0.01,
        qc_upper     = 0.99,
        min_gene_cells = 0.01,
        point_size   = 0.5,
        plot_dpi     = 150,
        seed         = 42,
        verbose      = TRUE,
        quiet        = FALSE,
        skip_plots   = FALSE,
        skip_umap    = FALSE
    )

    # Check for help flag
    if ("--help" %in% raw_args || "-h" %in% raw_args) {
        cat("
BANKSY Spatial Transcriptomics Clustering Analysis
===================================================

USAGE:
  Rscript banksy.R --input <FILE> --output <DIR> [OPTIONS]

REQUIRED ARGUMENTS:
  --input, -i <FILE>       Path to input h5ad file with spatial transcriptomics data
  --output, -o <DIR>       Output directory for results (created if needed)

OPTIONAL ARGUMENTS:
  Sample Identification:
    --sample-name <NAME>   Sample name for output files [default: input filename]

  BANKSY Parameters:
    --lambda <VALUES>      Lambda parameter(s) controlling spatial weight [default: 0.2,0.8]
                           0 = non-spatial, 0.2 = cell-typing, 0.8 = domain segmentation
                           Multiple values: comma-separated (e.g., '0.2,0.8')
    --k-geom <VALUES>      Neighbors for feature computation [default: 15,30]
    --resolution <NUM>     Leiden clustering resolution [default: 0.8]
    --n-pcs <INT>          Number of principal components [default: 20]

  Quality Control:
    --qc-lower <NUM>       Lower quantile for QC filtering [default: 0.01]
    --qc-upper <NUM>       Upper quantile for QC filtering [default: 0.99]
    --min-gene-cells <NUM> Min fraction of cells for gene retention [default: 0.01]

  Visualization:
    --point-size <NUM>     Point size for spatial plots [default: 0.5]
    --plot-dpi <INT>       DPI resolution for plots [default: 150]
    --skip-plots           Skip visualization generation
    --skip-umap            Skip UMAP computation

  Other:
    --seed <INT>           Random seed for reproducibility [default: 42]
    --verbose, -v          Enable verbose output [default: TRUE]
    --quiet, -q            Suppress most output
    --help, -h             Show this help message

EXAMPLES:
  # Basic analysis
  Rscript banksy.R --input data/sample.h5ad --output results/

  # Cell-typing focused (lambda=0.2)
  Rscript banksy.R --input data/sample.h5ad --output results/ --lambda 0.2

  # Domain segmentation (lambda=0.8)
  Rscript banksy.R --input data/sample.h5ad --output results/ --lambda 0.8

  # Multiple lambda values
  Rscript banksy.R --input data/sample.h5ad --output results/ --lambda 0.2,0.8

REFERENCES:
  Singhal et al. (2024) Nature Genetics
  https://prabhakarlab.github.io/Banksy/

")
        quit(save = "no", status = 0)
    }

    # Parse arguments
    args <- defaults
    i <- 1
    while (i <= length(raw_args)) {
        arg <- raw_args[i]

        # Handle flags (no value needed)
        if (arg %in% c("--verbose", "-v")) {
            args$verbose <- TRUE
            i <- i + 1
            next
        }
        if (arg %in% c("--quiet", "-q")) {
            args$quiet <- TRUE
            i <- i + 1
            next
        }
        if (arg == "--skip-plots") {
            args$skip_plots <- TRUE
            i <- i + 1
            next
        }
        if (arg == "--skip-umap") {
            args$skip_umap <- TRUE
            i <- i + 1
            next
        }

        # Handle key=value or key value pairs
        if (grepl("^--", arg) || grepl("^-[a-z]$", arg)) {
            # Check if value is provided
            if (i + 1 > length(raw_args)) {
                stop(sprintf("ERROR: Argument %s requires a value", arg))
            }
            value <- raw_args[i + 1]

            # Map argument to parameter name
            param <- switch(arg,
                "--input" = , "-i" = "input",
                "--output" = , "-o" = "output",
                "--sample-name" = "sample_name",
                "--lambda" = , "-l" = "lambda",
                "--k-geom" = "k_geom",
                "--resolution" = , "-r" = "resolution",
                "--n-pcs" = , "-p" = "n_pcs",
                "--qc-lower" = "qc_lower",
                "--qc-upper" = "qc_upper",
                "--min-gene-cells" = "min_gene_cells",
                "--point-size" = "point_size",
                "--plot-dpi" = "plot_dpi",
                "--seed" = "seed",
                NULL
            )

            if (is.null(param)) {
                stop(sprintf("ERROR: Unknown argument: %s", arg))
            }

            # Convert to appropriate type
            if (param %in% c("resolution", "qc_lower", "qc_upper", "min_gene_cells", "point_size")) {
                args[[param]] <- as.numeric(value)
            } else if (param %in% c("n_pcs", "plot_dpi", "seed")) {
                args[[param]] <- as.integer(value)
            } else {
                args[[param]] <- value
            }

            i <- i + 2
        } else {
            stop(sprintf("ERROR: Unexpected argument: %s", arg))
        }
    }

    return(args)
}

# Parse command-line arguments
args <- parse_arguments()

# =============================================================================
# 2. ARGUMENT VALIDATION
# =============================================================================

# Check required arguments
if (is.null(args$input)) {
    stop("ERROR: Input file (--input) is required. Use --help for usage information.")
}
if (is.null(args$output)) {
    stop("ERROR: Output directory (--output) is required. Use --help for usage information.")
}

# Validate input file exists
if (!file.exists(args$input)) {
    stop(sprintf("ERROR: Input file not found: %s", args$input))
}

#' Parse comma-separated numeric values
#'
#' @param x Character string with comma-separated values
#' @param param_name Name of parameter (for error messages)
#' @return Numeric vector
parse_numeric_list <- function(x, param_name) {
    tryCatch({
        values <- as.numeric(strsplit(x, ",")[[1]])
        if (any(is.na(values))) {
            stop(sprintf("Invalid numeric values in %s: %s", param_name, x))
        }
        return(values)
    }, error = function(e) {
        stop(sprintf("ERROR parsing %s: %s", param_name, e$message))
    })
}

# Parse lambda and k_geom values
lambda_values <- parse_numeric_list(args$lambda, "lambda")
k_geom_values <- as.integer(parse_numeric_list(args$k_geom, "k_geom"))

# Validate parameter ranges
if (any(lambda_values < 0) || any(lambda_values > 1)) {
    stop("ERROR: Lambda values must be between 0 and 1")
}
if (any(k_geom_values < 1)) {
    stop("ERROR: k-geom values must be positive integers")
}
if (args$resolution <= 0) {
    stop("ERROR: Resolution must be positive")
}
if (args$n_pcs < 2) {
    stop("ERROR: Number of PCs must be at least 2")
}

# Derive sample name from input if not provided
if (is.null(args$sample_name)) {
    args$sample_name <- tools::file_path_sans_ext(basename(args$input))
}

# Set verbosity
verbose <- args$verbose && !args$quiet

# =============================================================================
# 3. LOGGING UTILITIES
# =============================================================================

#' Log a message with timestamp
#'
#' @param msg Message to log
#' @param level Log level (INFO, WARNING, ERROR)
log_message <- function(msg, level = "INFO") {
    if (verbose || level == "ERROR") {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
    }
}

#' Log a section header
#'
#' @param title Section title
log_section <- function(title) {
    if (verbose) {
        cat("\n")
        cat(paste(rep("=", 60), collapse = ""), "\n")
        cat(title, "\n")
        cat(paste(rep("=", 60), collapse = ""), "\n")
    }
}

#' Log a subsection
#'
#' @param title Subsection title
log_subsection <- function(title) {
    if (verbose) {
        cat("\n")
        cat(paste(rep("-", 40), collapse = ""), "\n")
        cat(title, "\n")
        cat(paste(rep("-", 40), collapse = ""), "\n")
    }
}

# =============================================================================
# 4. MAIN ANALYSIS
# =============================================================================

# Record start time for duration tracking
start_time <- Sys.time()

# Print banner
log_section("BANKSY SPATIAL TRANSCRIPTOMICS ANALYSIS")
log_message(sprintf("Analysis started: %s", format(start_time, "%Y-%m-%d %H:%M:%S")))

# -----------------------------------------------------------------------------
# 4.1 Load Required Packages
# -----------------------------------------------------------------------------

log_subsection("Loading Required Packages")

# List of required packages with their purposes
required_packages <- c(
    "Banksy",              # Core BANKSY algorithm
    "SpatialExperiment",   # Spatial data structure
    "SummarizedExperiment",# Base experiment class
    "SingleCellExperiment",# Single-cell methods
    "scuttle",             # QC and normalization
    "scater",              # Visualization utilities
    "rhdf5",               # HDF5/h5ad file reading
    "Matrix",              # Sparse matrix support
    "ggplot2",             # Visualization
    "data.table"           # Efficient data handling
)

# Check all packages are available before proceeding
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
    stop(sprintf(
        "ERROR: Missing required packages: %s\nInstall with BiocManager::install()",
        paste(missing_packages, collapse = ", ")
    ))
}

# Load packages (suppress startup messages for cleaner output)
suppressPackageStartupMessages({
    library(Banksy)
    library(SpatialExperiment)
    library(SummarizedExperiment)
    library(SingleCellExperiment)
    library(scuttle)
    library(scater)
    library(rhdf5)
    library(Matrix)
    library(ggplot2)
    library(data.table)
})

log_message("All packages loaded successfully")

# -----------------------------------------------------------------------------
# 4.2 Configuration Summary
# -----------------------------------------------------------------------------

log_subsection("Analysis Configuration")

# Build configuration object
config <- list(
    input_file   = normalizePath(args$input),
    output_dir   = args$output,
    sample_name  = args$sample_name,
    lambda       = lambda_values,
    k_geom       = k_geom_values,
    resolution   = args$resolution,
    n_pcs        = args$n_pcs,
    qc_lower     = args$qc_lower,
    qc_upper     = args$qc_upper,
    min_gene_pct = args$min_gene_cells,
    point_size   = args$point_size,
    plot_dpi     = args$plot_dpi,
    seed         = args$seed,
    skip_plots   = args$skip_plots,
    skip_umap    = args$skip_umap
)

# Print configuration
log_message(sprintf("Input file:      %s", config$input_file))
log_message(sprintf("Output directory: %s", config$output_dir))
log_message(sprintf("Sample name:     %s", config$sample_name))
log_message(sprintf("Lambda values:   %s", paste(config$lambda, collapse = ", ")))
log_message(sprintf("k_geom values:   %s", paste(config$k_geom, collapse = ", ")))
log_message(sprintf("Resolution:      %.2f", config$resolution))
log_message(sprintf("Number of PCs:   %d", config$n_pcs))
log_message(sprintf("Random seed:     %d", config$seed))

# Create output directory
dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(config$output_dir)) {
    stop(sprintf("ERROR: Failed to create output directory: %s", config$output_dir))
}
log_message(sprintf("Output directory created: %s", config$output_dir))

# Set random seed for reproducibility
set.seed(config$seed)

# -----------------------------------------------------------------------------
# 4.3 Data Loading
# -----------------------------------------------------------------------------

log_subsection("Loading H5AD Data")

# Read h5ad file using rhdf5 (faster than anndata for large files)
tryCatch({
    # Inspect h5ad file structure
    h5_contents <- h5ls(config$input_file)
    log_message("H5AD file structure detected")

    # Determine if expression matrix is sparse or dense
    is_sparse <- any(h5_contents$name == "data" & h5_contents$group == "/X")

    if (is_sparse) {
        # Sparse matrix format (CSR/CSC) - common for large datasets
        log_message("Reading sparse expression matrix...")

        data_vec <- h5read(config$input_file, "/X/data")
        indices  <- h5read(config$input_file, "/X/indices")
        indptr   <- h5read(config$input_file, "/X/indptr")

        # Read observation (cell/spot) identifiers
        obs_names <- h5read(config$input_file, "/obs/_index")

        # Read variable (gene) identifiers - try multiple possible locations
        var_names <- tryCatch({
            h5read(config$input_file, "/var/_index")
        }, error = function(e) {
            tryCatch({
                h5read(config$input_file, "/var/Symbol")
            }, error = function(e2) {
                h5read(config$input_file, "/var/gene_ids")
            })
        })

        n_obs <- length(obs_names)
        n_var <- length(var_names)
        log_message(sprintf("Dimensions: %d cells x %d genes", n_obs, n_var))

        # Construct sparse matrix from CSR components
        # Note: Python uses 0-indexing, R uses 1-indexing
        expr_matrix <- sparseMatrix(
            i = as.integer(indices) + 1L,
            p = as.integer(indptr),
            x = as.numeric(data_vec),
            dims = c(n_var, n_obs)
        )
    } else {
        # Dense matrix format
        log_message("Reading dense expression matrix...")
        expr_matrix <- t(h5read(config$input_file, "/X"))

        obs_names <- h5read(config$input_file, "/obs/_index")
        var_names <- tryCatch({
            h5read(config$input_file, "/var/_index")
        }, error = function(e) {
            h5read(config$input_file, "/var/Symbol")
        })
    }

    # Assign row/column names
    rownames(expr_matrix) <- var_names
    colnames(expr_matrix) <- obs_names

    log_message(sprintf("Expression matrix: %d genes x %d cells/spots",
                        nrow(expr_matrix), ncol(expr_matrix)))

    # Extract spatial coordinates
    # Try multiple possible storage locations
    spatial_coords <- NULL

    if (any(h5_contents$group == "/obsm" & h5_contents$name == "spatial")) {
        spatial_coords <- h5read(config$input_file, "/obsm/spatial")
        log_message("Spatial coordinates found in obsm/spatial")
    } else if (any(h5_contents$group == "/obsm" & h5_contents$name == "X_spatial")) {
        spatial_coords <- h5read(config$input_file, "/obsm/X_spatial")
        log_message("Spatial coordinates found in obsm/X_spatial")
    } else {
        # Try to find coordinates in obs columns
        obs_cols <- h5_contents[h5_contents$group == "/obs", "name"]
        x_candidates <- c("x", "X", "x_coord", "spatial_x", "array_col", "imagecol", "x_centroid")
        y_candidates <- c("y", "Y", "y_coord", "spatial_y", "array_row", "imagerow", "y_centroid")

        x_col <- intersect(x_candidates, obs_cols)[1]
        y_col <- intersect(y_candidates, obs_cols)[1]

        if (!is.na(x_col) && !is.na(y_col)) {
            x_vals <- h5read(config$input_file, paste0("/obs/", x_col))
            y_vals <- h5read(config$input_file, paste0("/obs/", y_col))
            spatial_coords <- cbind(x_vals, y_vals)
            log_message(sprintf("Spatial coordinates found in obs: %s, %s", x_col, y_col))
        }
    }

    if (is.null(spatial_coords)) {
        h5closeAll()
        stop("ERROR: No spatial coordinates found in h5ad file")
    }

    # Ensure correct orientation (cells x coordinates)
    if (is.matrix(spatial_coords)) {
        if (nrow(spatial_coords) < ncol(spatial_coords)) {
            spatial_coords <- t(spatial_coords)
        }
        if (ncol(spatial_coords) > 2) {
            spatial_coords <- spatial_coords[, 1:2]
        }
    }

    # Verify dimensions match
    if (nrow(spatial_coords) != ncol(expr_matrix)) {
        if (ncol(spatial_coords) == ncol(expr_matrix)) {
            spatial_coords <- t(spatial_coords)
        } else {
            stop("ERROR: Spatial coordinate dimensions don't match expression matrix")
        }
    }

    log_message(sprintf("Spatial coordinates: %d points", nrow(spatial_coords)))

    # Create SpatialExperiment object
    log_message("Creating SpatialExperiment object...")
    spe <- SpatialExperiment(
        assays = list(counts = expr_matrix),
        spatialCoords = spatial_coords
    )

    # Close HDF5 connections
    h5closeAll()

    log_message("Data loading complete")

}, error = function(e) {
    h5closeAll()
    stop(sprintf("ERROR loading data: %s", e$message))
})

# -----------------------------------------------------------------------------
# 4.4 Quality Control
# -----------------------------------------------------------------------------

log_subsection("Quality Control")

# Calculate per-cell QC metrics
# These metrics help identify low-quality cells/spots
qc_metrics <- perCellQCMetrics(spe)
colData(spe) <- cbind(colData(spe), qc_metrics)

log_message(sprintf("Total UMI range: %d - %d",
                    min(qc_metrics$sum), max(qc_metrics$sum)))
log_message(sprintf("Genes detected range: %d - %d",
                    min(qc_metrics$detected), max(qc_metrics$detected)))

# Apply adaptive QC thresholds based on data distribution
# This approach is robust to dataset-specific characteristics
sum_threshold      <- quantile(qc_metrics$sum, c(config$qc_lower, config$qc_upper))
detected_threshold <- quantile(qc_metrics$detected, c(config$qc_lower, config$qc_upper))

# Identify cells passing QC
keep_cells <- (qc_metrics$sum >= sum_threshold[1]) &
              (qc_metrics$sum <= sum_threshold[2]) &
              (qc_metrics$detected >= detected_threshold[1]) &
              (qc_metrics$detected <= detected_threshold[2])

n_removed <- sum(!keep_cells)
pct_removed <- round(100 * n_removed / ncol(spe), 1)
log_message(sprintf("Cells removed by QC: %d (%.1f%%)", n_removed, pct_removed))

# Apply cell filter
spe <- spe[, keep_cells]
log_message(sprintf("Cells after QC: %d", ncol(spe)))

# Filter genes based on detection across cells
# Genes detected in too few cells provide unreliable signal
min_cells <- max(10, floor(ncol(spe) * config$min_gene_pct))
gene_counts <- rowSums(assay(spe, "counts") > 0)
keep_genes <- gene_counts >= min_cells

n_genes_removed <- sum(!keep_genes)
log_message(sprintf("Genes removed (detected in <%d cells): %d", min_cells, n_genes_removed))

spe <- spe[keep_genes, ]
log_message(sprintf("Genes after filtering: %d", nrow(spe)))

# -----------------------------------------------------------------------------
# 4.5 Normalization
# -----------------------------------------------------------------------------

log_subsection("Normalization")

# Library size normalization
# Accounts for differences in sequencing depth between cells
spe <- computeLibraryFactors(spe)

# Generate normalized counts (linear scale)
assay(spe, "normcounts") <- normalizeCounts(spe, log = FALSE)

# Generate log-transformed normalized counts
# Log transformation stabilizes variance for downstream analysis
assay(spe, "lognormcounts") <- log1p(assay(spe, "normcounts"))

log_message("Normalization complete (library size + log1p)")

# -----------------------------------------------------------------------------
# 4.6 BANKSY Analysis
# -----------------------------------------------------------------------------

log_subsection("BANKSY Analysis")

# Step 1: Compute spatial neighborhood features
# BANKSY incorporates information from neighboring cells to capture spatial context
log_message("Computing neighborhood features...")
log_message(sprintf("  k_geom values: %s", paste(config$k_geom, collapse = ", ")))

spe <- computeBanksy(
    spe,
    assay_name = "normcounts",
    compute_agf = TRUE,          # Compute azimuthal Gabor filter features
    k_geom = config$k_geom       # Number of neighbors for mean and Gabor
)

# Step 2: Dimensionality reduction with BANKSY
# PCA on combined own-cell and neighborhood features
log_message("Running BANKSY PCA...")
log_message(sprintf("  Lambda values: %s", paste(config$lambda, collapse = ", ")))
log_message(sprintf("  Number of PCs: %d", config$n_pcs))

spe <- runBanksyPCA(
    spe,
    use_agf = TRUE,
    lambda = config$lambda,
    npcs = config$n_pcs
)

# Step 3: UMAP embedding (optional but useful for visualization)
if (!config$skip_umap) {
    log_message("Computing UMAP embeddings...")
    spe <- runBanksyUMAP(
        spe,
        use_agf = TRUE,
        lambda = config$lambda
    )
}

# Step 4: Leiden clustering
# Identifies cell populations with similar expression + spatial patterns
log_message("Clustering with Leiden algorithm...")
log_message(sprintf("  Resolution: %.2f", config$resolution))

spe <- clusterBanksy(
    spe,
    use_agf = TRUE,
    lambda = config$lambda,
    resolution = config$resolution
)

# Step 5: Connect clusters across lambda values (if multiple)
# This helps track how clusters relate at different spatial weights
if (length(config$lambda) > 1) {
    spe <- connectClusters(spe)
    log_message("Clusters connected across lambda values")
}

log_message("BANKSY analysis complete")

# -----------------------------------------------------------------------------
# 4.7 Results Extraction
# -----------------------------------------------------------------------------

log_subsection("Extracting Results")

# Identify cluster columns in the results
cluster_cols <- grep("^clust_", colnames(colData(spe)), value = TRUE)
log_message(sprintf("Cluster columns found: %s", paste(cluster_cols, collapse = ", ")))

# Summarize cluster counts for each lambda
for (cc in cluster_cols) {
    n_clusters <- length(unique(colData(spe)[[cc]]))
    log_message(sprintf("  %s: %d clusters", cc, n_clusters))
}

# Create summary data frame for reporting
results_summary <- data.frame(
    sample          = config$sample_name,
    n_cells         = ncol(spe),
    n_genes         = nrow(spe),
    lambda_values   = paste(config$lambda, collapse = ","),
    k_geom_values   = paste(config$k_geom, collapse = ","),
    resolution      = config$resolution,
    n_pcs           = config$n_pcs,
    qc_lower        = config$qc_lower,
    qc_upper        = config$qc_upper,
    seed            = config$seed,
    analysis_date   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
)

# Add cluster counts to summary
for (cc in cluster_cols) {
    results_summary[[paste0("n_", cc)]] <- length(unique(colData(spe)[[cc]]))
}

# -----------------------------------------------------------------------------
# 4.8 Save Outputs
# -----------------------------------------------------------------------------

log_subsection("Saving Outputs")

# Save processed SpatialExperiment object (RDS format)
# This preserves all analysis results for downstream use in R
rds_file <- file.path(config$output_dir, paste0(config$sample_name, "_banksy_results.rds"))
saveRDS(spe, rds_file)
log_message(sprintf("RDS object saved: %s", rds_file))

# Save cluster assignments as CSV
# Useful for integration with other tools or clinical reports
cluster_data <- as.data.frame(colData(spe)[, c(cluster_cols, "sum", "detected")])
cluster_data$cell_id <- colnames(spe)
cluster_data$x <- spatialCoords(spe)[, 1]
cluster_data$y <- spatialCoords(spe)[, 2]

csv_file <- file.path(config$output_dir, paste0(config$sample_name, "_clusters.csv"))
write.csv(cluster_data, csv_file, row.names = FALSE)
log_message(sprintf("Cluster assignments saved: %s", csv_file))

# Save analysis summary
summary_file <- file.path(config$output_dir, paste0(config$sample_name, "_summary.csv"))
write.csv(results_summary, summary_file, row.names = FALSE)
log_message(sprintf("Summary saved: %s", summary_file))

# -----------------------------------------------------------------------------
# 4.9 Generate Visualizations
# -----------------------------------------------------------------------------

if (!config$skip_plots) {
    log_subsection("Generating Visualizations")

    # Spatial cluster plots for each lambda value
    for (cc in cluster_cols) {
        plot_data <- data.frame(
            x = spatialCoords(spe)[, 1],
            y = spatialCoords(spe)[, 2],
            cluster = factor(colData(spe)[[cc]])
        )

        # Extract lambda value from column name for title
        lambda_str <- gsub("clust_M1_lambda", "λ=", cc)
        lambda_str <- gsub("_res.*", "", lambda_str)

        p <- ggplot(plot_data, aes(x = x, y = y, color = cluster)) +
            geom_point(size = config$point_size, alpha = 0.7) +
            theme_minimal() +
            theme(
                legend.position = "right",
                plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                plot.subtitle = element_text(hjust = 0.5, size = 10),
                panel.grid = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank()
            ) +
            labs(
                title = sprintf("BANKSY Spatial Clustering (%s)", lambda_str),
                subtitle = sprintf("Sample: %s | Resolution: %.2f | n=%d",
                                   config$sample_name, config$resolution, nrow(plot_data)),
                x = NULL,
                y = NULL,
                color = "Cluster"
            ) +
            coord_fixed() +
            scale_color_discrete()

        plot_file <- file.path(config$output_dir, paste0(config$sample_name, "_", cc, "_spatial.png"))
        ggsave(plot_file, p, width = 10, height = 8, dpi = config$plot_dpi)
        log_message(sprintf("Spatial plot saved: %s", basename(plot_file)))
    }

    # UMAP plots (if computed)
    if (!config$skip_umap) {
        umap_names <- reducedDimNames(spe)[grep("^UMAP", reducedDimNames(spe))]

        for (umap_name in umap_names) {
            umap_coords <- reducedDim(spe, umap_name)

            # Match UMAP to corresponding cluster column
            # UMAP names: UMAP_M1_lam0.2_k50
            # Cluster names: clust_M1_lam0.2_k50_res0.8
            # Find cluster column that starts with the UMAP pattern
            umap_pattern <- gsub("UMAP_", "clust_", umap_name)
            cc <- cluster_cols[grep(paste0("^", umap_pattern), cluster_cols)]

            if (length(cc) > 0) {
                cc <- cc[1]  # Take first match if multiple
                plot_data <- data.frame(
                    UMAP1 = umap_coords[, 1],
                    UMAP2 = umap_coords[, 2],
                    cluster = factor(colData(spe)[[cc]])
                )

                # Extract lambda value for title (e.g., "lam0.2" -> "λ=0.2")
                lambda_match <- regmatches(umap_name, regexpr("lam[0-9.]+", umap_name))
                lambda_label <- if (length(lambda_match) > 0) {
                    gsub("lam", "λ=", lambda_match)
                } else {
                    umap_name
                }

                p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = cluster)) +
                    geom_point(size = config$point_size, alpha = 0.7) +
                    theme_minimal() +
                    theme(
                        legend.position = "right",
                        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                        plot.subtitle = element_text(hjust = 0.5, size = 10)
                    ) +
                    labs(
                        title = sprintf("BANKSY UMAP (%s)", lambda_label),
                        subtitle = sprintf("Sample: %s", config$sample_name),
                        color = "Cluster"
                    ) +
                    scale_color_discrete()

                plot_file <- file.path(config$output_dir, paste0(config$sample_name, "_", umap_name, ".png"))
                ggsave(plot_file, p, width = 8, height = 6, dpi = config$plot_dpi)
                log_message(sprintf("UMAP plot saved: %s", basename(plot_file)))
            }
        }
    }
}

# -----------------------------------------------------------------------------
# 4.10 Completion
# -----------------------------------------------------------------------------

# Calculate analysis duration
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

# Write completion marker file
completion_file <- file.path(config$output_dir, ".banksy_complete")
writeLines(c(
    sprintf("Completed: %s", format(end_time, "%Y-%m-%d %H:%M:%S")),
    sprintf("Duration: %.2f minutes", as.numeric(duration)),
    sprintf("Sample: %s", config$sample_name)
), completion_file)

# Final summary
log_section("ANALYSIS COMPLETE")
log_message(sprintf("Duration: %.2f minutes", as.numeric(duration)))
log_message(sprintf("Output directory: %s", config$output_dir))
log_message("Output files:")
log_message(sprintf("  - %s_banksy_results.rds (R object)", config$sample_name))
log_message(sprintf("  - %s_clusters.csv (cluster assignments)", config$sample_name))
log_message(sprintf("  - %s_summary.csv (analysis summary)", config$sample_name))
if (!config$skip_plots) {
    log_message(sprintf("  - %s_*_spatial.png (spatial plots)", config$sample_name))
    if (!config$skip_umap) {
        log_message(sprintf("  - %s_*_UMAP.png (UMAP plots)", config$sample_name))
    }
}

cat("\nSUCCESS: BANKSY analysis completed successfully\n")
