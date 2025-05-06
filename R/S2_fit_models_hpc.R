library(Hmsc)
library(parallel)
print("library loaded")

# --- Get Command-Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript S2_fit_models.R <output_dir> <unfitted_models_file> <nChains> <nParallel>")
}

output_dir <- args[1]
unfitted_models_file <- args[2]
nChains <- as.numeric(args[3])
nParallel <- as.numeric(args[4])
#samples_index <- as.integer(args[5])
#thin_index <- as.integer(args[6])

# --- Create Output Directory (if needed) ---
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load Unfitted Models ---
load(file = unfitted_models_file)

nm <- length(models)

# --- Get MCMC Parameters from Lists ---
# samples_list <- c(5, 250, 1000)
# thin_list <- c(1, 1, 1000)

samples <- as.numeric(1000)
thin <- as.numeric(1000)

print(paste0("samples = ", samples, ", thining = ", thin))

print("--- Script started ---")

# --- Main Loop for Fitting Models ---

for (mi in 1:nm) {
  model_name <- names(models)[mi]
  filename <- file.path(output_dir, paste0(
    model_name,
    "_thin_", as.character(thin),
    "_samples_", as.character(samples),
    "_chains_", as.character(nChains),
    ".Rdata"
  ))

  if (file.exists(filename)) {
    print(paste("Skipping (exists):", basename(filename)))
  } else {
    print(paste("Processing model:", model_name, " | Saving to:", basename(filename)))
    start_time <- Sys.time()
    print(start_time)

    m_unfitted <- models[[mi]]

    m_fitted <- tryCatch({
      sampleMcmc(
        m_unfitted,
        samples = samples,
        thin = thin,
        adaptNf = rep(ceiling(0.4 * samples * thin), m_unfitted$nr),
        transient = ceiling(0.5 * samples * thin),
        nChains = nChains,
        nParallel = nParallel,
        verbose = ceiling(samples / 10)
      )
    }, error = function(e) {
      print(paste("ERROR sampling model:", model_name))
      print(e$message)
      NULL
    })

    if (!is.null(m_fitted)) {
      save(m_fitted, file = filename)
      end_time <- Sys.time()
      print(paste("Finished model:", model_name, "| Time:", format(end_time - start_time)))
    } else {
      print(paste("Failed model:", model_name, "- File not saved."))
    }
  }
}

print("--- Script finished ---")