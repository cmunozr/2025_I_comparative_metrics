library(Hmsc)
library(parallel)
print("library loaded")

# --- Get Command-Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) { # Changed to expect 5 arguments
  stop("Usage: Rscript S2_fit_models.R <output_dir> <unfitted_models_file> <nChains> <nParallel>")# <model_index>")
}

output_dir <- args[1]
unfitted_models_file <- args[2]
nChains <- as.numeric(args[3])
nParallel <- as.numeric(args[4])
model_index <- as.integer(1) #args[5]) # Get the model index from the arguments

# --- Create Output Directory (if needed) ---
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load Unfitted Models ---
load(file = unfitted_models_file)

nm <- length(models)

# --- Get MCMC Parameters ---
samples <- as.numeric(1000)
thin <- as.numeric(1000)

print(paste0("samples = ", samples, ", thining = ", thin))

print("--- Script started ---")

# --- Main ---

if (model_index >= 1 && model_index <= nm) { #added this if
  model_name <- names(models)[model_index]
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

    m_unfitted <- models[[model_index]] # Use the model_index

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
} else {
  print(paste("Error: model_index out of bounds (", model_index, ").  Must be between 1 and ", nm, sep=""))
} #added close bracket

print("--- Script finished ---")