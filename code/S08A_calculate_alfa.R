# Calculate alfa
# MTP version

# Load presence/abscence Model
fitted_full_model_path_PA <- file.path("models/fbs_M016PA_thin_150_samples_1000_chains_4/fitted_fbs_M016PA_thin_150_samples_1000_chains_4.rds")
hM_PA <- readRDS(fitted_full_model_path_PA)
spatial_level_name_PA <- hM_PA$rLNames[1] 

# Calculate MTP per species

predY_full <- predict(hM_PA, expected = TRUE)
predY_full <- simplify2array(predY_full)
Y_obs <- (hM_PA$Y > 0) * 1L # Ensure binary indicator (1 = presence)

n_sp      <- dim(predY_full)[2]
n_samples <- dim(predY_full)[3]

alfa_matrix <- matrix(NA_real_, nrow = n_sp, ncol = n_samples,
                     dimnames = list(colnames(hM_PA$Y), paste0("sample_", seq_len(n_samples))))

for (j in seq_len(n_sp)) {
  # Indices where species j was actually observed in the observed data
  pres_idx <- which(Y_obs[, j] == 1L)
  
  if (length(pres_idx) > 0) {
    # Extract predicted probabilities at presence sites across all posterior samples
    # Resulting sub-matrix: [presence_sites, n_samples]
    probs_at_presences <- predY_full[pres_idx, j, , drop = FALSE]
    
    alfa_matrix[j, ] <- apply(probs_at_presences, 3, min)
  }
}

dir.create("results", showWarnings = F)
saveRDS(alfa_matrix, file.path("results", "alfa_matrix.rds"))
