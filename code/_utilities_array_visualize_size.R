# calculate size of array

calculate_array_size <- function(rows, columns, layers_sequence) {
  # Data type sizes in bytes
  logical_size <- 1
  integer_size <- 4
  double_size <- 8
  
  # Initialize the results table
  results <- data.frame(
    Layers = layers_sequence,
    Logical_GB = numeric(length(layers_sequence)),
    Integer_GB = numeric(length(layers_sequence)),
    Double_GB = numeric(length(layers_sequence))
  )
  
  # Calculate the array size and convert it to GB for each data type
  for (i in seq_along(layers_sequence)) {
    layers <- layers_sequence[i]
    array_size <- rows * columns * layers
    
    results$Logical_GB[i] <- array_size * logical_size / (1024^3)
    results$Integer_GB[i] <- array_size * integer_size / (1024^3)
    results$Double_GB[i] <- array_size * double_size / (1024^3)
  }
  
  return(results)
}

#----------------------------
rows_sequence <- 10824
columns_sequence <- seq(20, 586, 20)
layers_sequence <- seq(50, 4000, 100)

all_results <- lapply(rows_sequence, function(rows) {
  lapply(columns_sequence, function(cols) {
    calculate_array_size(rows, cols, layers_sequence)
  })
})

# Combine results into a single data frame
combined_table <- do.call(rbind, lapply(all_results, function(row_results) {
  do.call(rbind, row_results)
}))

# Add row and column information to the combined table
combined_table$Rows <- rep(rep(rows_sequence, each = length(columns_sequence) * length(layers_sequence)), 1)
combined_table$Columns <- rep(rep(columns_sequence, each = length(layers_sequence)), length(rows_sequence))

# Reorder columns for clarity
combined_table <- combined_table[, c("Rows", "Columns", "Layers", "Logical_GB", "Integer_GB", "Double_GB")]

#-------------------------------
library(ggplot2)
library(viridis)

combined_table |> 
  ggplot(aes(x = Layers, y = Columns, fill = Logical_GB )) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE)

