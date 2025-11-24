library(ggplot2)
library(dplyr)
library(tidyr)
library(outliers)
library(pheatmap)
library(viridis)
library(GenomicRanges)
library(gridExtra)
library(grid)
library(scales)

setwd("/Users/atrayeeray/Desktop/Sid/CUT&Tag_cellno. standardization/Direct CUT&Tag/new_analysis_07082025/saturation_analysis/")

##### CTCF

##Downsample IgG (merged) reads and no downsample M100K_CTCF (replicate wise) reads 
M100K_CTCF_no_vs_M100K_IgG_down <- read.csv("M100K_CTCF_downsample_based_constant_merged_IgG/downsample_summary.csv")


# Define colors
colors <- c("MACS2_Default" = "#1f77b4",  # Blue
            "MACS2_Large" = "#ff7f0e")     # Orange

# Function to perform Grubbs test and remove outliers
remove_outliers_grubbs <- function(x) {
  if (length(x) < 3) {
    return(list(values = x, outliers = NULL))
  }
  
  outliers <- NULL
  values <- x
  
  # Perform Grubbs test iteratively (max 1 outlier for 3 replicates)
  if (length(values) == 3) {
    test_result <- grubbs.test(values, type = 10, opposite = FALSE, two.sided = FALSE)
    if (!is.null(test_result) && !is.na(test_result$p.value) && test_result$p.value < 0.05) {
      # Identify the outlier
      mean_val <- mean(values)
      outlier_idx <- which.max(abs(values - mean_val))
      outliers <- values[outlier_idx]
      values <- values[-outlier_idx]
    }
  }
  
  return(list(values = values, outliers = outliers))
}

# Prepare data
df <- M100K_CTCF_no_vs_M100K_IgG_down

# Convert target_reads_for_CTCF to numeric for ordering
# Handle "original" case - we want Original first (leftmost), then descending
# So order should be: Original, 5.00M, 2.50M, 1000K, 0.50M, 100K, 0.05M, 10K
# Assign Original a value larger than 5M (e.g., 6M) so it sorts first when descending
df$target_reads_numeric <- ifelse(df$target_reads_for_CTCF == "original", 
                                  6000000,  # Use 6M so Original sorts first (leftmost)
                                  as.numeric(df$target_reads_for_CTCF))

# Sort by target_reads_numeric descending (Original first, then descending)
df <- df[order(-df$target_reads_numeric), ]

# Group by target_reads and process outliers
df_summary <- df %>%
  group_by(target_reads_for_CTCF, target_reads_numeric) %>%
  summarise(
    # Process num_peaks_default with outlier removal
    peaks_default_raw = list(num_peaks_default),
    peaks_default_outliers = list(remove_outliers_grubbs(num_peaks_default)$outliers),
    peaks_default_clean = list(remove_outliers_grubbs(num_peaks_default)$values),
    mean_peaks_default = mean(remove_outliers_grubbs(num_peaks_default)$values),
    sd_peaks_default = sd(remove_outliers_grubbs(num_peaks_default)$values),
    n_default = length(remove_outliers_grubbs(num_peaks_default)$values),
    
    # Process num_peaks_large with outlier removal
    peaks_large_raw = list(num_peaks_large),
    peaks_large_outliers = list(remove_outliers_grubbs(num_peaks_large)$outliers),
    peaks_large_clean = list(remove_outliers_grubbs(num_peaks_large)$values),
    mean_peaks_large = mean(remove_outliers_grubbs(num_peaks_large)$values),
    sd_peaks_large = sd(remove_outliers_grubbs(num_peaks_large)$values),
    n_large = length(remove_outliers_grubbs(num_peaks_large)$values),
    
    # Use mean of actual_reads_for_CTCF for x-axis position
    mean_actual_reads = mean(actual_reads_for_CTCF),
    .groups = 'drop'
  ) %>%
  ungroup()

# Print outlier information
cat("=== Outlier Detection Summary ===\n")
for (i in 1:nrow(df_summary)) {
  target <- df_summary$target_reads_for_CTCF[i]
  default_outliers <- unlist(df_summary$peaks_default_outliers[i])
  large_outliers <- unlist(df_summary$peaks_large_outliers[i])
  
  if (!is.null(default_outliers) && length(default_outliers) > 0) {
    cat(sprintf("Target %s - Default peaks outliers: %s\n", target, paste(default_outliers, collapse=", ")))
  }
  if (!is.null(large_outliers) && length(large_outliers) > 0) {
    cat(sprintf("Target %s - Large peaks outliers: %s\n", target, paste(large_outliers, collapse=", ")))
  }
}
cat("\n")

# Prepare data for plotting (long format)
df_plot <- bind_rows(
  df_summary %>%
    select(target_reads_for_CTCF, target_reads_numeric, mean_actual_reads, 
           mean_peaks = mean_peaks_default, sd_peaks = sd_peaks_default, n = n_default) %>%
    mutate(method = "MACS2_Default"),
  df_summary %>%
    select(target_reads_for_CTCF, target_reads_numeric, mean_actual_reads,
           mean_peaks = mean_peaks_large, sd_peaks = sd_peaks_large, n = n_large) %>%
    mutate(method = "MACS2_Large")
) %>%
  # Arrange: Original first (leftmost), then descending
  arrange(-target_reads_numeric)

# Create x-axis labels (use target_reads but position at mean_actual_reads)
# Convert to numeric first, then format: >= 1M as "M", < 1M as "K"
df_plot <- df_plot %>%
  mutate(
    target_numeric = ifelse(target_reads_for_CTCF == "original", 
                            NA, 
                            as.numeric(target_reads_for_CTCF)),
    x_label = ifelse(target_reads_for_CTCF == "original",
                     "Original",
                     ifelse(target_numeric >= 1e6,
                            paste0(format(target_numeric / 1e6, digits=2), "M"),
                            paste0(format(target_numeric / 1e3, digits=1), "K")))
  )

# Create a mapping for x-axis breaks and labels with fixed intervals
# Order: Original first (leftmost), then descending
x_labels_map <- df_plot %>%
  select(target_reads_for_CTCF, target_reads_numeric, mean_actual_reads, x_label) %>%
  distinct() %>%
  arrange(-target_reads_numeric) %>%  # Original (6M) first, then descending
  mutate(x_position = row_number())  # Create fixed interval positions (1, 2, 3, ...)

# Add x_position to df_plot based on target_reads_for_CTCF
df_plot <- df_plot %>%
  left_join(x_labels_map %>% select(target_reads_for_CTCF, x_position), 
            by = "target_reads_for_CTCF")

# Calculate error bars for log scale (additive in log space)
df_plot <- df_plot %>%
  mutate(
    mean_peaks_log = log10(mean_peaks + 1),  # +1 to avoid log(0)
    # For log scale, error bars need special handling
    # Use multiplicative factor for symmetric error bars in log space
    lower_error = pmax(1, mean_peaks - sd_peaks),  # Ensure positive
    upper_error = mean_peaks + sd_peaks,
    lower_error_log = log10(lower_error + 1),
    upper_error_log = log10(upper_error + 1)
  )

# Create the plot
# Note: Using 'size' for compatibility; use 'linewidth' in ggplot2 3.4.0+
# Use x_position for fixed interval spacing on x-axis
p <- ggplot(df_plot, aes(x = x_position, y = mean_peaks, color = method, group = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3.5) +
  # Error bars in linear space (will be transformed by scale_y_log10)
  geom_errorbar(aes(ymin = lower_error, ymax = upper_error), 
                width = 0.1,  # Fixed width for error bars
                size = 0.8, alpha = 0.7) +
  scale_color_manual(values = colors, 
                     labels = c("MACS2_Default" = "MACS2 Default", 
                                "MACS2_Large" = "MACS2 Large"),
                     name = "Method") +
  scale_x_continuous(
    breaks = x_labels_map$x_position,
    labels = x_labels_map$x_label
  ) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    breaks = scales::trans_breaks("log10", function(x) 10^x)
  ) +
  labs(
    x = "Target Reads for CTCF",
    y = "Number of Peaks (log10 scale)",
    title = "CTCF Peak Calling Saturation Curve",
    subtitle = "X-axis: fixed intervals (target reads) | Y-axis: log10 scale"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels if needed
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor.y = element_line(color = "grey90", size = 0.3),
    panel.grid.major = element_line(color = "grey80", size = 0.5)
  )

# Print the plot
print(p)

# Save the plot
ggsave("CTCF_saturation_curve_based_on_constant_merged_IgG.pdf", plot = p, width = 10, height = 6, units = "in")
ggsave("CTCF_saturation_curve_based_on_constant_merged_IgG.png", plot = p, width = 10, height = 6, units = "in", dpi = 300)

cat("Plot saved as 'CTCF_saturation_curve.pdf' and 'CTCF_saturation_curve.png'\n")


########### K27me3

##Downsample IgG (merged) reads and no downsample M100K_K27me3 (replicate wise) reads 
M100K_K27me3_no_vs_M100K_IgG_down <- read.csv("M100K_K27me3_downsample_based_constant_merged_IgG/downsample_summary.csv")


# Define colors
colors <- c("MACS2_Default" = "#1f77b4",  # Blue
            "MACS2_Large" = "#ff7f0e")     # Orange

# Function to perform Grubbs test and remove outliers
remove_outliers_grubbs <- function(x) {
  if (length(x) < 3) {
    return(list(values = x, outliers = NULL))
  }
  
  outliers <- NULL
  values <- x
  
  # Perform Grubbs test iteratively (max 1 outlier for 3 replicates)
  if (length(values) == 3) {
    test_result <- grubbs.test(values, type = 10, opposite = FALSE, two.sided = FALSE)
    if (!is.null(test_result) && !is.na(test_result$p.value) && test_result$p.value < 0.05) {
      # Identify the outlier
      mean_val <- mean(values)
      outlier_idx <- which.max(abs(values - mean_val))
      outliers <- values[outlier_idx]
      values <- values[-outlier_idx]
    }
  }
  
  return(list(values = values, outliers = outliers))
}

# Prepare data
df <- M100K_K27me3_no_vs_M100K_IgG_down

# Convert target_reads_for_K27me3 to numeric for ordering
# Handle "original" case - we want Original first (leftmost), then descending
# So order should be: Original, 5.00M, 2.50M, 1000K, 0.50M, 100K, 0.05M, 10K
# Assign Original a value larger than 5M (e.g., 6M) so it sorts first when descending
df$target_reads_numeric <- ifelse(df$target_reads_for_K27me3 == "original", 
                                  16000000,  # Use 6M so Original sorts first (leftmost)
                                  as.numeric(df$target_reads_for_K27me3))

# Sort by target_reads_numeric descending (Original first, then descending)
df <- df[order(-df$target_reads_numeric), ]

# Group by target_reads and process outliers
df_summary <- df %>%
  group_by(target_reads_for_K27me3, target_reads_numeric) %>%
  summarise(
    # Process num_peaks_default with outlier removal
    peaks_default_raw = list(num_peaks_default),
    peaks_default_outliers = list(remove_outliers_grubbs(num_peaks_default)$outliers),
    peaks_default_clean = list(remove_outliers_grubbs(num_peaks_default)$values),
    mean_peaks_default = mean(remove_outliers_grubbs(num_peaks_default)$values),
    sd_peaks_default = sd(remove_outliers_grubbs(num_peaks_default)$values),
    n_default = length(remove_outliers_grubbs(num_peaks_default)$values),
    
    # Process num_peaks_large with outlier removal
    peaks_large_raw = list(num_peaks_large),
    peaks_large_outliers = list(remove_outliers_grubbs(num_peaks_large)$outliers),
    peaks_large_clean = list(remove_outliers_grubbs(num_peaks_large)$values),
    mean_peaks_large = mean(remove_outliers_grubbs(num_peaks_large)$values),
    sd_peaks_large = sd(remove_outliers_grubbs(num_peaks_large)$values),
    n_large = length(remove_outliers_grubbs(num_peaks_large)$values),
    
    # Use mean of actual_reads_for_K27me3 for x-axis position
    mean_actual_reads = mean(actual_reads_for_K27me3),
    .groups = 'drop'
  ) %>%
  ungroup()

# Print outlier information
cat("=== Outlier Detection Summary ===\n")
for (i in 1:nrow(df_summary)) {
  target <- df_summary$target_reads_for_K27me3[i]
  default_outliers <- unlist(df_summary$peaks_default_outliers[i])
  large_outliers <- unlist(df_summary$peaks_large_outliers[i])
  
  if (!is.null(default_outliers) && length(default_outliers) > 0) {
    cat(sprintf("Target %s - Default peaks outliers: %s\n", target, paste(default_outliers, collapse=", ")))
  }
  if (!is.null(large_outliers) && length(large_outliers) > 0) {
    cat(sprintf("Target %s - Large peaks outliers: %s\n", target, paste(large_outliers, collapse=", ")))
  }
}
cat("\n")

# Prepare data for plotting (long format)
df_plot <- bind_rows(
  df_summary %>%
    select(target_reads_for_K27me3, target_reads_numeric, mean_actual_reads, 
           mean_peaks = mean_peaks_default, sd_peaks = sd_peaks_default, n = n_default) %>%
    mutate(method = "MACS2_Default"),
  df_summary %>%
    select(target_reads_for_K27me3, target_reads_numeric, mean_actual_reads,
           mean_peaks = mean_peaks_large, sd_peaks = sd_peaks_large, n = n_large) %>%
    mutate(method = "MACS2_Large")
) %>%
  # Arrange: Original first (leftmost), then descending
  arrange(-target_reads_numeric)

# Create x-axis labels (use target_reads but position at mean_actual_reads)
# Convert to numeric first, then format: >= 1M as "M", < 1M as "K"
df_plot <- df_plot %>%
  mutate(
    target_numeric = ifelse(target_reads_for_K27me3 == "original", 
                            NA, 
                            as.numeric(target_reads_for_K27me3)),
    x_label = ifelse(target_reads_for_K27me3 == "original",
                     "Original",
                     ifelse(target_numeric >= 1e6,
                            paste0(format(target_numeric / 1e6, digits=2), "M"),
                            paste0(format(target_numeric / 1e3, digits=1), "K")))
  )

# Create a mapping for x-axis breaks and labels with fixed intervals
# Order: Original first (leftmost), then descending
x_labels_map <- df_plot %>%
  select(target_reads_for_K27me3, target_reads_numeric, mean_actual_reads, x_label) %>%
  distinct() %>%
  arrange(-target_reads_numeric) %>%  # Original (16M) first, then descending
  mutate(x_position = row_number())  # Create fixed interval positions (1, 2, 3, ...)

# Add x_position to df_plot based on target_reads_for_K27me3
df_plot <- df_plot %>%
  left_join(x_labels_map %>% select(target_reads_for_K27me3, x_position), 
            by = "target_reads_for_K27me3")

# Calculate error bars for log scale (additive in log space)
df_plot <- df_plot %>%
  mutate(
    mean_peaks_log = log10(mean_peaks + 1),  # +1 to avoid log(0)
    # For log scale, error bars need special handling
    # Use multiplicative factor for symmetric error bars in log space
    lower_error = pmax(1, mean_peaks - sd_peaks),  # Ensure positive
    upper_error = mean_peaks + sd_peaks,
    lower_error_log = log10(lower_error + 1),
    upper_error_log = log10(upper_error + 1)
  )

# Create the plot
# Note: Using 'size' for compatibility; use 'linewidth' in ggplot2 3.4.0+
# Use x_position for fixed interval spacing on x-axis
p <- ggplot(df_plot, aes(x = x_position, y = mean_peaks, color = method, group = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3.5) +
  # Error bars in linear space (will be transformed by scale_y_log10)
  geom_errorbar(aes(ymin = lower_error, ymax = upper_error), 
                width = 0.1,  # Fixed width for error bars
                size = 0.8, alpha = 0.7) +
  scale_color_manual(values = colors, 
                     labels = c("MACS2_Default" = "MACS2 Default", 
                                "MACS2_Large" = "MACS2 Large"),
                     name = "Method") +
  scale_x_continuous(
    breaks = x_labels_map$x_position,
    labels = x_labels_map$x_label
  ) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    breaks = scales::trans_breaks("log10", function(x) 10^x)
  ) +
  labs(
    x = "Target Reads for K27me3",
    y = "Number of Peaks (log10 scale)",
    title = "K27me3 Peak Calling Saturation Curve",
    subtitle = "X-axis: fixed intervals (target reads) | Y-axis: log10 scale"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels if needed
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor.y = element_line(color = "grey90", size = 0.3),
    panel.grid.major = element_line(color = "grey80", size = 0.5)
  )

# Print the plot
print(p)

# Save the plot
ggsave("K27me3_saturation_curve_based_on_constant_merged_IgG.pdf", plot = p, width = 10, height = 6, units = "in")
ggsave("K27me3_saturation_curve_based_on_constant_merged_IgG.png", plot = p, width = 10, height = 6, units = "in", dpi = 300)

cat("Plot saved as 'K27me3_saturation_curve.pdf' and 'K27me3_saturation_curve.png'\n")


########## compare CTCF and K27me3 MACS2 default###


# Define colors - more distinct shades of blue (only Default methods)
colors <- c("CTCF_MACS2_Default" = "#3182bd",      # Darker blue
            "K27me3_MACS2_Default" = "#9ecae1")    # Lighter blue (more distinct)

# Function to perform Grubbs test and remove outliers
remove_outliers_grubbs <- function(x) {
  if (length(x) < 3) {
    return(list(values = x, outliers = NULL))
  }
  
  outliers <- NULL
  values <- x
  
  # Perform Grubbs test iteratively (max 1 outlier for 3 replicates)
  if (length(values) == 3) {
    test_result <- grubbs.test(values, type = 10, opposite = FALSE, two.sided = FALSE)
    if (!is.null(test_result) && !is.na(test_result$p.value) && test_result$p.value < 0.05) {
      # Identify the outlier
      mean_val <- mean(values)
      outlier_idx <- which.max(abs(values - mean_val))
      outliers <- values[outlier_idx]
      values <- values[-outlier_idx]
    }
  }
  
  return(list(values = values, outliers = outliers))
}

# Function to process a single dataset
process_dataset <- function(df, sample_type) {
  # Determine the target_reads column name based on sample_type
  target_col <- ifelse(sample_type == "CTCF", 
                       "target_reads_for_CTCF",
                       "target_reads_for_K27me3")
  actual_col <- ifelse(sample_type == "CTCF",
                       "actual_reads_for_CTCF",
                       "actual_reads_for_K27me3")
  
  # Convert target_reads to numeric for ordering
  # We'll handle ordering later in the combined dataframe
  df$target_reads_numeric <- ifelse(df[[target_col]] == "original", 
                                    NA,  # Will handle separately
                                    as.numeric(df[[target_col]]))
  
  # Sort by target_reads_numeric descending
  df <- df[order(-df$target_reads_numeric, na.last = TRUE), ]
  
  # Group by target_reads and process outliers
  df_summary <- df %>%
    group_by(!!sym(target_col), target_reads_numeric) %>%
    summarise(
      # Process num_peaks_default with outlier removal
      peaks_default_raw = list(num_peaks_default),
      peaks_default_outliers = list(remove_outliers_grubbs(num_peaks_default)$outliers),
      peaks_default_clean = list(remove_outliers_grubbs(num_peaks_default)$values),
      mean_peaks_default = mean(remove_outliers_grubbs(num_peaks_default)$values),
      sd_peaks_default = sd(remove_outliers_grubbs(num_peaks_default)$values),
      n_default = length(remove_outliers_grubbs(num_peaks_default)$values),
      
      # Use mean of actual_reads for x-axis position
      mean_actual_reads = mean(!!sym(actual_col)),
      .groups = 'drop'
    ) %>%
    ungroup()
  
  # Prepare data for plotting (only Default method)
  df_plot <- df_summary %>%
    select(!!sym(target_col), target_reads_numeric, mean_actual_reads, 
           mean_peaks = mean_peaks_default, sd_peaks = sd_peaks_default, n = n_default) %>%
    mutate(method = paste0(sample_type, "_MACS2_Default"),
           sample_type = sample_type,
           target_reads = !!sym(target_col))
  
  return(df_plot)
}

# Check if dataframes exist
if (!exists("M100K_CTCF_no_vs_M100K_IgG_down")) {
  stop("Dataframe 'M100K_CTCF_no_vs_M100K_IgG_down' not found. Please load it first.")
}
if (!exists("M100K_K27me3_no_vs_M100K_IgG_down")) {
  stop("Dataframe 'M100K_K27me3_no_vs_M100K_IgG_down' not found. Please load it first.")
}

# Process both datasets
cat("Processing CTCF data...\n")
df_ctcf <- process_dataset(M100K_CTCF_no_vs_M100K_IgG_down, "CTCF")

cat("Processing K27me3 data...\n")
df_k27me3 <- process_dataset(M100K_K27me3_no_vs_M100K_IgG_down, "K27me3")

# Combine both datasets
df_plot <- bind_rows(df_ctcf, df_k27me3)

# Create x-axis labels with sample-specific "original" names
# Use consistent formatting: 15M, 12M, 10M, 5M, 2.5M, 1M, 500K, 100K, 50K, 10K
df_plot <- df_plot %>%
  mutate(
    target_numeric = ifelse(target_reads == "original", 
                            NA, 
                            as.numeric(target_reads))
  ) %>%
  mutate(
    x_label = case_when(
      target_reads == "original" ~ paste0(sample_type, "_original"),
      target_numeric == 15000000 ~ "15M",
      target_numeric == 12000000 ~ "12M",
      target_numeric == 10000000 ~ "10M",
      target_numeric == 5000000 ~ "5M",
      target_numeric == 2500000 ~ "2.5M",
      target_numeric == 1000000 ~ "1M",
      target_numeric == 500000 ~ "500K",
      target_numeric == 100000 ~ "100K",
      target_numeric == 50000 ~ "50K",
      target_numeric == 10000 ~ "10K",
      target_numeric >= 1e6 ~ paste0(format(target_numeric / 1e6, digits=2), "M"),
      TRUE ~ paste0(format(target_numeric / 1e3, digits=1), "K")
    )
  )

# Create custom ordering: K27me3_original, 15M, 12M, 10M, CTCF_original, 5M, 2.5M, 1M, 500K, 100K, 50K, 10K
# Define the desired order with exact label matches
desired_order <- c("K27me3_original", "15M", "12M", "10M", "CTCF_original", 
                   "5M", "2.5M", "1M", "500K", "100K", "50K", "10K")

# Create a unified mapping for x-axis breaks and labels with fixed intervals
# Get all unique x_labels (for numeric values, both samples share the same label)
# For "original", each sample has its own label
x_labels_map <- df_plot %>%
  select(x_label) %>%
  distinct() %>%
  mutate(
    # Create a sort key based on desired order
    sort_key = case_when(
      x_label == "K27me3_original" ~ 1,
      x_label == "15M" ~ 2,
      x_label == "12M" ~ 3,
      x_label == "10M" ~ 4,
      x_label == "CTCF_original" ~ 5,
      x_label == "5M" ~ 6,
      x_label == "2.5M" ~ 7,
      x_label == "1M" ~ 8,
      x_label == "500K" ~ 9,
      x_label == "100K" ~ 10,
      x_label == "50K" ~ 11,
      x_label == "10K" ~ 12,
      TRUE ~ 999  # Default for any unexpected values
    )
  ) %>%
  arrange(sort_key) %>%
  mutate(x_position = row_number())

# Add x_position to df_plot based on x_label (not target_reads, since original has different labels)
df_plot <- df_plot %>%
  left_join(x_labels_map %>% select(x_label, x_position), 
            by = "x_label")

# Calculate error bars for log scale
df_plot <- df_plot %>%
  mutate(
    mean_peaks_log = log10(mean_peaks + 1),  # +1 to avoid log(0)
    lower_error = pmax(1, mean_peaks - sd_peaks),  # Ensure positive
    upper_error = mean_peaks + sd_peaks,
    lower_error_log = log10(lower_error + 1),
    upper_error_log = log10(upper_error + 1)
  )

# Create method labels for legend
df_plot <- df_plot %>%
  mutate(method_label = case_when(
    method == "CTCF_MACS2_Default" ~ "CTCF",
    method == "K27me3_MACS2_Default" ~ "K27me3",
    TRUE ~ method
  ))

# Create the plot - both CTCF and K27me3 on the same plot
p <- ggplot(df_plot, aes(x = x_position, y = mean_peaks, color = method, group = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower_error, ymax = upper_error), 
                width = 0.1,
                size = 0.8, alpha = 0.7) +
  scale_color_manual(values = colors, 
                     labels = c("CTCF_MACS2_Default" = "CTCF",
                                "K27me3_MACS2_Default" = "K27me3"),
                     name = "Sample") +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    breaks = scales::trans_breaks("log10", function(x) 10^x)
  ) +
  labs(
    x = "Target Reads",
    y = "Number of Peaks (log10 scale)",
    title = "Peak Calling Saturation Curves: CTCF vs K27me3 (MACS2 Default)",
    subtitle = "X-axis: fixed intervals (target reads) | Y-axis: log10 scale"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor.y = element_line(color = "grey90", size = 0.3),
    panel.grid.major = element_line(color = "grey80", size = 0.5)
  )

# Set x-axis breaks and labels (unified for both datasets)
p <- p + scale_x_continuous(
  breaks = x_labels_map$x_position,
  labels = x_labels_map$x_label
)

# Print the plot
print(p)

# Save the plot
ggsave("M100K_CTCF_K27me3_saturation_curves_based_on_constant_merged_IgG.pdf", plot = p, width = 14, height = 6, units = "in")
ggsave("M100K_CTCF_K27me3_saturation_curves_based_on_constant_merged_IgG.png", plot = p, width = 14, height = 6, units = "in", dpi = 300)


