library(ggplot2)
library(dplyr)
library(tidyr)
library(outliers)
library(pheatmap)
library(viridis)
library(GenomicRanges)
library(gridExtra)
library(grid)
library(gtable)


setwd("/Users/atrayeeray/Desktop/Sid/CUT&Tag_cellno. standardization/Direct CUT&Tag/new_analysis_07082025/Downsampling")


##Downsample IgG (merged) reads and no downsample M100K_CTCF (replicate wise) reads 
M100K_CTCF_no_vs_M100K_IgG_down <- read.csv("M100K_CTCF_no_downsample_M100K_IgG_down/11132025/downsample_summary.csv")

# Function to detect outliers using Grubbs test
detect_outliers_grubbs <- function(data, sample_name, peak_type) {
  replicates <- data$replicate
  
  # Get the correct peak column name
  if(peak_type == "Default") {
    peak_counts <- data$num_peaks_default
  } else if(peak_type == "Large") {
    peak_counts <- data$num_peaks_large
  } else if(peak_type == "SEACR_Stringent") {
    peak_counts <- data$num_peaks_seacr_stringent
  } else if(peak_type == "SEACR_Relaxed") {
    peak_counts <- data$num_peaks_seacr_relaxed
  }
  
  if(length(unique(peak_counts)) < 2) {
    return(NULL)
  }
  
  # Perform Grubbs test
  tryCatch({
    grubbs_result <- grubbs.test(peak_counts)
    
    if(grubbs_result$p.value < 0.05) {
      # Find the outlier
      outlier_idx <- which.max(abs(peak_counts - mean(peak_counts)))
      outlier_replicate <- replicates[outlier_idx]
      outlier_value <- peak_counts[outlier_idx]
      
      return(list(
        sample = sample_name,
        peak_type = peak_type,
        outlier_replicate = outlier_replicate,
        outlier_value = outlier_value,
        p_value = grubbs_result$p.value
      ))
    } else {
      return(NULL)
    }
  }, error = function(e) {
    cat("  Error in Grubbs test for", sample_name, peak_type, ":", e$message, "\n")
    return(NULL)
  })
}

# Function to identify and mark outliers for all peak types
identify_all_outliers <- function(data) {
  # Initialize outlier columns
  data$is_outlier_default <- FALSE
  data$is_outlier_large <- FALSE
  data$is_outlier_seacr_stringent <- FALSE
  data$is_outlier_seacr_relaxed <- FALSE
  
  # Get unique samples (group by target_reads_for_IgG to treat each downsampling level separately)
  unique_targets <- unique(data$target_reads_for_IgG)
  outliers_list <- list()
  
  for(target in unique_targets) {
    if(target == "original") {
      # For original data, group by all replicates
      target_data <- data[data$target_reads_for_IgG == target, ]
      
      if(nrow(target_data) == 3) {
        # Test all four peak calling methods
        peak_types <- c("Default", "Large", "SEACR_Stringent", "SEACR_Relaxed")
        
        for(peak_type in peak_types) {
          outlier_result <- detect_outliers_grubbs(target_data, paste0("original_", target), peak_type)
          
          if(!is.null(outlier_result)) {
            outliers_list[[length(outliers_list) + 1]] <- outlier_result
            
            # Mark outliers in the appropriate column
            outlier_rows <- data$target_reads_for_IgG == target & 
              data$replicate == outlier_result$outlier_replicate
            
            if(peak_type == "Default") {
              data$is_outlier_default[outlier_rows] <- TRUE
            } else if(peak_type == "Large") {
              data$is_outlier_large[outlier_rows] <- TRUE
            } else if(peak_type == "SEACR_Stringent") {
              data$is_outlier_seacr_stringent[outlier_rows] <- TRUE
            } else if(peak_type == "SEACR_Relaxed") {
              data$is_outlier_seacr_relaxed[outlier_rows] <- TRUE
            }
          }
        }
      }
    } else {
      # For downsampled data, group by target and replicate
      target_data <- data[data$target_reads_for_IgG == target, ]
      
      if(nrow(target_data) == 3) {
        # Test all four peak calling methods
        peak_types <- c("Default", "Large", "SEACR_Stringent", "SEACR_Relaxed")
        
        for(peak_type in peak_types) {
          outlier_result <- detect_outliers_grubbs(target_data, paste0("target_", target), peak_type)
          
          if(!is.null(outlier_result)) {
            outliers_list[[length(outliers_list) + 1]] <- outlier_result
            
            # Mark outliers in the appropriate column
            outlier_rows <- data$target_reads_for_IgG == target & 
              data$replicate == outlier_result$outlier_replicate
            
            if(peak_type == "Default") {
              data$is_outlier_default[outlier_rows] <- TRUE
            } else if(peak_type == "Large") {
              data$is_outlier_large[outlier_rows] <- TRUE
            } else if(peak_type == "SEACR_Stringent") {
              data$is_outlier_seacr_stringent[outlier_rows] <- TRUE
            } else if(peak_type == "SEACR_Relaxed") {
              data$is_outlier_seacr_relaxed[outlier_rows] <- TRUE
            }
          }
        }
      }
    }
  }
  
  return(list(data = data, outliers = outliers_list))
}

# Function to create the plot with X-axis from 10M to 10K (left to right) - FIXED VERSION
create_peak_vs_igg_plot_10M_to_10K_fixed <- function(data) {
  
  # Convert to long format for plotting (INCLUDING outliers)
  plot_data_long <- data %>%
    select(actual_reads_for_IgG, replicate, 
           num_peaks_default, num_peaks_large, 
           num_peaks_seacr_stringent, num_peaks_seacr_relaxed) %>%
    pivot_longer(
      cols = c(num_peaks_default, num_peaks_large, 
               num_peaks_seacr_stringent, num_peaks_seacr_relaxed),
      names_to = "Peak_Method",
      values_to = "Peak_Count"
    ) %>%
    mutate(
      Peak_Method = case_when(
        Peak_Method == "num_peaks_default" ~ "MACS2_Default",
        Peak_Method == "num_peaks_large" ~ "MACS2_Large",
        Peak_Method == "num_peaks_seacr_stringent" ~ "SEACR_Stringent",
        Peak_Method == "num_peaks_seacr_relaxed" ~ "SEACR_Relaxed"
      ),
      # Convert to factor for proper grouping
      Peak_Method = factor(Peak_Method, 
                           levels = c("MACS2_Default", "MACS2_Large", 
                                      "SEACR_Stringent", "SEACR_Relaxed")),
      # Handle zero values for log transformation by adding 0.1
      Peak_Count_adj = ifelse(Peak_Count == 0, 0.1, Peak_Count)
    )
  
  # Create fixed interval mapping for X-axis
  # Order: 10M, 5M, 2.5M, 1M, 500K, 100K, 50K, 10K (left to right)
  x_axis_values <- c(10000000, 5000000, 2500000, 1000000, 500000, 100000, 50000, 10000)
  x_axis_labels <- c("10M", "5M", "2.5M", "1M", "500K", "100K", "50K", "10K")
  
  # Create mapping dataframe
  x_mapping <- data.frame(
    actual_reads_for_IgG = x_axis_values,
    x_label = x_axis_labels,
    x_position = 1:length(x_axis_values)
  )
  
  # Add x_position to plot_data_long based on closest match to actual_reads_for_IgG
  plot_data_long <- plot_data_long %>%
    mutate(
      x_position = sapply(actual_reads_for_IgG, function(x) {
        # Find the closest x_axis_value
        closest_idx <- which.min(abs(x_axis_values - x))
        x_mapping$x_position[closest_idx]
      })
    )
  
  # Custom Y-axis breaks including 0.1 (representing actual 0 values)
  custom_y_breaks <- c(0.1, 10, 1000, 10000, 100000)
  custom_y_labels <- c("0", "10", "1K", "10K", "100K")
  
  # Create the plot WITH outliers and zero values properly handled
  p <- ggplot(plot_data_long, aes(x = x_position, y = Peak_Count_adj, 
                                  color = Peak_Method, shape = replicate)) +
    # Individual points for each replicate with increased size
    geom_point(size = 4, alpha = 0.8) +
    # SMOOTH LINES connecting the points for each method (enhanced visibility)
    geom_smooth(aes(group = Peak_Method, color = Peak_Method), 
                method = "loess", se = FALSE, 
                linewidth = 3, alpha = 0.9) +
    # Use fixed interval X-axis (linear scale with fixed positions)
    scale_x_continuous(
      breaks = x_mapping$x_position,
      labels = x_mapping$x_label
    ) +
    # Custom Y-axis breaks with proper zero handling
    scale_y_log10(breaks = custom_y_breaks, labels = custom_y_labels) +
    # Enhanced color scheme for better contrast
    scale_color_manual(
      values = c("MACS2_Default" = "#1f77b4",      # Blue
                 "MACS2_Large" = "#ff7f0e",        # Orange  
                 "SEACR_Stringent" = "#2ca02c",    # Green
                 "SEACR_Relaxed" = "#d62728"),     # Red
      labels = c("MACS2 Default", "MACS2 Large", 
                 "SEACR Stringent", "SEACR Relaxed")
    ) +
    scale_shape_manual(
      values = c("R1" = 16, "R2" = 17, "R3" = 18),
      labels = c("R1", "R2", "R3")
    ) +
    labs(
      title = "Downsampling of IgG Reads",
      subtitle = "Peak numbers vs IgG read counts across different downsampling levels\nZero values are shown at 0.1 on log scale for visibility",
      x = "Actual Reads for IgG (fixed intervals, 10M to 10K)",
      y = "Peak Count (log scale, 0.1 = 0 peaks)",
      color = "Peak Calling Method",
      shape = "Replicate"
    ) +
    theme_minimal(base_size = 16) +  # Increased base size
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      # Horizontal X-axis labels for better readability
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, 
                                 margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold", 
                                  margin = margin(t = 15, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
      # Add more space for X-axis labels
      plot.margin = margin(t = 20, r = 20, b = 50, l = 20)
    ) +
    # Set X-axis limits to match the fixed positions (1 to 8)
    coord_cartesian(xlim = c(1, length(x_axis_values)))
  
  return(p)
}


# Alternative approach: Create a transformed variable
create_peak_vs_igg_plot_transformed <- function(data) {
  
  # Convert to long format for plotting (INCLUDING outliers)
  plot_data_long <- data %>%
    select(actual_reads_for_IgG, replicate, 
           num_peaks_default, num_peaks_large, 
           num_peaks_seacr_stringent, num_peaks_seacr_relaxed) %>%
    pivot_longer(
      cols = c(num_peaks_default, num_peaks_large, 
               num_peaks_seacr_stringent, num_peaks_seacr_relaxed),
      names_to = "Peak_Method",
      values_to = "Peak_Count"
    ) %>%
    mutate(
      Peak_Method = case_when(
        Peak_Method == "num_peaks_default" ~ "MACS2_Default",
        Peak_Method == "num_peaks_large" ~ "MACS2_Large",
        Peak_Method == "num_peaks_seacr_stringent" ~ "SEACR_Stringent",
        Peak_Method == "num_peaks_seacr_relaxed" ~ "SEACR_Relaxed"
      ),
      # Convert to factor for proper grouping
      Peak_Method = factor(Peak_Method, 
                           levels = c("MACS2_Default", "MACS2_Large", 
                                      "SEACR_Stringent", "SEACR_Relaxed")),
      # Handle zero values for log transformation by adding 0.1
      Peak_Count_adj = ifelse(Peak_Count == 0, 0.1, Peak_Count),
      # Transform IgG reads so that 10M becomes 1 and 10K becomes 10
      IgG_reads_transformed = 10000000 / actual_reads_for_IgG
    )
  
  # Custom Y-axis breaks including 0.1 (representing actual 0 values)
  custom_y_breaks <- c(0.1, 10, 1000, 10000, 100000)
  custom_y_labels <- c("0", "10", "1K", "10K", "100K")
  
  # Create the plot WITH outliers and zero values properly handled
  p <- ggplot(plot_data_long, aes(x = IgG_reads_transformed, y = Peak_Count_adj, 
                                  color = Peak_Method, shape = replicate)) +
    # Individual points for each replicate with increased size
    geom_point(size = 4, alpha = 0.8) +
    # SMOOTH LINES connecting the points for each method (enhanced visibility)
    geom_smooth(aes(group = Peak_Method, color = Peak_Method), 
                method = "loess", se = FALSE, 
                linewidth = 3, alpha = 0.9) +
    # Use log scale for transformed values
    scale_x_log10(breaks = c(1, 2, 4, 10, 20, 100, 200, 1000),
                  labels = c("10M", "5M", "2.5M", "1M", "500K", "100K", "50K", "10K")) +
    # Custom Y-axis breaks with proper zero handling
    scale_y_log10(breaks = custom_y_breaks, labels = custom_y_labels) +
    # Enhanced color scheme for better contrast
    scale_color_manual(
      values = c("MACS2_Default" = "#1f77b4",      # Blue
                 "MACS2_Large" = "#ff7f0e",        # Orange  
                 "SEACR_Stringent" = "#2ca02c",    # Green
                 "SEACR_Relaxed" = "#d62728"),     # Red
      labels = c("MACS2 Default", "MACS2 Large", 
                 "SEACR Stringent", "SEACR Relaxed")
    ) +
    scale_shape_manual(
      values = c("R1" = 16, "R2" = 17, "R3" = 18),
      labels = c("R1", "R2", "R3")
    ) +
    labs(
      title = "Downsampling of IgG Reads",
      subtitle = "Peak numbers vs IgG read counts across different downsampling levels\nZero values are shown at 0.1 on log scale for visibility",
      x = "Actual Reads for IgG (log scale, 10M to 10K)",
      y = "Peak Count (log scale, 0.1 = 0 peaks)",
      color = "Peak Calling Method",
      shape = "Replicate"
    ) +
    theme_minimal(base_size = 16) +  # Increased base size
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      # Horizontal X-axis labels for better readability
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, 
                                 margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold", 
                                  margin = margin(t = 15, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
      # Add more space for X-axis labels
      plot.margin = margin(t = 20, r = 20, b = 50, l = 20)
    )
  
  return(p)
}

# Updated main execution function
main_analysis_M100K_CTCF_no_vs_M100K_IgG_down <- function() {
  cat("Starting enhanced peak analysis with X-axis from 10M to 10K...\n")
  
  # Step 1: Identify outliers
  cat("Step 1: Identifying outliers using Grubbs test...\n")
  outlier_results <- identify_all_outliers(M100K_CTCF_no_vs_M100K_IgG_down)
  M100K_CTCF_no_vs_M100K_IgG_down_with_outliers <- outlier_results$data
  
  # Count outliers for each peak type
  outlier_counts <- list(
    Default = sum(M100K_CTCF_no_vs_M100K_IgG_down_with_outliers$is_outlier_default),
    Large = sum(M100K_CTCF_no_vs_M100K_IgG_down_with_outliers$is_outlier_large),
    SEACR_Stringent = sum(M100K_CTCF_no_vs_M100K_IgG_down_with_outliers$is_outlier_seacr_stringent),
    SEACR_Relaxed = sum(M100K_CTCF_no_vs_M100K_IgG_down_with_outliers$is_outlier_seacr_relaxed)
  )
  
  cat("Outliers detected:\n")
  for(peak_type in names(outlier_counts)) {
    cat("  ", peak_type, ":", outlier_counts[[peak_type]], "\n")
  }
  
  # Step 2: Create plots with X-axis from 10M to 10K
  cat("\nStep 2: Creating plots with X-axis from 10M to 10K...\n")
  
  # Method 1: Using coord_cartesian
  p_coord_cartesian <- create_peak_vs_igg_plot_10M_to_10K_fixed(M100K_CTCF_no_vs_M100K_IgG_down_with_outliers)
  
  # Method 2: Using transformed variable
  p_transformed <- create_peak_vs_igg_plot_transformed(M100K_CTCF_no_vs_M100K_IgG_down_with_outliers)
  
  # Step 3: Save plots
  cat("\nStep 3: Saving plots...\n")
  
  ggsave("peak_numbers_vs_igg_reads_10M_to_10K_coord.pdf", p_coord_cartesian, 
         width = 16, height = 12, dpi = 300)
  ggsave("peak_numbers_vs_igg_reads_10M_to_10K_transformed.pdf", p_transformed, 
         width = 16, height = 12, dpi = 300)
  
  cat("Plots saved:\n")
  cat("  peak_numbers_vs_igg_reads_10M_to_10K_coord.pdf\n")
  cat("  peak_numbers_vs_igg_reads_10M_to_10K_transformed.pdf\n")
  
  # Step 4: Summary statistics
  cat("\nStep 4: Summary statistics...\n")
  
  # Summary by peak method and IgG read level
  summary_stats <- M100K_CTCF_no_vs_M100K_IgG_down_with_outliers %>%
    group_by(target_reads_for_IgG) %>%
    summarise(
      MACS2_Default_Mean = mean(num_peaks_default, na.rm = TRUE),
      MACS2_Large_Mean = mean(num_peaks_large, na.rm = TRUE),
      SEACR_Stringent_Mean = mean(num_peaks_seacr_stringent, na.rm = TRUE),
      SEACR_Relaxed_Mean = mean(num_peaks_seacr_relaxed, na.rm = TRUE),
      .groups = 'drop'
    )
  
  print(summary_stats)
  
  cat("\nAnalysis complete!\n")
  
  return(list(
    data = M100K_CTCF_no_vs_M100K_IgG_down_with_outliers,
    outliers = outlier_results$outliers,
    plot_coord = p_coord_cartesian,
    plot_transformed = p_transformed,
    summary = summary_stats
  ))
}

# Run the analysis
results_fixed <- main_analysis_M100K_CTCF_no_vs_M100K_IgG_down()

# Display the plots
print(results_fixed$plot_coord)
print(results_fixed$plot_transformed)

##### Frip score based heatmap ####

M100K_CTCF_no_vs_M100K_IgG_down_frip <- read.csv("frip_results_downsampled_peak_results_for_merge_M100K_IgG_no_downsample_M100K_CTCF/frip_scores_summary.csv")
M100K_CTCF_no_vs_M100K_IgG_down_frip <- M100K_CTCF_no_vs_M100K_IgG_down_frip[,c(1:4,6:9)]

#Function to calculate standard deviation across replicates for FRiP scores
calculate_frip_std_dev <- function(frip_data) {
  
  cat("Calculating standard deviation across replicates for FRiP scores...\n")
  
  # Check if we have replicates
  if("Replicate" %in% colnames(frip_data)) {
    cat("Found replicates, calculating mean ± SD across R1, R2, R3...\n")
    
    # Calculate mean and SD FRiP scores for each IgG level across replicates
    summary_data <- frip_data %>%
      group_by(IgG_Downsampled_Reads) %>%
      summarise(
        MACS2_Default_Mean = mean(FRiP_macs_default, na.rm = TRUE),
        MACS2_Default_SD = sd(FRiP_macs_default, na.rm = TRUE),
        MACS2_Large_Mean = mean(FRiP_macs_large, na.rm = TRUE),
        MACS2_Large_SD = sd(FRiP_macs_large, na.rm = TRUE),
        SEACR_Stringent_Mean = mean(FRiP_seacr_stringent, na.rm = TRUE),
        SEACR_Stringent_SD = sd(FRiP_seacr_stringent, na.rm = TRUE),
        SEACR_Relaxed_Mean = mean(FRiP_seacr_relaxed, na.rm = TRUE),
        SEACR_Relaxed_SD = sd(FRiP_seacr_relaxed, na.rm = TRUE),
        .groups = 'drop'
      )
    
    cat("✓ Mean and SD calculated for each method across replicates\n")
    
  } else {
    cat("No replicates found, using single measurements...\n")
    
    # For single measurements, set SD to 0
    summary_data <- frip_data %>%
      group_by(IgG_Downsampled_Reads) %>%
      summarise(
        MACS2_Default_Mean = mean(FRiP_macs_default, na.rm = TRUE),
        MACS2_Default_SD = 0,
        MACS2_Large_Mean = mean(FRiP_macs_large, na.rm = TRUE),
        MACS2_Large_SD = 0,
        SEACR_Stringent_Mean = mean(FRiP_seacr_stringent, na.rm = TRUE),
        SEACR_Stringent_SD = 0,
        SEACR_Relaxed_Mean = mean(FRiP_seacr_relaxed, na.rm = TRUE),
        SEACR_Relaxed_SD = 0,
        .groups = 'drop'
      )
    
    cat("⚠ Single measurements detected - SD set to 0\n")
  }
  
  return(summary_data)
}

# Function to create FRiP Score heatmap with standard deviation
create_frip_heatmap_with_sd <- function(frip_data) {
  
  cat("Creating FRiP Score heatmap with standard deviation (±)...\n")
  
  # Calculate mean and SD
  summary_data <- calculate_frip_std_dev(frip_data)
  
  # Convert to long format for plotting
  plot_data <- summary_data %>%
    pivot_longer(cols = c(MACS2_Default_Mean, MACS2_Large_Mean, 
                          SEACR_Stringent_Mean, SEACR_Relaxed_Mean), 
                 names_to = "Method", values_to = "FRiP_Score") %>%
    mutate(
      Method = case_when(
        Method == "MACS2_Default_Mean" ~ "MACS2 Default",
        Method == "MACS2_Large_Mean" ~ "MACS2 Large",
        Method == "SEACR_Stringent_Mean" ~ "SEACR Stringent",
        Method == "SEACR_Relaxed_Mean" ~ "SEACR Relaxed"
      ),
      Method = factor(Method, levels = c("SEACR Relaxed", "SEACR Stringent", "MACS2 Large", "MACS2 Default"))
    )
  
  # Add corresponding SD values
  plot_data <- plot_data %>%
    mutate(
      SD_Value = case_when(
        Method == "MACS2 Default" ~ summary_data$MACS2_Default_SD[match(IgG_Downsampled_Reads, summary_data$IgG_Downsampled_Reads)],
        Method == "MACS2 Large" ~ summary_data$MACS2_Large_SD[match(IgG_Downsampled_Reads, summary_data$IgG_Downsampled_Reads)],
        Method == "SEACR Stringent" ~ summary_data$SEACR_Stringent_SD[match(IgG_Downsampled_Reads, summary_data$IgG_Downsampled_Reads)],
        Method == "SEACR Relaxed" ~ summary_data$SEACR_Relaxed_SD[match(IgG_Downsampled_Reads, summary_data$IgG_Downsampled_Reads)]
      )
    )
  
  # Create IgG depth labels with EXACT string matching
  plot_data <- plot_data %>%
    mutate(
      IgG_Depth_Label = case_when(
        IgG_Downsampled_Reads == "original" ~ "Original",
        IgG_Downsampled_Reads == "5000000" ~ "5M",
        IgG_Downsampled_Reads == "2500000" ~ "2.5M",
        IgG_Downsampled_Reads == "1000000" ~ "1M",
        IgG_Downsampled_Reads == "500000" ~ "500K",
        IgG_Downsampled_Reads == "100000" ~ "100K",
        IgG_Downsampled_Reads == "50000" ~ "50K",
        IgG_Downsampled_Reads == "10000" ~ "10K",
        TRUE ~ paste0("Unknown_", IgG_Downsampled_Reads)
      )
    ) %>%
    # Set factor levels in the order we want
    mutate(IgG_Depth_Label = factor(IgG_Depth_Label, 
                                    levels = c("Original", "5M", "2.5M", "1M", "500K", "100K", "50K", "10K")))
  
  cat("Final IgG depth labels:", paste(levels(plot_data$IgG_Depth_Label), collapse = ", "), "\n")
  
  # Remove any rows with Unknown labels
  plot_data <- plot_data %>%
    filter(!grepl("Unknown_", IgG_Depth_Label))
  
  cat("Number of data points after filtering:", nrow(plot_data), "\n")
  
  # Get the range of FRiP scores for color scaling
  score_range <- range(plot_data$FRiP_Score, na.rm = TRUE)
  cat("FRiP score range:", score_range[1], "to", score_range[2], "\n")
  
  # Create the FRiP Score heatmap with standard deviation
  p <- ggplot(plot_data, aes(x = IgG_Depth_Label, y = Method, fill = FRiP_Score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    
    # Add readable text labels (FRiP scores ± SD, rounded to 4 decimal places)
    geom_text(aes(label = sprintf("%.4f", FRiP_Score)), 
              color = "black", size = 3.5, fontface = "bold",
              position = position_nudge(y = -0.15)) +
    
    # Add standard deviation labels below the main values
    geom_text(aes(label = sprintf("±%.4f", SD_Value)), 
              color = "darkblue", size = 2.8, fontface = "italic",
              position = position_nudge(y = 0.15)) +
    
    # Use the same color scheme as your efficiency plot (adapted for FRiP scores)
    scale_fill_gradient2(
      low = "#313695",      # Dark blue for low values
      mid = "#ffffcc",      # Pale yellow for middle values
      high = "#a50026",     # Dark red for high values
      midpoint = median(plot_data$FRiP_Score, na.rm = TRUE),
      na.value = "grey90",
      trans = "sqrt"
    ) +
    
    # Labels and theme
    labs(
      title = "FRiP Score Heatmap with Standard Deviation",
      subtitle = "FRiP = Fraction of Reads in Peaks\nValues shown as Mean ± SD across replicates (R1, R2, R3)",
      x = "IgG Read Depth (Background Noise Level)",
      y = "Peak Calling Method",
      fill = "FRiP Score"
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid = element_blank(),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    
    # Ensure tiles are square
    coord_fixed(ratio = 1)
  
  return(list(plot = p, summary_data = summary_data))
}

# Function to create summary statistics table
create_summary_table <- function(summary_data) {
  
  # Prepare data for table display
  table_data <- summary_data %>%
    mutate(
      IgG_Depth_Label = case_when(
        IgG_Downsampled_Reads == "original" ~ "Original",
        IgG_Downsampled_Reads == "5000000" ~ "5M",
        IgG_Downsampled_Reads == "2500000" ~ "2.5M",
        IgG_Downsampled_Reads == "1000000" ~ "1M",
        IgG_Downsampled_Reads == "500000" ~ "500K",
        IgG_Downsampled_Reads == "100000" ~ "100K",
        IgG_Downsampled_Reads == "50000" ~ "50K",
        IgG_Downsampled_Reads == "10000" ~ "10K",
        TRUE ~ IgG_Downsampled_Reads
      )
    ) %>%
    select(
      `IgG Level` = IgG_Depth_Label,
      `MACS2 Default (Mean ± SD)` = MACS2_Default_Mean,
      `MACS2 Large (Mean ± SD)` = MACS2_Large_Mean,
      `SEACR Stringent (Mean ± SD)` = SEACR_Stringent_Mean,
      `SEACR Relaxed (Mean ± SD)` = SEACR_Relaxed_Mean
    ) %>%
    arrange(factor(`IgG Level`, levels = c("Original", "5M", "2.5M", "1M", "500K", "100K", "50K", "10K")))
  
  # Format the values to show Mean ± SD
  formatted_table <- table_data %>%
    mutate(
      `MACS2 Default (Mean ± SD)` = sprintf("%.4f ± %.4f", 
                                            summary_data$MACS2_Default_Mean[match(`IgG Level`, table_data$`IgG Level`)],
                                            summary_data$MACS2_Default_SD[match(`IgG Level`, table_data$`IgG Level`)]),
      `MACS2 Large (Mean ± SD)` = sprintf("%.4f ± %.4f", 
                                          summary_data$MACS2_Large_Mean[match(`IgG Level`, table_data$`IgG Level`)],
                                          summary_data$MACS2_Large_SD[match(`IgG Level`, table_data$`IgG Level`)]),
      `SEACR Stringent (Mean ± SD)` = sprintf("%.4f ± %.4f", 
                                              summary_data$SEACR_Stringent_Mean[match(`IgG Level`, table_data$`IgG Level`)],
                                              summary_data$SEACR_Stringent_SD[match(`IgG Level`, table_data$`IgG Level`)]),
      `SEACR Relaxed (Mean ± SD)` = sprintf("%.4f ± %.4f", 
                                            summary_data$SEACR_Relaxed_Mean[match(`IgG Level`, table_data$`IgG Level`)],
                                            summary_data$SEACR_Relaxed_SD[match(`IgG Level`, table_data$`IgG Level`)])
    )
  
  # Create table grob
  table_grob <- tableGrob(formatted_table, 
                          rows = NULL,
                          theme = ttheme_minimal(
                            base_size = 9,
                            core = list(
                              fg_params = list(hjust = 0, x = 0.1),
                              bg_params = list(fill = c("white", "lightgray"))
                            ),
                            colhead = list(
                              fg_params = list(fontface = "bold"),
                              bg_params = list(fill = "lightblue")
                            )
                          ))
  
  # Add title to table
  table_title <- textGrob("FRiP Score Summary (Mean ± SD)", 
                          gp = gpar(fontsize = 14, fontface = "bold"),
                          just = "center")
  
  # Combine table elements
  table_with_title <- arrangeGrob(table_title, table_grob, 
                                  heights = c(0.1, 0.9),
                                  ncol = 1)
  
  return(table_with_title)
}

# Main function for FRiP Score analysis with standard deviation
main_frip_analysis_with_sd <- function(data) {
  
  cat("=== FRiP SCORE ANALYSIS WITH STANDARD DEVIATION ===\n")
  cat("1. ✓ Calculate FRiP scores per replicate then average\n")
  cat("2. ✓ Calculate standard deviation (±) across replicates\n")
  cat("3. ✓ Display Mean ± SD in each heatmap cell\n")
  cat("4. ✓ Create heatmap with working color gradient\n")
  cat("5. ✓ Include summary statistics table\n")
  cat("6. ✓ All IgG levels: Original, 5M, 2.5M, 1M, 500K, 100K, 50K, 10K\n\n")
  
  cat("FRiP = Fraction of Reads in Peaks\n")
  cat("Will calculate FRiP for EACH REPLICATE individually, then show Mean ± SD\n")
  cat("Standard deviation shows reproducibility across R1, R2, R3\n\n")
  
  # Step 1: Create FRiP heatmap with standard deviation
  cat("Step 1: Creating FRiP heatmap with standard deviation...\n")
  heatmap_results <- create_frip_heatmap_with_sd(data)
  frip_heatmap <- heatmap_results$plot
  summary_data <- heatmap_results$summary_data
  
  # Step 2: Create summary statistics table
  cat("Step 2: Creating summary statistics table...\n")
  summary_table <- create_summary_table(summary_data)
  
  # Step 3: Save plots
  cat("Step 3: Saving all plots...\n")
  
  # Save individual heatmap
  ggsave("FRiP_Score_Heatmap_with_SD.pdf", frip_heatmap, width = 14, height = 8, dpi = 300)
  
  # Save combined plot with table
  pdf("FRiP_Score_Complete_Analysis_with_SD.pdf", width = 18, height = 10)
  
  # Arrange plots side by side
  grid.arrange(
    frip_heatmap, summary_table,
    ncol = 2,
    widths = c(0.65, 0.35),
    top = textGrob("FRiP Score Analysis with Standard Deviation", 
                   gp = gpar(fontsize = 18, fontface = "bold"))
  )
  
  dev.off()
  
  # Step 4: Display results
  cat("\n=== FRiP SCORE ANALYSIS RESULTS ===\n")
  cat("✓ Standard deviation calculated across replicates (R1, R2, R3)\n")
  cat("✓ Mean ± SD displayed in each heatmap cell\n")
  cat("✓ Working color gradient applied\n")
  cat("✓ All IgG levels included: Original, 5M, 2.5M, 1M, 500K, 100K, 50K, 10K\n")
  cat("✓ FRiP calculated for EACH REPLICATE individually, then averaged\n")
  cat("✓ Summary statistics table created with Mean ± SD\n\n")
  
  cat("All plots saved!\n")
  cat("Analysis complete with standard deviation!\n")
  
  return(list(
    summary_data = summary_data,
    frip_heatmap = frip_heatmap,
    summary_table = summary_table
  ))
}

# Example usage:
results <- main_frip_analysis_with_sd(M100K_CTCF_no_vs_M100K_IgG_down_frip)


############ macs2_large overlap percentage with IgG downsampled reads ##########

setwd("/Users/atrayeeray/Desktop/Sid/CUT&Tag_cellno. standardization/Direct CUT&Tag/new_analysis_07082025/Downsampling/M100K_CTCF_no_downsample_M100K_IgG_down/11132025/macs2_peaks_large/")

######## Replicate level analysis 

# Function to read narrowPeak files
read_narrowPeak <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(GRanges())
  }
  
  # Read the narrowPeak file
  peaks <- read.table(file_path, header = FALSE, sep = "\t")
  
  # Create GRanges object
  gr <- GRanges(
    seqnames = peaks$V1,
    ranges = IRanges(start = peaks$V2 + 1, end = peaks$V3),
    strand = "*",
    name = peaks$V4,
    score = peaks$V5,
    signalValue = peaks$V7,
    pValue = peaks$V8,
    qValue = peaks$V9
  )
  
  return(gr)
}

# Function to calculate overlap percentage (ROW-WISE: relative to row peaks)
calculate_overlap_percentage_rowwise <- function(gr1, gr2, gr1_name, gr2_name) {
  if (length(gr1) == 0 || length(gr2) == 0) {
    return(0)
  }
  
  # Find overlapping peaks
  overlapping_peaks <- GenomicRanges::intersect(gr1, gr2)
  overlap_count <- length(overlapping_peaks)
  
  # Calculate overlap percentage relative to gr1 (ROW) - how much of row peaks overlap with column
  overlap_percentage <- ifelse(length(gr1) > 0, (overlap_count / length(gr1)) * 100, 0)
  
  return(overlap_percentage)
}

# Function to calculate overlap percentage (COLUMN-WISE: relative to column peaks)
calculate_overlap_percentage_columnwise <- function(gr1, gr2, gr1_name, gr2_name) {
  if (length(gr1) == 0 || length(gr2) == 0) {
    return(0)
  }
  
  # Find overlapping peaks
  overlapping_peaks <- GenomicRanges::intersect(gr1, gr2)
  overlap_count <- length(overlapping_peaks)
  
  # Calculate overlap percentage relative to gr2 (COLUMN) - how much of column peaks overlap with row
  overlap_percentage <- ifelse(length(gr2) > 0, (overlap_count / length(gr2)) * 100, 0)
  
  return(overlap_percentage)
}

# Function to create overlap matrix for a specific replicate
create_overlap_matrix <- function(file_list, replicate = "R1", overlap_type = "rowwise") {
  cat("Creating", overlap_type, "overlap matrix for", replicate, "...\n")
  
  # Read all peak files for this replicate
  peak_granges <- list()
  file_names <- c()
  
  for (file_path in file_list) {
    if (grepl(replicate, file_path)) {
      cat("Reading:", basename(file_path), "\n")
      gr <- read_narrowPeak(file_path)
      peak_granges[[length(peak_granges) + 1]] <- gr
      
      # Extract the IgG read count from filename for labeling
      if (grepl("original", file_path)) {
        file_names <- c(file_names, "Original")
      } else {
        number_match <- regexpr("\\d+_q0\\.00001_large_peaks", file_path)
        if (number_match != -1) {
          number_str <- substr(file_path, number_match, number_match + attr(number_match, "match.length") - 1)
          number_val <- as.numeric(gsub("_q0\\.00001_large_peaks", "", number_str))
          
          if (number_val >= 1000000) {
            label <- paste0(number_val / 1000000, "M")
          } else if (number_val >= 1000) {
            label <- paste0(number_val / 1000, "K")
          } else {
            label <- as.character(number_val)
          }
          file_names <- c(file_names, label)
        } else {
          file_names <- c(file_names, basename(file_path))
        }
      }
    }
  }
  
  if (length(peak_granges) == 0) {
    stop("No files found for replicate ", replicate)
  }
  
  cat("Found", length(peak_granges), "peak files for", replicate, "\n")
  cat("File labels:", paste(file_names, collapse = ", "), "\n")
  
  # Calculate pairwise overlaps
  n_files <- length(peak_granges)
  overlap_matrix <- matrix(0, nrow = n_files, ncol = n_files)
  
  # Choose overlap calculation function based on type
  overlap_func <- if(overlap_type == "rowwise") {
    calculate_overlap_percentage_rowwise
  } else {
    calculate_overlap_percentage_columnwise
  }
  
  # Fill the matrix
  for (i in 1:n_files) {
    for (j in 1:n_files) {
      if (i == j) {
        overlap_matrix[i, j] <- 100  # Self-overlap = 100%
      } else {
        overlap_matrix[i, j] <- overlap_func(peak_granges[[i]], peak_granges[[j]], 
                                             file_names[i], file_names[j])
      }
    }
  }
  
  # Set row and column names
  rownames(overlap_matrix) <- file_names
  colnames(overlap_matrix) <- file_names
  
  return(list(
    overlap_matrix = overlap_matrix,
    file_names = file_names,
    peak_counts = sapply(peak_granges, length),
    overlap_type = overlap_type
  ))
}

# Function to calculate mean overlap matrices across replicates with standard deviation
calculate_replicate_statistics <- function(replicate_results_list, use_mean = TRUE) {
  cat("Calculating", ifelse(use_mean, "mean", "median"), "overlap matrices across replicates...\n")
  
  # Extract matrices from each replicate
  n_replicates <- length(replicate_results_list)
  n_files <- nrow(replicate_results_list[[1]]$overlap_matrix)
  
  # Initialize array to store values for each quadrant across replicates
  overlap_array <- array(0, dim = c(n_files, n_files, n_replicates))
  
  # Fill array with values from each replicate
  for (i in 1:n_replicates) {
    overlap_array[,,i] <- replicate_results_list[[i]]$overlap_matrix
  }
  
  # Calculate mean or median for each quadrant
  if (use_mean) {
    central_overlap <- apply(overlap_array, c(1,2), mean, na.rm = TRUE)
    stat_name <- "mean"
  } else {
    central_overlap <- apply(overlap_array, c(1,2), median, na.rm = TRUE)
    stat_name <- "median"
  }
  
  # Calculate standard deviation for each quadrant
  sd_overlap <- apply(overlap_array, c(1,2), sd, na.rm = TRUE)
  
  # Set row and column names
  rownames(central_overlap) <- replicate_results_list[[1]]$file_names
  colnames(central_overlap) <- replicate_results_list[[1]]$file_names
  rownames(sd_overlap) <- replicate_results_list[[1]]$file_names
  colnames(sd_overlap) <- replicate_results_list[[1]]$file_names
  
  return(list(
    central_overlap = central_overlap,
    sd_overlap = sd_overlap,
    file_names = replicate_results_list[[1]]$file_names,
    n_replicates = n_replicates,
    stat_name = stat_name,
    overlap_type = replicate_results_list[[1]]$overlap_type
  ))
}

# Function to create ggplot heatmap with standard deviation (LOWER TRIANGULAR)
create_ggplot_heatmap_with_sd <- function(central_results, title = "Peak Overlap Analysis") {
  
  n_rows <- nrow(central_results$central_overlap)
  n_cols <- ncol(central_results$central_overlap)
  
  # Build plotting dataframe for LOWER TRIANGULAR only
  plot_data <- data.frame(
    Row = rep(1:n_rows, each = n_cols),
    Column = rep(1:n_cols, times = n_rows),
    Value = as.vector(central_results$central_overlap),
    SD_Value = as.vector(central_results$sd_overlap),
    RowName = rep(rownames(central_results$central_overlap), each = n_cols),
    ColName = rep(colnames(central_results$central_overlap), times = n_rows),
    stringsAsFactors = FALSE
  )
  
  # Filter to show only LOWER TRIANGULAR (diagonal and below)
  plot_data <- plot_data %>%
    filter(Row >= Column)  # This creates the lower triangular pattern
  
  # Set factor levels to preserve order
  plot_data$RowName <- factor(plot_data$RowName, levels = rownames(central_results$central_overlap))
  plot_data$ColName <- factor(plot_data$ColName, levels = colnames(central_results$central_overlap))
  
  # Build heatmap - CLEAN STYLE with NUMBERS like before
  p <- ggplot(plot_data, aes(x = ColName, y = RowName, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    
    # ADD BACK THE NUMERICAL VALUES - Main overlap percentage
    geom_text(aes(label = sprintf("%.1f", Value)), 
              color = "black", size = 3.5, fontface = "bold") +
    
    # ADD BACK THE STANDARD DEVIATION - below main values
    geom_text(aes(label = sprintf("±%.2f", SD_Value)), 
              color = "darkblue", size = 2.5, fontface = "italic",
              position = position_nudge(y = -0.15)) +
    
    # Use diverging color palette like the image (red-white-blue)
    scale_fill_gradient2(
      low = "#313695",      # Dark blue for low values
      mid = "white",        # White for middle values  
      high = "#a50026",     # Dark red for high values
      midpoint = 50,        # 50% overlap as neutral point
      limits = c(0, 100),
      na.value = "grey90"
    ) +
    
    labs(
      title = title,
      subtitle = sprintf("Lower triangular matrix showing Mean ± SD across %d replicates (R1, R2, R3)", central_results$n_replicates),
      fill = "Overlap %",
      x = "Column",
      y = "Row"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(hjust = 1, size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right",
      plot.margin = margin(20, 20, 20, 20),
      panel.grid = element_blank()  # Remove grid lines for cleaner look
    ) +
    
    coord_fixed()
  
  return(p)
}

# Function to create summary statistics table with Mean ± SD
create_summary_table <- function(central_results) {
  
  # Create a comprehensive table showing Mean ± SD for each cell
  n_files <- length(central_results$file_names)
  
  # Build table data for lower triangular matrix
  table_data <- data.frame(
    Row = rep(1:n_files, each = n_files),
    Column = rep(1:n_files, times = n_files),
    RowName = rep(central_results$file_names, each = n_files),
    ColName = rep(central_results$file_names, times = n_files),
    Mean_Value = as.vector(central_results$central_overlap),
    SD_Value = as.vector(central_results$sd_overlap),
    stringsAsFactors = FALSE
  )
  
  # Filter to lower triangular only
  table_data <- table_data %>%
    filter(Row >= Column) %>%
    mutate(
      `Comparison` = paste(RowName, "→", ColName),
      `Mean ± SD` = sprintf("%.1f ± %.2f", Mean_Value, SD_Value)
    ) %>%
    select(`Comparison`, `Mean ± SD`)
  
  # Create table grob
  table_grob <- tableGrob(table_data, 
                          rows = NULL,
                          theme = ttheme_minimal(
                            base_size = 8,
                            core = list(
                              fg_params = list(hjust = 0, x = 0.1),
                              bg_params = list(fill = c("white", "lightgray"))
                            ),
                            colhead = list(
                              fg_params = list(fontface = "bold"),
                              bg_params = list(fill = "lightblue")
                            )
                          ))
  
  # Add title to table
  table_title <- textGrob("Overlap Statistics (Mean ± SD)", 
                          gp = gpar(fontsize = 12, fontface = "bold"),
                          just = "center")
  
  # Combine table elements
  table_with_title <- arrangeGrob(table_title, table_grob, 
                                  heights = c(0.1, 0.9),
                                  ncol = 1)
  
  return(table_with_title)
}

# Main execution function for REPLICATE ANALYSIS with TWO SEPARATE HEATMAPS
main_replicate_analysis_two_heatmaps <- function() {
  cat("Starting REPLICATE ANALYSIS with TWO SEPARATE HEATMAPS (Row-wise vs Column-wise)...\n")
  
  # Define your file paths for all replicates
  file_list <- c(
    "M100K_CTCF_R1_based_on_IgG_merged_original_large_peaks.narrowPeak",
    "M100K_CTCF_R1_based_on_IgG_merged_5000000_large_peaks.narrowPeak",
    "M100K_CTCF_R1_based_on_IgG_merged_2500000_large_peaks.narrowPeak",
    "M100K_CTCF_R1_based_on_IgG_merged_1000000_large_peaks.narrowPeak",
    "M100K_CTCF_R1_based_on_IgG_merged_500000_large_peaks.narrowPeak",
    "M100K_CTCF_R1_based_on_IgG_merged_100000_large_peaks.narrowPeak",
    "M100K_CTCF_R1_based_on_IgG_merged_50000_large_peaks.narrowPeak",
    "M100K_CTCF_R1_based_on_IgG_merged_10000_large_peaks.narrowPeak",
    
    "M100K_CTCF_R2_based_on_IgG_merged_original_large_peaks.narrowPeak",
    "M100K_CTCF_R2_based_on_IgG_merged_5000000_large_peaks.narrowPeak",
    "M100K_CTCF_R2_based_on_IgG_merged_2500000_large_peaks.narrowPeak",
    "M100K_CTCF_R2_based_on_IgG_merged_1000000_large_peaks.narrowPeak",
    "M100K_CTCF_R2_based_on_IgG_merged_500000_large_peaks.narrowPeak",
    "M100K_CTCF_R2_based_on_IgG_merged_100000_large_peaks.narrowPeak",
    "M100K_CTCF_R2_based_on_IgG_merged_50000_large_peaks.narrowPeak",
    "M100K_CTCF_R2_based_on_IgG_merged_10000_large_peaks.narrowPeak",
    
    "M100K_CTCF_R3_based_on_IgG_merged_original_large_peaks.narrowPeak",
    "M100K_CTCF_R3_based_on_IgG_merged_5000000_large_peaks.narrowPeak",
    "M100K_CTCF_R3_based_on_IgG_merged_2500000_large_peaks.narrowPeak",
    "M100K_CTCF_R3_based_on_IgG_merged_1000000_large_peaks.narrowPeak",
    "M100K_CTCF_R3_based_on_IgG_merged_500000_large_peaks.narrowPeak",
    "M100K_CTCF_R3_based_on_IgG_merged_100000_large_peaks.narrowPeak",
    "M100K_CTCF_R3_based_on_IgG_merged_50000_large_peaks.narrowPeak",
    "M100K_CTCF_R3_based_on_IgG_merged_10000_large_peaks.narrowPeak"
  )
  
  # Check which files exist
  existing_files <- file_list[file.exists(file_list)]
  cat("Found", length(existing_files), "existing files\n")
  
  # Analyze each replicate for BOTH overlap types
  replicates <- c("R1", "R2", "R3")
  replicate_results_rowwise <- list()
  replicate_results_columnwise <- list()
  
  for (replicate in replicates) {
    cat("\n=== Analyzing replicate", replicate, "===\n")
    
    # Check if files exist for this replicate
    replicate_files <- existing_files[grepl(replicate, existing_files)]
    
    if (length(replicate_files) > 0) {
      # Perform ROW-WISE overlap analysis for this replicate
      cat("Creating ROW-WISE overlap matrix...\n")
      results_rowwise <- create_overlap_matrix(replicate_files, replicate, "rowwise")
      replicate_results_rowwise[[replicate]] <- results_rowwise
      
      # Perform COLUMN-WISE overlap analysis for this replicate
      cat("Creating COLUMN-WISE overlap matrix...\n")
      results_columnwise <- create_overlap_matrix(replicate_files, replicate, "columnwise")
      replicate_results_columnwise[[replicate]] <- results_columnwise
      
      cat("Completed analysis for", replicate, "\n")
      
      # Print individual replicate results
      cat("ROW-WISE Overlap matrix for", replicate, ":\n")
      print(round(results_rowwise$overlap_matrix, 1))
      
      cat("COLUMN-WISE Overlap matrix for", replicate, ":\n")
      print(round(results_columnwise$overlap_matrix, 1))
      
    } else {
      cat("No files found for replicate", replicate, "\n")
    }
  }
  
  # Check if we have at least 2 replicates for analysis
  if (length(replicate_results_rowwise) < 2) {
    cat("ERROR: Need at least 2 replicates for analysis!\n")
    return(NULL)
  }
  
  # Calculate mean overlap matrices across replicates for BOTH types
  cat("\n=== Calculating MEAN overlap matrices ===\n")
  
  cat("ROW-WISE analysis...\n")
  mean_results_rowwise <- calculate_replicate_statistics(replicate_results_rowwise, use_mean = TRUE)
  
  cat("COLUMN-WISE analysis...\n")
  mean_results_columnwise <- calculate_replicate_statistics(replicate_results_columnwise, use_mean = TRUE)
  
  # Create plots for BOTH overlap types
  cat("\n=== Creating TWO SEPARATE overlap plots ===\n")
  
  # Create ROW-WISE heatmap
  cat("Creating ROW-WISE heatmap...\n")
  p_rowwise <- create_ggplot_heatmap_with_sd(
    mean_results_rowwise, 
    title = "ROW-WISE Peak Overlap Analysis\n(Overlap % = Overlapping Peaks / Total Peaks in ROW)"
  )
  
  # Create COLUMN-WISE heatmap
  cat("Creating COLUMN-WISE heatmap...\n")
  p_columnwise <- create_ggplot_heatmap_with_sd(
    mean_results_columnwise, 
    title = "COLUMN-WISE Peak Overlap Analysis\n(Overlap % = Overlapping Peaks / Total Peaks in COLUMN)"
  )
  
  # Save individual plots
  ggsave("ROW_WISE_overlap_heatmap.pdf", p_rowwise, width = 12, height = 8, dpi = 300)
  ggsave("COLUMN_WISE_overlap_heatmap.pdf", p_columnwise, width = 12, height = 8, dpi = 300)
  
  # Save combined plot with both heatmaps side by side
  pdf("BOTH_overlap_heatmaps_combined.pdf", width = 20, height = 10)
  
  grid.arrange(
    p_rowwise, p_columnwise,
    ncol = 2,
    top = textGrob("Peak Overlap Analysis: Row-wise vs Column-wise", 
                   gp = gpar(fontsize = 18, fontface = "bold"))
  )
  
  dev.off()
  
  # Display results summary
  cat("\n=== REPLICATE ANALYSIS SUMMARY ===\n")
  cat("Number of replicates analyzed:", mean_results_rowwise$n_replicates, "\n")
  cat("File names:", paste(mean_results_rowwise$file_names, collapse = ", "), "\n")
  
  cat("\n=== ROW-WISE RESULTS ===\n")
  cat("Mean overlap percentage matrix (relative to ROW peaks):\n")
  print(round(mean_results_rowwise$central_overlap, 1))
  
  cat("\n=== COLUMN-WISE RESULTS ===\n")
  cat("Mean overlap percentage matrix (relative to COLUMN peaks):\n")
  print(round(mean_results_columnwise$central_overlap, 1))
  
  cat("\n=== STANDARD DEVIATION (ROW-WISE) ===\n")
  cat("Standard deviation matrix:\n")
  print(round(mean_results_rowwise$sd_overlap, 2))
  
  cat("\n=== STANDARD DEVIATION (COLUMN-WISE) ===\n")
  cat("Standard deviation matrix:\n")
  print(round(mean_results_columnwise$sd_overlap, 2))
  
  cat("\n=== INTERPRETATION ===\n")
  cat("ROW-WISE: Shows how much of each ROW's peaks overlap with each COLUMN\n")
  cat("COLUMN-WISE: Shows how much of each COLUMN's peaks overlap with each ROW\n")
  cat("Standard deviation shows reproducibility across R1, R2, R3\n")
  cat("Lower SD = more consistent results across replicates\n")
  
  return(list(
    replicate_results_rowwise = replicate_results_rowwise,
    replicate_results_columnwise = replicate_results_columnwise,
    mean_results_rowwise = mean_results_rowwise,
    mean_results_columnwise = mean_results_columnwise
  ))
}

# Run the TWO HEATMAP ANALYSIS
results_two_heatmaps <- main_replicate_analysis_two_heatmaps() 


################ for seacr peaks ###########################

setwd("/Users/atrayeeray/Desktop/Sid/CUT&Tag_cellno. standardization/Direct CUT&Tag/new_analysis_07082025/Downsampling/M100K_CTCF_no_downsample_M100K_IgG_down/11132025/seacr_peaks/")

# Function to check file status before analysis
check_file_status <- function(file_list) {
  cat("=== FILE STATUS CHECK ===\n")
  
  for (file_path in file_list) {
    if (file.exists(file_path)) {
      file_info <- file.info(file_path)
      if (file_info$size == 0) {
        cat("WARNING:", basename(file_path), "- File exists but is EMPTY\n")
      } else {
        cat("OK:", basename(file_path), "- Size:", file_info$size, "bytes\n")
      }
    } else {
      cat("MISSING:", basename(file_path), "\n")
    }
  }
  
  cat("\n")
}

# Function to read stringent.stringent.bed files with error handling
read_stringent_bed <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(GRanges())
  }
  
  # Check if file is empty
  file_info <- file.info(file_path)
  if (file_info$size == 0) {
    warning(paste("File is empty:", file_path))
    return(GRanges())
  }
  
  # Try to read the file
  tryCatch({
    peaks <- read.table(file_path, header = FALSE, sep = "\t")
    
    # Check if any data was read
    if (nrow(peaks) == 0) {
      warning(paste("File contains no data rows:", file_path))
      return(GRanges())
    }
    
    # Check if we have the expected number of columns (BED format typically has 6 columns)
    if (ncol(peaks) < 3) {
      warning(paste("File has insufficient columns (expected at least 3, got", ncol(peaks), "):", file_path))
      return(GRanges())
    }
    
    # Create GRanges object (BED format: chr, start, end, name, score, strand)
    gr <- GRanges(
      seqnames = peaks$V1,
      ranges = IRanges(start = peaks$V2 + 1, end = peaks$V3),  # Convert to 1-based coordinates
      strand = "*",
      name = ifelse(ncol(peaks) >= 4, peaks$V4, ""),
      score = ifelse(ncol(peaks) >= 5, peaks$V5, 0)
    )
    
    cat("Successfully read", length(gr), "peaks from", basename(file_path), "\n")
    return(gr)
    
  }, error = function(e) {
    warning(paste("Error reading file", file_path, ":", e$message))
    return(GRanges())
  })
}

# Function to calculate overlap percentage (ROW-WISE: relative to row peaks)
calculate_overlap_percentage_rowwise <- function(gr1, gr2, gr1_name, gr2_name) {
  if (length(gr1) == 0 || length(gr2) == 0) {
    return(0)
  }
  
  # Find overlapping peaks
  overlapping_peaks <- GenomicRanges::intersect(gr1, gr2)
  overlap_count <- length(overlapping_peaks)
  
  # Calculate overlap percentage relative to gr1 (ROW) - how much of row peaks overlap with column
  overlap_percentage <- ifelse(length(gr1) > 0, (overlap_count / length(gr1)) * 100, 0)
  
  return(overlap_percentage)
}

# Function to calculate overlap percentage (COLUMN-WISE: relative to column peaks)
calculate_overlap_percentage_columnwise <- function(gr1, gr2, gr1_name, gr2_name) {
  if (length(gr1) == 0 || length(gr2) == 0) {
    return(0)
  }
  
  # Find overlapping peaks
  overlapping_peaks <- GenomicRanges::intersect(gr1, gr2)
  overlap_count <- length(overlapping_peaks)
  
  # Calculate overlap percentage relative to gr2 (COLUMN) - how much of column peaks overlap with row
  overlap_percentage <- ifelse(length(gr2) > 0, (overlap_count / length(gr2)) * 100, 0)
  
  return(overlap_percentage)
}

# Function to extract label from filename
extract_label_from_filename <- function(file_path) {
  # Pattern: M100K_CTCF_R1_based_on_IgG_merged_10000000_stringent.stringent.bed
  # Extract the number before _stringent.stringent.bed
  # Use regex to extract the number
  number_match <- regmatches(file_path, regexpr("_([0-9]+)_stringent\\.stringent\\.bed", file_path))
  
  if (length(number_match) > 0) {
    # Extract just the number part (remove underscores and suffix)
    number_str <- gsub("_([0-9]+)_stringent\\.stringent\\.bed", "\\1", number_match)
    number_val <- as.numeric(number_str)
    
    # Convert to label format
    if (number_val == 10000000) {
      return("Original")  # 10M is the original
    } else if (number_val >= 1000000) {
      return(paste0(number_val / 1000000, "M"))
    } else if (number_val >= 1000) {
      return(paste0(number_val / 1000, "K"))
    } else {
      return(as.character(number_val))
    }
  } else {
    # Fallback to basename if pattern doesn't match
    return(basename(file_path))
  }
}

# Function to create overlap matrix for a specific replicate
create_overlap_matrix <- function(file_list, replicate = "R1", overlap_type = "rowwise") {
  cat("Creating", overlap_type, "overlap matrix for", replicate, "...\n")
  
  # Read all peak files for this replicate
  peak_granges <- list()
  file_names <- c()
  
  for (file_path in file_list) {
    if (grepl(replicate, file_path)) {
      cat("Reading:", basename(file_path), "\n")
      gr <- read_stringent_bed(file_path)
      peak_granges[[length(peak_granges) + 1]] <- gr
      
      # Extract label from filename
      label <- extract_label_from_filename(file_path)
      file_names <- c(file_names, label)
    }
  }
  
  if (length(peak_granges) == 0) {
    stop("No files found for replicate ", replicate)
  }
  
  cat("Found", length(peak_granges), "peak files for", replicate, "\n")
  cat("File labels:", paste(file_names, collapse = ", "), "\n")
  
  # Calculate pairwise overlaps
  n_files <- length(peak_granges)
  overlap_matrix <- matrix(0, nrow = n_files, ncol = n_files)
  
  # Choose overlap calculation function based on type
  overlap_func <- if(overlap_type == "rowwise") {
    calculate_overlap_percentage_rowwise
  } else {
    calculate_overlap_percentage_columnwise
  }
  
  # Fill the matrix
  for (i in 1:n_files) {
    for (j in 1:n_files) {
      if (i == j) {
        overlap_matrix[i, j] <- 100  # Self-overlap = 100%
      } else {
        overlap_matrix[i, j] <- overlap_func(peak_granges[[i]], peak_granges[[j]], 
                                             file_names[i], file_names[j])
      }
    }
  }
  
  # Set row and column names
  rownames(overlap_matrix) <- file_names
  colnames(overlap_matrix) <- file_names
  
  return(list(
    overlap_matrix = overlap_matrix,
    file_names = file_names,
    peak_counts = sapply(peak_granges, length),
    overlap_type = overlap_type
  ))
}

# Function to calculate mean overlap matrices across replicates with standard deviation
calculate_replicate_statistics <- function(replicate_results_list, use_mean = TRUE) {
  cat("Calculating", ifelse(use_mean, "mean", "median"), "overlap matrices across replicates...\n")
  
  # Extract matrices from each replicate
  n_replicates <- length(replicate_results_list)
  n_files <- nrow(replicate_results_list[[1]]$overlap_matrix)
  
  # Initialize array to store values for each quadrant across replicates
  overlap_array <- array(0, dim = c(n_files, n_files, n_replicates))
  
  # Fill array with values from each replicate
  for (i in 1:n_replicates) {
    overlap_array[,,i] <- replicate_results_list[[i]]$overlap_matrix
  }
  
  # Calculate mean or median for each quadrant
  if (use_mean) {
    central_overlap <- apply(overlap_array, c(1,2), mean, na.rm = TRUE)
    stat_name <- "mean"
  } else {
    central_overlap <- apply(overlap_array, c(1,2), median, na.rm = TRUE)
    stat_name <- "median"
  }
  
  # Calculate standard deviation for each quadrant
  sd_overlap <- apply(overlap_array, c(1,2), sd, na.rm = TRUE)
  
  # Set row and column names
  rownames(central_overlap) <- replicate_results_list[[1]]$file_names
  colnames(central_overlap) <- replicate_results_list[[1]]$file_names
  rownames(sd_overlap) <- replicate_results_list[[1]]$file_names
  colnames(sd_overlap) <- replicate_results_list[[1]]$file_names
  
  return(list(
    central_overlap = central_overlap,
    sd_overlap = sd_overlap,
    file_names = replicate_results_list[[1]]$file_names,
    n_replicates = n_replicates,
    stat_name = stat_name,
    overlap_type = replicate_results_list[[1]]$overlap_type
  ))
}

# Function to create ggplot heatmap with standard deviation (LOWER TRIANGULAR)
create_ggplot_heatmap_with_sd <- function(central_results, title = "Peak Overlap Analysis") {
  
  n_rows <- nrow(central_results$central_overlap)
  n_cols <- ncol(central_results$central_overlap)
  
  # Build plotting dataframe for LOWER TRIANGULAR only
  plot_data <- data.frame(
    Row = rep(1:n_rows, each = n_cols),
    Column = rep(1:n_cols, times = n_rows),
    Value = as.vector(central_results$central_overlap),
    SD_Value = as.vector(central_results$sd_overlap),
    RowName = rep(rownames(central_results$central_overlap), each = n_cols),
    ColName = rep(colnames(central_results$central_overlap), times = n_rows),
    stringsAsFactors = FALSE
  )
  
  # Filter to show only LOWER TRIANGULAR (diagonal and below)
  plot_data <- plot_data %>%
    filter(Row >= Column)  # This creates the lower triangular pattern
  
  # Set factor levels to preserve order
  plot_data$RowName <- factor(plot_data$RowName, levels = rownames(central_results$central_overlap))
  plot_data$ColName <- factor(plot_data$ColName, levels = colnames(central_results$central_overlap))
  
  # Build heatmap - CLEAN STYLE with NUMBERS like before
  p <- ggplot(plot_data, aes(x = ColName, y = RowName, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    
    # ADD BACK THE NUMERICAL VALUES - Main overlap percentage
    geom_text(aes(label = sprintf("%.1f", Value)), 
              color = "black", size = 3.5, fontface = "bold") +
    
    # ADD BACK THE STANDARD DEVIATION - below main values
    geom_text(aes(label = sprintf("±%.2f", SD_Value)), 
              color = "darkblue", size = 2.5, fontface = "italic",
              position = position_nudge(y = -0.15)) +
    
    # Use diverging color palette like the image (red-white-blue)
    scale_fill_gradient2(
      low = "#313695",      # Dark blue for low values
      mid = "white",        # White for middle values  
      high = "#a50026",     # Dark red for high values
      midpoint = 50,        # 50% overlap as neutral point
      limits = c(0, 100),
      na.value = "grey90"
    ) +
    
    labs(
      title = title,
      subtitle = sprintf("Lower triangular matrix showing Mean ± SD across %d replicates (R1, R2, R3)", central_results$n_replicates),
      fill = "Overlap %",
      x = "Column",
      y = "Row"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(hjust = 1, size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right",
      plot.margin = margin(20, 20, 20, 20),
      panel.grid = element_blank()  # Remove grid lines for cleaner look
    ) +
    
    coord_fixed()
  
  return(p)
}

# Function to create summary statistics table with Mean ± SD
create_summary_table <- function(central_results) {
  
  # Create a comprehensive table showing Mean ± SD for each cell
  n_files <- length(central_results$file_names)
  
  # Build table data for lower triangular matrix
  table_data <- data.frame(
    Row = rep(1:n_files, each = n_files),
    Column = rep(1:n_files, times = n_files),
    RowName = rep(central_results$file_names, each = n_files),
    ColName = rep(central_results$file_names, times = n_files),
    Mean_Value = as.vector(central_results$central_overlap),
    SD_Value = as.vector(central_results$sd_overlap),
    stringsAsFactors = FALSE
  )
  
  # Filter to lower triangular only
  table_data <- table_data %>%
    filter(Row >= Column) %>%
    mutate(
      `Comparison` = paste(RowName, "→", ColName),
      `Mean ± SD` = sprintf("%.1f ± %.2f", Mean_Value, SD_Value)
    ) %>%
    select(`Comparison`, `Mean ± SD`)
  
  # Create table grob
  table_grob <- tableGrob(table_data, 
                          rows = NULL,
                          theme = ttheme_minimal(
                            base_size = 8,
                            core = list(
                              fg_params = list(hjust = 0, x = 0.1),
                              bg_params = list(fill = c("white", "lightgray"))
                            ),
                            colhead = list(
                              fg_params = list(fontface = "bold"),
                              bg_params = list(fill = "lightblue")
                            )
                          ))
  
  # Add title to table
  table_title <- textGrob("Overlap Statistics (Mean ± SD)", 
                          gp = gpar(fontsize = 12, fontface = "bold"),
                          just = "center")
  
  # Combine table elements
  table_with_title <- arrangeGrob(table_title, table_grob, 
                                  heights = c(0.1, 0.9),
                                  ncol = 1)
  
  return(table_with_title)
}

# Main execution function for REPLICATE ANALYSIS with TWO SEPARATE HEATMAPS
main_replicate_analysis_two_heatmaps <- function() {
  cat("Starting REPLICATE ANALYSIS with TWO SEPARATE HEATMAPS (Row-wise vs Column-wise)...\n")
  
  # Define your file paths for all replicates - UPDATED for stringent.stringent.bed files
  file_list <- c(
    "M100K_CTCF_R1_based_on_IgG_merged_original_stringent.stringent.bed",  # Original (was 10M)
    "M100K_CTCF_R1_based_on_IgG_merged_5000000_stringent.stringent.bed",   # 5M
    "M100K_CTCF_R1_based_on_IgG_merged_2500000_stringent.stringent.bed",   # 2.5M
    "M100K_CTCF_R1_based_on_IgG_merged_1000000_stringent.stringent.bed",   # 1M
    "M100K_CTCF_R1_based_on_IgG_merged_500000_stringent.stringent.bed",    # 500K
    "M100K_CTCF_R1_based_on_IgG_merged_100000_stringent.stringent.bed",    # 100K
    "M100K_CTCF_R1_based_on_IgG_merged_50000_stringent.stringent.bed",     # 50K
    "M100K_CTCF_R1_based_on_IgG_merged_10000_stringent.stringent.bed",     # 10K
    
    "M100K_CTCF_R2_based_on_IgG_merged_original_stringent.stringent.bed",  # Original (was 10M)
    "M100K_CTCF_R2_based_on_IgG_merged_5000000_stringent.stringent.bed",   # 5M
    "M100K_CTCF_R2_based_on_IgG_merged_2500000_stringent.stringent.bed",   # 2.5M
    "M100K_CTCF_R2_based_on_IgG_merged_1000000_stringent.stringent.bed",   # 1M
    "M100K_CTCF_R2_based_on_IgG_merged_500000_stringent.stringent.bed",    # 500K
    "M100K_CTCF_R2_based_on_IgG_merged_100000_stringent.stringent.bed",    # 100K
    "M100K_CTCF_R2_based_on_IgG_merged_50000_stringent.stringent.bed",     # 50K
    "M100K_CTCF_R2_based_on_IgG_merged_10000_stringent.stringent.bed",     # 10K
    
    "M100K_CTCF_R3_based_on_IgG_merged_original_stringent.stringent.bed",  # Original (was 10M)
    "M100K_CTCF_R3_based_on_IgG_merged_5000000_stringent.stringent.bed",   # 5M
    "M100K_CTCF_R3_based_on_IgG_merged_2500000_stringent.stringent.bed",   # 2.5M
    "M100K_CTCF_R3_based_on_IgG_merged_1000000_stringent.stringent.bed",   # 1M
    "M100K_CTCF_R3_based_on_IgG_merged_500000_stringent.stringent.bed",    # 500K
    "M100K_CTCF_R3_based_on_IgG_merged_100000_stringent.stringent.bed",    # 100K
    "M100K_CTCF_R3_based_on_IgG_merged_50000_stringent.stringent.bed",     # 50K
    "M100K_CTCF_R3_based_on_IgG_merged_10000_stringent.stringent.bed"      # 10K
  )
  
  # Check file status
  check_file_status(file_list)
  
  # Check which files exist
  existing_files <- file_list[file.exists(file_list)]
  cat("Found", length(existing_files), "existing files\n")
  
  # Analyze each replicate for BOTH overlap types
  replicates <- c("R1", "R2", "R3")
  replicate_results_rowwise <- list()
  replicate_results_columnwise <- list()
  
  for (replicate in replicates) {
    cat("\n=== Analyzing replicate", replicate, "===\n")
    
    # Check if files exist for this replicate
    replicate_files <- existing_files[grepl(replicate, existing_files)]
    
    if (length(replicate_files) > 0) {
      # Perform ROW-WISE overlap analysis for this replicate
      cat("Creating ROW-WISE overlap matrix...\n")
      results_rowwise <- create_overlap_matrix(replicate_files, replicate, "rowwise")
      replicate_results_rowwise[[replicate]] <- results_rowwise
      
      # Perform COLUMN-WISE overlap analysis for this replicate
      cat("Creating COLUMN-WISE overlap matrix...\n")
      results_columnwise <- create_overlap_matrix(replicate_files, replicate, "columnwise")
      replicate_results_columnwise[[replicate]] <- results_columnwise
      
      cat("Completed analysis for", replicate, "\n")
      
      # Print individual replicate results
      cat("ROW-WISE Overlap matrix for", replicate, ":\n")
      print(round(results_rowwise$overlap_matrix, 1))
      
      cat("COLUMN-WISE Overlap matrix for", replicate, ":\n")
      print(round(results_columnwise$overlap_matrix, 1))
      
    } else {
      cat("No files found for replicate", replicate, "\n")
    }
  }
  
  # Check if we have at least 2 replicates for analysis
  if (length(replicate_results_rowwise) < 2) {
    cat("ERROR: Need at least 2 replicates for analysis!\n")
    return(NULL)
  }
  
  # Calculate mean overlap matrices across replicates for BOTH types
  cat("\n=== Calculating MEAN overlap matrices ===\n")
  
  cat("ROW-WISE analysis...\n")
  mean_results_rowwise <- calculate_replicate_statistics(replicate_results_rowwise, use_mean = TRUE)
  
  cat("COLUMN-WISE analysis...\n")
  mean_results_columnwise <- calculate_replicate_statistics(replicate_results_columnwise, use_mean = TRUE)
  
  # Create plots for BOTH overlap types
  cat("\n=== Creating TWO SEPARATE overlap plots ===\n")
  
  # Create ROW-WISE heatmap
  cat("Creating ROW-WISE heatmap...\n")
  p_rowwise <- create_ggplot_heatmap_with_sd(
    mean_results_rowwise, 
    title = "ROW-WISE Peak Overlap Analysis\n(Overlap % = Overlapping Peaks / Total Peaks in ROW)"
  )
  
  # Create COLUMN-WISE heatmap
  cat("Creating COLUMN-WISE heatmap...\n")
  p_columnwise <- create_ggplot_heatmap_with_sd(
    mean_results_columnwise, 
    title = "COLUMN-WISE Peak Overlap Analysis\n(Overlap % = Overlapping Peaks / Total Peaks in COLUMN)"
  )
  
  # Save individual plots
  ggsave("ROW_WISE_overlap_heatmap.pdf", p_rowwise, width = 12, height = 8, dpi = 300)
  ggsave("COLUMN_WISE_overlap_heatmap.pdf", p_columnwise, width = 12, height = 8, dpi = 300)
  
  # Save combined plot with both heatmaps side by side
  pdf("BOTH_overlap_heatmaps_combined.pdf", width = 20, height = 10)
  
  grid.arrange(
    p_rowwise, p_columnwise,
    ncol = 2,
    top = textGrob("Peak Overlap Analysis: Row-wise vs Column-wise", 
                   gp = gpar(fontsize = 18, fontface = "bold"))
  )
  
  dev.off()
  
  # Display results summary
  cat("\n=== REPLICATE ANALYSIS SUMMARY ===\n")
  cat("Number of replicates analyzed:", mean_results_rowwise$n_replicates, "\n")
  cat("File names:", paste(mean_results_rowwise$file_names, collapse = ", "), "\n")
  
  cat("\n=== ROW-WISE RESULTS ===\n")
  cat("Mean overlap percentage matrix (relative to ROW peaks):\n")
  print(round(mean_results_rowwise$central_overlap, 1))
  
  cat("\n=== COLUMN-WISE RESULTS ===\n")
  cat("Mean overlap percentage matrix (relative to COLUMN peaks):\n")
  print(round(mean_results_columnwise$central_overlap, 1))
  
  cat("\n=== STANDARD DEVIATION (ROW-WISE) ===\n")
  cat("Standard deviation matrix:\n")
  print(round(mean_results_rowwise$sd_overlap, 2))
  
  cat("\n=== STANDARD DEVIATION (COLUMN-WISE) ===\n")
  cat("Standard deviation matrix:\n")
  print(round(mean_results_columnwise$sd_overlap, 2))
  
  cat("\n=== INTERPRETATION ===\n")
  cat("ROW-WISE: Shows how much of each ROW's peaks overlap with each COLUMN\n")
  cat("COLUMN-WISE: Shows how much of each COLUMN's peaks overlap with each ROW\n")
  cat("Standard deviation shows reproducibility across R1, R2, R3\n")
  cat("Lower SD = more consistent results across replicates\n")
  
  return(list(
    replicate_results_rowwise = replicate_results_rowwise,
    replicate_results_columnwise = replicate_results_columnwise,
    mean_results_rowwise = mean_results_rowwise,
    mean_results_columnwise = mean_results_columnwise
  ))
}

# Run the TWO HEATMAP ANALYSIS
results_two_heatmaps <- main_replicate_analysis_two_heatmaps()

