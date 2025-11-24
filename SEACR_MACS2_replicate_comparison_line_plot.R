library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
library(outliers)
library(stats)



setwd("/Users/atrayeeray/Desktop/Sid/CUT&Tag_cellno. standardization/Direct CUT&Tag/new_analysis_07082025/11072025/")


Peaksummary <- read.csv("aligned_reads_based_on_merged_IgG.csv")
Peaksummary <- Peaksummary[, !grepl("^X", colnames(Peaksummary))]
Peaksummary <- Peaksummary[,c(1:10)]

cat("Total.Reads.in.Million class:", class(Peaksummary$Total.Reads.in.Million), "\n")
#Total.Reads.in.Million class: numeric
cat("Mapped_Reads.in.Million class:", class(Peaksummary$Mapped_Reads.in.Million), "\n")
#Mapped_Reads.in.Million: numeric
cat("Corresponding.IgG.reads class:", class(Peaksummary$Corresponding.IgG.reads), "\n")
#Corresponding.IgG.reads class: integer 
cat("MACS2_Default_Peak:", class(Peaksummary$Corresponding.IgG.reads), "\n")
cat("MACS2_Large_Peak:", class(Peaksummary$Corresponding.IgG.reads), "\n")


# Step 1: Filter samples with less than 0.2M mapped reads
cat("Step 1: Filtering samples with less than 0.2M mapped reads...\n")
samples_before <- nrow(Peaksummary)
Peaksummary <- Peaksummary %>% 
  filter(Mapped_Reads.in.Million >= 0.2)  # Keep only samples with >= 0.2M reads (0.2 million)
samples_after <- nrow(Peaksummary)
cat("Samples before filtering:", samples_before, "\n")
cat("Samples after filtering:", samples_after, "\n")
cat("Samples removed:", samples_before - samples_after, "\n")

# Step 2: Function to detect outliers using Grubbs test
detect_outliers_grubbs <- function(data, sample_name, peak_type) {
  replicates <- data$Replicate
  
  # Get the correct peak column name
  if(peak_type == "Default") {
    peak_counts <- data$MACS2_Default_Peak
  } else if(peak_type == "Large") {
    peak_counts <- data$MACS2_Large_Peak
  } else if(peak_type == "SEACR_Stringent") {
    peak_counts <- data$SEACR_Stringent_Peak
  } else if(peak_type == "SEACR_Relaxed") {
    peak_counts <- data$SEACR_Relaxed_Peak
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

# Step 3: Function to identify and mark outliers for all peak types
identify_all_outliers <- function(data) {
  # Initialize outlier columns
  data$is_outlier_default <- FALSE
  data$is_outlier_large <- FALSE
  data$is_outlier_seacr_stringent <- FALSE
  data$is_outlier_seacr_relaxed <- FALSE
  
  # Get unique samples
  unique_samples <- unique(data$Sample)
  outliers_list <- list()
  
  for(sample in unique_samples) {
    sample_data <- data[data$Sample == sample, ]
    
    if(nrow(sample_data) == 3) {
      # Test all four peak calling methods
      peak_types <- c("Default", "Large", "SEACR_Stringent", "SEACR_Relaxed")
      
      for(peak_type in peak_types) {
        outlier_result <- detect_outliers_grubbs(sample_data, sample, peak_type)
        
        if(!is.null(outlier_result)) {
          outliers_list[[length(outliers_list) + 1]] <- outlier_result
          
          # Mark outliers in the appropriate column
          outlier_rows <- data$Sample == sample & 
            data$Replicate == outlier_result$outlier_replicate
          
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
  
  return(list(data = data, outliers = outliers_list))
}

# Step 4: Function to filter data (remove 5K and 25K, add species and mark info)
filter_and_prepare_data <- function(data) {
  data_filtered <- data %>%
    # Remove 5K and 25K cell numbers
    filter(!Cell_numbers %in% c("5K", "25K")) %>%
    # Add species information
    mutate(
      Species = ifelse(grepl("^M", Sample), "Mouse", "Rat"),
      Mark = ifelse(grepl("CTCF", Sample), "CTCF", "K27me3")
    ) %>%
    # Convert cell numbers to numeric for proper ordering
    mutate(
      Cell_numbers_numeric = as.numeric(gsub("K", "000", Cell_numbers))
    )
  
  return(data_filtered)
}

# Step 5: Function to create summary statistics for plotting
create_summary_stats <- function(data, peak_type, outlier_col) {
  # Get the correct peak column name
  if(peak_type == "Default") {
    peak_col <- "MACS2_Default_Peak"
  } else if(peak_type == "Large") {
    peak_col <- "MACS2_Large_Peak"
  } else if(peak_type == "SEACR_Stringent") {
    peak_col <- "SEACR_Stringent_Peak"
  } else if(peak_type == "SEACR_Relaxed") {
    peak_col <- "SEACR_Relaxed_Peak"
  }
  
  # Create summary statistics without outliers
  summary_data <- data[!data[[outlier_col]], ] %>%
    group_by(Cell_numbers_numeric, Species, Mark) %>%
    summarise(
      mean_peaks = mean(get(peak_col), na.rm = TRUE),
      sd_peaks = sd(get(peak_col), na.rm = TRUE),
      n_replicates = n(),
      .groups = 'drop'
    )
  
  return(summary_data)
}

#### this is for the line plot

# Step 6: Function to create combined plot for a specific species and mark
create_combined_plot <- function(data, species, mark) {
  # Filter data for specific species and mark
  plot_data <- data %>%
    filter(Species == species, Mark == mark)
  
  # Add small constant to zero values to avoid log(0) issues
  plot_data <- plot_data %>%
    mutate(mean_peaks_adj = ifelse(mean_peaks == 0, 0.1, mean_peaks))
  
  # Create factor for cell numbers with proper ordering (1K, 10K, 50K, 100K)
  plot_data <- plot_data %>%
    mutate(
      Cell_numbers_factor = factor(
        paste0(Cell_numbers_numeric/1000, "K"),
        levels = c("1K", "10K", "50K", "100K"),
        ordered = TRUE
      )
    )
  
  # Create plot with all four peak types
  p <- ggplot(plot_data, aes(x = Cell_numbers_factor, y = mean_peaks_adj, color = Peak_Type, shape = Peak_Type)) +
    geom_point(size = 3) +
    geom_line(aes(group = Peak_Type), alpha = 0.7, linewidth = 1) +
    geom_errorbar(aes(ymin = mean_peaks - sd_peaks, ymax = mean_peaks + sd_peaks), 
                  width = 0.1, alpha = 0.5) +
    # Use discrete scale for X-axis (factors) and log scale for Y-axis
    scale_x_discrete(name = "Cell Number") +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(
      values = c("MACS2_Default" = "#1f77b4", "MACS2_Large" = "#ff7f0e", 
                 "SEACR_Stringent" = "#2ca02c", "SEACR_Relaxed" = "#d62728"),
      labels = c("MACS2_Default" = "MACS2 Default", "MACS2_Large" = "MACS2 Large",
                 "SEACR_Stringent" = "SEACR Stringent", "SEACR_Relaxed" = "SEACR Relaxed")
    ) +
    scale_shape_manual(
      values = c("MACS2_Default" = 16, "MACS2_Large" = 17, 
                 "SEACR_Stringent" = 15, "SEACR_Relaxed" = 18),
      labels = c("MACS2_Default" = "MACS2 Default", "MACS2_Large" = "MACS2 Large",
                 "SEACR_Stringent" = "SEACR Stringent", "SEACR_Relaxed" = "SEACR Relaxed")
    ) +
    labs(
      title = paste(species, mark, "Peak Counts by Method"),
      subtitle = "Cell Number vs Peak Count (Outliers Removed, Single Replicates Dropped)",
      x = "Cell Number",
      y = "Mean Peak Count (log scale)",
      color = "Peak Calling Method",
      shape = "Peak Calling Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  
  return(p)
}

# Step 7: Main execution
main_analysis <- function() {
  cat("Starting peak analysis...\n")
  
  # Step 1: Identify outliers (after read filtering)
  cat("\nStep 2: Identifying outliers using Grubbs test...\n")
  outlier_results <- identify_all_outliers(Peaksummary)
  Peaksummary_with_outliers <- outlier_results$data
  
  # Count outliers for each peak type
  outlier_counts <- list(
    Default = sum(Peaksummary_with_outliers$is_outlier_default),
    Large = sum(Peaksummary_with_outliers$is_outlier_large),
    SEACR_Stringent = sum(Peaksummary_with_outliers$is_outlier_seacr_stringent),
    SEACR_Relaxed = sum(Peaksummary_with_outliers$is_outlier_seacr_relaxed)
  )
  
  cat("Outliers detected:\n")
  for(peak_type in names(outlier_counts)) {
    cat("  ", peak_type, ":", outlier_counts[[peak_type]], "\n")
  }
  
  # Step 2: Filter and prepare data
  cat("\nStep 3: Filtering and preparing data...\n")
  Peaksummary_filtered <- filter_and_prepare_data(Peaksummary_with_outliers)
  
  cat("Data points after filtering:", nrow(Peaksummary_filtered), "\n")
  cat("Unique cell numbers:", unique(Peaksummary_filtered$Cell_numbers), "\n")
  
  # DIAGNOSTIC: Check what peak columns exist and their content
  cat("\nDiagnostic: Checking peak columns...\n")
  cat("Available columns containing 'Peak':\n")
  peak_cols <- grep("Peak", names(Peaksummary_filtered), value = TRUE)
  print(peak_cols)
  
  if("MACS2_Default_Peak" %in% names(Peaksummary_filtered)) {
    cat("\nMACS2_Default_Peak column:\n")
    cat("  Class:", class(Peaksummary_filtered$MACS2_Default_Peak), "\n")
    cat("  Non-NA count:", sum(!is.na(Peaksummary_filtered$MACS2_Default_Peak)), "\n")
    cat("  Sample values:", head(as.character(Peaksummary_filtered$MACS2_Default_Peak[!is.na(Peaksummary_filtered$MACS2_Default_Peak)]), 5), "\n")
  } else {
    cat("\nWARNING: MACS2_Default_Peak column NOT FOUND!\n")
  }
  
  if("MACS2_Large_Peak" %in% names(Peaksummary_filtered)) {
    cat("\nMACS2_Large_Peak column:\n")
    cat("  Class:", class(Peaksummary_filtered$MACS2_Large_Peak), "\n")
    cat("  Non-NA count:", sum(!is.na(Peaksummary_filtered$MACS2_Large_Peak)), "\n")
    cat("  Sample values:", head(as.character(Peaksummary_filtered$MACS2_Large_Peak[!is.na(Peaksummary_filtered$MACS2_Large_Peak)]), 5), "\n")
  } else {
    cat("\nWARNING: MACS2_Large_Peak column NOT FOUND!\n")
  }
  
  # Step 3: Create combined data for plotting
  cat("\nStep 4: Creating combined data for plotting...\n")
  peak_types <- c("Default", "Large", "SEACR_Stringent", "SEACR_Relaxed")
  
  combined_data <- data.frame()
  
  for(peak_type in peak_types) {
    outlier_col <- paste0("is_outlier_", tolower(gsub("SEACR_", "seacr_", peak_type)))
    summary_data <- create_summary_stats(Peaksummary_filtered, peak_type, outlier_col)
    
    # FIXED: Use consistent naming for peak types
    if(peak_type == "Default") {
      summary_data$Peak_Type <- "MACS2_Default"
    } else if(peak_type == "Large") {
      summary_data$Peak_Type <- "MACS2_Large"
    } else if(peak_type == "SEACR_Stringent") {
      summary_data$Peak_Type <- "SEACR_Stringent"
    } else if(peak_type == "SEACR_Relaxed") {
      summary_data$Peak_Type <- "SEACR_Relaxed"
    }
    
    combined_data <- rbind(combined_data, summary_data)
  }
  
  # Debug: Check what peak types are actually in the data
  cat("Peak types in combined data:", unique(combined_data$Peak_Type), "\n")
  cat("Number of rows per peak type:\n")
  if(nrow(combined_data) > 0) {
    print(table(combined_data$Peak_Type))
  }
  
  # Step 4: Create the four plots
  cat("\nStep 5: Creating plots...\n")
  species_marks <- list(
    c("Mouse", "CTCF"),
    c("Mouse", "K27me3"),
    c("Rat", "CTCF"),
    c("Rat", "K27me3")
  )
  
  all_plots <- list()
  
  for(i in 1:length(species_marks)) {
    species <- species_marks[[i]][1]
    mark <- species_marks[[i]][2]
    
    plot_name <- paste0(tolower(species), "_", tolower(mark))
    
    # Debug: Check data availability for this plot
    plot_subset <- combined_data %>%
      filter(Species == species, Mark == mark)
    cat("  Data for", plot_name, "- Rows:", nrow(plot_subset), "\n")
    if(nrow(plot_subset) > 0) {
      cat("    Peak types available:", unique(plot_subset$Peak_Type), "\n")
    }
    
    p <- create_combined_plot(combined_data, species, mark)
    all_plots[[plot_name]] <- p
    
    cat("  Created plot:", plot_name, "\n")
  }
  
  # Step 5: Save all plots
  cat("\nStep 6: Saving plots...\n")
  for(plot_name in names(all_plots)) {
    filename <- paste0(plot_name, "_combined_peaks.pdf")
    ggsave(filename, all_plots[[plot_name]], width = 12, height = 8, dpi = 300)
    cat("  Saved:", filename, "\n")
  }
  
  cat("\nAnalysis complete! 4 plots have been saved.\n")
  
  # Return processed data and plots for further use
  return(list(
    data = Peaksummary_filtered,
    combined_data = combined_data,
    outliers = outlier_results$outliers,
    plots = all_plots
  ))
}

# Run the analysis
results <- main_analysis()


##### this is for the box and line plot

# Step 5: Function to prepare raw data for box plots (INCLUDING ALL DATA - NO OUTLIER EXCLUSION)
create_raw_data_for_plots <- function(data, peak_type) {
  # Get the correct peak column name
  if(peak_type == "Default") {
    peak_col <- "MACS2_Default_Peak"
  } else if(peak_type == "Large") {
    peak_col <- "MACS2_Large_Peak"
  } else if(peak_type == "SEACR_Stringent") {
    peak_col <- "SEACR_Stringent_Peak"
  } else if(peak_type == "SEACR_Relaxed") {
    peak_col <- "SEACR_Relaxed_Peak"
  }
  
  # FIX: Convert peak column to numeric
  data[[peak_col]] <- as.numeric(as.character(data[[peak_col]]))
  
  # Remove rows where peak values are NA after conversion
  data <- data[!is.na(data[[peak_col]]), ]
  
  if(nrow(data) == 0) {
    return(data.frame())
  }
  
  # Prepare raw data with Peak_Type column
  raw_data <- data %>%
    select(Cell_numbers_numeric, Species, Mark, all_of(peak_col), Sample, Replicate) %>%
    rename(peak_value = !!peak_col) %>%
    mutate(Peak_Type = case_when(
      peak_type == "Default" ~ "MACS2_Default",
      peak_type == "Large" ~ "MACS2_Large",
      peak_type == "SEACR_Stringent" ~ "SEACR_Stringent",
      peak_type == "SEACR_Relaxed" ~ "SEACR_Relaxed"
    ))
  
  return(raw_data)
}

# Step 6: Function to create plot with box plots only (all methods side by side)
create_combined_plot <- function(raw_data, species, mark) {
  # Filter data for specific species and mark
  plot_data <- raw_data %>%
    filter(Species == species, Mark == mark) %>%
    # Filter out groups with only 1 replicate (for box plots)
    group_by(Cell_numbers_numeric, Peak_Type) %>%
    filter(n() >= 2) %>%
    ungroup()
  
  # Add small constant to zero values to avoid log(0) issues
  plot_data <- plot_data %>%
    mutate(peak_value_adj = ifelse(peak_value == 0, 0.1, peak_value))
  
  # Create factor for cell numbers with proper ordering (1K, 10K, 50K, 100K)
  plot_data <- plot_data %>%
    mutate(
      Cell_numbers_factor = factor(
        paste0(Cell_numbers_numeric/1000, "K"),
        levels = c("1K", "10K", "50K", "100K"),
        ordered = TRUE
      )
    )
  
  # Ensure Replicate is a factor for proper shape mapping
  plot_data <- plot_data %>%
    mutate(Replicate = as.factor(Replicate))
  
  # Get unique replicate values for shape mapping
  unique_replicates <- sort(unique(as.character(plot_data$Replicate)))
  n_replicates <- length(unique_replicates)
  # Use standard shapes: circle (16), triangle (17), diamond (18), square (15), plus (3)
  shape_values <- c(16, 17, 18, 15, 3, 4, 8)[1:min(n_replicates, 7)]
  names(shape_values) <- unique_replicates
  
  # Calculate median values for trend lines (connect medians across cell numbers for each method)
  trend_data <- plot_data %>%
    group_by(Cell_numbers_factor, Peak_Type) %>%
    summarise(
      median_peaks = median(peak_value_adj, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(Cell_numbers_factor)
  
  # Create plot with box plots, individual points, and trend lines
  p <- ggplot(plot_data, aes(x = Cell_numbers_factor, y = peak_value_adj, 
                             color = Peak_Type, fill = Peak_Type, group = interaction(Cell_numbers_factor, Peak_Type))) +
    # Box plots
    geom_boxplot(
      alpha = 0.4,
      outlier.size = 0,  # Hide default outliers, we'll show all points manually
      width = 0.15,
      position = position_dodge(width = 0.2)
    ) +
    # Trend lines connecting medians across cell numbers for each method
    geom_line(
      data = trend_data,
      aes(x = Cell_numbers_factor, y = median_peaks, group = Peak_Type, color = Peak_Type),
      linewidth = 1.2,
      alpha = 0.8,
      linetype = "solid"
    ) +
    # Individual replicate points
    geom_point(
      aes(shape = Replicate),
      size = 2.5,
      alpha = 0.8,
      position = position_jitterdodge(
        jitter.width = 0.05,
        dodge.width = 0.2,
        seed = 123
      )
    ) +
    # Use discrete scale for X-axis (factors) and log scale for Y-axis
    scale_x_discrete(name = "Cell Number") +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(
      values = c("MACS2_Default" = "#1f77b4", "MACS2_Large" = "#ff7f0e", 
                 "SEACR_Stringent" = "#2ca02c", "SEACR_Relaxed" = "#d62728"),
      labels = c("MACS2_Default" = "MACS2 Default", "MACS2_Large" = "MACS2 Large",
                 "SEACR_Stringent" = "SEACR Stringent", "SEACR_Relaxed" = "SEACR Relaxed")
    ) +
    scale_fill_manual(
      values = c("MACS2_Default" = "#1f77b4", "MACS2_Large" = "#ff7f0e", 
                 "SEACR_Stringent" = "#2ca02c", "SEACR_Relaxed" = "#d62728"),
      labels = c("MACS2_Default" = "MACS2 Default", "MACS2_Large" = "MACS2 Large",
                 "SEACR_Stringent" = "SEACR Stringent", "SEACR_Relaxed" = "SEACR Relaxed")
    ) +
    scale_shape_manual(
      values = shape_values
    ) +
    labs(
      title = paste(species, mark, "Peak Counts by Method"),
      subtitle = "Boxplot by cell number, all methods side by side",
      x = "Cell Number",
      y = "Number of Peaks",
      color = "Peak Calling Method",
      fill = "Peak Calling Method",
      shape = "Replicate"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "horizontal"
    )
  
  return(p)
}

# Step 7: Main execution
main_analysis <- function() {
  cat("Starting peak analysis...\n")
  
  # Step 1: Identify outliers (for information only, not for exclusion)
  cat("\nStep 2: Identifying outliers using Grubbs test...\n")
  outlier_results <- identify_all_outliers(Peaksummary)
  Peaksummary_with_outliers <- outlier_results$data
  
  # Count outliers for each peak type
  outlier_counts <- list(
    Default = sum(Peaksummary_with_outliers$is_outlier_default),
    Large = sum(Peaksummary_with_outliers$is_outlier_large),
    SEACR_Stringent = sum(Peaksummary_with_outliers$is_outlier_seacr_stringent),
    SEACR_Relaxed = sum(Peaksummary_with_outliers$is_outlier_seacr_relaxed)
  )
  
  cat("Outliers detected (not excluded from plots):\n")
  for(peak_type in names(outlier_counts)) {
    cat("  ", peak_type, ":", outlier_counts[[peak_type]], "\n")
  }
  
  # Step 2: Filter and prepare data
  cat("\nStep 3: Filtering and preparing data...\n")
  Peaksummary_filtered <- filter_and_prepare_data(Peaksummary_with_outliers)
  
  cat("Data points after filtering:", nrow(Peaksummary_filtered), "\n")
  cat("Unique cell numbers:", unique(Peaksummary_filtered$Cell_numbers), "\n")
  
  # Step 3: Create raw data for box plots (ALL DATA INCLUDED)
  cat("\nStep 4: Creating raw data for box plots (all data included)...\n")
  peak_types <- c("Default", "Large", "SEACR_Stringent", "SEACR_Relaxed")
  
  combined_raw <- data.frame()
  
  for(peak_type in peak_types) {
    # Create raw data for box plots (no outlier exclusion)
    raw_data <- create_raw_data_for_plots(Peaksummary_filtered, peak_type)
    combined_raw <- rbind(combined_raw, raw_data)
  }
  
  # Debug: Check what peak types are actually in the data
  cat("Peak types in raw data:", unique(combined_raw$Peak_Type), "\n")
  cat("Number of rows per peak type:\n")
  if(nrow(combined_raw) > 0) {
    print(table(combined_raw$Peak_Type))
  }
  
  # Step 4: Create the four plots
  cat("\nStep 5: Creating box plots...\n")
  species_marks <- list(
    c("Mouse", "CTCF"),
    c("Mouse", "K27me3"),
    c("Rat", "CTCF"),
    c("Rat", "K27me3")
  )
  
  all_plots <- list()
  
  for(i in 1:length(species_marks)) {
    species <- species_marks[[i]][1]
    mark <- species_marks[[i]][2]
    
    plot_name <- paste0(tolower(species), "_", tolower(mark))
    
    p <- create_combined_plot(combined_raw, species, mark)
    all_plots[[plot_name]] <- p
    
    cat("  Created plot:", plot_name, "\n")
  }
  
  # Step 5: Save all plots
  cat("\nStep 6: Saving plots...\n")
  for(plot_name in names(all_plots)) {
    filename <- paste0(plot_name, "_boxplots.pdf")
    ggsave(filename, all_plots[[plot_name]], width = 12, height = 8, dpi = 300)
    cat("  Saved:", filename, "\n")
  }
  
  cat("\nAnalysis complete! 4 plots have been saved.\n")
  
  # Return processed data and plots for further use
  return(list(
    data = Peaksummary_filtered,
    raw_data = combined_raw,
    outliers = outlier_results$outliers,
    plots = all_plots
  ))
}

# Run the analysis
results <- main_analysis()

########## this is for only box plot

# Step 6: Function to create plot with box plots only (all methods side by side)

create_combined_plot <- function(raw_data, species, mark) {
  # Filter data for specific species and mark
  plot_data <- raw_data %>%
    filter(Species == species, Mark == mark) %>%
    # Filter out groups with only 1 replicate (for box plots)
    group_by(Cell_numbers_numeric, Peak_Type) %>%
    filter(n() >= 2) %>%
    ungroup()
  
  # Add small constant to zero values to avoid log(0) issues
  plot_data <- plot_data %>%
    mutate(peak_value_adj = ifelse(peak_value == 0, 0.1, peak_value))
  
  # Create factor for cell numbers with proper ordering (1K, 10K, 50K, 100K)
  plot_data <- plot_data %>%
    mutate(
      Cell_numbers_factor = factor(
        paste0(Cell_numbers_numeric/1000, "K"),
        levels = c("1K", "10K", "50K", "100K"),
        ordered = TRUE
      )
    )
  
  # Ensure Replicate is a factor for proper shape mapping
  plot_data <- plot_data %>%
    mutate(Replicate = as.factor(Replicate))
  
  # Get unique replicate values for shape mapping
  unique_replicates <- sort(unique(as.character(plot_data$Replicate)))
  n_replicates <- length(unique_replicates)
  # Use standard shapes: circle (16), triangle (17), diamond (18), square (15), plus (3)
  shape_values <- c(16, 17, 18, 15, 3, 4, 8)[1:min(n_replicates, 7)]
  names(shape_values) <- unique_replicates
  
  # Create plot with box plots and individual points
  p <- ggplot(plot_data, aes(x = Cell_numbers_factor, y = peak_value_adj, 
                             color = Peak_Type, fill = Peak_Type, group = interaction(Cell_numbers_factor, Peak_Type))) +
    # Box plots - INCREASED WIDTH from 0.15 to 0.4
    geom_boxplot(
      alpha = 0.4,
      outlier.size = 0,  # Hide default outliers, we'll show all points manually
      width = 0.4,  # Increased from 0.15 to make boxes wider
      position = position_dodge(width = 0.5)  # Increased from 0.2 to accommodate wider boxes
    ) +
    # Individual replicate points
    geom_point(
      aes(shape = Replicate),
      size = 2.5,
      alpha = 0.8,
      position = position_jitterdodge(
        jitter.width = 0.05,
        dodge.width = 0.5,  # Increased from 0.2 to match box plot dodge width
        seed = 123
      )
    ) +
    # Use discrete scale for X-axis (factors) and log scale for Y-axis
    scale_x_discrete(name = "Cell Number") +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(
      values = c("MACS2_Default" = "#1f77b4", "MACS2_Large" = "#ff7f0e", 
                 "SEACR_Stringent" = "#2ca02c", "SEACR_Relaxed" = "#d62728"),
      labels = c("MACS2_Default" = "MACS2 Default", "MACS2_Large" = "MACS2 Large",
                 "SEACR_Stringent" = "SEACR Stringent", "SEACR_Relaxed" = "SEACR Relaxed")
    ) +
    scale_fill_manual(
      values = c("MACS2_Default" = "#1f77b4", "MACS2_Large" = "#ff7f0e", 
                 "SEACR_Stringent" = "#2ca02c", "SEACR_Relaxed" = "#d62728"),
      labels = c("MACS2_Default" = "MACS2 Default", "MACS2_Large" = "MACS2 Large",
                 "SEACR_Stringent" = "SEACR Stringent", "SEACR_Relaxed" = "SEACR Relaxed")
    ) +
    scale_shape_manual(
      values = shape_values
    ) +
    labs(
      title = paste(species, mark, "Peak Counts by Method"),
      subtitle = "Boxplot by cell number, all methods side by side",
      x = "Cell Number",
      y = "Number of Peaks",
      color = "Peak Calling Method",
      fill = "Peak Calling Method",
      shape = "Replicate"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "horizontal"
    )
  
  return(p)
}



# Step 7: Main execution
main_analysis <- function() {
  cat("Starting peak analysis...\n")
  
  # Step 1: Identify outliers (for information only, not for exclusion)
  cat("\nStep 2: Identifying outliers using Grubbs test...\n")
  outlier_results <- identify_all_outliers(Peaksummary)
  Peaksummary_with_outliers <- outlier_results$data
  
  # Count outliers for each peak type
  outlier_counts <- list(
    Default = sum(Peaksummary_with_outliers$is_outlier_default),
    Large = sum(Peaksummary_with_outliers$is_outlier_large),
    SEACR_Stringent = sum(Peaksummary_with_outliers$is_outlier_seacr_stringent),
    SEACR_Relaxed = sum(Peaksummary_with_outliers$is_outlier_seacr_relaxed)
  )
  
  cat("Outliers detected (not excluded from plots):\n")
  for(peak_type in names(outlier_counts)) {
    cat("  ", peak_type, ":", outlier_counts[[peak_type]], "\n")
  }
  
  # Step 2: Filter and prepare data
  cat("\nStep 3: Filtering and preparing data...\n")
  Peaksummary_filtered <- filter_and_prepare_data(Peaksummary_with_outliers)
  
  cat("Data points after filtering:", nrow(Peaksummary_filtered), "\n")
  cat("Unique cell numbers:", unique(Peaksummary_filtered$Cell_numbers), "\n")
  
  # Step 3: Create raw data for box plots (ALL DATA INCLUDED)
  cat("\nStep 4: Creating raw data for box plots (all data included)...\n")
  peak_types <- c("Default", "Large", "SEACR_Stringent", "SEACR_Relaxed")
  
  combined_raw <- data.frame()
  
  for(peak_type in peak_types) {
    # Create raw data for box plots (no outlier exclusion)
    raw_data <- create_raw_data_for_plots(Peaksummary_filtered, peak_type)
    combined_raw <- rbind(combined_raw, raw_data)
  }
  
  # Debug: Check what peak types are actually in the data
  cat("Peak types in raw data:", unique(combined_raw$Peak_Type), "\n")
  cat("Number of rows per peak type:\n")
  if(nrow(combined_raw) > 0) {
    print(table(combined_raw$Peak_Type))
  }
  
  # Step 4: Create the four plots
  cat("\nStep 5: Creating box plots...\n")
  species_marks <- list(
    c("Mouse", "CTCF"),
    c("Mouse", "K27me3"),
    c("Rat", "CTCF"),
    c("Rat", "K27me3")
  )
  
  all_plots <- list()
  
  for(i in 1:length(species_marks)) {
    species <- species_marks[[i]][1]
    mark <- species_marks[[i]][2]
    
    plot_name <- paste0(tolower(species), "_", tolower(mark))
    
    p <- create_combined_plot(combined_raw, species, mark)
    all_plots[[plot_name]] <- p
    
    cat("  Created plot:", plot_name, "\n")
  }
  
  # Step 5: Save all plots
  cat("\nStep 6: Saving plots...\n")
  for(plot_name in names(all_plots)) {
    filename <- paste0(plot_name, "_boxplots.pdf")
    ggsave(filename, all_plots[[plot_name]], width = 12, height = 8, dpi = 300)
    cat("  Saved:", filename, "\n")
  }
  
  cat("\nAnalysis complete! 4 plots have been saved.\n")
  
  # Return processed data and plots for further use
  return(list(
    data = Peaksummary_filtered,
    raw_data = combined_raw,
    outliers = outlier_results$outliers,
    plots = all_plots
  ))
}

# Run the analysis
results <- main_analysis()


########### Henikoff data analysis ########

setwd("/Users/atrayeeray/Desktop/Sid/CUT&Tag_cellno. standardization/Direct CUT&Tag/new_analysis_07082025/Henikoff_cell_no_standardization/final_peakfiles_cutoff_0.01_basedon_IgG_R1/")

Henikoff_peak_summary <- read.csv("peak_summary.csv")
# Remove NA values
Henikoff_peak_summary <- Henikoff_peak_summary[rowSums(!is.na(Henikoff_peak_summary[, -1])) > 0, ][, colSums(!is.na(Henikoff_peak_summary[rowSums(!is.na(Henikoff_peak_summary[, -1])) > 0, ])) > 0]

###remove 0.2K, and 0.06K, and 0.6K data
Henikoff_peak_summary <- Henikoff_peak_summary %>% filter(!Sample %in% c("H0.2K_K27me3", "H0.06K_K27me3", "H0.6K_K27me3") )

# Step 1: Function to filter and prepare data
filter_and_prepare_data <- function(data) {
  data_filtered <- data %>%
    # Add species information
    mutate(
      Species = "Henikoff",
      Mark = "K27me3"
    ) %>%
    # Convert cell numbers to numeric for proper ordering - FIXED ORDER
    mutate(
      Cell_numbers_numeric = case_when(
        grepl("60K", Sample) ~ 60000,       # Check first
        grepl("20K", Sample) ~ 20000,       # Check second
        grepl("6K", Sample) ~ 6000,         # Check third
        grepl("2K", Sample) ~ 2000          # Check LAST - least specific
      )
    )
  
  # Debug: Print what we have
  cat("Data after filtering:\n")
  print(data_filtered[, c("Sample", "Cell_numbers_numeric")])
  
  return(data_filtered)
}

# Step 2: Function to create summary statistics for plotting
create_summary_stats <- function(data, peak_type, igg_type) {
  # Get the correct peak column name
  if(igg_type == "IgG_R1") {
    if(peak_type == "Default") {
      peak_col <- "MACS2_peak_default_based.on.IgG_R1"
    } else if(peak_type == "Large") {
      peak_col <- "MACS2_peak_large_based.on.IgG_R1"
    } else if(peak_type == "SEACR_Stringent") {
      peak_col <- "Seacr_Stringent_based.on.IgG_R1"
    } else if(peak_type == "SEACR_Relaxed") {
      peak_col <- "Seacr_relaxed_based.on.IgG_R1"
    }
  }

  
  # Create summary statistics (no outliers to remove)
  summary_data <- data %>%
    group_by(Cell_numbers_numeric, Species, Mark) %>%
    summarise(
      mean_peaks = mean(get(peak_col), na.rm = TRUE),
      sd_peaks = 0,  # No replicates in this dataset
      n_replicates = n(),
      .groups = 'drop'
    )
  
  return(summary_data)
}

# Step 3: Function to create combined plot for a specific IgG type (to exclude 0.2K, and 0.06K data)
create_combined_plot <- function(data, igg_type) {
  # Create plot with all four peak types
  p <- ggplot(data, aes(x = Cell_numbers_numeric, y = mean_peaks, color = Peak_Type, shape = Peak_Type)) +
    geom_point(size = 3) +
    geom_line(aes(group = Peak_Type), alpha = 0.7, linewidth = 1) +
    # Use log scale for BOTH axes to spread out the smaller cell numbers
    scale_x_log10(breaks = c(2000, 6000, 20000, 60000),
                  labels = c("2K", "6K", "20K", "60K")) +
    scale_y_log10(
      labels = scales::comma,
      limits = c(1000, NA)  # Start Y-axis from 100
    ) +
    scale_color_manual(
      values = c("MACS2_Default" = "#1f77b4", "MACS2_Large" = "#ff7f0e", 
                 "SEACR_Stringent" = "#2ca02c", "SEACR_Relaxed" = "#d62728"),
      labels = c("MACS2_Default" = "MACS2 Default", "MACS2_Large" = "MACS2 Large",
                 "SEACR_Stringent" = "SEACR Stringent", "SEACR_Relaxed" = "SEACR Relaxed")
    ) +
    scale_shape_manual(
      values = c("MACS2_Default" = 16, "MACS2_Large" = 17, 
                 "SEACR_Stringent" = 15, "SEACR_Relaxed" = 18),
      labels = c("MACS2_Default" = "MACS2 Default", "MACS2_Large" = "MACS2 Large",
                 "SEACR_Stringent" = "SEACR Stringent", "SEACR_Relaxed" = "SEACR Relaxed")
    ) +
    labs(
      title = paste("Henikoff K27me3 Peak Counts by Method -", igg_type),
      subtitle = "Cell Number vs Peak Count (Log Scale)",
      x = "Cell Number",
      y = "Mean Peak Count (log scale)",
      color = "Peak Calling Method",
      shape = "Peak Calling Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal",
      # Add minor grid lines for better readability
      panel.grid.minor = element_line(color = "gray90"),
      panel.grid.major = element_line(color = "gray70"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
    )
  
  return(p)
}

# Step 4: Main execution (to exclude 0.2K, and 0.06K data)
main_analysis <- function() {
  cat("Starting Henikoff peak analysis (NO outlier removal)...\n")
  
  # Step 1: Filter and prepare data
  cat("Step 1: Filtering and preparing data...\n")
  Henikoff_filtered <- filter_and_prepare_data(Henikoff_peak_summary)
  
  cat("Data points after filtering:", nrow(Henikoff_filtered), "\n")
  cat("Unique cell numbers:", unique(Henikoff_filtered$Cell_numbers_numeric), "\n")
  
  # Step 2: Create combined data for plotting for each IgG type
  cat("\nStep 2: Creating combined data for plotting...\n")
  peak_types <- c("Default", "Large", "SEACR_Stringent", "SEACR_Relaxed")
  igg_types <- c("IgG_R1")
  
  all_plots <- list()
  
  for(igg_type in igg_types) {
    cat("Processing", igg_type, "...\n")
    
    combined_data <- data.frame()
    
    for(peak_type in peak_types) {
      summary_data <- create_summary_stats(Henikoff_filtered, peak_type, igg_type)
      
      # Use consistent naming for peak types
      if(peak_type == "Default") {
        summary_data$Peak_Type <- "MACS2_Default"
      } else if(peak_type == "Large") {
        summary_data$Peak_Type <- "MACS2_Large"
      } else if(peak_type == "SEACR_Stringent") {
        summary_data$Peak_Type <- "SEACR_Stringent"
      } else if(peak_type == "SEACR_Relaxed") {
        summary_data$Peak_Type <- "SEACR_Relaxed"
      }
      
      combined_data <- rbind(combined_data, summary_data)
    }
    
    # Debug: Check what peak types are actually in the data
    cat("  Peak types in combined data for", igg_type, ":", unique(combined_data$Peak_Type), "\n")
    cat("  Total data points for plotting:", nrow(combined_data), "\n")
    
    # Create plot for this IgG type
    plot_name <- paste0("henikoff_k27me3_", tolower(gsub("_", "", igg_type)))
    
    p <- create_combined_plot(combined_data, igg_type)
    all_plots[[plot_name]] <- p
    
    cat("  Created plot:", plot_name, "\n")
  }
  
  # Step 3: Save all plots
  cat("\nStep 3: Saving plots...\n")
  for(plot_name in names(all_plots)) {
    filename <- paste0(plot_name, "_combined_peaks.pdf")
    ggsave(filename, all_plots[[plot_name]], width = 12, height = 8, dpi = 300)
    cat("  Saved:", filename, "\n")
  }
  
  cat("\nAnalysis complete! 2 plots have been saved.\n")
  
  # Return processed data and plots for further use
  return(list(
    data = Henikoff_filtered,
    plots = all_plots
  ))
}

# Run the analysis
results <- main_analysis()


###remove 0.2K, and 0.06K data
#Henikoff_peak_summary <- Henikoff_peak_summary %>% filter(!Sample %in% c("H0.2K_K27me3", "H0.06K_K27me3") )


# Step 1: Function to filter and prepare data
#filter_and_prepare_data <- function(data) {
  data_filtered <- data %>%
    # Add species information
    mutate(
      Species = "Henikoff",
      Mark = "K27me3"
    ) %>%
    # Convert cell numbers to numeric for proper ordering - FIXED ORDER
    mutate(
      Cell_numbers_numeric = case_when(
        grepl("0\\.06K", Sample) ~ 60,      # Check FIRST - most specific
        grepl("0\\.2K", Sample) ~ 200,      # Check SECOND - more specific
        grepl("0\\.6K", Sample) ~ 600,      # Check THIRD - more specific
        grepl("60K", Sample) ~ 60000,       # Check FOURTH
        grepl("20K", Sample) ~ 20000,       # Check FIFTH
        grepl("6K", Sample) ~ 6000,         # Check SIXTH
        grepl("2K", Sample) ~ 2000          # Check LAST - least specific
      )
    )
  
  # Debug: Print what we have
  cat("Data after filtering:\n")
  print(data_filtered[, c("Sample", "Cell_numbers_numeric")])
  
  return(data_filtered)
}

# Step 2: Function to create summary statistics for plotting
#create_summary_stats <- function(data, peak_type, igg_type) {
  # Get the correct peak column name
  if(igg_type == "IgG_R1") {
    if(peak_type == "Default") {
      peak_col <- "MACS2_peak_default_based.on.IgG_R1"
    } else if(peak_type == "Large") {
      peak_col <- "MACS2_peak_large_based.on.IgG_R1"
    } else if(peak_type == "SEACR_Stringent") {
      peak_col <- "Seacr_Stringent_based.on.IgG_R1"
    } else if(peak_type == "SEACR_Relaxed") {
      peak_col <- "Seacr_relaxed_based.on.IgG_R1"
    }
  } else if(igg_type == "IgG_R2") {
    if(peak_type == "Default") {
      peak_col <- "MACS2_peak_default_based.on.IgG_R2"
    } else if(peak_type == "Large") {
      peak_col <- "MACS2_peak_large_based.on.IgG_R2"
    } else if(peak_type == "SEACR_Stringent") {
      peak_col <- "Seacr_Stringent_based.on.IgG_R2"
    } else if(peak_type == "SEACR_Relaxed") {
      peak_col <- "Seacr_relaxed_based.on.IgG_R2"
    }
  }
  
  # Create summary statistics (no outliers to remove)
  summary_data <- data %>%
    group_by(Cell_numbers_numeric, Species, Mark) %>%
    summarise(
      mean_peaks = mean(get(peak_col), na.rm = TRUE),
      sd_peaks = 0,  # No replicates in this dataset
      n_replicates = n(),
      .groups = 'drop'
    )
  
  return(summary_data)
}

# Step 3: Function to create combined plot for a specific IgG type (to exclude 0.2K, and 0.06K data)
#create_combined_plot <- function(data, igg_type) {
  # Create plot with all four peak types
  p <- ggplot(data, aes(x = Cell_numbers_numeric, y = mean_peaks, color = Peak_Type, shape = Peak_Type)) +
    geom_point(size = 3) +
    geom_line(aes(group = Peak_Type), alpha = 0.7, linewidth = 1) +
    # Use log scale for BOTH axes to spread out the smaller cell numbers
    scale_x_log10(breaks = c(600, 2000, 6000, 20000, 60000),
                  labels = c("0.6K", "2K", "6K", "20K", "60K")) +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(
      values = c("MACS2_Default" = "#1f77b4", "MACS2_Large" = "#ff7f0e", 
                 "SEACR_Stringent" = "#2ca02c", "SEACR_Relaxed" = "#d62728"),
      labels = c("MACS2_Default" = "MACS2 Default", "MACS2_Large" = "MACS2 Large",
                 "SEACR_Stringent" = "SEACR Stringent", "SEACR_Relaxed" = "SEACR Relaxed")
    ) +
    scale_shape_manual(
      values = c("MACS2_Default" = 16, "MACS2_Large" = 17, 
                 "SEACR_Stringent" = 15, "SEACR_Relaxed" = 18),
      labels = c("MACS2_Default" = "MACS2 Default", "MACS2_Large" = "MACS2 Large",
                 "SEACR_Stringent" = "SEACR Stringent", "SEACR_Relaxed" = "SEACR Relaxed")
    ) +
    labs(
      title = paste("Henikoff K27me3 Peak Counts by Method -", igg_type),
      subtitle = "Cell Number vs Peak Count (Log Scale)",
      x = "Cell Number",
      y = "Mean Peak Count (log scale)",
      color = "Peak Calling Method",
      shape = "Peak Calling Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal",
      # Add minor grid lines for better readability
      panel.grid.minor = element_line(color = "gray90"),
      panel.grid.major = element_line(color = "gray70"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
    )
  
  return(p)
}

# Step 4: Main execution (to exclude 0.2K, and 0.06K data)
#main_analysis <- function() {
  cat("Starting Henikoff peak analysis (NO outlier removal)...\n")
  
  # Step 1: Filter and prepare data
  cat("Step 1: Filtering and preparing data...\n")
  Henikoff_filtered <- filter_and_prepare_data(Henikoff_peak_summary)
  
  cat("Data points after filtering:", nrow(Henikoff_filtered), "\n")
  cat("Unique cell numbers:", unique(Henikoff_filtered$Cell_numbers_numeric), "\n")
  
  # Step 2: Create combined data for plotting for each IgG type
  cat("\nStep 2: Creating combined data for plotting...\n")
  peak_types <- c("Default", "Large", "SEACR_Stringent", "SEACR_Relaxed")
  igg_types <- c("IgG_R1", "IgG_R2")
  
  all_plots <- list()
  
  for(igg_type in igg_types) {
    cat("Processing", igg_type, "...\n")
    
    combined_data <- data.frame()
    
    for(peak_type in peak_types) {
      summary_data <- create_summary_stats(Henikoff_filtered, peak_type, igg_type)
      
      # Use consistent naming for peak types
      if(peak_type == "Default") {
        summary_data$Peak_Type <- "MACS2_Default"
      } else if(peak_type == "Large") {
        summary_data$Peak_Type <- "MACS2_Large"
      } else if(peak_type == "SEACR_Stringent") {
        summary_data$Peak_Type <- "SEACR_Stringent"
      } else if(peak_type == "SEACR_Relaxed") {
        summary_data$Peak_Type <- "SEACR_Relaxed"
      }
      
      combined_data <- rbind(combined_data, summary_data)
    }
    
    # Debug: Check what peak types are actually in the data
    cat("  Peak types in combined data for", igg_type, ":", unique(combined_data$Peak_Type), "\n")
    cat("  Total data points for plotting:", nrow(combined_data), "\n")
    
    # Create plot for this IgG type
    plot_name <- paste0("henikoff_k27me3_", tolower(gsub("_", "", igg_type)))
    
    p <- create_combined_plot(combined_data, igg_type)
    all_plots[[plot_name]] <- p
    
    cat("  Created plot:", plot_name, "\n")
  }
  
  # Step 3: Save all plots
  cat("\nStep 3: Saving plots...\n")
  for(plot_name in names(all_plots)) {
    filename <- paste0(plot_name, "_combined_peaks.pdf")
    ggsave(filename, all_plots[[plot_name]], width = 12, height = 8, dpi = 300)
    cat("  Saved:", filename, "\n")
  }
  
  cat("\nAnalysis complete! 2 plots have been saved.\n")
  
  # Return processed data and plots for further use
  return(list(
    data = Henikoff_filtered,
    plots = all_plots
  ))
}

# Run the analysis
#results <- main_analysis()

### correlation analysis
# Create comprehensive ratio impact plot
create_ratio_impact_plot <- function(data) {
  
  # Calculate sample/IgG ratios
  analysis_data <- data %>%
    mutate(
      Ratio_to_IgG_R1 = round(Mapped.Reads.in.Million / Mapped.IgG_R1.reads, 2),
      Ratio_to_IgG_R2 = round(Mapped.Reads.in.Million / Mapped.IgG_R2.reads, 2)
    )
  
  # Create long format data for plotting
  plot_data_r1 <- analysis_data %>%
    select(Sample, Ratio_to_IgG_R1, 
           MACS2_Default = MACS2_peak_default_based.on.IgG_R1,
           MACS2_Large = MACS2_peak_large_based.on.IgG_R1,
           SEACR_Stringent = Seacr_Stringent_based.on.IgG_R1,
           SEACR_Relaxed = Seacr_relaxed_based.on.IgG_R1) %>%
    pivot_longer(cols = c(MACS2_Default, MACS2_Large, SEACR_Stringent, SEACR_Relaxed),
                 names_to = "Method", values_to = "Peak_Count")
  
  plot_data_r2 <- analysis_data %>%
    select(Sample, Ratio_to_IgG_R2, 
           MACS2_Default = MACS2_peak_default_based.on.IgG_R2,
           MACS2_Large = MACS2_peak_large_based.on.IgG_R2,
           SEACR_Stringent = Seacr_Stringent_based.on.IgG_R2,
           SEACR_Relaxed = Seacr_relaxed_based.on.IgG_R2) %>%
    pivot_longer(cols = c(MACS2_Default, MACS2_Large, SEACR_Stringent, SEACR_Relaxed),
                 names_to = "Method", values_to = "Peak_Count")
  
  # Create plots
  p1 <- ggplot(plot_data_r1, aes(x = Ratio_to_IgG_R1, y = Peak_Count, color = Method, shape = Method)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
    scale_color_manual(values = c("MACS2_Default" = "#1f77b4", "MACS2_Large" = "#ff7f0e", 
                                  "SEACR_Stringent" = "#2ca02c", "SEACR_Relaxed" = "#d62728")) +
    scale_shape_manual(values = c("MACS2_Default" = 16, "MACS2_Large" = 17, 
                                  "SEACR_Stringent" = 15, "SEACR_Relaxed" = 18)) +
    labs(title = "Sample/IgG_R1 Ratio vs Peak Counts",
         subtitle = "How background normalization depth affects peak detection",
         x = "Sample Mapped Reads / IgG_R1 Reads", y = "Peak Count", 
         color = "Peak Calling Method", shape = "Peak Calling Method") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  p2 <- ggplot(plot_data_r2, aes(x = Ratio_to_IgG_R2, y = Peak_Count, color = Method, shape = Method)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
    scale_color_manual(values = c("MACS2_Default" = "#1f77b4", "MACS2_Large" = "#ff7f0e", 
                                  "SEACR_Stringent" = "#2ca02c", "SEACR_Relaxed" = "#d62728")) +
    scale_shape_manual(values = c("MACS2_Default" = 16, "MACS2_Large" = 17, 
                                  "SEACR_Stringent" = 15, "SEACR_Relaxed" = 18)) +
    labs(title = "Sample/IgG_R2 Ratio vs Peak Counts",
         subtitle = "How background normalization depth affects peak detection",
         x = "Sample Mapped Reads / IgG_R2 Reads", y = "Peak Count", 
         color = "Peak Calling Method", shape = "Peak Calling Method") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  # Save plots
  ggsave("ratio_impact_igg_r1.pdf", p1, width = 12, height = 8, dpi = 300)
  ggsave("ratio_impact_igg_r2.pdf", p2, width = 12, height = 8, dpi = 300)
  
  return(list(plot_r1 = p1, plot_r2 = p2))
}

# Create faceted plot for better comparison
create_faceted_ratio_plot <- function(data) {
  
  # Calculate ratios
  analysis_data <- data %>%
    mutate(
      Ratio_to_IgG_R1 = round(Mapped.Reads.in.Million / Mapped.IgG_R1.reads, 2),
      Ratio_to_IgG_R2 = round(Mapped.Reads.in.Million / Mapped.IgG_R2.reads, 2)
    )
  
  # Create long format data for both IgG types
  plot_data_combined <- bind_rows(
    # IgG_R1 data
    analysis_data %>%
      select(Sample, Ratio = Ratio_to_IgG_R1, 
             MACS2_Default = MACS2_peak_default_based.on.IgG_R1,
             MACS2_Large = MACS2_peak_large_based.on.IgG_R1,
             SEACR_Stringent = Seacr_Stringent_based.on.IgG_R1,
             SEACR_Relaxed = Seacr_relaxed_based.on.IgG_R1) %>%
      pivot_longer(cols = c(MACS2_Default, MACS2_Large, SEACR_Stringent, SEACR_Relaxed),
                   names_to = "Method", values_to = "Peak_Count") %>%
      mutate(IgG_Type = "IgG_R1"),
    
    # IgG_R2 data
    analysis_data %>%
      select(Sample, Ratio = Ratio_to_IgG_R2, 
             MACS2_Default = MACS2_peak_default_based.on.IgG_R2,
             MACS2_Large = MACS2_peak_large_based.on.IgG_R2,
             SEACR_Stringent = Seacr_Stringent_based.on.IgG_R2,
             SEACR_Relaxed = Seacr_relaxed_based.on.IgG_R2) %>%
      pivot_longer(cols = c(MACS2_Default, MACS2_Large, SEACR_Stringent, SEACR_Relaxed),
                   names_to = "Method", values_to = "Peak_Count") %>%
      mutate(IgG_Type = "IgG_R2")
  )
  
  # Create faceted plot
  p <- ggplot(plot_data_combined, aes(x = Ratio, y = Peak_Count, color = Method)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    facet_grid(Method ~ IgG_Type, scales = "free_y") +
    scale_color_manual(values = c("MACS2_Default" = "#1f77b4", "MACS2_Large" = "#ff7f0e", 
                                  "SEACR_Stringent" = "#2ca02c", "SEACR_Relaxed" = "#d62728")) +
    labs(title = "Impact of Sample/IgG Ratio on Peak Detection",
         subtitle = "Each panel shows how a specific method responds to ratio changes",
         x = "Sample Mapped Reads / IgG Reads", y = "Peak Count", 
         color = "Peak Calling Method") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          strip.text = element_text(face = "bold"))
  
  # Save plot
  ggsave("faceted_ratio_impact.pdf", p, width = 14, height = 10, dpi = 300)
  
  return(p)
}

# Create correlation heatmap
create_ratio_correlation_heatmap <- function(data) {
  
  # Calculate ratios
  analysis_data <- data %>%
    mutate(
      Ratio_to_IgG_R1 = Mapped.Reads.in.Million / Mapped.IgG_R1.reads,
      Ratio_to_IgG_R2 = Mapped.Reads.in.Million / Mapped.IgG_R2.reads
    )
  
  # Create correlation matrix
  cor_data <- analysis_data %>%
    select(Ratio_to_IgG_R1, Ratio_to_IgG_R2,
           MACS2_Default_R1 = MACS2_peak_default_based.on.IgG_R1,
           MACS2_Large_R1 = MACS2_peak_large_based.on.IgG_R1,
           SEACR_Stringent_R1 = Seacr_Stringent_based.on.IgG_R1,
           SEACR_Relaxed_R1 = Seacr_relaxed_based.on.IgG_R1,
           MACS2_Default_R2 = MACS2_peak_default_based.on.IgG_R2,
           MACS2_Large_R2 = MACS2_peak_large_based.on.IgG_R2,
           SEACR_Stringent_R2 = Seacr_Stringent_based.on.IgG_R2,
           SEACR_Relaxed_R2 = Seacr_relaxed_based.on.IgG_R2)
  
  # Calculate correlation matrix
  cor_matrix <- cor(cor_data, use = "complete.obs", method = "pearson")
  
  # Create heatmap
  p <- ggplot(as.data.frame(as.table(cor_matrix)), aes(Var1, Var2, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.3f", Freq)), size = 3) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Correlation: Sample/IgG Ratio vs Peak Counts",
         subtitle = "Red = Negative correlation, Blue = Positive correlation, White = No correlation",
         x = "", y = "", fill = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Save plot
  ggsave("ratio_correlation_heatmap.pdf", p, width = 12, height = 10, dpi = 300)
  
  return(p)
}

# Create all three types of plots
ratio_plots <- create_ratio_impact_plot(Henikoff_peak_summary)
faceted_plot <- create_faceted_ratio_plot(Henikoff_peak_summary)
correlation_heatmap <- create_ratio_correlation_heatmap(Henikoff_peak_summary)

### peak quality check

# Load required libraries
library(GenomicRanges)
library(ggplot2)
library(dplyr)

# Function to extract cell number from filename
extract_cell_number <- function(filename) {
  cell_match <- regexpr("H([0-9.]+K)", filename)
  if(cell_match > 0) {
    cell_str <- substr(filename, cell_match + 1, cell_match + attr(cell_match, "match.length") - 1)
    return(cell_str)
  }
  return(NA)
}

# Function to read broadPeak files (for MACS2)
read_broadPeak <- function(file_path) {
  if(!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(GRanges())
  }
  
  peaks <- read.table(file_path, header=FALSE, sep="\t")
  if(nrow(peaks) == 0) return(GRanges())
  
  gr <- GRanges(
    seqnames = peaks$V1,
    ranges = IRanges(start = peaks$V2, end = peaks$V3),
    score = peaks$V7
  )
  
  return(gr)
}

# Function to read BED files (for SEACR)
read_bed <- function(file_path) {
  if(!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(GRanges())
  }
  
  peaks <- read.table(file_path, header=FALSE, sep="\t")
  if(nrow(peaks) == 0) return(GRanges())
  
  gr <- GRanges(
    seqnames = peaks$V1,
    ranges = IRanges(start = peaks$V2, end = peaks$V3)
  )
  
  return(gr)
}

# Function to analyze peak overlaps
analyze_peak_overlaps <- function(base_dir) {
  
  cat("=== PEAK OVERLAP ANALYSIS AMONG CELL NUMBERS ===\n\n")
  
  # Define the EXACT cell number order you want
  cell_order <- c("0.06K", "0.2K", "0.6K", "2K", "6K", "20K", "60K")
  methods <- c("macs2_default", "macs2_large", "seacr_stringent", "seacr_relaxed")
  igg_backgrounds <- c("R1", "R2")
  
  # Initialize results
  overlap_results <- data.frame()
  
  for(method in methods) {
    for(igg in igg_backgrounds) {
      cat("Processing:", method, "with", igg, "\n")
      
      # Get all files for this method and IgG background
      if(method == "macs2_default") {
        folder <- file.path(base_dir, "macs2_default")
        pattern <- paste0("H.*K27me3\\.dedup_vs_H100K_IgG_", igg, "\\.dedup\\.macs2\\.peak\\.q0\\.05-broad_peaks\\.broadPeak")
      } else if(method == "macs2_large") {
        folder <- file.path(base_dir, "macs2_large")
        pattern <- paste0("H.*K27me3\\.dedup_vs_H100K_IgG_", igg, "\\.dedup\\.macs2\\.peak\\.q0\\.05-broad_peaks\\.broadPeak")
      } else if(method == "seacr_stringent") {
        folder <- file.path(base_dir, "seacr_stringent")
        pattern <- paste0("H.*K27me3_vs_H100K_IgG_", igg, "_stringent\\.stringent\\.bed")
      } else if(method == "seacr_relaxed") {
        folder <- file.path(base_dir, "seacr_relaxed")
        pattern <- paste0("H.*K27me3_vs_H100K_IgG_", igg, "_relaxed\\.relaxed\\.bed")
      }
      
      # Check if folder exists
      if(!dir.exists(folder)) {
        cat("  Folder does not exist:", folder, "\n")
        next
      }
      
      files <- list.files(path = folder, pattern = pattern, full.names = TRUE)
      cat("  Found", length(files), "files\n")
      
      if(length(files) < 2) {
        cat("  Less than 2 files found, skipping...\n")
        next
      }
      
      # Read all peak files and organize by the specified cell order
      peaks_list <- list()
      for(file in files) {
        cell_num <- extract_cell_number(file)
        if(!is.na(cell_num)) {
          cat("    Reading", cell_num, "from", basename(file), "\n")
          if(method %in% c("macs2_default", "macs2_large")) {
            peaks_list[[cell_num]] <- read_broadPeak(file)
          } else if(method %in% c("seacr_stringent", "seacr_relaxed")) {
            peaks_list[[cell_num]] <- read_bed(file)
          }
        }
      }
      
      cat("  Successfully read peaks for", length(peaks_list), "cell numbers\n")
      
      # Reorder peaks_list according to the specified cell order
      peaks_list <- peaks_list[intersect(cell_order, names(peaks_list))]
      
      # Calculate overlaps between cell numbers in the specified order
      cell_nums <- names(peaks_list)
      for(i in 1:(length(cell_nums)-1)) {
        for(j in (i+1):length(cell_nums)) {
          cell1 <- cell_nums[i]
          cell2 <- cell_nums[j]
          
          # Calculate overlap statistics
          peaks1 <- peaks_list[[cell1]]
          peaks2 <- peaks_list[[cell2]]
          
          # Find overlaps
          overlaps <- findOverlaps(peaks1, peaks2)
          
          # Calculate metrics
          total_peaks1 <- length(peaks1)
          total_peaks2 <- length(peaks2)
          overlapping_peaks <- length(unique(queryHits(overlaps)))
          
          # Calculate overlap percentages
          overlap_percent_1to2 <- (overlapping_peaks / total_peaks1) * 100
          overlap_percent_2to1 <- (overlapping_peaks / total_peaks2) * 100
          jaccard_index <- overlapping_peaks / (total_peaks1 + total_peaks2 - overlapping_peaks)
          
          # Store results
          overlap_results <- rbind(overlap_results, data.frame(
            Method = method,
            IgG_Background = igg,
            Cell1 = cell1,
            Cell2 = cell2,
            Cell_Pair = paste(cell1, "vs", cell2),
            Total_Peaks_Cell1 = total_peaks1,
            Total_Peaks_Cell2 = total_peaks2,
            Overlapping_Peaks = overlapping_peaks,
            Overlap_Percent_1to2 = overlap_percent_1to2,
            Overlap_Percent_2to1 = overlap_percent_2to1,
            Jaccard_Index = jaccard_index
          ))
        }
      }
    }
  }
  
  # Convert Cell1 and Cell2 to factors with the specified order
  overlap_results$Cell1 <- factor(overlap_results$Cell1, levels = cell_order)
  overlap_results$Cell2 <- factor(overlap_results$Cell2, levels = cell_order)
  
  return(overlap_results)
}

# Function to create overlap visualization with better label visibility
create_overlap_plots <- function(overlap_results) {
  
  if(nrow(overlap_results) == 0) {
    cat("No overlap results to plot\n")
    return(NULL)
  }
  
  # Plot 1: Overlap percentage heatmap
  p1 <- ggplot(overlap_results, aes(x = Cell1, y = Cell2, fill = Overlap_Percent_1to2)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1f%%", Overlap_Percent_1to2)), size = 2.5) +
    scale_fill_gradient2(low = "white", mid = "yellow", high = "red", 
                         midpoint = 50, limits = c(0, 100)) +
    facet_grid(Method ~ IgG_Background) +
    labs(title = "Peak Overlap Percentages Among Cell Numbers",
         subtitle = "Shows how many peaks from Cell1 overlap with Cell2",
         x = "Cell Number 1", y = "Cell Number 2", 
         fill = "Overlap %") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10, margin = margin(t = 10)),
      axis.title.y = element_text(size = 10, margin = margin(r = 10)),
      strip.text = element_text(size = 9),
      plot.title = element_text(size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    )
  
  # Plot 2: Jaccard Index heatmap
  p2 <- ggplot(overlap_results, aes(x = Cell1, y = Cell2, fill = Jaccard_Index)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.3f", Jaccard_Index)), size = 2.5) +
    scale_fill_gradient2(low = "white", mid = "yellow", high = "red", 
                         midpoint = 0.5, limits = c(0, 1)) +
    facet_grid(Method ~ IgG_Background) +
    labs(title = "Peak Similarity (Jaccard Index) Among Cell Numbers",
         subtitle = "Higher values = more similar peak sets",
         x = "Cell Number 1", y = "Cell Number 2", 
         fill = "Jaccard Index") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10, margin = margin(t = 10)),
      axis.title.y = element_text(size = 10, margin = margin(r = 10)),
      strip.text = element_text(size = 9),
      plot.title = element_text(size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    )
  
  # Save plots with larger dimensions and better margins
  ggsave("peak_overlap_percentages.pdf", p1, 
         width = 24, height = 20, dpi = 300, 
         limitsize = FALSE)
  
  ggsave("peak_jaccard_similarity.pdf", p2, 
         width = 24, height = 20, dpi = 300, 
         limitsize = FALSE)
  
  return(list(overlap_plot = p1, jaccard_plot = p2))
}

# Function to identify conserved peaks
identify_conserved_peaks <- function(base_dir, method, igg_background) {
  
  cat("=== IDENTIFYING CONSERVED PEAKS ===\n")
  cat("Method:", method, "IgG Background:", igg_background, "\n\n")
  
  # Define the EXACT cell number order
  cell_order <- c("0.06K", "0.2K", "0.6K", "2K", "6K", "20K", "60K")
  
  # Get all files for this method and IgG background
  if(method == "macs2_default") {
    folder <- file.path(base_dir, "macs2_default")
    pattern <- paste0("H.*K27me3\\.dedup_vs_H100K_IgG_", igg_background, "\\.dedup\\.macs2\\.peak\\.q0\\.05-broad_peaks\\.broadPeak")
  } else if(method == "macs2_large") {
    folder <- file.path(base_dir, "macs2_large")
    pattern <- paste0("H.*K27me3\\.dedup_vs_H100K_IgG_", igg_background, "\\.dedup\\.macs2\\.peak\\.q0\\.05-broad_peaks\\.broadPeak")
  } else if(method == "seacr_stringent") {
    folder <- file.path(base_dir, "seacr_stringent")
    pattern <- paste0("H.*K27me3_vs_H100K_IgG_", igg_background, "_stringent\\.stringent\\.bed")
  } else if(method == "seacr_relaxed") {
    folder <- file.path(base_dir, "seacr_relaxed")
    pattern <- paste0("H.*K27me3_vs_H100K_IgG_", igg_background, "_relaxed\\.relaxed\\.bed")
  }
  
  if(!dir.exists(folder)) {
    cat("Folder does not exist:", folder, "\n")
    return(NULL)
  }
  
  files <- list.files(path = folder, pattern = pattern, full.names = TRUE)
  
  if(length(files) < 2) {
    cat("Less than 2 files found\n")
    return(NULL)
  }
  
  # Read all peak files
  peaks_list <- list()
  for(file in files) {
    cell_num <- extract_cell_number(file)
    if(!is.na(cell_num)) {
      if(method %in% c("macs2_default", "macs2_large")) {
        peaks_list[[cell_num]] <- read_broadPeak(file)
      } else if(method %in% c("seacr_stringent", "seacr_relaxed")) {
        peaks_list[[cell_num]] <- read_bed(file)
      }
    }
  }
  
  if(length(peaks_list) < 2) {
    cat("Could not read enough peak files\n")
    return(NULL)
  }
  
  # Reorder peaks_list according to the specified cell order
  peaks_list <- peaks_list[intersect(cell_order, names(peaks_list))]
  
  # Find peaks that appear in multiple cell numbers
  cell_nums <- names(peaks_list)
  all_peaks <- peaks_list[[1]]  # Start with first cell number
  
  for(cell_num in cell_nums[-1]) {
    all_peaks <- intersect(all_peaks, peaks_list[[cell_num]])
  }
  
  cat("Conserved peaks found in ALL cell numbers:", length(all_peaks), "\n")
  
  # Find peaks that appear in at least 50% of cell numbers
  min_cells <- ceiling(length(cell_nums) * 0.5)
  frequent_peaks <- GRanges()
  
  for(cell_num in cell_nums) {
    peaks <- peaks_list[[cell_num]]
    if(length(peaks) > 0) {
      for(peak_idx in 1:length(peaks)) {
        peak <- peaks[peak_idx]
        # Count how many cell numbers this peak appears in
        count <- 0
        for(other_cell in cell_nums) {
          if(other_cell != cell_num) {
            overlaps <- findOverlaps(peak, peaks_list[[other_cell]])
            if(length(overlaps) > 0) count <- count + 1
          }
        }
        if(count >= min_cells - 1) {
          frequent_peaks <- c(frequent_peaks, peak)
        }
      }
    }
  }
  
  frequent_peaks <- unique(frequent_peaks)
  cat("Peaks found in at least", min_cells, "cell numbers:", length(frequent_peaks), "\n")
  
  return(list(
    conserved_peaks = all_peaks,
    frequent_peaks = frequent_peaks,
    total_cell_numbers = length(cell_nums)
  ))
}

# Run the analysis
cat("Starting peak overlap analysis...\n")
overlap_results <- analyze_peak_overlaps(base_dir)

# Create overlap plots
overlap_plots <- create_overlap_plots(overlap_results)

# Identify conserved peaks for each method
for(method in c("macs2_default", "macs2_large", "seacr_stringent", "seacr_relaxed")) {
  for(igg in c("R1", "R2")) {
    cat("\n", paste(rep("=", 50), collapse = ""), "\n")
    conserved_peaks <- identify_conserved_peaks(base_dir, method, igg)
  }
}

# Function to create individual plots for each method
create_individual_overlap_plots <- function(overlap_results) {
  
  methods <- unique(overlap_results$Method)
  plots_list <- list()
  
  for(method in methods) {
    method_data <- overlap_results %>% filter(Method == method)
    
    p <- ggplot(method_data, aes(x = Cell1, y = Cell2, fill = Jaccard_Index)) +
      geom_tile() +
      geom_text(aes(label = sprintf("%.3f", Jaccard_Index)), size = 4) +
      scale_fill_gradient2(low = "white", mid = "yellow", high = "red", 
                           midpoint = 0.5, limits = c(0, 1)) +
      facet_wrap(~IgG_Background, ncol = 2) +
      labs(title = paste("Peak Similarity (Jaccard Index) -", gsub("_", " ", method)),
           subtitle = "Higher values = more similar peak sets",
           x = "Cell Number 1", y = "Cell Number 2", 
           fill = "Jaccard Index") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 11, hjust = 0.5)
      )
    
    # Save individual plot
    filename <- paste0("jaccard_", gsub("_", "_", method), ".pdf")
    ggsave(filename, p, width = 12, height = 8, dpi = 300)
    
    plots_list[[method]] <- p
    cat("Saved:", filename, "\n")
  }
  
  return(plots_list)
}

# Create individual plots
individual_plots <- create_individual_overlap_plots(overlap_results)

# Check what warnings occurred
cat("\nChecking warnings...\n")
warnings()

