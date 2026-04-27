generate_deltar_plots <- function(DDA1, 
                                  sample_sizes = NULL, 
                                  divergence_filter = NULL,
                                  save_plot = TRUE, 
                                  output_dir = "Plots_Deltar_without_BHC_notitle") {
  all_plots <- list()
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Optional filter by sample size
  if (!is.null(sample_sizes)) {
    DDA1 <- DDA1 %>% filter(sample_size %in% sample_sizes)
  }
  
  # ---- Remove unwanted panels and BHC ----
  DDA1 <- DDA1 %>%
    filter(!method %in% c("bhc_ref", "full_ref"), alg != "bhc")  # remove BHC completely
  
  # Unique parameter combinations (now uses q)
  unique_combos <- unique(DDA1[, c("k", "q", "sample_size")])
  
  for (i in seq_len(nrow(unique_combos))) {
    k_val <- unique_combos$k[i]
    q_val <- unique_combos$q[i]
    sample_val <- unique_combos$sample_size[i]
    
    # Filter for this combination
    df_sub <- DDA1 %>%
      filter(
        ((is.na(k_val) & is.na(k)) | k == k_val),
        ((is.na(q_val) & is.na(q)) | q == q_val),
        sample_size == sample_val
      )
    
    if (!is.null(divergence_filter)) {
      df_sub <- df_sub %>% filter(divergence %in% divergence_filter)
    }
    
    if (nrow(df_sub) == 0) {
      message(sprintf("Skipping: k=%s q=%s sample_size=%s", k_val, q_val, sample_val))
      next
    }
    
    # Define p values
    all_p <- sort(unique(df_sub$p))
    
    df_sub <- df_sub %>%
      mutate(
        group_label = paste(alg, divergence),
        p = factor(p, levels = all_p)
      )
    
    # Define color map using consistent palette
    legend_labels <- unique(df_sub$group_label)
    dark_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(legend_labels))
    color_map <- setNames(dark_colors, legend_labels)
    
    # ---- Plot ----
    plot <- ggplot(df_sub, aes(
      x = p,
      y = deltar,
      color = group_label,
      group = group_label
    )) +
      geom_line(linewidth = 1.0, alpha = 1) +
      geom_point(size = 1.6, alpha = 0.9) +
      
      scale_color_manual(values = color_map) +
      
      facet_grid(
        rows = vars(variable),
        cols = vars(method),
        scales = "free",
        labeller = labeller(
          variable = c(
            BIC = expression(Delta[BIC]),
            hamming = expression(Delta[HD])
          ),
          .default = label_parsed
        )
      ) +
      
      theme_bw(base_size = 13) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(size = 0.3, linetype = "dotted"),
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_blank()  # <-- removed title
      ) +
      labs(color = "", linetype = "")
    
    # ---- Save ----
    filename <- paste0(
      "_k_", ifelse(is.na(k_val), "NA", k_val),
      "_q_", ifelse(is.na(q_val), "NA", q_val),
      "_sample_", sample_val, ".pdf"
    )
    filepath <- file.path(output_dir, filename)
    
    if (save_plot) {
      ggsave(filepath, plot = plot, width = 8, height = 6, device = "pdf")
      message("Saved: ", filepath)
    }
    
    all_plots[[length(all_plots) + 1]] <- plot
  }
  
  return(all_plots)
}



### with bhc
generate_deltar_plots <- function(DDA1, 
                                  p_values = NULL, 
                                  divergence_filter = NULL,
                                  save_plot = TRUE, 
                                  output_dir = "Plots_Deltar_with_BHC_inside_methods_p") {
  all_plots <- list()
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Optional filter by p
  if (!is.null(p_values)) {
    DDA1 <- DDA1 %>% filter(p %in% p_values)
  }
  
  # ---- Remove unwanted panels ----
  DDA1 <- DDA1 %>%
    filter(!method %in% c("bhc_ref", "full_ref"))
  
  # Split data
  bhc_data <- DDA1 %>% filter(alg == "bhc")
  DDA1 <- DDA1 %>% filter(alg != "bhc")
  
  # Unique parameter combinations (NOW USES p)
  unique_combos <- unique(DDA1[, c("k", "q", "p")])
  
  for (i in seq_len(nrow(unique_combos))) {
    k_val <- unique_combos$k[i]
    q_val <- unique_combos$q[i]
    p_val <- unique_combos$p[i]
    
    # Filter for this combination
    df_sub <- DDA1 %>%
      filter(
        ((is.na(k_val) & is.na(k)) | k == k_val),
        ((is.na(q_val) & is.na(q)) | q == q_val),
        p == p_val
      )
    
    bhc_sub <- bhc_data %>%
      filter(p == p_val, q == q_val)
    
    if (!is.null(divergence_filter)) {
      df_sub <- df_sub %>% filter(divergence %in% divergence_filter)
    }
    
    if (nrow(df_sub) == 0 | nrow(bhc_sub) == 0) {
      message(sprintf("Skipping: k=%s q=%s p=%s", k_val, q_val, p_val))
      next
    }
    
    # Define sample sizes (NOW ON X-AXIS)
    all_sizes <- sort(unique(df_sub$sample_size))
    
    df_sub <- df_sub %>%
      mutate(
        group_label = paste(alg, divergence),
        sample_size = factor(sample_size, levels = all_sizes)
      )
    
    bhc_sub <- bhc_sub %>%
      mutate(
        group_label = "BHC (reference)",
        sample_size = factor(sample_size, levels = all_sizes)
      )
    
    # Define color map
    legend_labels <- unique(c(df_sub$group_label, "BHC (reference)"))
    dark_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(legend_labels) - 1)
    color_map <- c(
      setNames(dark_colors, legend_labels[legend_labels != "BHC (reference)"]),
      "BHC (reference)" = "black"
    )
    
    # ---- Plot ----
    plot <- ggplot(df_sub, aes(
      x = sample_size,
      y = deltar,
      color = group_label,
      group = group_label
    )) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 1.6) +
      
      # BHC reference inside each method
      geom_line(
        data = bhc_sub,
        aes(y = deltar, linetype = group_label),
        color = "black",
        linewidth = 1
      ) +
      geom_point(
        data = bhc_sub,
        aes(y = deltar),
        color = "black",
        size = 1.6
      ) +
      
      scale_color_manual(values = color_map) +
      scale_linetype_manual(values = c("BHC (reference)" = "dashed")) +
      
      facet_grid(
        rows = vars(variable),
        cols = vars(method),
        scales = "free",
        labeller = labeller(
          variable = c(
            BIC = "Delta[BIC]",
            hamming = "Delta[HD]"
          ),
          .default = label_parsed
        )
      ) +
      
      theme_bw(base_size = 13) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text  = element_text(size = 9),
        strip.text   = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(size = 0.3, linetype = "dotted"),
        axis.text.x  = element_text(angle = 35, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title   = element_blank()
      ) +
      labs(color = "", linetype = "")
    
    # ---- Save ----
    filename <- paste0(
      "_k_", ifelse(is.na(k_val), "NA", k_val),
      "_q_", ifelse(is.na(q_val), "NA", q_val),
      "_p_", p_val, ".pdf"
    )
    
    filepath <- file.path(output_dir, filename)
    
    if (save_plot) {
      ggsave(filepath, plot = plot, width = 10, height = 6)
      message("Saved: ", filepath)
    }
    
    all_plots[[length(all_plots) + 1]] <- plot
  }
  
  return(all_plots)
  
}
plots_list <- generate_deltar_plots(DDA_fixed,
                                    divergence_filter = c("jensenshannon", "hellinger", "jeffreys",
                                                          "kaniadakis","totvar"),
                                    save_plot = TRUE)

generate_deltar_plots <- function(DDA1, 
                                  sample_sizes = NULL, 
                                  divergence_filter = NULL,
                                  save_plot = TRUE, 
                                  output_dir = "Plots_Deltar_with_BHC_inside_methods_sample") {
  all_plots <- list()
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Optional filter by sample_size
  if (!is.null(sample_sizes)) {
    DDA1 <- DDA1 %>% dplyr::filter(sample_size %in% sample_sizes)
  }
  
  # ---- Remove unwanted panels ----
  DDA1 <- DDA1 %>%
    dplyr::filter(!method %in% c("bhc_ref", "full_ref"))
  
  # Split data
  bhc_data <- DDA1 %>% dplyr::filter(alg == "bhc")
  DDA1_nonbhc <- DDA1 %>% dplyr::filter(alg != "bhc")
  
  # Unique parameter combinations (NOW USES sample_size instead of p)
  unique_combos <- unique(DDA1_nonbhc[, c("k", "q", "sample_size")])
  
  for (i in seq_len(nrow(unique_combos))) {
    k_val      <- unique_combos$k[i]
    q_val      <- unique_combos$q[i]
    sample_val <- unique_combos$sample_size[i]
    
    # Filter for this combination
    df_sub <- DDA1_nonbhc %>%
      dplyr::filter(
        ((is.na(k_val) & is.na(k)) | k == k_val),
        ((is.na(q_val) & is.na(q)) | q == q_val),
        sample_size == sample_val
      )
    
    bhc_sub <- bhc_data %>%
      dplyr::filter(
        ((is.na(q_val) & is.na(q)) | q == q_val),
        sample_size == sample_val
      )
    
    if (!is.null(divergence_filter)) {
      df_sub <- df_sub %>% dplyr::filter(divergence %in% divergence_filter)
    }
    
    if (nrow(df_sub) == 0 | nrow(bhc_sub) == 0) {
      message(sprintf("Skipping: k=%s q=%s sample_size=%s", k_val, q_val, sample_val))
      next
    }
    
    # Put sample_size on x-axis (as ordered factor for equal spacing)
    all_sizes <- sort(unique(df_sub$sample_size))
    
    df_sub <- df_sub %>%
      dplyr::mutate(
        group_label = paste(alg, divergence),
        sample_size = factor(sample_size, levels = all_sizes, ordered = TRUE)
      )
    
    bhc_sub <- bhc_sub %>%
      dplyr::mutate(
        group_label = "BHC (reference)",
        sample_size = factor(sample_size, levels = all_sizes, ordered = TRUE)
      )
    
    # Color map
    legend_labels <- unique(c(df_sub$group_label, "BHC (reference)"))
    dark_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(legend_labels) - 1)
    color_map <- c(
      setNames(dark_colors, legend_labels[legend_labels != "BHC (reference)"]),
      "BHC (reference)" = "black"
    )
    
    # ---- Plot ----
    plot <- ggplot2::ggplot(df_sub, ggplot2::aes(
      x = p,
      y = deltar,
      color = group_label,
      group = group_label
    )) +
      ggplot2::geom_line(linewidth = 1.0) +
      ggplot2::geom_point(size = 1.6) +
      
      # Add BHC reference
      ggplot2::geom_line(
        data = bhc_sub,
        ggplot2::aes(y = deltar, linetype = group_label, group = 1),
        color = "black",
        linewidth = 1
      ) +
      ggplot2::geom_point(
        data = bhc_sub,
        ggplot2::aes(y = deltar),
        color = "black",
        size = 1.6
      ) +
      
      ggplot2::scale_color_manual(values = color_map) +
      ggplot2::scale_linetype_manual(values = c("BHC (reference)" = "dashed")) +
      
      ggplot2::facet_grid(
        rows = ggplot2::vars(variable),
        cols = ggplot2::vars(method),
        scales = "free",
        labeller = ggplot2::labeller(
          variable = c(
            BIC     = expression(Delta[BIC]),
            hamming = expression(Delta[HD])
          ),
          .default = ggplot2::label_parsed
        )
      ) +
      
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.title = ggplot2::element_blank(),
        legend.text  = ggplot2::element_text(size = 9),
        strip.text   = ggplot2::element_text(size = 13, face = "bold"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_line(size = 0.3, linetype = "dotted"),
        axis.text.x  = ggplot2::element_text(angle = 35, hjust = 1),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        plot.title   = ggplot2::element_blank()
      ) +
      ggplot2::labs(color = "", linetype = "")
    
    # ---- Save ----
    filename <- paste0(
      "_k_", ifelse(is.na(k_val), "NA", k_val),
      "_q_", ifelse(is.na(q_val), "NA", q_val),
      "_sample_", sample_val, ".pdf"
    )
    
    filepath <- file.path(output_dir, filename)
    
    if (save_plot) {
      ggplot2::ggsave(filepath, plot = plot, width = 10, height = 6)
      message("Saved: ", filepath)
    }
    
    all_plots[[length(all_plots) + 1]] <- plot
  }
  
  return(all_plots)
}

plots_selected <- generate_deltar_plots(
  DDA1 = DDA_fixed,
  sample_sizes = c(128, 512, 2048, 8192),   # choose sample sizes you want
  divergence_filter = c("jensenshannon", "hellinger", "jeffreys",
                        "kaniadakis","totvar"),  # divergences to plot
  save_plot = TRUE,
  output_dir = "Plots_Deltar_with_BHC_example_run"
)

###TIME
generate_deltar_plots <- function(DDA1, 
                                  sample_sizes = NULL,
                                  p_values = NULL,
                                  divergence_filter = NULL,
                                  save_plot = TRUE, 
                                  output_dir = "Plots_Deltar_p_on_x_by_sample") {
  
  all_plots <- list()
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Optional filters
  if (!is.null(sample_sizes)) {
    DDA1 <- DDA1 %>% dplyr::filter(sample_size %in% sample_sizes)
  }
  if (!is.null(p_values)) {
    DDA1 <- DDA1 %>% dplyr::filter(p %in% p_values)
  }
  
  # Remove unwanted panels
  DDA1 <- DDA1 %>% dplyr::filter(!method %in% c("bhc_ref", "full_ref"))
  
  # Split BHC and non-BHC
  bhc_data <- DDA1 %>% dplyr::filter(alg == "bhc")
  nonbhc   <- DDA1 %>% dplyr::filter(alg != "bhc")
  
  # ✅ Loop over (k, q, sample_size) so each plot has multiple p values
  unique_combos <- unique(nonbhc[, c("k", "q", "sample_size")])
  
  for (i in seq_len(nrow(unique_combos))) {
    k_val      <- unique_combos$k[i]
    q_val      <- unique_combos$q[i]
    sample_val <- unique_combos$sample_size[i]
    
    df_sub <- nonbhc %>%
      dplyr::filter(
        ((is.na(k_val) & is.na(k)) | k == k_val),
        ((is.na(q_val) & is.na(q)) | q == q_val),
        sample_size == sample_val
      )
    
    bhc_sub <- bhc_data %>%
      dplyr::filter(
        ((is.na(q_val) & is.na(q)) | q == q_val),
        sample_size == sample_val
      )
    
    # Apply divergence filter to non-BHC only
    if (!is.null(divergence_filter)) {
      df_sub <- df_sub %>% dplyr::filter(divergence %in% divergence_filter)
    }
    
    if (nrow(df_sub) == 0 | nrow(bhc_sub) == 0) {
      message(sprintf("Skipping: k=%s q=%s sample_size=%s", k_val, q_val, sample_val))
      next
    }
    
    # ✅ x-axis = p with equal spacing
    p_levels <- sort(unique(df_sub$p))
    
    df_sub <- df_sub %>%
      dplyr::mutate(
        p = factor(p, levels = p_levels, ordered = TRUE),
        group_label = paste(alg, divergence),
        # ✅ Use raw median for time, deltar for others
        y_plot = dplyr::if_else(variable == "time", median, deltar)
      )
    
    bhc_sub <- bhc_sub %>%
      dplyr::mutate(
        p = factor(p, levels = p_levels, ordered = TRUE),
        group_label = "BHC (reference)",
        y_plot = dplyr::if_else(variable == "time", median, deltar)
      )
    
    # Color map
    legend_labels <- unique(c(df_sub$group_label, "BHC (reference)"))
    dark_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(legend_labels) - 1)
    color_map <- c(
      setNames(dark_colors, legend_labels[legend_labels != "BHC (reference)"]),
      "BHC (reference)" = "black"
    )
    
    # ---- Plot ----
    plot <- ggplot2::ggplot(df_sub, ggplot2::aes(
      x = p, y = y_plot, color = group_label, group = group_label
    )) +
      ggplot2::geom_line(linewidth = 1.0) +
      ggplot2::geom_point(size = 1.6) +
      
      # BHC reference line
      ggplot2::geom_line(
        data = bhc_sub,
        ggplot2::aes(x = p, y = y_plot, linetype = group_label, group = 1),
        color = "black", linewidth = 1
      ) +
      ggplot2::geom_point(
        data = bhc_sub,
        ggplot2::aes(x = p, y = y_plot),
        color = "black", size = 1.6
      ) +
      
      ggplot2::scale_color_manual(values = color_map) +
      ggplot2::scale_linetype_manual(values = c("BHC (reference)" = "dashed")) +
      
      ggplot2::facet_grid(
        rows = ggplot2::vars(variable),
        cols = ggplot2::vars(method),
        scales = "free_y",
        labeller = ggplot2::labeller(
          variable = c(
            BIC = "Delta[BIC]",
            hamming = "Delta[HD]",
            time = "Time"
          ),
          .default = ggplot2::label_parsed
        )
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.title = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_blank()
      ) +
      ggplot2::labs(color = "", linetype = "")
    
    # Save
    if (save_plot) {
      filename <- paste0(
        "_k_", ifelse(is.na(k_val), "NA", k_val),
        "_q_", ifelse(is.na(q_val), "NA", q_val),
        "_sample_", sample_val, ".pdf"
      )
      ggplot2::ggsave(file.path(output_dir, filename), plot = plot, width = 10, height = 7)
    }
    
    all_plots[[length(all_plots) + 1]] <- plot
  }
  
  return(all_plots)
}

plots_selected <- generate_deltar_plots(
  DDA1 = DDA_fixed,
  sample_sizes = c(128, 512, 2048, 8192),   # choose sample sizes you want
  divergence_filter = c("jensenshannon", "hellinger", "jeffreys",
                        "kaniadakis","totvar"),  # divergences to plot
  save_plot = TRUE,
  output_dir = "Plots_Deltar_with_BHC_example_run"
)
#Ktrue
generate_deltar_plots <- function(DDA1, 
                                  sample_sizes = NULL,
                                  p_values = NULL,
                                  divergence_filter = NULL,
                                  save_plot = TRUE, 
                                  output_dir = "Plots_Deltar_p_on_x_by_sample_ktrue") {
  
  all_plots <- list()
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Optional filters
  if (!is.null(sample_sizes)) {
    DDA1 <- DDA1 %>% dplyr::filter(sample_size %in% sample_sizes)
  }
  if (!is.null(p_values)) {
    DDA1 <- DDA1 %>% dplyr::filter(p %in% p_values)
  }
  
  # Remove unwanted panels
  DDA1 <- DDA1 %>% dplyr::filter(!method %in% c("bhc_ref", "full_ref"))
  
  # Split BHC and non-BHC
  bhc_data <- DDA1 %>% dplyr::filter(alg == "bhc")
  nonbhc   <- DDA1 %>% dplyr::filter(alg != "bhc")
  
  # ✅ Loop over (k, ktrue, sample_size)
  unique_combos <- unique(nonbhc[, c("k", "ktrue", "sample_size")])
  
  for (i in seq_len(nrow(unique_combos))) {
    k_val      <- unique_combos$k[i]
    ktrue_val  <- unique_combos$ktrue[i]
    sample_val <- unique_combos$sample_size[i]
    
    df_sub <- nonbhc %>%
      dplyr::filter(
        ((is.na(k_val) & is.na(k)) | k == k_val),
        ((is.na(ktrue_val) & is.na(ktrue)) | ktrue == ktrue_val),
        sample_size == sample_val
      )
    
    bhc_sub <- bhc_data %>%
      dplyr::filter(
        ((is.na(ktrue_val) & is.na(ktrue)) | ktrue == ktrue_val),
        sample_size == sample_val
      )
    
    # Apply divergence filter to non-BHC only
    if (!is.null(divergence_filter)) {
      df_sub <- df_sub %>% dplyr::filter(divergence %in% divergence_filter)
    }
    
    if (nrow(df_sub) == 0 | nrow(bhc_sub) == 0) {
      message(sprintf("Skipping: k=%s ktrue=%s sample_size=%s", k_val, ktrue_val, sample_val))
      next
    }
    
    # ✅ x-axis = p with equal spacing
    p_levels <- sort(unique(df_sub$p))
    
    df_sub <- df_sub %>%
      dplyr::mutate(
        p = factor(p, levels = p_levels, ordered = TRUE),
        group_label = paste(alg, divergence),
        # ✅ raw median for time, deltar for others
        y_plot = dplyr::if_else(variable == "time", median, deltarf)
      )
    
    bhc_sub <- bhc_sub %>%
      dplyr::mutate(
        p = factor(p, levels = p_levels, ordered = TRUE),
        group_label = "BHC (reference)",
        y_plot = dplyr::if_else(variable == "time", median, deltarf)
      )
    
    # Color map
    legend_labels <- unique(c(df_sub$group_label, "BHC (reference)"))
    dark_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(legend_labels) - 1)
    color_map <- c(
      setNames(dark_colors, legend_labels[legend_labels != "BHC (reference)"]),
      "BHC (reference)" = "black"
    )
    
    # ---- Plot ----
    plot <- ggplot2::ggplot(df_sub, ggplot2::aes(
      x = p, y = y_plot, color = group_label, group = group_label
    )) +
      ggplot2::geom_line(linewidth = 1.0) +
      ggplot2::geom_point(size = 1.6) +
      
      # BHC reference line
      ggplot2::geom_line(
        data = bhc_sub,
        ggplot2::aes(x = p, y = y_plot, linetype = group_label, group = 1),
        color = "black", linewidth = 1
      ) +
      ggplot2::geom_point(
        data = bhc_sub,
        ggplot2::aes(x = p, y = y_plot),
        color = "black", size = 1.6
      ) +
      
      ggplot2::scale_color_manual(values = color_map) +
      ggplot2::scale_linetype_manual(values = c("BHC (reference)" = "dashed")) +
      
      ggplot2::facet_grid(
        rows = ggplot2::vars(variable),
        cols = ggplot2::vars(method),
        scales = "free_y",
        labeller = ggplot2::labeller(
          variable = c(
            BIC = "Delta[BIC]",
            hamming = "Delta[HD]",
            time = "Time"
          ),
          .default = ggplot2::label_parsed
        )
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.title = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_blank()
      ) +
      ggplot2::labs(color = "", linetype = "")
    
    # Save
    if (save_plot) {
      filename <- paste0(
        "_k_", ifelse(is.na(k_val), "NA", k_val),
        "_ktrue_", ifelse(is.na(ktrue_val), "NA", ktrue_val),
        "_sample_", sample_val, ".pdf"
      )
      ggplot2::ggsave(file.path(output_dir, filename), plot = plot, width = 10, height = 7)
    }
    
    all_plots[[length(all_plots) + 1]] <- plot
  }
  
  return(all_plots)
}


plots_selected[1] <- generate_deltar_plots(
  DDA1 = DDA_fixed,
  sample_sizes = c(128, 512, 2048, 8192),   # choose sample sizes you want
  divergence_filter = c("jensenshannon", "hellinger", "jeffreys",
                        "kaniadakis","totvar"),  # divergences to plot
  save_plot = TRUE,
  output_dir = "Plots_Deltar_with_BHC_example_run_1"
)













#with ktrue
generate_deltar_plots <- function(DDA1, 
                                  p_values = NULL, 
                                  divergence_filter = NULL,
                                  save_plot = TRUE, 
                                  output_dir = "Plots_Deltar_with_BHC_inside_methods_p") {
  all_plots <- list()
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Optional filter by p
  if (!is.null(p_values)) {
    DDA1 <- DDA1 %>% filter(p %in% p_values)
  }
  
  # ---- Remove unwanted panels ----
  DDA1 <- DDA1 %>%
    filter(!method %in% c("bhc_ref", "full_ref"))
  
  # Split data
  bhc_data <- DDA1 %>% filter(alg == "bhc")
  DDA1 <- DDA1 %>% filter(alg != "bhc")
  
  # Unique parameter combinations (NOW USES ktrue instead of q)
  unique_combos <- unique(DDA1[, c("k", "ktrue", "p")])
  
  for (i in seq_len(nrow(unique_combos))) {
    k_val <- unique_combos$k[i]
    ktrue_val <- unique_combos$ktrue[i]
    p_val <- unique_combos$p[i]
    
    # Filter for this combination
    df_sub <- DDA1 %>%
      filter(
        ((is.na(k_val) & is.na(k)) | k == k_val),
        ((is.na(ktrue_val) & is.na(ktrue)) | ktrue == ktrue_val),
        p == p_val
      )
    
    bhc_sub <- bhc_data %>%
      filter(p == p_val, ktrue == ktrue_val)
    
    if (!is.null(divergence_filter)) {
      df_sub <- df_sub %>% filter(divergence %in% divergence_filter)
    }
    
    if (nrow(df_sub) == 0 | nrow(bhc_sub) == 0) {
      message(sprintf("Skipping: k=%s ktrue=%s p=%s", k_val, ktrue_val, p_val))
      next
    }
    
    # Define sample sizes (NOW ON X-AXIS)
    all_sizes <- sort(unique(df_sub$sample_size))
    
    df_sub <- df_sub %>%
      mutate(
        group_label = paste(alg, divergence),
        sample_size = factor(sample_size, levels = all_sizes)
      )
    
    bhc_sub <- bhc_sub %>%
      mutate(
        group_label = "BHC (reference)",
        sample_size = factor(sample_size, levels = all_sizes)
      )
    
    # Define color map
    legend_labels <- unique(c(df_sub$group_label, "BHC (reference)"))
    dark_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(legend_labels) - 1)
    color_map <- c(
      setNames(dark_colors, legend_labels[legend_labels != "BHC (reference)"]),
      "BHC (reference)" = "black"
    )
    
    # ---- Plot ----
    plot <- ggplot(df_sub, aes(
      x = sample_size,
      y = deltar,
      color = group_label,
      group = group_label
    )) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 1.6) +
      
      # BHC reference inside each method
      geom_line(
        data = bhc_sub,
        aes(y = deltar, linetype = group_label),
        color = "black",
        linewidth = 1
      ) +
      geom_point(
        data = bhc_sub,
        aes(y = deltar),
        color = "black",
        size = 1.6
      ) +
      
      scale_color_manual(values = color_map) +
      scale_linetype_manual(values = c("BHC (reference)" = "dashed")) +
      
      facet_grid(
        rows = vars(variable),
        cols = vars(method),
        scales = "free",
        labeller = labeller(
          variable = c(
            BIC = "Delta[BIC]",
            hamming = "Delta[HD]"
          ),
          .default = label_parsed
        )
      ) +
      
      theme_bw(base_size = 13) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text  = element_text(size = 9),
        strip.text   = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(size = 0.3, linetype = "dotted"),
        axis.text.x  = element_text(angle = 35, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title   = element_blank()
      ) +
      labs(color = "", linetype = "")
    
    # ---- Save ----
    filename <- paste0(
      "_k_", ifelse(is.na(k_val), "NA", k_val),
      "_ktrue_", ifelse(is.na(ktrue_val), "NA", ktrue_val),
      "_p_", p_val, ".pdf"
    )
    
    filepath <- file.path(output_dir, filename)
    
    if (save_plot) {
      ggsave(filepath, plot = plot, width = 10, height = 6)
      message("Saved: ", filepath)
    }
    
    all_plots[[length(all_plots) + 1]] <- plot
  }
  
  return(all_plots)
}





#without bhc with q
generate_deltar_plots <- function(DDA1, 
                                  p_values = NULL, 
                                  divergence_filter = NULL,
                                  save_plot = TRUE, 
                                  output_dir = "Plots_Deltar_with_BHC_ktrue") {
  all_plots <- list()
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Optional filter by p
  if (!is.null(p_values)) {
    DDA1 <- DDA1 %>% filter(p %in% p_values)
  }
  
  # ---- Remove unwanted panels and BHC completely ----
  DDA1 <- DDA1 %>%
    filter(
      alg != "bhc",
      !method %in% c("bhc_ref", "full_ref")
    )
  
  # Unique parameter combinations
  unique_combos <- unique(DDA1[, c("k", "q", "p")])
  
  for (i in seq_len(nrow(unique_combos))) {
    k_val <- unique_combos$k[i]
    q_val <- unique_combos$q[i]
    p_val <- unique_combos$p[i]
    
    # Filter for this combination
    df_sub <- DDA1 %>%
      filter(
        ((is.na(k_val) & is.na(k)) | k == k_val),
        ((is.na(q_val) & is.na(q)) | q == q_val),
        p == p_val
      )
    
    if (!is.null(divergence_filter)) {
      df_sub <- df_sub %>% filter(divergence %in% divergence_filter)
    }
    
    if (nrow(df_sub) == 0) {
      message(sprintf("Skipping: k=%s q=%s p=%s", k_val, q_val, p_val))
      next
    }
    
    # Sample sizes on x-axis
    all_sizes <- sort(unique(df_sub$sample_size))
    
    df_sub <- df_sub %>%
      mutate(
        group_label = paste(alg, divergence),
        sample_size = factor(sample_size, levels = all_sizes)
      )
    
    # Color map (divergences only)
    legend_labels <- unique(df_sub$group_label)
    dark_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(legend_labels))
    color_map <- setNames(dark_colors, legend_labels)
    
    # ---- Plot ----
    plot <- ggplot(df_sub, aes(
      x = sample_size,
      y = deltar,
      color = group_label,
      group = group_label
    )) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 1.6) +
      
      scale_color_manual(values = color_map) +
      
      facet_grid(
        rows = vars(variable),
        cols = vars(method),
        scales = "free",
        labeller = labeller(
          variable = c(
            BIC = "Delta[BIC]",
            hamming = "Delta[HD]"
          ),
          .default = label_parsed
        )
      ) +
      
      theme_bw(base_size = 13) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text  = element_text(size = 9),
        strip.text   = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(size = 0.3, linetype = "dotted"),
        axis.text.x  = element_text(angle = 35, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title   = element_blank()
      ) +
      labs(color = "")
    
    # ---- Save ----
    filename <- paste0(
      "_k_", ifelse(is.na(k_val), "NA", k_val),
      "_q_", ifelse(is.na(q_val), "NA", q_val),
      "_p_", p_val, ".pdf"
    )
    
    filepath <- file.path(output_dir, filename)
    
    if (save_plot) {
      ggsave(filepath, plot = plot, width = 8, height = 6)
      message("Saved: ", filepath)
    }
    
    all_plots[[length(all_plots) + 1]] <- plot
  }
  
  return(all_plots)
}


# without bhc with ktrue
# without bhc, using ktrue instead of q
generate_deltar_plots <- function(DDA1, 
                                  p_values = NULL, 
                                  divergence_filter = NULL,
                                  save_plot = TRUE, 
                                  output_dir = "Plots_Deltar_with_BHC_ktrue") {
  all_plots <- list()
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Optional filter by p
  if (!is.null(p_values)) {
    DDA1 <- DDA1 %>% filter(p %in% p_values)
  }
  
  # ---- Remove unwanted panels and BHC completely ----
  DDA1 <- DDA1 %>%
    filter(
      alg != "bhc",
      !method %in% c("bhc_ref", "full_ref")
    )
  
  # Unique parameter combinations (USES ktrue)
  unique_combos <- unique(DDA1[, c("k", "ktrue", "p")])
  
  for (i in seq_len(nrow(unique_combos))) {
    k_val     <- unique_combos$k[i]
    ktrue_val <- unique_combos$ktrue[i]
    p_val     <- unique_combos$p[i]
    
    # Filter for this combination
    df_sub <- DDA1 %>%
      filter(
        ((is.na(k_val) & is.na(k)) | k == k_val),
        ((is.na(ktrue_val) & is.na(ktrue)) | ktrue == ktrue_val),
        p == p_val
      )
    
    if (!is.null(divergence_filter)) {
      df_sub <- df_sub %>% filter(divergence %in% divergence_filter)
    }
    
    if (nrow(df_sub) == 0) {
      message(sprintf("Skipping: k=%s ktrue=%s p=%s",
                      k_val, ktrue_val, p_val))
      next
    }
    
    # Sample sizes on x-axis
    all_sizes <- sort(unique(df_sub$sample_size))
    
    df_sub <- df_sub %>%
      mutate(
        group_label = paste(alg, divergence),
        sample_size = factor(sample_size, levels = all_sizes)
      )
    
    # Color map
    legend_labels <- unique(df_sub$group_label)
    dark_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(legend_labels))
    color_map <- setNames(dark_colors, legend_labels)
    
    # ---- Plot ----
    plot <- ggplot(df_sub, aes(
      x = sample_size,
      y = deltar,
      color = group_label,
      group = group_label
    )) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 1.6) +
      
      scale_color_manual(values = color_map) +
      
      facet_grid(
        rows = vars(variable),
        cols = vars(method),
        scales = "free",
        labeller = labeller(
          variable = c(
            BIC     = "Delta[BIC]",
            hamming = "Delta[HD]"
          ),
          .default = label_parsed
        )
      ) +
      
      theme_bw(base_size = 13) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text  = element_text(size = 9),
        strip.text   = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(size = 0.3, linetype = "dotted"),
        axis.text.x  = element_text(angle = 35, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title   = element_blank()
      ) +
      labs(color = "")
    
    # ---- Save ----
    filename <- paste0(
      "_k_", ifelse(is.na(k_val), "NA", k_val),
      "_ktrue_", ifelse(is.na(ktrue_val), "NA", ktrue_val),
      "_p_", p_val, ".pdf"
    )
    
    filepath <- file.path(output_dir, filename)
    
    if (save_plot) {
      ggsave(filepath, plot = plot, width = 8, height = 6)
      message("Saved: ", filepath)
    }
    
    all_plots[[length(all_plots) + 1]] <- plot
  }
  
  return(all_plots)
}


# Exact Time
generate_deltar_plots <- function(DDA1, 
                                  sample_sizes = NULL, 
                                  divergence_filter = NULL,
                                  save_plot = TRUE, 
                                  output_dir = "Plots_Deltar_p_axis_noBHC") {
  
  all_plots <- list()
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ---- Filter sample sizes ----
  if (!is.null(sample_sizes)) {
    DDA1 <- DDA1 %>% 
      dplyr::filter(sample_size %in% sample_sizes)
  }
  
  # ---- Keep only p = 5,7,9,11 ----
  DDA1 <- DDA1 %>%
    dplyr::filter(p %in% c(5, 7, 9, 11))
  
  # ---- Remove BHC and reference panels ----
  DDA1 <- DDA1 %>%
    dplyr::filter(
      alg != "bhc",
      !method %in% c("bhc_ref", "full_ref")
    )
  
  # ---- Unique parameter combinations ----
  unique_combos <- unique(DDA1[, c("k", "q", "sample_size")])
  
  for (i in seq_len(nrow(unique_combos))) {
    
    k_val      <- unique_combos$k[i]
    q_val      <- unique_combos$q[i]
    sample_val <- unique_combos$sample_size[i]
    
    # ---- Subset data ----
    df_sub <- DDA1 %>%
      dplyr::filter(
        ((is.na(k_val) & is.na(k)) | k == k_val),
        ((is.na(q_val) & is.na(q)) | q == q_val),
        sample_size == sample_val
      )
    
    if (!is.null(divergence_filter)) {
      df_sub <- df_sub %>%
        dplyr::filter(divergence %in% divergence_filter)
    }
    
    if (nrow(df_sub) == 0) {
      message(sprintf(
        "Skipping: k=%s q=%s sample_size=%s",
        k_val, q_val, sample_val
      ))
      next
    }
    
    # ---- X-axis: p ----
    p_levels <- sort(unique(df_sub$p))
    
    # ---- Define plotted y-values ----
    df_sub <- df_sub %>%
      dplyr::mutate(
        y_plot = dplyr::case_when(
          variable == "time"    ~ median,   # absolute time
          TRUE                  ~ deltar    # ΔBIC, ΔHD
        ),
        group_label = paste(alg, divergence),
        p = factor(p, levels = p_levels)
      )
    
    # ---- Colors ----
    legend_labels <- unique(df_sub$group_label)
    colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(
      length(legend_labels)
    )
    color_map <- setNames(colors, legend_labels)
    
    # ---- Plot ----
    plot <- ggplot(df_sub, aes(
      x = p,
      y = y_plot,
      color = group_label,
      group = group_label
    )) +
      geom_line(linewidth = 1) +
      geom_point(size = 1.6) +
      scale_color_manual(values = color_map) +
      facet_grid(
        rows = vars(variable),
        cols = vars(method),
        scales = "free",
        labeller = labeller(
          variable = c(
            BIC     = "Delta[BIC]",
            hamming = "Delta[HD]",
            time    = "Time"
          ),
          .default = label_parsed
        )
      ) +
      theme_bw(base_size = 13) +
      theme(
        legend.position   = "bottom",
        legend.title      = element_blank(),
        legend.text       = element_text(size = 9),
        strip.text        = element_text(size = 13, face = "bold"),
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_line(
          size = 0.3, linetype = "dotted"
        ),
        axis.text.x       = element_text(angle = 35, hjust = 1),
        axis.title.x      = element_blank(),
        axis.title.y      = element_blank(),
        plot.title        = element_blank()
      ) +
      labs(color = "")
    
    # ---- Save ----
    filename <- paste0(
      "_k_", ifelse(is.na(k_val), "NA", k_val),
      "_q_", ifelse(is.na(q_val), "NA", q_val),
      "_n_", sample_val,
      ".pdf"
    )
    
    filepath <- file.path(output_dir, filename)
    
    if (save_plot) {
      ggsave(filepath, plot = plot, width = 8, height = 6)
      message("Saved: ", filepath)
    }
    
    all_plots[[length(all_plots) + 1]] <- plot
  }
  
  return(all_plots)
}

plot<-generate_deltar_plots (DDA, 
                                  sample_sizes = c(128,512,2048,2192), 
                                  divergence_filter = c("jensenshannon", "hellinger", "jeffreys",
                                                        "kaniadakis","totvar"),  # divergences to plot
                                  save_plot = FALSE, 
                                  output_dir = "Plots_Deltar_p_axis_noBHC")




# Time with ktrue
generate_deltar_plots <- function(DDA1, 
                                  sample_sizes = NULL, 
                                  divergence_filter = NULL,
                                  save_plot = TRUE, 
                                  output_dir = "Plots_Deltar_p_axis_noBHC_ktrue") {
  
  all_plots <- list()
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ---- Filter sample sizes ----
  if (!is.null(sample_sizes)) {
    DDA1 <- DDA1 %>% 
      dplyr::filter(sample_size %in% sample_sizes)
  }
  
  # ---- Keep only p = 5,7,9,11 ----
  DDA1 <- DDA1 %>%
    dplyr::filter(p %in% c(5, 7, 9, 11,13,15))
  
  # ---- Remove BHC and reference panels ----
  DDA1 <- DDA1 %>%
    dplyr::filter(
      alg != "bhc",
      !method %in% c("bhc_ref", "full_ref")
    )
  
  # ---- Unique parameter combinations (USES ktrue) ----
  unique_combos <- unique(DDA1[, c("k", "ktrue", "sample_size")])
  
  for (i in seq_len(nrow(unique_combos))) {
    
    k_val      <- unique_combos$k[i]
    ktrue_val  <- unique_combos$ktrue[i]
    sample_val <- unique_combos$sample_size[i]
    
    # ---- Subset data ----
    df_sub <- DDA1 %>%
      dplyr::filter(
        ((is.na(k_val) & is.na(k)) | k == k_val),
        ((is.na(ktrue_val) & is.na(ktrue)) | ktrue == ktrue_val),
        sample_size == sample_val
      )
    
    if (!is.null(divergence_filter)) {
      df_sub <- df_sub %>%
        dplyr::filter(divergence %in% divergence_filter)
    }
    
    if (nrow(df_sub) == 0) {
      message(sprintf(
        "Skipping: k=%s ktrue=%s sample_size=%s",
        k_val, ktrue_val, sample_val
      ))
      next
    }
    
    # ---- X-axis: p ----
    p_levels <- sort(unique(df_sub$p))
    
    # ---- Define plotted y-values ----
    df_sub <- df_sub %>%
      dplyr::mutate(
        y_plot = dplyr::case_when(
          variable == "time" ~ median,   # absolute time
          TRUE               ~ deltarf    # ΔBIC, ΔHD
        ),
        group_label = paste(alg, divergence),
        p = factor(p, levels = p_levels)
      )
    
    # ---- Colors ----
    legend_labels <- unique(df_sub$group_label)
    colors <- colorRampPalette(
      RColorBrewer::brewer.pal(8, "Dark2")
    )(length(legend_labels))
    
    color_map <- setNames(colors, legend_labels)
    
    # ---- Plot ----
    plot <- ggplot(df_sub, aes(
      x = p,
      y = y_plot,
      color = group_label,
      group = group_label
    )) +
      geom_line(linewidth = 1) +
      geom_point(size = 1.6) +
      scale_color_manual(values = color_map) +
      facet_grid(
        rows = vars(variable),
        cols = vars(method),
        scales = "free",
        labeller = labeller(
          variable =  c(
            BIC     = "Delta[BIC]",
            hamming = "Delta[HD]",
            time    = "Time"
          ),
          .default = label_parsed
        )
      ) +
      theme_bw(base_size = 13) +
      theme(
        legend.position    = "bottom",
        legend.title       = element_blank(),
        legend.text        = element_text(size = 9),
        strip.text         = element_text(size = 13, face = "bold"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(size = 0.3, linetype = "dotted"),
        axis.text.x        = element_text(angle = 35, hjust = 1),
        axis.title.x       = element_blank(),
        axis.title.y       = element_blank(),
        plot.title         = element_blank()
      ) +
      labs(color = "")
    
    # ---- Save ----
    filename <- paste0(
      "_k_", ifelse(is.na(k_val), "NA", k_val),
      "_ktrue_", ifelse(is.na(ktrue_val), "NA", ktrue_val),
      "_n_", sample_val,
      ".pdf"
    )
    
    filepath <- file.path(output_dir, filename)
    
    if (save_plot) {
      ggsave(filepath, plot = plot, width = 8, height = 6)
      message("Saved: ", filepath)
    }
    
    all_plots[[length(all_plots) + 1]] <- plot
  }
  
  return(all_plots)
}

p<-generate_deltar_plots(DDA_fixed, 
                       sample_sizes = c(128,512,2048,8192), 
                      divergence_filter = c("jensenshannon", "hellinger", "jeffreys",
                                            "kaniadakis","totvar"),
                      save_plot = TRUE, 
                      output_dir = "Plots_Deltar_p_axis_noBHC_ktru_!") 




# Time with q
generate_deltar_plots <- function(DDA1, 
                                  sample_sizes = NULL, 
                                  divergence_filter = NULL,
                                  save_plot = TRUE, 
                                  output_dir = "Plots_Deltar_p_axis_noBHC_q") {
  
  all_plots <- list()
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ---- Filter sample sizes ----
  if (!is.null(sample_sizes)) {
    DDA1 <- DDA1 %>% 
      dplyr::filter(sample_size %in% sample_sizes)
  }
  
  # ---- Keep only p = 5,7,9,11 ----
  DDA1 <- DDA1 %>%
    dplyr::filter(p %in% c(5, 7, 9, 11, 13, 15))
  
  # ---- Remove BHC and reference panels ----
  DDA1 <- DDA1 %>%
    dplyr::filter(
      alg != "bhc",
      !method %in% c("bhc_ref", "full_ref")
    )
  
  # ---- Unique parameter combinations (USES q) ----
  unique_combos <- unique(DDA1[, c("k", "q", "sample_size")])
  
  for (i in seq_len(nrow(unique_combos))) {
    
    k_val      <- unique_combos$k[i]
    q_val      <- unique_combos$q[i]
    sample_val <- unique_combos$sample_size[i]
    
    # ---- Subset data ----
    df_sub <- DDA1 %>%
      dplyr::filter(
        ((is.na(k_val) & is.na(k)) | k == k_val),
        ((is.na(q_val) & is.na(q)) | q == q_val),
        sample_size == sample_val
      )
    
    if (!is.null(divergence_filter)) {
      df_sub <- df_sub %>%
        dplyr::filter(divergence %in% divergence_filter)
    }
    
    if (nrow(df_sub) == 0) {
      message(sprintf(
        "Skipping: k=%s q=%s sample_size=%s",
        k_val, q_val, sample_val
      ))
      next
    }
    
    # ---- X-axis: p ----
    p_levels <- sort(unique(df_sub$p))
    
    # ---- Define plotted y-values ----
    df_sub <- df_sub %>%
      dplyr::mutate(
        y_plot = dplyr::case_when(
          variable == "time" ~ median,   # absolute time
          TRUE               ~ deltarf   # ΔBIC, ΔHD
        ),
        group_label = paste(alg, divergence),
        p = factor(p, levels = p_levels)
      )
    
    # ---- Colors ----
    legend_labels <- unique(df_sub$group_label)
    colors <- colorRampPalette(
      RColorBrewer::brewer.pal(8, "Dark2")
    )(length(legend_labels))
    
    color_map <- setNames(colors, legend_labels)
    
    # ---- Plot ----
    plot <- ggplot(df_sub, aes(
      x = p,
      y = y_plot,
      color = group_label,
      group = group_label
    )) +
      geom_line(linewidth = 1) +
      geom_point(size = 1.6) +
      scale_color_manual(values = color_map) +
      facet_grid(
        rows = vars(variable),
        cols = vars(method),
        scales = "free",
        labeller = labeller(
          variable =  c(
            BIC     = "Delta[BIC]",
            hamming = "Delta[HD]",
            time    = "Time"
          ),
          .default = label_parsed
        )
      ) +
      theme_bw(base_size = 13) +
      theme(
        legend.position    = "bottom",
        legend.title       = element_blank(),
        legend.text        = element_text(size = 9),
        strip.text         = element_text(size = 13, face = "bold"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(size = 0.3, linetype = "dotted"),
        axis.text.x        = element_text(angle = 35, hjust = 1),
        axis.title.x       = element_blank(),
        axis.title.y       = element_blank(),
        plot.title         = element_blank()
      ) +
      labs(color = "")
    
    # ---- Save ----
    filename <- paste0(
      "_k_", ifelse(is.na(k_val), "NA", k_val),
      "_q_", ifelse(is.na(q_val), "NA", q_val),
      "_n_", sample_val,
      ".pdf"
    )
    
    filepath <- file.path(output_dir, filename)
    
    if (save_plot) {
      ggsave(filepath, plot = plot, width = 8, height = 6)
      message("Saved: ", filepath)
    }
    
    all_plots[[length(all_plots) + 1]] <- plot
  }
  
  return(all_plots)
}
p<-generate_deltar_plots(DDA, 
                         sample_sizes = c(128,512,2048,8192), 
                         divergence_filter = c("jensenshannon", "hellinger", "jeffreys",
                                               "kaniadakis","totvar"),
                         save_plot = TRUE, 
                         output_dir = "Plots_Deltar_p_axis_noBHC_a_!") 

  