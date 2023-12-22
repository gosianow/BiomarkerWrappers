

# biomarker_var = "biomarker"; 
# hr_prefix = "HR"; hr_ci_lower_prefix = "HR.CI95.lower"; hr_ci_upper_prefix = "HR.CI95.upper";
# pval_prefix = "P.Value"; adjp_prefix = "adj.P.Val";  
# sep = "_"; pval = 0.05; title = ""; 
# color_low = '#42399B'; color_mid = "white"; color_high = '#D70131'; 
# trim_values = c(0.25, 4); trim_prop = NULL; trim_range = NULL; ceiling = FALSE; 
# radius_range = c(3, 10); legend_position = "right"; 
# axis_text_y_size = NULL; axis_text_y_width = 80; title_size = NULL



#' Dot plot with HR and 95% CIs and p-values for multiple contrasts
#' 
#' @param x TopTable with Cox regression results 
#' @export
wrapper_HR_dotplot <- function(x, biomarker_var = "biomarker", 
  hr_prefix = "HR", hr_ci_lower_prefix = "HR.CI95.lower", hr_ci_upper_prefix = "HR.CI95.upper",
  pval_prefix = "P.Value", adjp_prefix = "adj.P.Val",  
  sep = "_", pval = 0.05, title = "", 
  color_low = '#42399B', color_mid = "white", color_high = '#D70131', 
  trim_values = c(0.25, 4), trim_prop = NULL, trim_range = NULL, ceiling = FALSE, 
  radius_range = c(3, 10), legend_position = "right", 
  axis_text_y_size = NULL, axis_text_y_width = 80, title_size = NULL){
  
  
  
  stopifnot(length(biomarker_var) == 1)
  
  
  data_hr <- wrapper_extract_from_topTable(x, extract_prefix = hr_prefix, sep = sep)
  
  data_hr_ci_lower <- wrapper_extract_from_topTable(x, extract_prefix = hr_ci_lower_prefix, sep = sep)
  
  data_hr_ci_upper <- wrapper_extract_from_topTable(x, extract_prefix = hr_ci_upper_prefix, sep = sep)
  
  data_pval <- wrapper_extract_from_topTable(x, extract_prefix = pval_prefix, sep = sep)
  
  data_adjp <- wrapper_extract_from_topTable(x, extract_prefix = adjp_prefix, sep = sep)
  
  
  contrasts <- colnames(data_hr)
  contrasts
  
  
  data_hr <- pivot_longer(data.frame(x[, biomarker_var, drop = FALSE], data_hr, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = all_of(contrasts), names_to = "contrast", values_to = hr_prefix)
  
  data_hr_ci_lower <- pivot_longer(data.frame(x[, biomarker_var, drop = FALSE], data_hr_ci_lower, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = all_of(contrasts), names_to = "contrast", values_to = hr_ci_lower_prefix)
  
  data_hr_ci_upper <- pivot_longer(data.frame(x[, biomarker_var, drop = FALSE], data_hr_ci_upper, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = all_of(contrasts), names_to = "contrast", values_to = hr_ci_upper_prefix)
  
  data_pval <- pivot_longer(data.frame(x[, biomarker_var, drop = FALSE], data_pval, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = all_of(contrasts), names_to = "contrast", values_to = pval_prefix)
  
  data_adjp <- pivot_longer(data.frame(x[, biomarker_var, drop = FALSE], data_adjp, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = all_of(contrasts), names_to = "contrast", values_to = adjp_prefix)
  
  
  if(pval_prefix == adjp_prefix){
    data <- data_hr %>% 
      dplyr::left_join(data_hr_ci_lower, by = c(biomarker_var, "contrast")) %>% 
      dplyr::left_join(data_hr_ci_upper, by = c(biomarker_var, "contrast")) %>% 
      dplyr::left_join(data_pval, by = c(biomarker_var, "contrast")) %>% 
      data.frame()
  }else{
    data <- data_hr %>% 
      dplyr::left_join(data_hr_ci_lower, by = c(biomarker_var, "contrast")) %>% 
      dplyr::left_join(data_hr_ci_upper, by = c(biomarker_var, "contrast")) %>% 
      dplyr::left_join(data_pval, by = c(biomarker_var, "contrast")) %>% 
      dplyr::left_join(data_adjp, by = c(biomarker_var, "contrast")) %>% 
      data.frame()
  }
  
  
  data$significance <- factor(ifelse(data[, adjp_prefix] <= pval, paste0("<=", pval), paste0(">", pval)), levels = paste0(c("<=", ">"), pval))
  
  values_shape <- c(4, 32)
  names(values_shape) <- levels(data$significance)
  
  
  data$contrast <- factor(data$contrast, levels = contrasts)
  
  ## Shorten the gene set names so they can be nicely displayed in the plots
  rownames_wrap <- stringr::str_wrap(data[, biomarker_var], width = axis_text_y_width)
  
  rownames_wrap <- limma::strsplit2(rownames_wrap, split = "\\\n")[, 1]
  
  data[, biomarker_var] <- factor(rownames_wrap, levels = rev(unique(rownames_wrap)))
  
  
  ### Use -log10(p-value) to define the size of the bubbles 
  
  data$log.P.Val <- -log10(data[, pval_prefix])
  
  radius_labels = rev(c(1, 0.1, 0.05, 0.01, 0.001))
  radius_breaks = -log10(radius_labels)
  radius_limits = range(radius_breaks)
  radius_labels <- formatC(radius_labels, format = "f", drop0trailing = TRUE, digits = 10)
  
  ### Squish the limits because oob = scales::squish is not possible for scale_radius
  data$log.P.Val[data$log.P.Val > max(radius_limits)] <- max(radius_limits)
  
  
  if(is.null(trim_values)){
    trim_values <- compute_trim_values(x = data[, hr_prefix], centered = FALSE, trim_prop = trim_prop, trim_range = trim_range, ceiling = ceiling)
  }else{
    max_abs_value <- max(abs(log2(trim_values)))
    trim_values <- c(2^(-max_abs_value), 2^max_abs_value)
  }
  
  limits <- trim_values

  
  # ---------------------------------------------------------------------------
  # ggplot
  # ---------------------------------------------------------------------------
  
  # This is great. squish in this context converts clamps all values to be within the min and max of the limits argument. i.e., if value < min(limits) then value = min(limits) else if value > max(limits) then value = max(limits).
  # scale_colour_gradient2(limits = c(-1.5, 1.5), oob = scales::squish)
  
  
  
  ggp <- ggplot(data, aes(x = .data[[hr_prefix]], y = .data[[biomarker_var]], size = .data[["log.P.Val"]], fill = .data[[hr_prefix]])) +
    geom_vline(xintercept = 1, color = "gray", linetype = "dashed") +
    geom_point(shape = 21) +
    geom_point(aes(size = .data[["log.P.Val"]] - 0.5, shape = .data[["significance"]]), color = "black", show.legend = TRUE) +
    geom_errorbar(aes(xmin = .data[[hr_ci_lower_prefix]], xmax = .data[[hr_ci_upper_prefix]]), size = 0.4, width = 0.2, color = "black") +
    ggtitle(title) +
    theme(plot.title = element_text(size = title_size),
      axis.line = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.y = element_text(size = axis_text_y_size),
      legend.position = legend_position,
      strip.background = element_rect(colour = "white", fill = "white", linewidth = 0.8)) +
    panel_border(colour = "black", linetype = 1, size = 0.5, remove = FALSE) +
    background_grid(major = "xy", minor = "none", size.major = 0.25) +
    scale_shape_manual(name = adjp_prefix, values = values_shape, drop = FALSE) +
    scale_fill_gradient2(name = hr_prefix, trans = "log2", breaks = scales::log_breaks(n = 7, base = 2), low = color_low, mid = color_mid, high = color_high, limits = limits, oob = scales::squish) +
    scale_radius(name = pval_prefix, range = radius_range, breaks = radius_breaks, labels = radius_labels, limits = radius_limits) + 
    scale_x_continuous(trans = "log2", breaks = scales::log_breaks(n = 7, base = 2)) +
    facet_grid(~contrast)
  
  
  ggp
  
  
  
}










