




#' Forest plot
#' 
#' 
#' @param data Data frame.
wrapper_core_forest_plot <- function(out, title = NULL, hrzl_lines_var = "Biomarker", clip = c(0, 4), lineheight = "auto"){
  
  # lineheight = unit(1, "cm")
  
  
  out$cox_out[, "HR"] <- out$cox_out$`Hazard ratio`
  out$cox_out[, "95% CI"] <- paste0("(", out$cox_out$`95% CI lower`, "-", out$cox_out$`95% CI upper`, ")")
  
  labeltext <- out$cox_out
  
  ### Place P-value at the end 
  if("P-value" %in% colnames(labeltext)){
    labeltext <- labeltext %>% select(-`P-value`, everything())
  }
  
  ### Place Adj. P-value at the end
  if("Adj. P-value" %in% colnames(labeltext)){
    labeltext <- labeltext %>% select(-`Adj. P-value`, everything())
  }
  
  
  labeltext[, c("Hazard ratio", "95% CI lower", "95% CI upper")] <- NULL
  labeltext <- rbind(colnames(labeltext), labeltext)
  
  
  out$cox_out$`95% CI upper`[out$cox_out$`95% CI upper` == Inf] <- NA
  
  
  ### To separate Biomarkers with a horizontal line
  
  ## Find numbers of rows where the biomarker blocks end. This value has to be shifted by 2
  line_row <- which(!duplicated(out$cox_out[, hrzl_lines_var], fromLast = TRUE)) + 2
  line_row <- line_row[-length(line_row)]
  
  hrzl_lines <- list()
  
  if(length(line_row) >= 1){
    hrzl_lines <- lapply(line_row, function(x){
      gpar(col = "#b4b4b4", lwd = 0.5)
    })
    names(hrzl_lines) <- line_row
  }
  
  hrzl_lines[["2"]] <- gpar(col = "#444444", lwd = 1)
  
  
  forestplot(labeltext, mean = c(NA, out$cox_out$`Hazard ratio`), lower = c(NA, out$cox_out$`95% CI lower`), upper = c(NA, out$cox_out$`95% CI upper`), 
    is.summary = c(TRUE, rep(FALSE, nrow(out$cox_out))), xlog = FALSE, xlab = "Hazard Rate", zero = 1,
    title = title,
    col = fpColors(box = "darkblue", line = "darkblue"), 
    boxsize = 0.3,
    hrzl_lines = hrzl_lines, 
    graphwidth = unit(10, "cm"), colgap = unit(6, "mm"),
    lineheight = lineheight,
    lwd.ci = 2, lwd.xaxis = 2, lwd.zero = 2, 
    txt_gp = fpTxtGp(xlab = gpar(fontsize = 18), ticks = gpar(fontsize = 16)), 
    mar = unit(c(20, rep(5, times = 3)), "mm"), 
    clip = clip, 
    ci.vertices = TRUE)
  
  footnote <- stringr::str_wrap(out$cox_footnote, width = 120)
  
  x <- unit(30, 'mm')
  y <- unit(5, 'mm')
  
  grid.text(footnote, x, y, gp = gpar(fontsize = 10, font = 3), just = c("left", "bottom"))
  
  
}
















