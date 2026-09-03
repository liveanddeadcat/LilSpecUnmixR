#' Build and visualise a few pre-unmixing diagnostics
#'
#' Take the spectral signature matrix, preview covariance matrix, cosine similarity matrix, and the condition number.
#'
#' @param M Matrix, the spectral signature matrix
#' @param plot Logical. Default TRUE. Whether plots should be generated.
#' @returns return the covariance matrix, cosine similarity matrix, and condition number
#' @export
panel_diagnose <- function(M, plot=TRUE){
  # Calculate M^T
  M <- as.data.table(M)
  M_t <- transpose(M)
  M <- as.matrix(M)
  M_t <- as.matrix(M_t)
  rownames(M_t)<- colnames(M)

  # (M^T M)^(-1) (Kmet, 2025, the endmember-dependent component of the spectra)
  mtm <- solve(M_t %*% M)
  rownames(mtm) <- colnames(mtm)

  # Pearson's cor
  cor_M <- cor(M)

  # cosine sim
  cos_M <- proxy::simil(M_t, method = "cosine")
  cos_M <- as.matrix(cos_M)
  colnames(cos_M)<- colnames(M)
  rownames(cos_M) <- colnames(M)
  diag(cos_M)<- 1

  kappa <- kappa(M)

  if (plot==FALSE){
   stop("The matrices and the condition number have been generated. No plotting.")
  }else{
    # plotting
    mtm[upper.tri(mtm)] <- NA
    mtm %>%
      reshape2::melt() %>%
      ggplot(aes(x = Var2, y = Var1, fill = value)) +
      geom_tile(color = "white", linewidth = 0.5) + # white grid lines
      scale_fill_gradient2(
        low = "white",
        mid = "grey",
        high = "red",
        midpoint = 0,
        limits = c(0, max(mtm)),
        na.value = "transparent",
        oob = scales::squish
      ) +
      coord_fixed() + # square tiles
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank()
      ) +
      labs(fill = "Value", x = NULL, y = NULL)

    ggsave(filename="./plots/covariance.png")

    png(filename = "./plots/correlation_plot.png", width = 2000, height = 2000, res = 300)
    p_cor <- corrplot(cor_M, method = "shade",tl.col = "black", type = c("low"), addCoef.col = "black")
    dev.off()

    png(filename = "./plots/cosine_simil_plot.png", width = 2000, height = 2000, res = 300)
    p_cos <- corrplot(cos_M, method = "shade",tl.col = "black", type = c("low"),  na.label = " ",cl.cex = 1,
                      na.label.col = "white", addCoef.col = "black")
    dev.off()
  }

  return(list(Covariance=mtm,PearsonsCor=cor_M, CosineSimil=cos_M, ConditionNumber=kappa))
}
