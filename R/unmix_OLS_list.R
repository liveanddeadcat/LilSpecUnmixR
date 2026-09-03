#' Unmixing OLS, standard solution for a list of samples
#'
#' Take the spectral signature M and unmix a list of full-stained samples via closed-form solution of OLS.
#'
#' @param M Matrix, the signature matrix
#' @param fs_raw_dir a string of the path that stores a list of fcs files to unmix.
#' @param aurora_channel_change Logical, Whether you would like to change Aurora configuration. Default to WIMR, Australia.
#' @param channels If aurora_channel_change is TRUE, enter a combined string of names of channels
#' @param plot Plotting the unmixing matrix. Default TRUE.
#' @returns return a list of matrices of unmixed data
#' @export
unmix_OLS_list <- function(M,fs_raw_dir, aurora_channel_change = FALSE, channels = NULL, plot = TRUE){

  # Calculate M^T
  M <- as.data.table(M)
  M_t <- transpose(M)
  M <- as.matrix(M)
  M_t <- as.matrix(M_t)
  rownames(M_t)<- colnames(M)

  # (M^T M)^(-1) (Kmet, 2025, the endmember-dependent component of the spectra)
  mtm_1 <- solve(M_t %*% M)

  # unmixing matrix
  pseudo_inver_M <- mtm_1 %*% M_t

  if (aurora_channel_change==TRUE){
    if (is.null(channels)){
      stop("If aurora_channel_change = TRUE, supply strings of your Aurora configuration.")
    }
    rownames(pseudo_inver_M) <- colnames(M)
    colnames(pseudo_inver_M) <- aurora_channels

  }else{
    aurora_channels <- c('UV1-A',  'UV2-A',  'UV3-A',  'UV4-A',  'UV5-A',  'UV6-A',  'UV7-A',  'UV8-A',  'UV9-A',  'UV10-A','UV11-A', 'UV12-A', 'UV13-A', 'UV14-A', 'UV15-A', 'UV16-A', 'V1-A',   'V2-A',
                         'V3-A',   'V4-A',   'V5-A',   'V6-A',   'V7-A',   'V8-A',   'V9-A',   'V10-A',  'V11-A',  'V12-A',  'V13-A',  'V14-A',  'V15-A',  'V16-A',  'B1-A',   'B2-A',   'B3-A',   'B4-A',  'B5-A',   'B6-A',   'B7-A',   'B8-A',   'B9-A',   'B10-A',  'B11-A',  'B12-A',  'B13-A',  'B14-A',  'YG1-A',  'YG2-A',  'YG3-A',  'YG4-A',  'YG5-A',  'YG6-A',  'YG7-A',  'YG8-A', 'YG9-A',  'YG10-A', 'R1-A',   'R2-A',   'R3-A',   'R4-A',   'R5-A',   'R6-A',   'R7-A',   'R8-A')
    rownames(pseudo_inver_M) <- colnames(M)
    colnames(pseudo_inver_M) <- aurora_channels
  }

  if (plot==FALSE){
    stop("No plotting. The pseudoinverse has been generated.")

  } else {
    pheatmap(pseudo_inver_M,
             main = "Unmixing matrix by OLS",
             color = colorRampPalette(c("blue", "white", "red"))(100),
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             angle_col = c("45"),
             filename = paste0("./plots/OLS_unmixing_matrix.png"),
             #  cellwidth = 5.5,
             width = 12,
             height = 4)

  }

  fs_raw <- list.files(fs_raw_dir, pattern = "\\.fcs$", full.names=T)

# Unmixing
  unmixed_sample_list <- lapply(fs_raw, FUN = function(fs){

    fcs <- read.FCS(fs, transformation = "linearize-with-PnG-scaling")
    sample_fluo<- exprs(fcs[,fcs@parameters@data$name[fcs@parameters@data$name %in% aurora_channels]])
    colnames(sample_fluo) <- NULL

    # transpose it into the cells as columns vs channels as rows format
    sample_fluo<- as.data.table(sample_fluo)
    sample_fluo<- transpose(sample_fluo)
    sample_fluo <- as.matrix(sample_fluo)

    # Applying unmixing matrix
    unmixed_fluo <- t(pseudo_inver_M %*% sample_fluo)
    colnames(unmixed_fluo) <- gsub("_","_unmixed_", colnames(unmixed_fluo))

    # Attach back to the time and forward-side scatters
   cbind(exprs(fcs[,fcs@parameters@data$name[!fcs@parameters@data$name %in% aurora_channels]]),
          unmixed_fluo)

  })
  names(unmixed_sample_list)<- list.files(fs_raw_dir, pattern = "\\.fcs$", full.names=F)
  return(unmixed_sample_list)

}
