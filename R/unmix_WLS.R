#' Unmixing Weighted Least-Squares, standard solution
#'
#' Take the spectral signature M, calculate the weights by 1/detector MFI, and unmix the full-stained sample via closed-form solution.
#'
#' @param M Matrix, the signature matrix
#' @param fs_raw string, the name of the raw full-stained .fcs file
#' @param aurora_channel_change Logical, Whether you would like to change Aurora configuration. Default to WIMR, Australia.
#' @param channels If aurora_channel_change is TRUE, enter a combined string of names of channels
#' @returns return a matrix of unmixed data
#' @export
unmix_WLS <- function(M,fs_raw, aurora_channel_change = FALSE, channels = NULL){

  # weights - 1/mean
  Weights <- 1/rowMeans(M)
  W <- diag(Weights)

  # Calculate M^T
  M <- as.data.table(M)
  M_t <- transpose(M)
  M <- as.matrix(M)
  M_t <- as.matrix(M_t)
  rownames(M_t)<- colnames(M)

  # the weighted M inverse matrix
  wmtm_1 <- solve(M_t %*% W %*% M)

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

  sample <- read.FCS(fs_raw, transformation = "linearize-with-PnG-scaling")

  sample_fluo<-exprs(sample[,sample@parameters@data$name[sample@parameters@data$name %in% aurora_channels]])

  colnames(sample_fluo) <- NULL

  # transpose it into the cells as columns vs channels as rows format
  sample_fluo<- as.data.table(sample_fluo)
  sample_fluo<- transpose(sample_fluo)
  sample_fluo <- as.matrix(sample_fluo)

  # Solution
  unmixed_matrix <- wmtm_1 %*% M_t %*% W %*% sample_fluo

  # Transpose to make the cells on rows again to match the original data
  t_unmixed_matrix <- t(unmixed_matrix)
  rm(unmixed_matrix, sample_fluo)

  colnames(t_unmixed_matrix) <- gsub("_","_unmixed_", colnames(t_unmixed_matrix))

  #Add these columns to the Time and FSC, SSC channels
  unmixed_data <- cbind(exprs(sample[,sample@parameters@data$name[!sample@parameters@data$name %in% aurora_channels]]),
                        t_unmixed_matrix)
  return(unmixed_data)

}
