#' Generate spectral signatures
#'
#' Load reference controls and unstained controls and generate signatures
#'
#' @param ref_dir Directory containing reference control .fcs, string
#' @param af Whether to do autofluorescence extraction. Default = TRUE, recommended.
#' @param af_dir Directory containing AF .fcs, string
#' @param change_order Logical, only use when exporting gated .fcs files from flowjo. Default order is as in Aurora in WIMR, Australia.
#' @param custom_order Change order to another Aurora's channels.
#' @param channel_order When custom_order = TRUE, enter a string vector of the Aurora channels.
#' @returns return a data.frame of channels, normalised signatures, marker names
#' @export
load_reference <- function(ref_dir, af=TRUE, af_dir, change_order=FALSE, custom_order=FALSE, channel_order=NULL){
  files <- dir(ref_dir, pattern = "fcs$", full.names=T)
  as.matrix(files)
  # read in fcs files
  dat <- lapply(files, flowCore::read.FCS)
  names(dat) <-gsub("./ref controls/","",files)
  as.matrix(names(dat))
  p1 <- lapply(seq(2,length(dat), by = 2), function(i){flowSpectrum::spectralplot(dat[[i]], normalize = TRUE, guessPop = TRUE, unstained = dat[[i-1]])$data})
  #Name the list of new dfs
  names(p1) <- gsub("./ref control/", "",names(dat[seq(2,length(dat),2)]))
  names(p1) <- gsub("_pos.fcs", "",names(dat[seq(2,length(dat),2)]))
  print(as.matrix(names(p1)))

  if (change_order==TRUE){
    if (custom_order==FALSE){
      aurora_channels <- c('UV1-A',  'UV2-A',  'UV3-A',  'UV4-A',  'UV5-A',  'UV6-A',  'UV7-A',  'UV8-A',  'UV9-A',  'UV10-A','UV11-A', 'UV12-A', 'UV13-A', 'UV14-A', 'UV15-A', 'UV16-A', 'V1-A',   'V2-A',
                           'V3-A',   'V4-A',   'V5-A',   'V6-A',   'V7-A',   'V8-A',   'V9-A',   'V10-A',  'V11-A',  'V12-A',  'V13-A',  'V14-A',  'V15-A',  'V16-A',  'B1-A',   'B2-A',   'B3-A',   'B4-A',  'B5-A',   'B6-A',   'B7-A',   'B8-A',   'B9-A',   'B10-A',  'B11-A',  'B12-A',  'B13-A',  'B14-A',  'YG1-A',  'YG2-A',  'YG3-A',  'YG4-A',  'YG5-A',  'YG6-A',  'YG7-A',  'YG8-A', 'YG9-A',  'YG10-A', 'R1-A',   'R2-A',   'R3-A',   'R4-A',   'R5-A',   'R6-A',   'R7-A',   'R8-A' )
      # Reorder the channels
      for (i in 1:length(p1)){
        p1[[i]] <- p1[[i]][match(aurora_channels, p1[[i]]$`rownames(df)`),]
      }
    }else{
      aurora_channels <- channel_order
      for (i in 1:length(p1)){
        p1[[i]] <- p1[[i]][match(aurora_channels, p1[[i]]$`rownames(df)`),]
      }
        }

# Add a sample name to the dfs
    for(i in 1:length(p1)){p1[[i]]$sample = names(p1)[i]}
    cvars<-gsub("\ ","",apply(expand.grid(c("UV", 'V', 'B', "YG", "R"),1:16,"-A") %>%
                                dplyr::mutate(Var1=factor(Var1,levels=unique(Var1))) %>%
                                arrange(Var1,Var2),1,function(X) paste0(X,collapse="")))

    # combine the df into a big df
    df <- df <- dplyr::bind_rows(c(p1))
  } else {
    # add a sample column for plotting
    for(i in 1:length(p1)){p1[[i]]$sample = names(p1)[i]}

    cvars<-gsub("\ ","",apply(expand.grid(c("UV", 'V', 'B', "YG", "R"),1:16,"-A") %>%
                                dplyr::mutate(Var1=factor(Var1,levels=unique(Var1))) %>%
                                arrange(Var1,Var2),1,function(X) paste0(X,collapse="")))
    # combine the df into a big df
    df <- df <- dplyr::bind_rows(c(p1))
  }

if (af==FALSE){
  stop("No AF extraction.")
}else{
  # Import AF profile
  files <- dir(af_dir, pattern = "fcs$", full.names=T)
  as.matrix(files)
  # read in fcs files
  dat <- lapply(files, read.FCS)
  # tidy up the file names as you wish, if you need to
  names(dat) <-gsub("./af/","",files)
  names(dat) <-gsub(".fcs.*","",names(dat))
  # Extract normalised spectra
  af <- lapply(1:length(dat), function(i){flowSpectrum::spectralplot(dat[[i]], normalize = TRUE)$data})
  names(af) <- names(dat)
  af[[1]]$sample <- "AF"
  af <- af[[1]]
  # attaching AF profile to the spectral profile data frame and list
  df <- rbind(df,af)
  p1$AF <- af
}

  return(list(df=df, refs=p1))
}
