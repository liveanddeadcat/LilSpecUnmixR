#' A function to generate spectral signature matrix
#'
#' @param df data.frame, the signature data frame, for plotting
#' @param plot Logical, whether a plot of spectral signatures should be generated
#' @param refs A list of reference controls & unstained control
#' @returns return a spectral signature matrix M
#' @export
M_generate <- function(df, refs, plot=TRUE){
  cvars<-gsub("\ ","",apply(expand.grid(c("UV", 'V', 'B', "YG", "R"),1:16,"-A") %>%
                              dplyr::mutate(Var1=factor(Var1,levels=unique(Var1))) %>%
                              arrange(Var1,Var2),1,function(X) paste0(X,collapse="")))

  if (plot==TRUE){

    markers <- unique(df$sample)
    # Visualise the spectral profiles
    df %>%
      dplyr::mutate(Detector=factor(`rownames(df)`,levels=cvars[cvars %in% `rownames(df)`])) %>%
      ggplot2::ggplot(aes(Detector,medians.max.medians.,group=sample,color=sample)) +
      geom_line(linewidth = 0.7) + geom_point() +
      guides(x=guide_axis(angle=45)) + theme_bw(base_size=8) +
      theme(legend.text=element_text(size=5)) +
      theme(legend.position = "right")+
      theme(legend.title = element_text(size=0.1))+
      ylab("")+
      xlab("")
    #save the plot
    ggplot2::ggsave("./plots/bead_and_cell_ref control spectral profiles.png", width = 40, height = 10, units = c("cm"))

    # construct the spectra matrix
    Mlist <- list()
    for (i in 1:length(table(df$sample))){
      Mlist[[i]] <- as.matrix(refs[[i]]$medians.max.medians.)
    }
    spectra_matrix <- do.call(cbind, Mlist)
    rm(Mlist)
    colnames(spectra_matrix) <- names(refs)

  } else {
    # construct the spectra matrix
    Mlist <- list()
    for (i in 1:length(table(df$sample))){
      Mlist[[i]] <- as.matrix(refs[[i]]$medians.max.medians.)
    }
    spectra_matrix <- do.call(cbind, Mlist)
    rm(Mlist)
    colnames(spectra_matrix) <- names(refs)
  }

  return(spectra_matrix)
}
