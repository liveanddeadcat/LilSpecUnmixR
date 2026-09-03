#' Write out the unmixed data and save as RDS
#'
#' Take the matrix of the unmixed data, write to FCS files
#'
#' @param unmixed_data Matrix, the unmixed full-stained data
#' @param sample_name string, name of your sample
#' @returns return the fcs files and saved RDS
#' @export
fcs_write_save <- function(unmixed_data, sample_name){

    if(!dir.exists("./rds")){
      dir.create("./rds")
    }
    if(!dir.exists("./unmixed")){
      dir.create("./unmixed")
    }

    if (!is.character(sample_name)){
      stop("Sample name must be a character.")
    }
    # write the fcs file
    ff <- flowFrame(as.matrix(unmixed_data))
    write.FCS(ff, filename = paste0("./unmixed/ols_unmixed_", sample_name, ".fcs"))
    # Save the unmixing results
    saveRDS(unmixed_data,paste0("./rds/OLS_unmixed_",sample_name, ".rds"))
  }



