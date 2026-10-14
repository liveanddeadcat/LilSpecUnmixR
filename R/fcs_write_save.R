#' Write out the unmixed data and save as RDS
#'
#' Take the matrix of the unmixed data, write to FCS files
#'
#' @param unmixed_data Matrix, the unmixed full-stained data
#' @param sample_name string, name of your sample
#' @param model String, "OLS", "WLS" or something
#' @returns return the fcs files and saved RDS
#' @export
fcs_write_save <- function(unmixed_data, sample_name, model=NULL){

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
    write.FCS(ff, filename = paste0("./unmixed/", model, "_unmixed_", name, ".fcs"))
    # Save the unmixing results
    saveRDS(unmixed_data,paste0("./rds/", model, "_unmixed_",rds_name, ".rds"))
  }



