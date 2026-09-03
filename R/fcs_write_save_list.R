#' Write out to fcs for a list of unmixed data and save as RDS
#'
#' Take the list of matrices of the unmixed data, write to several FCS files, and save RDS.
#'
#' @param unmixed_data_list List of matrices, the unmixed data
#' @param rds_name String, a name to save the list of unmixed data matrices as an .rds
#' @returns return the fcs files and saved RDS
#' @export
fcs_write_save_list <- function(unmixed_data_list, rds_name=NULL){

  if(!dir.exists("./rds")){
    dir.create("./rds")
  }
  if(!dir.exists("./unmixed")){
    dir.create("./unmixed")
  }

 lapply(names(unmixed_data_list), FUN=function(name){
    ff <- flowFrame(as.matrix(unmixed_data_list[[name]]))
    # write the fcs file
    write.FCS(ff, filename = paste0("./unmixed/ols_unmixed_", name, ".fcs"))
  })

 # Save the unmixing results
 saveRDS(unmixed_data_list,paste0("./rds/OLS_unmixed_",rds_name, ".rds"))








}
