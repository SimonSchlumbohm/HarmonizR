#' Rebuilding
#'
#' The rebuild function rebuilds the sub-dataframes to one big output
#' data.frame.
#'
#' @param cured_subdfs a list of data.frames, which are the result from
#' splitting().
#' @param output_file a string with the desired file name or path. By default
#' this parameter will be set equal to "cured_data".
#' @return The rebuild() function returns the adjusted data.frame and writes
#' out cured_data.tsv
#' @export


# ----- REBUILDING THE DATAFRAME ----------------------------------------------
# Contents: rebuild()


# ----- rebuild(cured_subdf, output_file) -------------------------------------
# This function takes the result of the splitting() function and binds all of
# the sub-dataframes together to one big dataframe, which is written out and
# therefore returned to the user (cured_data.tsv)
rebuild <- function(cured_subdfs, output_file) {
  prot <- unlist(lapply(cured_subdfs, row.names))
  cured_subdfs <- do.call(plyr::rbind.fill, cured_subdfs)
  row.names(cured_subdfs) <- prot
  # The cured_data.tsv file is written
  outfilename <- paste(output_file, "tsv", sep = ".")
  unlink(outfilename)
  write.table(cured_subdfs, outfilename, sep = "\t", col.names = NA)

  return(cured_subdfs)
}
