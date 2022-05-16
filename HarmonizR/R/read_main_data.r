#' Reading main data
#'
#' The read_main_data function reads in a file via its file path and converts
#' it to a for the rest of the workflow readable format.
#'
#' @param data_source Usually the path to the input data. It can also be a
#' correctly formatted data.frame.
#' @return Main Data as data.frame
#' @export


# ----- READING FILES ---------------------------------------------------------
# Contents: read_main_data()


# ----- read_main_data(data_source) -------------------------------------------
# The data is read and then every column is turned into either a numeric or NA
read_main_data <- function(data_source) {
  main_data <- read.table(data_source,
    sep = "\t", header = TRUE,
    row.names = 1, check.names = FALSE,
    comment.char = "", stringsAsFactors = FALSE
  )

  # Janitor package: remove all empty rows and columns from the data.frame
  # rows aren't really a problem, but extra empty columns screw everything up
  main_data <- janitor::remove_empty(main_data,
    which = c("rows", "cols"), quiet = TRUE
  )

  # Make numerical (I do suppress the 'introduced NAs by coercion' warnings)
  cols.num <- names(main_data)
  suppressWarnings(main_data[cols.num] <- sapply(main_data[cols.num], as.numeric))

  return(main_data)
}
