#' Cropping current row
#'
#' This function is repeadiately used by the splitting() function to make edits
#' on the current row from the given sub-data.frame.
#'
#' @param main_data This is the input data.frame read in by the HarmonizR.
#' @param batch_data This is the description data.frame read in by the
#' HarmonizR.
#' @param batch_list An overview of the batch groupings in list form.
#' @param each_vector A single vector extracted from affiliation_list. For
#' further information see the return value from the spotting_missing_values()
#' function.
#' @param current_row This is the data the function is editing.
#' @return A correctly cropped row for the splitting() function
#' @export


# ----- crop_row(main_data, batch_data, batch_list, each_vector, current_row) -
# This function gets called within the splitting() function and removes
# all columns which are not interesting for the protein currently looked
# at and returns that row
crop_row <- function(main_data,
                     batch_data,
                     batch_list,
                     each_vector,
                     current_row) {

  # This line gets a list of sample names/IDs
  batch_list_id <- c(batch_data[1]) # 1 corresponds to the column number
  temp_id_list <- list()
  my_row <- main_data[current_row, ]
  position_counter <- 1
  position_counter_2 <- 1

  # This loop determines which columns to use in the current row using
  # The given vector as a template and cropping the row accordingly
  for (element in batch_list) {
    if (element %in% each_vector == TRUE) {
      temp_id_list[[position_counter_2]] <- unlist(batch_list_id)[position_counter]
      position_counter_2 <- position_counter_2 + 1
    }
    position_counter <- position_counter + 1
  }
  # Only use the rows which appear in temp_id_list
  my_row <- subset(my_row, select = unlist(temp_id_list))
  temp_id_list <- list()

  return(my_row)
}
