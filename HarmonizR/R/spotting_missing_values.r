#' Spotting
#'
#' This function spots missing values within the given data.frame.
#'
#' @param main_data This is the input data.frame read in by the HarmonizR.
#' @param batch_list An overview of the batch groupings in list form.
#' @param needed_values The number of values needed to be present in a batch
#' in order to be valid.
#' @return A list of vectors to pass to the upcoming splitting() function.
#' @export


# ----- CHECK FOR MISSING VALUES (A LIST OF VECTORS) --------------------------
# Contents: spotting_missing_values()


# ----- spotting_missing_values(main_data, batch_list, needed_values) ---------
# This function builds the "affiliation_list".
# "affiliation list" is a list of vectors and is the return of this function.
spotting_missing_values <- function(main_data, batch_list, needed_values) {

  # The initialization
  missing_value_found <- FALSE
  value_existence_counter <- 0
  temp_vec <- vector()
  affiliation_list <- list()

  # Iterate through the entire dataset "main_data"
  for (i in 1:nrow(main_data)) {
    for (j in 1:ncol(main_data)) {
      # If a value is found, increase "value_existence_counter"
      if (is.na(main_data[i, j]) == FALSE) {
        value_existence_counter <- value_existence_counter + 1
      }
      # The following if-statement means:
      # The next iteration will look into a new batch
      if (is.na(batch_list[j + 1]) || batch_list[j + 1] != batch_list[j]) {
        # Needed_values is either 2 or 1 based on the chosen ComBat_mode
        if (value_existence_counter < needed_values) {
          missing_value_found <- TRUE
        }
        # We found a batch with enough valid values
        if (missing_value_found != TRUE) {
          # Since we can use the protein data for the batch, we add it to the
          # corresponding vector, which later ends up in the "affiliation_list"
          temp_vec[batch_list[j]] <- batch_list[j]
        }
        # Reset "missing_value_found" to FALSE
        # and "value_existence_counter" to 0
        missing_value_found <- FALSE
        value_existence_counter <- 0
      }
    }
    # Append the vector to "affiliation_list"
    temp_vec <- temp_vec[!is.na(temp_vec)]
    affiliation_list[[i]] <- temp_vec
    # Reset "temp_vec"
    temp_vec <- vector()
  }

  # Returns a list of vectors ("affiliation_list").
  # Each element/vector in the list is correspondent to a row in "main_data"
  return(affiliation_list)
}
