#' Visualization
#'
#' The visual functions turn their input dataframes into easily plottable
#' results.
#'
#' @param input_dataframe A data.frame object as input.
#' @param batch_list A list object giving information about which column
#' corresponds to which batch.
#' @return A data.frame object, which is ready to be plotted
#' @rdname visual
#' @export


# ----- Visualization ---------------------------------------------------------
# Contents: visual(), visual2(), visual3()


# ----- visual(input_dataframe, batch_list) -----------------------------------
# Visual calculates the mean over rows for each protein
# and then makes a boxplot-ready data.frame with all feature means
# for each batch
visual <- function(input_dataframe, batch_list) {

  # Create the plottable data.frame
  which_batch <- 1
  plotting_df <- data.frame()
  batch_list_index <- 1
  starting <- 1
  for (element in batch_list) {

    # This if() is true for the last iteration
    if (batch_list_index == length(batch_list)) {
      ending <- batch_list_index
      batch_name <- paste0("Batch_", toString(which_batch))
      plotting_df[batch_name] <- rowMeans(input_dataframe[, starting:ending], na.rm = TRUE)

      break
    }

    if (batch_list_index != 1) {
      if (element != batch_list[batch_list_index - 1]) {
        ending <- batch_list_index - 1

        # This is always the first batch
        if (starting == 1) {
          plotting_df <- data.frame(Batch_1 = rowMeans(input_dataframe[, starting:ending], na.rm = TRUE))
        } else {
          batch_name <- paste0("Batch_", toString(which_batch))
          plotting_df[batch_name] <- rowMeans(input_dataframe[, starting:ending], na.rm = TRUE)
        }

        which_batch <- which_batch + 1
        starting <- batch_list_index
      }
    }

    batch_list_index <- batch_list_index + 1
  }

  return(plotting_df)
}


#' Visualization
#'
#' @rdname visual
#' @export


# ----- visual2(input_dataframe, batch_list) ----------------------------------
# Visual 2 calculates the means of each sample
# and then makes a boxplot-ready data.frame of all means for each batch
visual2 <- function(input_dataframe, batch_list) {
  showing_df <- data.frame(colMeans(input_dataframe, na.rm = TRUE))

  list_of_dfs <- list()

  which_batch <- 1
  batch_list_index <- 1
  starting <- 1

  for (element in batch_list) {

    # This if() is true for the last iteration
    if (batch_list_index == length(batch_list)) {
      ending <- batch_list_index

      # This is the data.frame we want to add to the list_of_dfs
      to_add <- data.frame(showing_df[starting:ending, ])
      colnames(to_add) <- paste0("Batch_", toString(which_batch))

      list_of_dfs <- c(list_of_dfs, to_add)

      break
    }

    if (batch_list_index != 1) {
      if (element != batch_list[batch_list_index - 1]) {
        ending <- batch_list_index - 1

        # This is the data.frame we want to add to the list_of_dfs
        to_add <- data.frame(showing_df[starting:ending, ])
        colnames(to_add) <- paste0("Batch_", toString(which_batch))

        list_of_dfs <- c(list_of_dfs, to_add)

        which_batch <- which_batch + 1
        starting <- batch_list_index
      }
    }

    batch_list_index <- batch_list_index + 1
  }

  return(list_of_dfs)
}


#' Visualization
#'
#' @rdname visual
#' @export


# ----- visual3(input_dataframe, batch_list) ----------------------------------
# Visual 3 calculates the CV of each sample
# and then makes a boxplot-ready data.frame of all CVs for each batch
visual3 <- function(input_dataframe, batch_list) {
  showing_df <- data.frame(sapply(input_dataframe, function(x) sd(x, na.rm = T) / mean(x, na.rm = T)))

  list_of_dfs <- list()

  which_batch <- 1
  batch_list_index <- 1
  starting <- 1

  for (element in batch_list) {

    # This if() is true for the last iteration
    if (batch_list_index == length(batch_list)) {
      ending <- batch_list_index

      # This is the data.frame we want to add to the list_of_dfs
      to_add <- data.frame(showing_df[starting:ending, ])
      colnames(to_add) <- paste0("Batch_", toString(which_batch))

      list_of_dfs <- c(list_of_dfs, to_add)

      break
    }

    if (batch_list_index != 1) {
      if (element != batch_list[batch_list_index - 1]) {
        ending <- batch_list_index - 1

        # This is the data.frame we want to add to the list_of_dfs
        to_add <- data.frame(showing_df[starting:ending, ])
        colnames(to_add) <- paste0("Batch_", toString(which_batch))

        list_of_dfs <- c(list_of_dfs, to_add)

        which_batch <- which_batch + 1
        starting <- batch_list_index
      }
    }

    batch_list_index <- batch_list_index + 1
  }

  return(list_of_dfs)
}
