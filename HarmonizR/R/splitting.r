#' Splitting
#'
#' This function splits the data.frame. The data is very sensitive to its
#' specific input. Only to be called via harmonizR()
#'
#' @param affiliation_list An overview of which protein has which missing value
#' distribution.
#' @param main_data This is the input data.frame read in by the HarmonizR.
#' @param batch_data This is the description data.frame read in by the
#' HarmonizR.
#' @param batch_list An overview of the batch groupings in list form.
#' @param algorithm Either "ComBat" or "limma". Based on the selected algorithm
#' for the harmonizR() function.
#' @param ComBat_mode The chosen ComBat mode influences the parameters the
#' ComBat algorithm is using. Based on the ComBat_mode parameter given to the
#' harmonizR() function. Not active during limma execution.
#' @return Returns a list of 'chopped up' data.frames
#' @export


# ----- Splitting the dataframe and performing ComBat -------------------------
# Contents: splitting()


# -- splitting(affiliation_list, main_data, batch_data, batch_list, ComBat_mode)
# The first part of the function function splits the dataframe (main_data)
# based on 'affiliation_list' and stores them in corresponding variables
# in the second part of the function, the ComBat algorithm is called for all
# suitable sub-dataframes, which are afterwards combined and written into a
# .csv file, yielding the result and terminating the HarmonizR algorithm
splitting <- function(affiliation_list,
                      main_data,
                      batch_data,
                      batch_list,
                      algorithm,
                      ComBat_mode) {

  # ----- THE FIRST PART OF THE splitting() FUNCTION STARTS -------------------
  storage_list <- list()
  zero_names <- list()
  one_names <- list()
  two_or_more_names <- list()
  current_row <- 1
  zero_counter <- 1
  one_counter <- 1
  two_or_more_counter <- 1

  # In this loop, the entire input dataframe (main_data) gets split into
  # sub-dataframes based on the template provided by 'affiliation list'
  for (each_vector in affiliation_list) {
    # The function crop_row() prepares the row
    some_row <- crop_row(
      main_data,
      batch_data,
      batch_list,
      each_vector,
      current_row
    )

    # This if-statement is used only the first time a specific vector was used
    if ((toString(each_vector) %in% storage_list) == FALSE) {
      # The variable gets declared and it's value assigned (the df row)
      assign(paste0("batches", toString(each_vector)), some_row)
    }
    # Every time a familiar vector is found, this else-statement is executed
    else {
      # This is how to call the variable with only the current vector as
      # information
      current_name <- eval(parse(text = paste0("`batches", toString(each_vector), "`")))
      current_name <- rbind(current_name, some_row)
      assign(paste0("batches", toString(each_vector)), current_name)
    }

    # Here the sub-dataframe's names are appended to the correct list of names,
    # effectively storing which sub-dataframe belongs to which category, which
    # are defined by data appearing in either 0, 1 or 2+ batches
    if (length(each_vector) >= 2) {
      two_or_more_names[[two_or_more_counter]] <- paste0("`batches", toString(each_vector), "`")
      two_or_more_counter <- two_or_more_counter + 1
    } else if (length(each_vector) == 1) {
      one_names[[one_counter]] <- paste0("`batches", toString(each_vector), "`")
      one_counter <- one_counter + 1
    } else {
      zero_names[[zero_counter]] <- paste0("`batches", toString(each_vector), "`")
      zero_counter <- zero_counter + 1
    }

    # The current vector derived from 'affiliation_list' is remembered in
    # 'storage_list', so no existing sub-dataframe gets overwritten by accident
    storage_list[[current_row]] <- toString(each_vector)
    current_row <- current_row + 1
  }

  two_or_more_names <- unique(two_or_more_names)
  one_names <- unique(one_names)
  zero_names <- unique(zero_names)

  # Here the sub-dataframes are appended to the list they belong to, based
  # on whether the protein appears in 0, 1 or 2+ batches
  two_or_more_subdf <- list()
  one_subdf <- list()
  zero_subdf <- list()

  # The lists of names filled up in the previous loop are used to append the
  # actual sub-dataframe's contents to the correct list below
  j <- 1
  for (each_name in two_or_more_names) {
    two_or_more_subdf[[j]] <- eval(parse(text = each_name))
    j <- j + 1
  }
  j <- 1
  for (each_name in one_names) {
    one_subdf[[j]] <- eval(parse(text = each_name))
    j <- j + 1
  }
  j <- 1
  for (each_name in zero_names) {
    zero_subdf[[j]] <- eval(parse(text = each_name))
    j <- j + 1
  }

  # ----- THE SECOND PART OF THE splitting() FUNCTION STARTS ------------------

  cured_subdfs <- list()

  # This line checks how many cores are available and saves them in 'numCores'
  numCores <- parallel::detectCores()
  # Here, parallelization is activated and set to 'numCores'
  doParallel::registerDoParallel(cores = numCores)

  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`

  # In this parallelized loop the ComBat algorithm is used for batch adjustment
  cured_subdfs <- foreach::foreach(current_subdf = two_or_more_subdf, .combine = c) %dopar% {
    description <- data.frame()
    # The description for each sub-dataframe is fetched in this loop
    for (row in 1:nrow(batch_data)) {
      # The respective description file is build with the current
      # sub-dataframe (current_subdf) as a template
      if (batch_data[row, 1] %in% names(current_subdf) == TRUE) {
        description <- rbind(description, batch_data[row, ])
      }
    }
    # The sub-dataframe and the description are changed into type 'double'
    # in order to feed the ComBat algorithm properly
    # it is now type 'double'
    ComBat_input_dat <- as.matrix(current_subdf)
    # it is now type 'double'
    ComBat_input_batch <- c(description[3])
    ComBat_input_batch <- as.numeric(ComBat_input_batch[[1]])

    # All 'one protein' files are skipped due to a ComBat limitation
    if (nrow(ComBat_input_dat) > 1) {

      # limma part
      if (algorithm == "limma") {
        ComBat_result <- as.data.frame(limma::removeBatchEffect(
          ComBat_input_dat,
          ComBat_input_batch
        ))
      } else {
        # Here the actual ComBat call takes place
        # an if-else-statement is used to execute the chosen ComBat_mode
        if (ComBat_mode == 1) {
          # This is the default ComBat mode
          ComBat_result <- as.data.frame(sva::ComBat(
            dat = ComBat_input_dat,
            batch = ComBat_input_batch
          ))
        } else if (ComBat_mode == 2) {
          ComBat_result <- as.data.frame(sva::ComBat(
            dat = ComBat_input_dat,
            batch = ComBat_input_batch,
            par.prior = TRUE,
            mean.only = TRUE
          ))
        } else if (ComBat_mode == 3) {
          ComBat_result <- as.data.frame(sva::ComBat(
            dat = ComBat_input_dat,
            batch = ComBat_input_batch,
            par.prior = FALSE,
            mean.only = FALSE
          ))
        } else if (ComBat_mode == 4) {
          ComBat_result <- as.data.frame(sva::ComBat(
            dat = ComBat_input_dat,
            batch = ComBat_input_batch,
            par.prior = FALSE,
            mean.only = TRUE
          ))
        }
      }

      # In the list 'cured_subdfs' the ComBat-adjusted dataframes are stored
      # this return defines what will be gotten back from the foreach-loop
      return(list(ComBat_result))
    } else {
      return(list())
    }
  }

  # This line combines all the ComBat-adjusted dataframes with all
  # sub-dataframes containing proteins holding data in exactly 1 batch
  cured_subdfs <- c(cured_subdfs, one_subdf)

  return(cured_subdfs)
}
