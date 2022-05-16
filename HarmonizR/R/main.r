#' Main function
#'
#' This function executes the entire harmonizR program and executes all other
#' functions found in this package. Therefore, this is the only function in
#' need of calling.
#'
#' @param data_as_input Path to input data. Please read the SOP for the correct
#' input format found on https://github.com/SimonSchlumbohm/HarmonizR.
#' Additionally, the input can be a data.frame with proper row- and column
#' names.
#' @param description_as_input Path to input description. Please read the SOP
#' for the correct input format found on
#' https://github.com/SimonSchlumbohm/HarmonizR. Additionally, the input can be
#' a data.frame with three columns total.
#' @param algorithm Optional. Pass either "ComBat" or "limma" to select the
#' preferred adjustment method. Defaults to ComBat.
#' @param ComBat_mode Optional. Pass a number between 1 and 4 to select the
#' desired ComBat parameters. Can only be set when ComBat is used. For
#' information on the meaning of the numbers, please view the SOP found on
#' https://github.com/SimonSchlumbohm/HarmonizR. Defaults to 1.
#' @param plot Optional. Takes either "samplemeans" for sample specific means,
#' "featuremeans" for feature specific means or "CV" for the coeffivient of
#' variation as input and creates before/after plots for the given data.
#' Defaults to FALSE -> Turned off.
#' @param output_file Optional. Takes a string as input for the .tsv file name.
#' This can also be a path. Defaults to "cured_data", hence yielding a
#' "cured_data.tsv" file.
#' @return A .tsv file by default called cured_data.tsv will be written out as
#' a result.
#' As a return, the harmonizR function will yield the batch effect adjusted
#' data.frame.
#' @export


# ----- MAIN (FUNCTION CALLS) -------------------------------------------------
# Contents: harmonizR()


# ----- harmonizR(
#         data_as_input,
#         description_as_input,
#         algorithm,
#         ComBat_mode,
#         plot,
#         output_file
#       )
# This function will be called by the User and will execute the entire program
harmonizR <- function(data_as_input = NULL,
                      description_as_input = NULL,
                      ...,
                      algorithm = "ComBat",
                      ComBat_mode = 1,
                      plot = FALSE,
                      output_file = "cured_data") {

  # Check whether data and description are present
  if (is.null(data_as_input) && is.null(description_as_input)) {
    stop("No parameters given. Usage: harmonizR(\"path/to/data\", \"path/to/description\")")
  } else if (is.null(description_as_input)) {
    stop("Not enough parameters. Usage: main(\"path/to/data\", \"path/to/description\")")
  }

  # Check algorithm input
  if (algorithm != "ComBat" && algorithm != "limma") {
    print("Please set the algorithm parameter to either ComBat or limma to choose the prefered adjustment. Parameter is now being set to default (ComBat).")
    algorithm <- "ComBat"
  }

  # Check ComBat_mode input
  if (ComBat_mode < 1 || ComBat_mode > 4) {
    print("Invalid ComBat Mode chosen. Select 1, 2, 3 or 4. Parameter is now being set to default (1).")
    ComBat_mode <- 1
  }

  # Check plot input
  if (plot != FALSE &&
    plot != "samplemeans" &&
    plot != "featuremeans" &&
    plot != "CV") {
    print("Please set the plot parameter to either samplemeans, featuremeans, CV or FALSE.")
    plot <- FALSE
  }

  # Check output_file input
  if (is.character(output_file) != TRUE) {
    print("Please only pass a string via the output_file parameter.")
    output_file <- "cured_data"
  }

  # This line potentially checks the OS
  # print(paste0("Current operating system: ", toupper(.Platform$OS.type)))


  # ----- READING SECTION -----------------------------------------------------
  print("Reading the files...")

  # Logic for the Perseus plugin
  if (is.character(data_as_input)) {
    # Read in the data
    main_data <- read_main_data(data_as_input)
  } else {
    # Read in the data
    main_data <- data_as_input
  }

  # Logic for the Perseus plugin
  if (is.character(description_as_input)) {
    # Read in the batch-descriptions
    batch_data <- read_description(description_as_input)
  } else {
    # Read in the batch-descriptions
    batch_data <- description_as_input
  }

  # Information about which row of batch_data is in which batch
  batch_list <- fetch_batch_overview(batch_data)

  # Removing duplicates
  main_data <- unique(main_data)

  # ----- SPOTTING SECTION ----------------------------------------------------
  print("Preparing...")

  # Create the "affiliation_list", the list of vectors
  if (ComBat_mode == 1 ||
    ComBat_mode == 3 ||
    algorithm == "limma") {
    affiliation_list <- spotting_missing_values(main_data, batch_list, 2)
  }
  if (ComBat_mode == 2 || ComBat_mode == 4) {
    affiliation_list <- spotting_missing_values(main_data, batch_list, 1)
  }


  # ----- SPLITTING SECTION ---------------------------------------------------
  print(paste("Splitting the data using", algorithm, "adjustment..."))

  # Split up the "main_data" dataframe into the sub-dataframes
  cured_subdfs <- splitting(
    affiliation_list,
    main_data,
    batch_data,
    batch_list,
    algorithm,
    ComBat_mode
  )


  # ----- REBUILD SECTION -----------------------------------------------------
  print("Rebuilding...")

  # Build the result file
  cured <- rebuild(cured_subdfs, output_file)
  # Save a version of the logged cured for return
  original_cured <- cured


  # ----- VISUALIZATION SECTION -----------------------------------------------

  main_data <- 2^main_data
  cured <- 2^cured

  if (plot == "featuremeans") {
    print("Visualizing feature means...")

    original <- visual(main_data, batch_list)
    corrected <- visual(cured, batch_list)

    par(mfrow = c(1, 2))
    boxplot(original, main = "Original", las = 2)
    boxplot(corrected, main = "Corrected", las = 2)
  } else if (plot == "samplemeans") {
    print("Visualizing sample means...")

    original <- visual2(main_data, batch_list)
    corrected <- visual2(cured, batch_list)

    lmts <- range(original, corrected)

    par(mfrow = c(1, 2))
    boxplot(original, main = "Original", las = 2, ylim = lmts)
    boxplot(corrected, main = "Corrected", las = 2, ylim = lmts)
  } else if (plot == "CV") {
    print("Visualizing CV...")

    original <- visual3(main_data, batch_list)
    corrected <- visual3(cured, batch_list)

    lmts <- range(original, corrected)

    par(mfrow = c(1, 2))
    boxplot(original, main = "Original", las = 2, ylim = lmts)
    boxplot(corrected, main = "Corrected", las = 2, ylim = lmts)
  }


  # ----- END -----------------------------------------------------------------
  print("Termination.")

  return(original_cured)
}
