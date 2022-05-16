# All plugins can be split into 3 parts
# 1. Reading the command line arguments provided by Perseus and parsing the data.
# 2. Perform the desired functionality on the data.
# 3. Write the results to the expected locations in the Perseus formats.


# 1. Parse command line arguments passed in from Perseus with 
# optional, additional arguments
argv = commandArgs(trailingOnly=TRUE)
library("argparser")
p <- arg_parser(description = "Head processing")
p <- add_argument(p, 'input', help="path of the input file")
p <- add_argument(p, 'output', help="path of the output file")
p <- add_argument(p, '--algorithm', type="character", default="ComBat", help="chosen algorithm")
p <- add_argument(p, '--ComBat_mode', type="numeric", default=1, help="chosen ComBat_mode")
argp <- parse_args(p, argv)

inFile <- argp$input
outFile <- argp$output
algorithm_chosen <- argp$algorithm
ComBat_mode_chosen <- argp$ComBat_mode

# 1.2Use PerseusR to read and write the data in Perseus text format.
library(PerseusR)
mdata <- read.perseus(inFile)
# The mdata object can be easily deconstructed into a number of different
# data frames. Check reference manual or help() for full list.

# Build input_data START
mainMatrix <- main(mdata)

# get correct rownames for input matrix
rowstuff <- annotCols(mdata)
colnames(rowstuff) <- "theprotIDs"
rownames(mainMatrix) <- rowstuff$theprotIDs
# DEBUGGING
#unlink("whats.txt")
#write.table(mainMatrix, "whats.txt", sep="\t", col.names=NA)

# Build description START
desc_perseus <- annotRows(mdata)

# rename the column name to Grouping
colnames(desc_perseus) <- "Grouping"

# make the current rownames (protein IDs) to an actual first column
desc_perseus <- cbind(rownames(desc_perseus), data.frame(desc_perseus, row.names=NULL))

# Add samples
count <- c(1:nrow(desc_perseus))
desc_perseus$sample <- count

# Get batch groupings
vec_grouping <- vector()
for(entry in desc_perseus$Grouping){
  new_entry <- strtoi(substr(entry, nchar(entry), nchar(entry)))
  vec_grouping <- c(vec_grouping, new_entry)
}
desc_perseus$Grouping <- vec_grouping

desc_new <- desc_perseus[,c("rownames(desc_perseus)","sample","Grouping")]


# 2. Run any kind of analysis on the extracted data. --------------------------
library(HarmonizR)

outharm <- harmonizR(mainMatrix, desc_new, algorithm=algorithm_chosen, ComBat_mode=ComBat_mode_chosen)


# 2.1 Can we add the rownames? 
outharm_mod <- cbind(rownames(outharm), data.frame(outharm, row.names=NULL))
colnames(outharm_mod)[1] <- "rowsuniquenamejustincase"
names_rows <- data.frame(rownames = outharm_mod$rowsuniquenamejustincase)


# 3. Create a matrixData object which can be conveniently written to file
# in the Perseus txt format. --------------------------------------------------
outMdata <- matrixData(main=outharm, annotCols=names_rows, annotRows=annotRows(mdata))
write.perseus(outMdata, outFile)