library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript 01-process_getMean.R <input_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]

process_file <- function(input_file, output_file) {

  data <- read.table(input_file, header = FALSE, sep = "\t", fill = TRUE, skip = 1)

  signal_cols <- data[, 7:ncol(data)]

  data$signal_sum <- rowSums(signal_cols, na.rm = TRUE)

  result <- data[, c("V4", "signal_sum")]

  colnames(result) <- c("hyperDHMR", "signal_sum")
  
  result$signal_sum_log10 <- log10(result$signal_sum + 1)

  write.table(result, file = output_file, row.names = FALSE, sep = "\t", quote = FALSE)
}

process_file(input_file, output_file)

cat("Processed:", input_file, "->", output_file, "\n")
