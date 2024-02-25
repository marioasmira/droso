# combine_logs.R

#' Combine log files into a single file and delete old files
#'
#' @param log_directory directory where to find the log files
#'
#' @return Nothing
#' @export
combine_logs <- function(log_directory) {
  # Get the names of the log files
  log_files <- list.files(path = log_directory, pattern = "log_")
  # Combine the log files
  log <- unlist(lapply(paste0(log_directory, log_files), readLines))
  # Write the combined log to a file
  writeLines(log,
             con = paste0(
               log_directory,
               as.character(Sys.Date()),
               "_combined_log.txt"
             ))

  # Delete the separate log files
  file.remove(paste0(log_directory, log_files))
}