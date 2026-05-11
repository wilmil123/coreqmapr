# Function which prints a message using shell echo; useful for printing messages from inside mclapply when running in Rstudio
message_parallel <- function(...) {
  system(sprintf('echo "%s"', paste0(..., collapse = "")))
}
