# Open a connection to a log file
logfile <- file("$LOGFILE", open = "a")
# Redirect both output and messages to the file and console
sink(logfile, append = TRUE, split = TRUE)

