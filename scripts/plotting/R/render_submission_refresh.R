render_submission_refresh_script_path <- function() {
  frame_files <- vapply(sys.frames(), function(frame) frame$ofile %||% "", character(1))
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0) {
    return(normalizePath(frame_files[[length(frame_files)]], winslash = "/", mustWork = TRUE))
  }
  cmd_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(cmd_file) == 0) {
    stop("Could not resolve render_submission_refresh.R path.", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", cmd_file[[1]]), winslash = "/", mustWork = TRUE)
}

`%||%` <- function(left, right) {
  if (is.null(left)) right else left
}

render_submission_refresh_root <- function() {
  normalizePath(file.path(dirname(render_submission_refresh_script_path()), "..", ".."), winslash = "/", mustWork = TRUE)
}

repo_root <- render_submission_refresh_root()
old_wd <- getwd()
setwd(repo_root)
on.exit(setwd(old_wd), add = TRUE)
Sys.setenv(M2M_REPO_ROOT = repo_root)

suppressPackageStartupMessages({
  source(file.path(repo_root, "plotting", "R", "submission_refresh", "figure1Refresh.R"))
  source(file.path(repo_root, "plotting", "R", "submission_refresh", "figure2Refresh.R"))
  source(file.path(repo_root, "plotting", "R", "submission_refresh", "figure3Refresh.R"))
})

rendered_paths <- c(
  render_submission_refresh_figure1(),
  render_submission_refresh_figure2(),
  render_submission_refresh_figure3()
)

for (path in rendered_paths) {
  cat(normalizePath(path, winslash = "/", mustWork = FALSE), "\n", sep = "")
}
