library(jsonlite)

AI_SNAPSHOT_TTL_DAYS <- 14
AI_SNAPSHOT_MAX_COUNT <- 100

AI_SNAPSHOT_DIR <- normalizePath(
  file.path(getwd(), "ai_snapshots"),
  winslash = "/",
  mustWork = FALSE
)

dir.create(AI_SNAPSHOT_DIR, recursive = TRUE, showWarnings = FALSE)

shiny::addResourcePath("ai_snapshots", AI_SNAPSHOT_DIR)

cleanup_ai_snapshots <- function() {
  files <- list.files(AI_SNAPSHOT_DIR, "\\.json$", full.names = TRUE)
  if (!length(files)) return()

  info <- file.info(files)
  age <- as.numeric(difftime(Sys.time(), info$mtime, units = "days"))

  old <- files[age > AI_SNAPSHOT_TTL_DAYS]
  if (length(old)) file.remove(old)

  files <- list.files(AI_SNAPSHOT_DIR, "\\.json$", full.names = TRUE)
  if (length(files) > AI_SNAPSHOT_MAX_COUNT) {
    info <- file.info(files)
    files <- files[order(info$mtime)]
    file.remove(head(files, length(files) - AI_SNAPSHOT_MAX_COUNT))
  }
}

save_ai_snapshot <- function(payload) {
  cleanup_ai_snapshots()

  id <- paste0(sample(c(letters, LETTERS, 0:9), 32, TRUE), collapse = "")

  payload$snapshot <- list(
    id = id,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    expires_after_days = AI_SNAPSHOT_TTL_DAYS
  )

  jsonlite::write_json(
    payload,
    file.path(AI_SNAPSHOT_DIR, paste0(id, ".json")),
    auto_unbox = TRUE,
    pretty = FALSE,
    na = "null",
    null = "null",
    digits = 10
  )

  id
}

records <- function(df) {
  if (is.null(df) || !nrow(df)) return(list())

  df[] <- lapply(df, function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.numeric(x)) x[!is.finite(x)] <- NA
    x
  })

  unname(split(df, seq_len(nrow(df))))
}

cleanup_ai_snapshots()