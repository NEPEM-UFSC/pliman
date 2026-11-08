#' Cryptographic SHA-256 and Base64 Utilities (Native C++)
#'
#' Fast, zero-dependency implementations of SHA-256 hashing and Base64
#' encoding/decoding implemented natively in C++ via Rcpp.
#'
#' @name utils_crypto
#' @rdname utils_crypto
NULL

#' Compute Cryptographic SHA-256 Checksum
#'
#' Computes the cryptographic SHA-256 hash of a file, raw bytes, or character
#' string. Files are processed in 64 KB binary streaming chunks to efficiently
#' handle large models (e.g., 500+ MB) with negligible memory footprint.
#'
#' @param x Character string (file path or text) or a `raw` vector of bytes.
#' @param file Logical. If `TRUE`, treats `x` as a file path. If `FALSE`, hashes
#'   the contents of `x`. Defaults to `file.exists(x[1])` when `x` is a single string.
#' @return A character string containing the 64-character lowercase hexadecimal hash.
#' @export
#' @examples
#' \dontrun{
#' # Hash a character string
#' sha256("Hello, pliman!")
#'
#' # Hash a local model file
#' sha256("model.onnx", file = TRUE)
#' }
sha256 <- function(x, file = NULL) {
  if (is.raw(x)) {
    return(rcpp_sha256(x))
  }
  if (!is.character(x) || length(x) == 0L) {
    cli::cli_abort("{.arg x} must be a character vector or raw bytes.")
  }
  is_file <- if (!is.null(file)) isTRUE(file) else (length(x) == 1L && file.exists(x))
  if (is_file) {
    rcpp_sha256_file(normalizePath(x[1], winslash = "/", mustWork = TRUE))
  } else {
    rcpp_sha256(charToRaw(paste(x, collapse = "\n")))
  }
}

#' Encode to Base64
#'
#' Encodes raw bytes, text, or a binary file into a Base64 string.
#'
#' @param x Character string (file path or text) or a `raw` vector of bytes.
#' @param file Logical. If `TRUE`, treats `x` as a file path and encodes its binary
#'   contents directly. Defaults to `file.exists(x[1])` when `x` is a single string.
#' @return A Base64-encoded character string.
#' @export
#' @examples
#' \dontrun{
#' # Encode raw bytes
#' base64_encode(as.raw(c(72, 101, 108, 108, 111)))
#'
#' # Encode an image file
#' b64 <- base64_encode("image.jpg", file = TRUE)
#' }
base64_encode <- function(x, file = NULL) {
  if (is.raw(x)) {
    return(rcpp_base64_encode(x))
  }
  if (!is.character(x) || length(x) == 0L) {
    cli::cli_abort("{.arg x} must be a character vector or raw bytes.")
  }
  is_file <- if (!is.null(file)) isTRUE(file) else (length(x) == 1L && file.exists(x))
  if (is_file) {
    rcpp_base64_encode_file(normalizePath(x[1], winslash = "/", mustWork = TRUE))
  } else {
    rcpp_base64_encode(charToRaw(paste(x, collapse = "\n")))
  }
}

#' Decode from Base64
#'
#' Decodes a Base64 character string back into a `raw` vector of bytes.
#'
#' @param text Base64-encoded character string.
#' @param as_char Logical. If `TRUE`, converts decoded bytes to a UTF-8 character string.
#'   Defaults to `FALSE` (returns a `raw` vector).
#' @return A `raw` vector of bytes, or a character string if `as_char = TRUE`.
#' @export
#' @examples
#' \dontrun{
#' raw_bytes <- base64_decode("SGVsbG8sIHBsZW1hbiE=")
#' rawToChar(raw_bytes)
#' }
base64_decode <- function(text, as_char = FALSE) {
  if (!is.character(text) || length(text) == 0L) {
    cli::cli_abort("{.arg text} must be a character string.")
  }
  out <- rcpp_base64_decode(paste(text, collapse = ""))
  if (isTRUE(as_char)) {
    rawToChar(out)
  } else {
    out
  }
}
