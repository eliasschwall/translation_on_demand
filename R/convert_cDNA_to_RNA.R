#' convert_cDNA_to_RNA
#' @author Elias Schwall
#' 
#' @importFrom stringr str_replace_all
#' @param cDNA \code{string} of your cDNA sequence
#'
#' @return RNA \code{string} RNA sequence
#' @export
#'
#' @examples convert_cDNA_to_RNA("AACTG")
convert_cDNA_to_RNA <- function(cDNA) {
  return(stringr::str_replace_all(cDNA, "T", "U"))
}
