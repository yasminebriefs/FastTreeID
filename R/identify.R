#' Identify parameters in tree-shaped SCMs in polynomial time
#' @export
fasttreeid_identify <- function(bidirected, directed, seed = NULL, prime = NULL) {
	.fasttreeid_identify_bridge(bidirected, directed, seed, prime)
}
