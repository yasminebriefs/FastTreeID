test_that ("prime generation using gmp works", {
	skip_if_not_installed("gmp")

	directed <- c(1, 2, 3, 4)

	bidirected <- matrix(c(1,2 , 1,3 , 1,4 , 1,5), byrow = TRUE, ncol = 2)

	result1 <- fasttreeid_identify(bidirected, directed)
	result2 <- fasttreeid_identify(bidirected, directed, seed=result1$seed)
	result3 <- fasttreeid_identify(bidirected, directed, seed=result1$seed)
	result4 <- fasttreeid_identify(bidirected, directed, seed=result1$seed, prime=result1$prime)

	expect_equal(result1, result2)
	expect_equal(result2, result3)
	expect_equal(result3, result4)
})
