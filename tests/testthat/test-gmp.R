test_that ("prime generation using gmp works", {
	skip_if_not_installed("gmp")

	directed <- c(1, 2, 3, 4)

	bidirected <- matrix(c(1,2 , 1,3, 1,4, 1,5), byrow = TRUE, ncol = 2)

	result <- fasttreeid_identify(bidirected, directed, "738687587898544965")

	expect_equal(result$prime, "335270452938866563")

	expect_no_error(fasttreeid_identify(bidirected, directed))
})
