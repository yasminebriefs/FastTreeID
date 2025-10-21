test_that ("identification by missing edges to the root works", {
	directed <- c(1, 2, 3)

	bidirected <- matrix(c(1,3 , 1,4 , 2,4), byrow = TRUE, ncol = 2)

	result <- fasttreeid_identify(bidirected, directed, "16630300901090146909", "533993101691002181")

	expect_equal(result$seed, "16630300901090146909")
	expect_equal(result$prime, "533993101691002181")

	expect_equal(result$identification[[1]]$identifiability, 0)
	expect_equal(result$identification[[2]]$identifiability, 1)
	expect_equal(result$identification[[3]]$identifiability, 1)
	expect_equal(result$identification[[4]]$identifiability, 1)

	expect_equal(result$identification[[2]]$type, "fraction")
	expect_equal(result$identification[[2]]$numerator$what, "sigma")
	expect_equal(result$identification[[2]]$numerator$i, 1)
	expect_equal(result$identification[[2]]$numerator$j, 2)
	expect_equal(result$identification[[2]]$denominator$i, 1)
	expect_equal(result$identification[[2]]$denominator$j, 1)

	expect_equal(result$identification[[3]]$type, "path")
	expect_equal(result$identification[[3]]$nodes, c(3, 2))

	expect_equal(result$identification[[4]]$type, "path")
	expect_equal(result$identification[[4]]$nodes, c(4, 3, 2))
})
