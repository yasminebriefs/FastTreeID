test_that ("identification with multiple components works", {
	directed <- c(1, 2, 2, 2, 2, 5, 4, 3, 2)

	bidirected <- matrix(c(4,9 , 3,6 , 2,3 , 1,7 , 6,10 , 1,4 , 9,10 , 5,6 , 5,10 , 2,7 , 1,2 , 8,10 , 1,9 , 1,8 , 1,5 , 2,8 , 8,9 , 2,10 , 5,8 , 1,3 , 7,10 , 4,7 , 3,9 , 4,10 , 6,9 , 3,4 , 6,8 , 3,7 , 4,5 , 7,9 , 1,6 , 5,7 , 2,4 , 4,6 , 7,8 , 3,5), byrow = TRUE, ncol = 2)

	result <- fasttreeid_identify(bidirected, directed, "16630300901090146909", "533993101691002181")

	expect_equal(result$seed, "16630300901090146909")
	expect_equal(result$prime, "533993101691002181")

	expect_equal(result$identification[[1]]$identifiability, 0)
	expect_equal(result$identification[[2]]$identifiability, 2)
	expect_equal(result$identification[[3]]$identifiability, 1)
	expect_equal(result$identification[[4]]$identifiability, 1)
	expect_equal(result$identification[[5]]$identifiability, 2)
	expect_equal(result$identification[[6]]$identifiability, 2)
	expect_equal(result$identification[[7]]$identifiability, 2)
	expect_equal(result$identification[[8]]$identifiability, 1)
	expect_equal(result$identification[[9]]$identifiability, 2)
	expect_equal(result$identification[[10]]$identifiability, 1)

	expect_equal(result$identification[[2]]$type, "cycle")
	expect_equal(result$identification[[2]]$nodes, c(2, 5, 9, 2))

	expect_equal(result$identification[[3]]$type, "path")
	expect_equal(result$identification[[3]]$nodes, c(3, 10))

	expect_equal(result$identification[[4]]$type, "path")
	expect_equal(result$identification[[4]]$nodes, c(4, 8, 3, 10))

	expect_equal(result$identification[[5]]$type, "path")
	expect_equal(result$identification[[5]]$nodes, c(5, 2))

	expect_equal(result$identification[[6]]$type, "path")
	expect_equal(result$identification[[6]]$nodes, c(6, 2))

	expect_equal(result$identification[[7]]$type, "path")
	expect_equal(result$identification[[7]]$nodes, c(7, 6, 2))

	expect_equal(result$identification[[8]]$type, "path")
	expect_equal(result$identification[[8]]$nodes, c(8, 3, 10))

	expect_equal(result$identification[[9]]$type, "path")
	expect_equal(result$identification[[9]]$nodes, c(9, 2))

	expect_equal(result$identification[[10]]$type, "fraction")
	expect_equal(result$identification[[10]]$numerator$what, "sigma")
	expect_equal(result$identification[[10]]$numerator$i, 1)
	expect_equal(result$identification[[10]]$numerator$j, 10)
	expect_equal(result$identification[[10]]$denominator$i, 1)
	expect_equal(result$identification[[10]]$denominator$j, 2)
})
