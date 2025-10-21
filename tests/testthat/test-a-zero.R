test_that ("identification works when a is zero", {
	directed <- c(1, 2, 3, 4, 5, 6, 7)

	bidirected <- matrix(c(1,2 , 1,3 , 1,4 , 1,5 , 1,6 , 1,7 , 1,8 , 2,3 , 2,4 , 2,5 , 2,8 , 3,4 , 3,5 , 3,7 , 3,8 , 4,5 , 4,6 , 4,8 , 5,6 , 5,7 , 6,7), byrow = TRUE, ncol = 2)

	result <- fasttreeid_identify(bidirected, directed, "16630300901090146909", "533993101691002181")

	expect_equal(result$seed, "16630300901090146909")
	expect_equal(result$prime, "533993101691002181")

	expect_equal(result$identification[[1]]$identifiability, 0)
	expect_equal(result$identification[[2]]$identifiability, 1)
	expect_equal(result$identification[[3]]$identifiability, 1)
	expect_equal(result$identification[[4]]$identifiability, 1)
	expect_equal(result$identification[[5]]$identifiability, 1)
	expect_equal(result$identification[[6]]$identifiability, 1)
	expect_equal(result$identification[[7]]$identifiability, 1)
	expect_equal(result$identification[[8]]$identifiability, 1)

	expect_equal(result$identification[[2]]$type, "cycle")
	expect_equal(result$identification[[2]]$nodes, c(2, 6, 8, 7, 2))
	expect_equal(result$identification[[2]]$reason, "a_is_zero")
	expect_equal(result$identification[[2]]$reason_edge$what, "lambda")
	expect_equal(result$identification[[2]]$reason_edge$i, 6)
	expect_equal(result$identification[[2]]$reason_edge$j, 7)

	expect_equal(result$identification[[3]]$type, "path")
	expect_equal(result$identification[[3]]$nodes, c(3, 6, 2))

	expect_equal(result$identification[[4]]$type, "path")
	expect_equal(result$identification[[4]]$nodes, c(4, 7, 2))

	expect_equal(result$identification[[5]]$type, "path")
	expect_equal(result$identification[[5]]$nodes, c(5, 8, 6, 2))

	expect_equal(result$identification[[6]]$type, "path")
	expect_equal(result$identification[[6]]$nodes, c(6, 2))

	expect_equal(result$identification[[7]]$type, "path")
	expect_equal(result$identification[[7]]$nodes, c(7, 2))

	expect_equal(result$identification[[8]]$type, "path")
	expect_equal(result$identification[[8]]$nodes, c(8, 6, 2))
})
