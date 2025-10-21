test_that ("identification by cycle with only one option works", {
	directed <- c(1, 2, 2, 1, 2, 4, 7, 8, 5)

	bidirected <- matrix(c(1,2 , 1,3 , 1,4 , 1,5 , 1,6 , 1,7 , 1,8 , 1,9 , 1,10 , 3,6 , 2,5 , 9,10 , 8,9 , 5,7 , 2,6 , 3,5), byrow = TRUE, ncol = 2)

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
	expect_equal(result$identification[[9]]$identifiability, 1)
	expect_equal(result$identification[[10]]$identifiability, 1)

	expect_equal(result$identification[[2]]$type, "cycle")
	expect_equal(result$identification[[2]]$nodes, c(2, 3, 4, 2))
	expect_equal(result$identification[[2]]$reason, "only_one_option")
	expect_equal(result$identification[[2]]$reason_edge$what, "missing_bidirected")
	expect_equal(result$identification[[2]]$reason_edge$i, 3)
	expect_equal(result$identification[[2]]$reason_edge$j, 4)

	expect_equal(result$identification[[3]]$type, "path")
	expect_equal(result$identification[[3]]$nodes, c(3, 2))

	expect_equal(result$identification[[4]]$type, "path")
	expect_equal(result$identification[[4]]$nodes, c(4, 2))

	expect_equal(result$identification[[5]]$type, "path")
	expect_equal(result$identification[[5]]$nodes, c(5, 4, 2))

	expect_equal(result$identification[[6]]$type, "path")
	expect_equal(result$identification[[6]]$nodes, c(6, 4, 2))

	expect_equal(result$identification[[7]]$type, "path")
	expect_equal(result$identification[[7]]$nodes, c(7, 2))

	expect_equal(result$identification[[8]]$type, "path")
	expect_equal(result$identification[[8]]$nodes, c(8, 2))

	expect_equal(result$identification[[9]]$type, "path")
	expect_equal(result$identification[[9]]$nodes, c(9, 2))

	expect_equal(result$identification[[10]]$type, "path")
	expect_equal(result$identification[[10]]$nodes, c(10, 2))
})
