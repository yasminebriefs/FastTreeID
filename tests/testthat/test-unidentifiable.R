test_that ("unidentifiable examples work", {
	directed <- c(1, 2, 3, 4, 5, 6, 7)

	bidirected <- matrix(c(1,2 , 1,3 , 1,4 , 1,5 , 1,6 , 1,7 , 1,8 , 2,3 , 2,4 , 2,7 , 2,8 , 3,4 , 3,7 , 3,8 , 4,5 , 4,6 , 4,7 , 4,8 , 5,6 , 5,7 , 5,8 , 6,7 , 6,8 , 7,8), byrow = TRUE, ncol = 2)

	result <- fasttreeid_identify(bidirected, directed, "16630300901090146909", "533993101691002181")

	expect_equal(result$seed, "16630300901090146909")
	expect_equal(result$prime, "533993101691002181")

	expect_equal(result$identification[[1]]$identifiability, 0)
	expect_equal(result$identification[[2]]$identifiability, 0)
	expect_equal(result$identification[[3]]$identifiability, 0)
	expect_equal(result$identification[[4]]$identifiability, 0)
	expect_equal(result$identification[[5]]$identifiability, 0)
	expect_equal(result$identification[[6]]$identifiability, 0)
	expect_equal(result$identification[[7]]$identifiability, 0)
	expect_equal(result$identification[[8]]$identifiability, 0)
})
